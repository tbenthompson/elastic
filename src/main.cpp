#include "spec.h"
#include "load.h"
#include "compute.h"
#include "filenames.h"
#include "data.h"
#include "3bem/3bem.h"
#include "3bem/function.h"
#include "3bem/armadillo_facade.h"

using namespace tbem;

template <size_t dim>
std::vector<std::vector<ConstraintMatrix>>
form_constraints(const std::vector<IntegralEquationSpec<dim>>& int_eqtns,
    const MeshMap<dim>& meshes)
{
    std::vector<std::vector<ConstraintMatrix>> out(int_eqtns.size());
    for (size_t i = 0; i < int_eqtns.size(); i++) {
        for (size_t d = 0; d < dim; d++) {
            out[i].push_back(int_eqtns[i].constraint_builder(meshes, d));
        }
    }
    return out;
}

ConcatenatedFunction
concatenate_condense(const std::vector<std::vector<ConstraintMatrix>>& cms,
    const std::vector<BlockFunction>& fncs) 
{
    std::vector<Function> condensed;

    for (size_t eqtn_idx = 0; eqtn_idx < fncs.size(); eqtn_idx++) {
        const auto& f = fncs[eqtn_idx];
        for (size_t d = 0; d < f.size(); d++) {
            condensed.push_back(condense_vector(cms[eqtn_idx][d], f[d]));
        }
    }
    return concatenate(condensed);
}


double condition_number(const LinearSystem& disp_system, 
    const LinearSystem& trac_system,
    const std::vector<ConstraintMatrix>& constraint_matrices)
{
    BlockOperator lhs{
        4, 4,
        {
            disp_system.lhs[1].op.ops[0],
            disp_system.lhs[1].op.ops[1],
            disp_system.lhs[0].op.ops[0],
            disp_system.lhs[0].op.ops[1],

            disp_system.lhs[1].op.ops[2],
            disp_system.lhs[1].op.ops[3],
            disp_system.lhs[0].op.ops[2],
            disp_system.lhs[0].op.ops[3],

            trac_system.lhs[1].op.ops[0],
            trac_system.lhs[1].op.ops[1],
            trac_system.lhs[0].op.ops[0],
            trac_system.lhs[0].op.ops[1],

            trac_system.lhs[1].op.ops[2],
            trac_system.lhs[1].op.ops[3],
            trac_system.lhs[0].op.ops[2],
            trac_system.lhs[0].op.ops[3],
        }
    };

    auto condensed_lhs = condense_block_operator(
        constraint_matrices, constraint_matrices, lhs);
    auto combined_lhs = combine_components(condensed_lhs);
    return arma_cond(combined_lhs.ops[0]);
}

/* A jacobi diagonal preconditioner. Super simple.
*/
void precondition(BlockFunction& f, const BlockOperator& block_op) 
{
    for (size_t d = 0; d < 2; d++) {
        const auto& op = block_op.ops[d * 2 + d];
        assert(op.n_rows == op.n_cols);
        for (size_t i = 0; i < f[d].size(); i++) {
            f[d][i] /= op.data[i * op.n_cols + i];
        }
    }
}

size_t find_diagonal_index(const LinearSystem& system, const std::string& obs_mesh) {
    for (size_t op_idx = 0; op_idx < system.lhs.size(); op_idx++) {
        if (system.lhs[op_idx].src_mesh == obs_mesh) {
            return op_idx;
        }
    }
    return -1;
}

template <size_t dim>
void solve(const std::string& filename) 
{
    auto bem_input = parse_into_bem<dim>(filename);
    auto constraint_matrices = form_constraints(bem_input.eqtn_specs, bem_input.meshes);

    size_t n_eqtns = bem_input.eqtn_specs.size();
    std::vector<LinearSystem> systems(n_eqtns);
    std::vector<BlockFunction> rhs(n_eqtns);
    std::vector<int> diag_idx(n_eqtns);
    for (size_t i = 0; i < n_eqtns; i++) {
        auto computed = compute_integral_equation(bem_input, bem_input.eqtn_specs[i]);
        systems[i] = separate(computed, bem_input.bcs),
        rhs[i] = systems[i].rhs;
        diag_idx[i] = find_diagonal_index(systems[i],
            bem_input.eqtn_specs[i].obs_mesh());
        precondition(rhs[i], systems[i].lhs[diag_idx[i]].op);
    }

    //reduce using constraints
    //stack rows
    auto stacked_rhs = concatenate_condense(constraint_matrices, rhs);

    //solve:
    int count = 0;
    auto reduced_soln = solve_system(stacked_rhs.data, bem_input.params.solver_tol,
        [&] (std::vector<double>& x, std::vector<double>& y) {
            std::cout << "iteration " << count << std::endl;
            count++;

            //unstack_rows
            auto unstacked_estimate = expand(stacked_rhs, x);

            //expand using constraints
            //TODO: better name than BCMap, FunctionMap?
            BCMap unknowns; 
            for (size_t i = 0; i < n_eqtns; i++) {
                BlockFunction unknown_field(dim);
                for (size_t d = 0; d < dim; d++) {
                    const auto& cm = constraint_matrices[i][d];
                    const auto& reduced_data = unstacked_estimate[i * dim + d];
                    auto n_dofs = systems[i].rhs[d].size();
                    unknown_field[d] = distribute_vector(cm, reduced_data, n_dofs);
                }
                auto src_mesh = bem_input.eqtn_specs[i].obs_mesh();
                auto field_name = systems[i].lhs[diag_idx[i]].function;
                unknowns[FieldDescriptor{src_mesh, field_name}] = unknown_field;
            }

            //apply operators
            std::vector<BlockFunction> eval(n_eqtns);
            for (size_t i = 0; i < n_eqtns; i++) {
                //TODO: If separate is split into a "evaluate_possible" function, these
                //two lines would be unnecessary
                auto fully_evaluated_system = separate(systems[i].lhs, unknowns);
                assert(fully_evaluated_system.lhs.size() == 0);
                auto eval_vec = -fully_evaluated_system.rhs;
                precondition(eval_vec, systems[i].lhs[diag_idx[i]].op);
                eval[i] = eval_vec;
            }
            //scale rows

            //reduce using constraints
            //stack rows
            auto evaluated_lhs = concatenate_condense(constraint_matrices, eval);

            for (size_t comp = 0; comp < evaluated_lhs.components; comp++) {
                assert(evaluated_lhs.component_lengths[comp] == 
                    stacked_rhs.component_lengths[comp]);
            }
            assert(y.size() == evaluated_lhs.data.size());

            for (size_t i = 0; i < y.size(); i++) {
                y[i] = evaluated_lhs.data[i];
            }
        }
    );
    
    //handle soln:
    //unstack rows
    auto unstacked_soln = expand(stacked_rhs, reduced_soln);

    //expand using constraints
    //output
    for (size_t i = 0; i < n_eqtns; i++) {
        BlockFunction soln(dim);
        for (size_t d = 0; d < dim; d++) {
            const auto& cm = constraint_matrices[i][d];
            const auto& reduced_data = unstacked_soln[i * dim + d];
            auto n_dofs = systems[i].rhs[d].size();
            soln[d] = distribute_vector(cm, reduced_data, n_dofs);
        }
        HDFOutputter file(bem_input.eqtn_specs[i].get_output_filename(filename));
        auto obs_mesh = bem_input.eqtn_specs[i].obs_mesh();
        if (soln[0].size() > 0) {
            out_surface(file, bem_input.meshes.at(obs_mesh), soln);
        }
    }
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage is 'solve filename'" << std::endl;
        return 1;
    }
    
    auto filename = argv[1];
    solve<2>(filename);
}
