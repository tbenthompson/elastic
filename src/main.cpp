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
std::vector<ConstraintMatrix> 
form_constraints(const std::vector<IntegralEquationSpec<dim>>& int_eqtns,
    const MeshMap<dim>& meshes)
{
    std::vector<ConstraintMatrix> out;
    for (const auto& ie: int_eqtns) {
        for (size_t d = 0; d < dim; d++) {
            out.push_back(ie.constraint_builder(meshes, d));
        }
    }
    return out;
}

ConcatenatedFunction concatenate_condense(const std::vector<ConstraintMatrix>& cms,
    const std::vector<BlockFunction>& fncs) 
{
    std::vector<Function> condensed;

    size_t constraint_matrix_idx = 0;
    for (size_t eqtn_idx = 0; eqtn_idx < fncs.size(); eqtn_idx++) {
        const auto& f = fncs[eqtn_idx];
        for (size_t d = 0; d < f.size(); d++) {
            condensed.push_back(condense_vector(cms[constraint_matrix_idx], f[d]));
            constraint_matrix_idx++;
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

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage is 'solve filename'" << std::endl;
        return 1;
    }
    
    auto filename = argv[1];
    auto bem_input = parse_into_bem<2>(filename);
    auto constraint_matrices = form_constraints(bem_input.eqtn_specs, bem_input.meshes);

    auto disp_system = separate(compute_integral_equation(bem_input, bem_input.eqtn_specs[0]), bem_input.bcs);
    auto trac_system = separate(compute_integral_equation(bem_input, bem_input.eqtn_specs[1]), bem_input.bcs);

    auto disp_rhs = disp_system.rhs;
    auto trac_rhs = trac_system.rhs;
    precondition(disp_rhs, disp_system.lhs[1].op);
    precondition(trac_rhs, trac_system.lhs[0].op);

    //reduce using constraints
    //stack rows
    auto stacked_rhs = concatenate_condense(
        constraint_matrices, {disp_rhs, trac_rhs}
    );

    auto n_unknown_trac_dofs = disp_system.rhs[0].size();
    auto n_unknown_disp_dofs = trac_system.rhs[0].size();

    //solve:
    int count = 0;
    auto reduced_soln = solve_system(stacked_rhs.data, bem_input.params.solver_tol,
        [&] (std::vector<double>& x, std::vector<double>& y) {
            std::cout << "iteration " << count << std::endl;
            count++;

            //unstack_rows
            auto unstacked_estimate = expand(stacked_rhs, x);

            //expand using constraints
            BlockFunction unknown_trac = {
                distribute_vector(constraint_matrices[0], unstacked_estimate[0],
                                                 n_unknown_trac_dofs),
                distribute_vector(constraint_matrices[1], unstacked_estimate[1],
                                                 n_unknown_trac_dofs)
            };
            BlockFunction unknown_disp = {
                distribute_vector(constraint_matrices[2], unstacked_estimate[2],
                                                 n_unknown_disp_dofs),
                distribute_vector(constraint_matrices[3], unstacked_estimate[3],
                                                 n_unknown_disp_dofs)
            };

            //apply operators
            //TODO: better name than BCMap, FunctionMap?
            BCMap unknowns; 
            unknowns[FieldDescriptor{"displacement", "traction"}] = unknown_trac;
            unknowns[FieldDescriptor{"traction", "displacement"}] = unknown_disp;
            auto disp_eval = separate(disp_system.lhs, unknowns);
            auto trac_eval = separate(trac_system.lhs, unknowns);
            assert(disp_eval.lhs.size() == 0);
            assert(trac_eval.lhs.size() == 0);
            //TODO: If separate is split into a "evaluate_possible" function, these
            //two lines would be unnecessary
            //scale rows
            auto disp_eval_vec = -disp_eval.rhs;
            auto trac_eval_vec = -trac_eval.rhs;
            precondition(disp_eval_vec, disp_system.lhs[1].op);
            precondition(trac_eval_vec, trac_system.lhs[0].op);

            //reduce using constraints
            //stack rows
            auto evaluated_lhs = concatenate({
                condense_vector(constraint_matrices[0], disp_eval_vec[0]),
                condense_vector(constraint_matrices[1], disp_eval_vec[1]),
                condense_vector(constraint_matrices[2], trac_eval_vec[0]),
                condense_vector(constraint_matrices[3], trac_eval_vec[1])
            });

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
    BlockFunction soln_trac = {
        distribute_vector(constraint_matrices[0], unstacked_soln[0],
                n_unknown_trac_dofs),
        distribute_vector(constraint_matrices[1], unstacked_soln[1],
                n_unknown_trac_dofs)
    };
    BlockFunction soln_disp = {
        distribute_vector(constraint_matrices[2], unstacked_soln[2],
                n_unknown_disp_dofs),
        distribute_vector(constraint_matrices[3], unstacked_soln[3],
                n_unknown_disp_dofs)
    };


    //output
    if (soln_disp[0].size() > 0) {
        HDFOutputter disp_file(disp_out_filename(filename));
        out_surface(disp_file, bem_input.meshes.at("traction"), soln_disp);
    }
    if (soln_trac[0].size() > 0) {
        HDFOutputter trac_file(trac_out_filename(filename));
        out_surface(trac_file, bem_input.meshes.at("displacement"), soln_trac);
    }
}
