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
ConstraintMatrix
form_constraints(const BlockDOFMap& dof_map,
    const std::vector<IntegralEquationSpec<dim>>& int_eqtns,
    const MeshMap<dim>& meshes)
{
    std::map<std::string,size_t> eqtn_component_map;
    for (size_t i = 0; i < int_eqtns.size(); i++) {
        eqtn_component_map[int_eqtns[i].obs_mesh()] = i * dim;
    }

    std::vector<ConstraintEQ> out;
    for (size_t i = 0; i < int_eqtns.size(); i++) {
        for (size_t d = 0; d < dim; d++) {
            auto constraints = 
                int_eqtns[i].constraint_builder(eqtn_component_map, dof_map, meshes, d);
            for (size_t c_idx = 0; c_idx < constraints.size(); c_idx++) {
                out.push_back(constraints[c_idx]);
            }
        }
    }
    return from_constraints(out);
}

template <size_t dim>
BlockDOFMap form_dof_map(const std::vector<IntegralEquationSpec<dim>>& int_eqtns,
    const MeshMap<dim>& meshes) 
{
    std::vector<size_t> components;
    for (size_t i = 0; i < int_eqtns.size(); i++) {
        const auto& obs_mesh = meshes.at(int_eqtns[i].obs_mesh());
        for (size_t d = 0; d < dim; d++) {
            components.push_back(obs_mesh.n_dofs());
        }
    }

    return build_block_dof_map(components);
}

std::vector<double>
concatenate_condense(const BlockDOFMap& dof_map,
    const ConstraintMatrix& cm,
    const std::vector<BlockFunction>& fncs) 
{
    BlockFunction flattened;
    //TODO: This should be avoidable.
    for (size_t i = 0; i < fncs.size(); i++) {
        for (size_t d = 0; d < fncs[i].size(); d++) {
            flattened.push_back(fncs[i][d]);
        }
    }

    return condense_vector(cm, concatenate(dof_map, flattened));
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
    auto dof_map = form_dof_map(bem_input.eqtn_specs, bem_input.meshes);
    auto constraint_matrix = form_constraints(dof_map, bem_input.eqtn_specs, bem_input.meshes);

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
        assert(diag_idx[i] != -1);
        precondition(rhs[i], systems[i].lhs[diag_idx[i]].op);
    }

    //reduce using constraints
    //stack rows
    auto stacked_rhs = concatenate_condense(dof_map, constraint_matrix, rhs);

    //solve:
    int count = 0;
    auto reduced_soln = solve_system(stacked_rhs, bem_input.params.solver_tol,
        [&] (std::vector<double>& x, std::vector<double>& y) {
            std::cout << "iteration " << count << std::endl;
            count++;

            //unstack_rows
            auto distributed = distribute_vector(constraint_matrix, x, dof_map.n_dofs);
            auto unstacked_estimate = expand(dof_map, distributed);

            //expand using constraints
            //TODO: better name than BCMap, FunctionMap?
            BCMap unknowns; 
            for (size_t i = 0; i < n_eqtns; i++) {
                BlockFunction unknown_field(dim);
                for (size_t d = 0; d < dim; d++) {
                    unknown_field[d] = unstacked_estimate[i * dim + d];
                }
                auto src_mesh = bem_input.eqtn_specs[i].obs_mesh();
                auto field_name = bem_input.eqtn_specs[i].unknown_field;
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
            auto evaluated_lhs = concatenate_condense(dof_map, constraint_matrix, eval);

            for (size_t i = 0; i < y.size(); i++) {
                y[i] = evaluated_lhs[i];
            }
        }
    );
    
    //handle soln:
    //unstack rows
    auto full_soln = distribute_vector(constraint_matrix, reduced_soln, dof_map.n_dofs);
    auto unstacked_soln = expand(dof_map, full_soln);

    //expand using constraints
    //output
    for (size_t i = 0; i < n_eqtns; i++) {
        BlockFunction soln(dim);
        for (size_t d = 0; d < dim; d++) {
            soln[d] = unstacked_soln[i * dim + d];
        }
        auto obs_mesh = bem_input.eqtn_specs[i].obs_mesh();
        if (soln[0].size() > 0) {
            HDFOutputter file(bem_input.eqtn_specs[i].get_output_filename(filename));
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
