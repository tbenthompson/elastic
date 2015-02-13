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


template <size_t dim>
std::vector<LinearSystem> compute_linear_systems(const BEM<dim>& bem_input) 
{
    std::vector<LinearSystem> systems;
    size_t n_eqtns = bem_input.eqtn_specs.size();
    for (size_t i = 0; i < n_eqtns; i++) {
        auto computed = compute_integral_equation(bem_input, bem_input.eqtn_specs[i]);
        systems.push_back(evaluate_computable_terms(computed, bem_input.bcs));
    }
    return systems;
}

/* A jacobi diagonal preconditioner. Super simple.
*/
void precondition(BlockFunction& f, const BlockOperator& block_op, size_t dim) 
{
    for (size_t d = 0; d < dim; d++) {
        const auto& op = block_op.ops[d * dim + d];
        assert(op.n_rows == op.n_cols);
        for (size_t i = 0; i < f[d].size(); i++) {
            f[d][i] /= op.data[i * op.n_cols + i];
        }
    }
}

std::vector<double>
preconditioned_rhs(const std::vector<LinearSystem>& systems,
    const ConstraintMatrix& constraint_matrix,
    const BlockDOFMap& dof_map, size_t dim)
{
    std::vector<BlockFunction> rhs(systems.size());
    for (size_t i = 0; i < systems.size(); i++) {
        rhs[i] = -systems[i].evaluated_terms;
        precondition(rhs[i], systems[i].get_diag_block().op, dim);
    }
    return concatenate_condense(dof_map, constraint_matrix, rhs);
}

template <size_t dim>
FunctionMap extract_solution_fields(const std::vector<Function>& concatenated_fields,
    const std::vector<IntegralEquationSpec<dim>>& eqtn_specs) 
{
    FunctionMap field_map;
    for (size_t i = 0; i < eqtn_specs.size(); i++) {
        BlockFunction field(dim);
        for (size_t d = 0; d < dim; d++) {
            field[d] = concatenated_fields[i * dim + d];
        }
        auto src_mesh = eqtn_specs[i].obs_mesh();
        auto field_name = eqtn_specs[i].unknown_field;
        field_map[FieldDescriptor{src_mesh, field_name}] = field;
    }
    return field_map;
}

std::vector<BlockFunction>
mat_vec_prod(const std::vector<LinearSystem>& systems, const FunctionMap& x, size_t dim)
{
    std::vector<BlockFunction> eval(systems.size());
    for (size_t i = 0; i < systems.size(); i++) {
        auto evaluated_system = evaluate_computable_terms(systems[i].lhs, x);
        assert(evaluated_system.lhs.size() == 0);
        auto eval_vec = evaluated_system.evaluated_terms;
        precondition(eval_vec, systems[i].get_diag_block().op, dim);
        eval[i] = eval_vec;
    }
    return eval;
}

template <size_t dim>
void output_solution(const std::vector<IntegralEquationSpec<dim>>& eqtn_specs,
                     std::string input_filename,
                     const MeshMap<dim>& meshes,
                     const std::vector<Function>& functions) {
    for (size_t i = 0; i < eqtn_specs.size(); i++) {
        BlockFunction soln(dim);
        for (size_t d = 0; d < dim; d++) {
            soln[d] = functions[i * dim + d];
        }
        auto obs_mesh = eqtn_specs[i].obs_mesh();
        if (soln[0].size() > 0) {
            auto out_filename = eqtn_specs[i].get_output_filename(input_filename);
            HDFOutputter file(out_filename);
            out_surface(file, meshes.at(obs_mesh), soln);
        }
    }
}

template <size_t dim>
void solve(const std::string& filename) 
{
    auto bem_input = parse_into_bem<dim>(filename);
    auto dof_map = form_dof_map(bem_input.eqtn_specs, bem_input.meshes);
    auto constraint_matrix = form_constraints(dof_map, bem_input.eqtn_specs, bem_input.meshes);

    auto systems = compute_linear_systems(bem_input);
    auto rhs = preconditioned_rhs(systems, constraint_matrix, dof_map, dim);

    int count = 0;
    auto reduced_soln = solve_system(rhs, bem_input.params.solver_tol,
        [&] (std::vector<double>& x, std::vector<double>& y) {
            std::cout << "iteration " << count << std::endl;
            count++;

            auto distributed = distribute_vector(constraint_matrix, x, dof_map.n_dofs);
            auto unstacked_soln = expand(dof_map, distributed);
            auto unknowns = extract_solution_fields(unstacked_soln, bem_input.eqtn_specs);

            auto eval = mat_vec_prod(systems, unknowns, dim);

            auto out = concatenate_condense(dof_map, constraint_matrix, eval);
            std::copy(out.begin(), out.end(), y.begin());
        }
    );
    
    auto full_soln = distribute_vector(constraint_matrix, reduced_soln, dof_map.n_dofs);
    auto unstacked_soln = expand(dof_map, full_soln);
    output_solution(bem_input.eqtn_specs, filename, bem_input.meshes, unstacked_soln); 
}
