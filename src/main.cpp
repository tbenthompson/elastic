#include "spec.h"
#include "load.h"
#include "compute.h"

using namespace tbem;

template <size_t dim>
ConstraintMatrix form_displacement_constraints(const MeshMap<dim>& meshes) 
{
    auto continuity = mesh_continuity(meshes.at("traction").begin());
    auto cut_continuity = cut_at_intersection(
        continuity, meshes.at("traction").begin(), meshes.at("slip").begin()
    );
    auto constraints = convert_to_constraints(cut_continuity);
    auto constraint_matrix = ConstraintMatrix::from_constraints(constraints);
    return constraint_matrix;
}

template <size_t dim>
BEM<dim> parse_into_bem(const std::string& filename)
{
    auto doc = parse_json(load_file(filename));
    auto params = get_parameters(doc);
    auto elements = get_elements<dim>(doc);
    auto meshes = get_meshes(elements);
    auto bcs = get_bcs(elements);

    return {
        params,
        meshes, 
        bcs,
        get_elastic_kernels<dim>(params.shear_modulus, params.poisson_ratio),
        QuadStrategy<dim>(params.obs_quad_order, params.src_far_quad_order,
            params.n_singular_steps, params.far_threshold, params.near_tol),
        {
            get_displacement_BIE(),
            get_traction_BIE()
        },
        {
            ConstraintMatrix::from_constraints({}),
            ConstraintMatrix::from_constraints({}),
            form_displacement_constraints(meshes),
            form_displacement_constraints(meshes)
        }
    };
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage is 'elastic_process filename'" << std::endl;
        return 1;
    }
    
    auto filename = argv[1];
    auto bem_input = parse_into_bem<2>(filename);

    auto disp_BIE_ops = compute_integral_equation(bem_input, bem_input.eqtn_specs[0]);
    auto trac_BIE_ops = compute_integral_equation(bem_input, bem_input.eqtn_specs[1]);
    assert(disp_BIE_ops.size() == 6);
    assert(trac_BIE_ops.size() == 6);

    auto disp_system = separate(disp_BIE_ops, bem_input.bcs);
    auto trac_system = separate(trac_BIE_ops, bem_input.bcs);
    assert(disp_system.lhs.size() == 2);
    assert(trac_system.lhs.size() == 2);

    // //prep:
    // scale rows (TODO: FIX)
    auto disp_rhs = disp_system.rhs;
    auto trac_rhs = trac_system.rhs * (1.0 / bem_input.params.shear_modulus);

    //reduce using constraints
    //stack rows
    auto stacked_rhs = concatenate({
        bem_input.constraints[0].get_reduced(disp_rhs[0]),
        bem_input.constraints[1].get_reduced(disp_rhs[1]),
        bem_input.constraints[2].get_reduced(trac_rhs[0]),
        bem_input.constraints[3].get_reduced(trac_rhs[1])
    });

    int n_unknown_trac_dofs = disp_system.rhs[0].size();
    int n_unknown_disp_dofs = trac_system.rhs[0].size();

    //solve:
    int count = 0;
    //TODO: Linear system tolerance should be a file parameter
    auto reduced_soln = solve_system(stacked_rhs.data, 1e-5,
        [&] (std::vector<double>& x, std::vector<double>& y) {
            std::cout << "iteration " << count << std::endl;
            count++;

            //unstack_rows
            auto unstacked_estimate = expand(stacked_rhs, x);

            //expand using constraints
            Function unknown_trac = {
                bem_input.constraints[0].get_all(unstacked_estimate[0],
                                                 n_unknown_trac_dofs),
                bem_input.constraints[1].get_all(unstacked_estimate[1],
                                                 n_unknown_trac_dofs)
            };
            Function unknown_disp = {
                bem_input.constraints[2].get_all(unstacked_estimate[2],
                                                 n_unknown_disp_dofs),
                bem_input.constraints[3].get_all(unstacked_estimate[3],
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
            auto disp_eval_vec = disp_eval.rhs;
            auto trac_eval_vec = trac_eval.rhs * (1.0 / bem_input.params.shear_modulus);

            //reduce using constraints
            //stack rows
            auto evaluated_lhs = concatenate({
                bem_input.constraints[0].get_reduced(disp_eval_vec[0]),
                bem_input.constraints[1].get_reduced(disp_eval_vec[1]),
                bem_input.constraints[2].get_reduced(trac_eval_vec[0]),
                bem_input.constraints[3].get_reduced(trac_eval_vec[1])
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
    Function soln_trac = {
        bem_input.constraints[0].get_all(unstacked_soln[0], n_unknown_trac_dofs),
        bem_input.constraints[1].get_all(unstacked_soln[1], n_unknown_trac_dofs)
    };
    Function soln_disp = {
        bem_input.constraints[2].get_all(unstacked_soln[2], n_unknown_disp_dofs),
        bem_input.constraints[3].get_all(unstacked_soln[3], n_unknown_disp_dofs)
    };


    //output
    auto in_filename_root = remove_extension(filename);
    auto out_filename_disp = in_filename_root + ".disp_out";
    auto out_filename_trac = in_filename_root + ".trac_out";
    HDFOutputter dispx_file(out_filename_disp + "x");
    HDFOutputter dispy_file(out_filename_disp + "y");
    HDFOutputter tracx_file(out_filename_trac + "x");
    HDFOutputter tracy_file(out_filename_trac + "y");
    if (soln_disp[0].size() > 0) {
        out_surface(dispx_file, bem_input.meshes.at("traction"), soln_disp[0], 1);
        out_surface(dispy_file, bem_input.meshes.at("traction"), soln_disp[1], 1);
    }
    if (soln_trac[0].size() > 0) {
        out_surface(tracx_file, bem_input.meshes.at("displacement"), soln_trac[0], 1);
        out_surface(tracy_file, bem_input.meshes.at("displacement"), soln_trac[1], 1);
    }

    //interior eval
}
