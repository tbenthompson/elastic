#include "spec.h"
#include "load.h"
#include "compute.h"
#include "3bem/function.h"
#include "3bem/armadillo_facade.h"

using namespace tbem;

template <size_t dim>
ConstraintMatrix form_traction_constraints(const MeshMap<dim>& meshes,
    const BCMap& bcs) 
{
    return from_constraints({});
}

template <size_t dim>
ConstraintMatrix form_displacement_constraints(const MeshMap<dim>& meshes, 
    const BCMap& bcs, size_t which_component) 
{
    auto continuity = mesh_continuity(meshes.at("traction").begin());
    auto cut_continuity = cut_at_intersection(
        continuity, meshes.at("traction").begin(), meshes.at("slip").begin()
    );
    auto constraints = convert_to_constraints(cut_continuity);
    auto constraint_matrix = from_constraints(constraints);
    return constraint_matrix;
}

template <size_t dim>
std::vector<ConstraintMatrix> form_constraints(const MeshMap<dim>& meshes,
    const BCMap& bcs)
{
    return {
        form_traction_constraints(meshes, bcs),
        form_traction_constraints(meshes, bcs),
        form_displacement_constraints(meshes, bcs, 0),
        form_displacement_constraints(meshes, bcs, 1)
    };
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
            get_displacement_BIE("displacement"),
            get_traction_BIE("traction")
        },
        form_constraints(meshes, bcs)
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

    // //prep:
    // scale rows
    auto disp_rhs = disp_system.rhs;
    auto trac_rhs = trac_system.rhs * (1.0 / bem_input.params.shear_modulus);

    //reduce using constraints
    //stack rows
    auto stacked_rhs = concatenate({
        condense_vector(bem_input.constraints[0], disp_rhs[0]),
        condense_vector(bem_input.constraints[1], disp_rhs[1]),
        condense_vector(bem_input.constraints[2], trac_rhs[0]),
        condense_vector(bem_input.constraints[3], trac_rhs[1])
    });

    int n_unknown_trac_dofs = disp_system.rhs[0].size();
    int n_unknown_disp_dofs = trac_system.rhs[0].size();

    std::cout << disp_system.lhs[0].src_mesh << std::endl;
    std::cout << disp_system.lhs[1].src_mesh << std::endl;
    std::cout << trac_system.lhs[0].src_mesh << std::endl;
    std::cout << trac_system.lhs[1].src_mesh << std::endl;
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
        bem_input.constraints, bem_input.constraints, lhs);
    auto combined_lhs = combine_components(condensed_lhs);
    std::cout << "Condition number: " << arma_cond(combined_lhs.ops[0]) << std::endl;
    return 0;

    //solve:
    int count = 0;
    //TODO: Linear system tolerance should be a file parameter
    auto reduced_soln = solve_system(stacked_rhs.data, 1e-10,
        [&] (std::vector<double>& x, std::vector<double>& y) {
            std::cout << "iteration " << count << std::endl;
            count++;

            //unstack_rows
            auto unstacked_estimate = expand(stacked_rhs, x);

            //expand using constraints
            BlockFunction unknown_trac = {
                distribute_vector(bem_input.constraints[0], unstacked_estimate[0],
                                                 n_unknown_trac_dofs),
                distribute_vector(bem_input.constraints[1], unstacked_estimate[1],
                                                 n_unknown_trac_dofs)
            };
            BlockFunction unknown_disp = {
                distribute_vector(bem_input.constraints[2], unstacked_estimate[2],
                                                 n_unknown_disp_dofs),
                distribute_vector(bem_input.constraints[3], unstacked_estimate[3],
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
            auto trac_eval_vec = -trac_eval.rhs * (1.0 / bem_input.params.shear_modulus);

            //reduce using constraints
            //stack rows
            auto evaluated_lhs = concatenate({
                condense_vector(bem_input.constraints[0], disp_eval_vec[0]),
                condense_vector(bem_input.constraints[1], disp_eval_vec[1]),
                condense_vector(bem_input.constraints[2], trac_eval_vec[0]),
                condense_vector(bem_input.constraints[3], trac_eval_vec[1])
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
        distribute_vector(bem_input.constraints[0], unstacked_soln[0],
                n_unknown_trac_dofs),
        distribute_vector(bem_input.constraints[1], unstacked_soln[1],
                n_unknown_trac_dofs)
    };
    BlockFunction soln_disp = {
        distribute_vector(bem_input.constraints[2], unstacked_soln[2],
                n_unknown_disp_dofs),
        distribute_vector(bem_input.constraints[3], unstacked_soln[3],
                n_unknown_disp_dofs)
    };


    //output
    auto in_filename_root = remove_extension(filename);
    auto out_filename_disp = in_filename_root + ".disp_out";
    auto out_filename_trac = in_filename_root + ".trac_out";
    if (soln_disp[0].size() > 0) {
        HDFOutputter disp_file(out_filename_disp);
        out_surface(disp_file, bem_input.meshes.at("traction"), soln_disp);
    }
    if (soln_trac[0].size() > 0) {
        HDFOutputter trac_file(out_filename_trac);
        out_surface(trac_file, bem_input.meshes.at("displacement"), soln_trac);
    }

    //interior eval
    // long-term, use an input file containing the list of points
    // find the nearest boundary point and make that the path to follow away
    // for nearfield integration --> nearest neighbor problems are a real theme!
    BCMap fields = bem_input.bcs; 
    fields[FieldDescriptor{"displacement", "traction"}] = soln_trac;
    fields[FieldDescriptor{"traction", "displacement"}] = soln_disp;

    auto x_vals = linspace(-1, 1, 20);
    auto y_vals = linspace(-1, 1, 20);
    std::vector<Vec<double,2>> locs;
    std::vector<ObsPt<2>> obs_pts;
    for (size_t i = 0; i < x_vals.size(); i++) {
        for (size_t j = 0; j < y_vals.size(); j++) {
            locs.push_back({x_vals[i], y_vals[j]});
        }
    }
    for (size_t i = 0; i < locs.size(); i++) {
        obs_pts.push_back({0.001, locs[i], {0, 1}, 
            Vec<double,2>{0.5, 0.5} - locs[i]});
    }
    auto interior = compute_interior(obs_pts, bem_input,
        get_displacement_BIE("displacement"), fields);
    HDFOutputter interior_disp_file(out_filename_disp + "int");
    out_volume(interior_disp_file, locs, interior);
}
