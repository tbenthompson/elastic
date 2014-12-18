#include "elastic.h"
#include "kernels.h"
#include "quadrature.h"
#include "bem.h"
#include "petsc_interface.h"
#include "util.h"

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage is 'elastic_process filename'" << std::endl;
        return 1;
    }

    auto doc = parse_json(load_file(argv[1]));
    auto elements = collect_elements(doc);
    auto elastic_prob = build_problem(doc, elements);

    Mesh<2> traction_mesh = elastic_prob.traction_mesh.refine_repeatedly(8);
    Mesh<2> slip_mesh = elastic_prob.slip_mesh.refine_repeatedly(0);

    // Gather the imposed boundary conditions
    auto displacement_continuity = ConstraintMatrix::from_constraints(
        mesh_continuity(traction_mesh)
    );

    // Remove continuity constraints at the intersection of the fault and the 
    // surface mesh.
    // TODO: Here, the problem is more complex than in the past because I 
    //       should allow the vertices not to match up.
    auto constraints = apply_discontinuities(
        traction_mesh, slip_mesh, displacement_continuity
    );
    
    // Setup the kernels that are necessary.
    //TODO: Make these a parameter in the input file/ElasticProblem
    double shear_modulus = 30e9;
    double poisson_ratio = 0.25;

    //TODO: template on dimension for the elastic kernels.
    //TODO: finish the plane strain elastic kernels
    ElasticKernels<2> ek(shear_modulus, poisson_ratio);

    // Setup the quadrature
    // TODO: Parameters for the file/ElasticProblem!
    QuadStrategy<2> qs(2);

    int n_trac_dofs = 2 * traction_mesh.facets.size();
    int n_slip_dofs = 2 * slip_mesh.facets.size();

    std::array<std::vector<double>,2> slip_bcs = {
        std::vector<double>(n_slip_dofs),   
        std::vector<double>(n_slip_dofs)
    };
    for (std::size_t i = 0; i < slip_mesh.facets.size(); i++) {
        for (int c = 0; c < 2; c++) {
            slip_bcs[0][i * 2 + c] = elastic_prob.slip_bcs[i][c][0];
            slip_bcs[1][i * 2 + c] = elastic_prob.slip_bcs[i][c][1];
        }
    }
    
    // Build the RHS for the traction_mesh DOFs 
    std::array<std::vector<double>,2> all_trac_rhs = {
        std::vector<double>(n_trac_dofs, 0.0),   
        std::vector<double>(n_trac_dofs, 0.0)
    };

    for (int k = 0; k < 2; k++) {
        for (int j = 0; j < 2; j++) {
            Problem<2> p = {slip_mesh, traction_mesh,
                            ek.hypersingular_mat[k][j], slip_bcs[j]};
            auto res = direct_interact(p, qs);
            for (unsigned int i = 0; i < res.size(); i++) {
                all_trac_rhs[k][i] += res[i];
            }
        }
    }

    std::array<std::vector<double>,2> trac_rhs = {
        constraints.get_reduced(all_trac_rhs[0]),
        constraints.get_reduced(all_trac_rhs[1]),
    };

    int n_reduced_trac_dofs = trac_rhs[0].size();

    std::vector<double> vector_trac_rhs(2 * n_reduced_trac_dofs);
    for (int d = 0; d < 2; d++) {
        std::copy(trac_rhs[d].begin(), trac_rhs[d].end(), vector_trac_rhs.begin() + 
                                                d * n_reduced_trac_dofs);
    }

    // Build the RHS for the displacment_mesh DOFs
    
    // Build the LHS matrix for the traction_mesh DOFs
    std::array<std::array<std::vector<double>,2>,2> trac_trac_mats;
    for (int k = 0; k < 2; k++) {
        for (int j = 0; j < 2; j++) {
            Problem<2> p = {traction_mesh, traction_mesh, 
                            ek.hypersingular_mat[k][j], {}};
            trac_trac_mats[k][j] = interact_matrix(p, qs);
        }
    }

    // Build the LHS matrix for the displacement_mesh DOFs.

    // Solve the linear system.
    int count = 0;
    auto surface_disp = solve_system(vector_trac_rhs, 1e-5,
        [&] (std::vector<double>& x, std::vector<double>& y) {
            std::cout << "iteration " << count << std::endl;
            count++;
            std::array<std::vector<double>,2> x_temp;
            for (int i = 0; i < 2; i++) {
                auto x_temp_reduced = std::vector<double>(n_trac_dofs);
                std::copy(x.begin() + i * n_reduced_trac_dofs,
                          x.begin() + (i + 1) * n_reduced_trac_dofs,
                          x_temp_reduced.begin());
                x_temp[i] = constraints.get_all(x_temp_reduced, n_trac_dofs);
            }

            std::array<std::vector<double>,2> y_temp;
            for (int k = 0; k < 2; k++) {
                auto y_temp_full = std::vector<double>(n_trac_dofs, 0.0);
                for (int j = 0; j < 2; j++) {
                    for (int mi = 0; mi < n_trac_dofs; mi++) {
                        for (int ni = 0; ni < n_trac_dofs; ni++) {
                            y_temp_full[mi] += 
                                trac_trac_mats[k][j][mi * n_trac_dofs + ni] 
                                * x_temp[j][ni];
                        }
                    }
                }
                auto y_temp_reduced = constraints.get_reduced(y_temp_full);
                std::copy(y_temp_reduced.begin(), y_temp_reduced.end(),
                          y.begin() + k * n_reduced_trac_dofs);
            }
        });

    // Calculate the tractions on the slip_mesh.

    std::array<std::vector<double>,2> soln;
    for (int i = 0; i < 2; i++) {
        auto reduced_soln = std::vector<double>(n_reduced_trac_dofs);
        std::copy(surface_disp.begin() + i * n_reduced_trac_dofs,
                  surface_disp.begin() + (i + 1) * n_reduced_trac_dofs,
                  reduced_soln.begin());
        soln[i] = constraints.get_all(reduced_soln, n_trac_dofs);
    }
    hdf_out_surface("2dthrust0.hdf5", traction_mesh, {soln[0], soln[1]});
}
