#include "elastic.h"
#include "3bem/kernels.h"
#include "3bem/quadrature.h"
#include "3bem/bem.h"
#include "3bem/petsc_interface.h"
#include "3bem/util.h"

using namespace tbem;

std::string remove_extension(const std::string& filename) {
    size_t lastdot = filename.find_last_of(".");
    if (lastdot == std::string::npos) {
        return filename;
    }
    return filename.substr(0, lastdot); 
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage is 'elastic_process filename'" << std::endl;
        return 1;
    }
    
    auto filename = argv[1];
    auto doc = parse_json(load_file(filename));
    auto elements = collect_elements(doc);
    auto params = parse_parameters(doc);
    auto elastic_prob = build_problem(elements);

    auto trac_mesh = elastic_prob.traction_mesh;
    auto slip_mesh = elastic_prob.slip_mesh;
    auto disp_mesh = elastic_prob.displacement_mesh;

    int n_disp_dofs = 2 * disp_mesh.facets.size();
    int n_trac_dofs = 2 * trac_mesh.facets.size();
    int n_slip_dofs = 2 * slip_mesh.facets.size();
    std::cout << "Number of displacement DOFs: " << n_disp_dofs << std::endl;
    std::cout << "Number of traction DOFs: " << n_trac_dofs << std::endl;
    std::cout << "Number of slip DOFs: " << n_slip_dofs << std::endl;

    // Gather the imposed displacement continuity constraints.
    auto displacement_continuity = 
        ConstraintMatrix::from_constraints(mesh_continuity<2>(trac_mesh));

    // Remove continuity constraints at the intersection of the fault and the 
    // surface mesh.
    // TODO: Here, the problem is more complex than in the past because I 
    //       should allow the vertices not to match up.
    auto constraints = apply_discontinuities<2>(
        trac_mesh, slip_mesh, displacement_continuity
    );

    
    // Setup the kernels that are necessary.
    ElasticKernels<2> ek(params.shear_modulus, params.poisson_ratio);

    
    // Setup the quadrature
    QuadStrategy<2> qs(
        params.obs_quad_order,
        params.src_far_quad_order,
        params.n_singular_steps,
        params.far_threshold,
        params.near_tol
    );

    // Gather the imposed boundary conditions
    //TODO: This will be unnecessary when Vec2 can be passed into problem
    // BEGIN UNNCESSARY
    std::array<std::vector<double>,2> slip_bcs = {
        std::vector<double>(n_slip_dofs),   
        std::vector<double>(n_slip_dofs)
    };
    for (std::size_t i = 0; i < slip_mesh.facets.size(); i++) {
        for (int c = 0; c < 2; c++) {
            slip_bcs[0][i * 2 + c] = elastic_prob.slip_bcs.facets[i].vertices[c][0];
            slip_bcs[1][i * 2 + c] = elastic_prob.slip_bcs.facets[i].vertices[c][1];
        }
    }

    std::array<std::vector<double>,2> trac_bcs = {
        std::vector<double>(n_trac_dofs),   
        std::vector<double>(n_trac_dofs)
    };
    for (std::size_t i = 0; i < trac_mesh.facets.size(); i++) {
        for (int c = 0; c < 2; c++) {
            trac_bcs[0][i * 2 + c] = elastic_prob.traction_bcs.facets[i].vertices[c][0];
            trac_bcs[1][i * 2 + c] = elastic_prob.traction_bcs.facets[i].vertices[c][1];
        }
    }

    std::array<std::vector<double>,2> disp_bcs = {
        std::vector<double>(n_disp_dofs),   
        std::vector<double>(n_disp_dofs)
    };
    for (std::size_t i = 0; i < disp_mesh.facets.size(); i++) {
        for (int c = 0; c < 2; c++) {
            disp_bcs[0][i * 2 + c] =
                elastic_prob.displacement_bcs.facets[i].vertices[c][0];
            disp_bcs[1][i * 2 + c] =
                elastic_prob.displacement_bcs.facets[i].vertices[c][1];
        }
    }
    //END UNNECESSARY
    
    // Build the RHS for the traction_mesh DOFs 
    std::array<std::vector<double>,2> all_trac_rhs = {
        std::vector<double>(n_trac_dofs, 0.0),   
        std::vector<double>(n_trac_dofs, 0.0)
    };

    for (int k = 0; k < 2; k++) {
        Problem<2> p_mass = {trac_mesh, trac_mesh,
                             one<2>, trac_bcs[k]};
        auto rhs_mass = mass_term(p_mass, qs);
        for (unsigned int i = 0; i < rhs_mass.size(); i++) {
            all_trac_rhs[k][i] -= rhs_mass[i];
        }

        for (int j = 0; j < 2; j++) {
            Problem<2> p_slip_trac = {slip_mesh, trac_mesh,
                                      ek.hypersingular_mat[k][j], slip_bcs[j]};
            auto int_slip_trac = direct_interact(p_slip_trac, qs);

            Problem<2> p_trac_trac = {trac_mesh, trac_mesh,
                                      ek.adjoint_traction_mat[k][j], trac_bcs[j]};
            auto int_trac_trac = direct_interact(p_trac_trac, qs);

            Problem<2> p_disp_trac = {disp_mesh, trac_mesh,
                                      ek.hypersingular_mat[k][j], disp_bcs[j]};
            auto int_disp_trac = direct_interact(p_disp_trac, qs);
                                        
            for (unsigned int i = 0; i < int_slip_trac.size(); i++) {
                all_trac_rhs[k][i] += 
                    int_slip_trac[i] - int_trac_trac[i] + int_disp_trac[i];
            }
        }
    }

    std::array<std::vector<double>,2> trac_rhs = {
        constraints.get_reduced(all_trac_rhs[0]),
        constraints.get_reduced(all_trac_rhs[1]),
    };

    int n_reduced_trac_dofs = trac_rhs[0].size();

    // Build the RHS for the displacment_mesh DOFs
    std::array<std::vector<double>,2> all_disp_rhs = {
        std::vector<double>(n_disp_dofs, 0.0),   
        std::vector<double>(n_disp_dofs, 0.0)
    };

    for (int k = 0; k < 2; k++) {
        Problem<2> p_mass = {disp_mesh, disp_mesh,
                             one<2>, disp_bcs[k]};
        auto rhs_mass = mass_term(p_mass, qs);
        for (unsigned int i = 0; i < rhs_mass.size(); i++) {
            all_disp_rhs[k][i] -= rhs_mass[i];
        }

        for (int j = 0; j < 2; j++) {
            Problem<2> p_trac_disp = {trac_mesh, disp_mesh,
                                      ek.displacement_mat[k][j], trac_bcs[j]};
            auto int_trac_disp = direct_interact(p_trac_disp, qs);

            Problem<2> p_disp_disp = {disp_mesh, disp_mesh,
                                      ek.traction_mat[k][j], disp_bcs[j]};
            auto int_disp_disp = direct_interact(p_disp_disp, qs);
                                        
            for (unsigned int i = 0; i < int_trac_disp.size(); i++) {
                all_disp_rhs[k][i] += int_trac_disp[i] - int_disp_disp[i];
            }
        }
    }

    std::array<std::vector<double>,2> disp_rhs = {
        constraints.get_reduced(all_disp_rhs[0]),
        constraints.get_reduced(all_disp_rhs[1]),
    };

    int n_reduced_disp_dofs = disp_rhs[0].size();

    
    // Build the LHS matrix for the traction_mesh DOFs
    std::array<std::array<std::vector<double>,2>,2> trac_trac_mats;
    std::array<std::array<std::vector<double>,2>,2> disp_trac_mats;
    for (int k = 0; k < 2; k++) {
        for (int j = 0; j < 2; j++) {
            Problem<2> p_trac_trac = {trac_mesh, trac_mesh, 
                                      ek.hypersingular_mat[k][j], {}};
            trac_trac_mats[k][j] = interact_matrix(p_trac_trac, qs);
            // TODO: For situations like this, it would be useful to slightly
            // encapsulate the std::vector storage mechanism
            for (std::size_t el = 0; el < trac_trac_mats[k][j].size(); el++) {
                trac_trac_mats[k][j][el] = -trac_trac_mats[k][j][el];
            }

            Problem<2> p_disp_trac = {disp_mesh, trac_mesh, 
                                      ek.adjoint_traction_mat[k][j], {}};
            disp_trac_mats[k][j] = interact_matrix(p_disp_trac, qs);
        }
    }

    // Build the LHS matrix for the displacement_mesh DOFs.
    std::array<std::array<std::vector<double>,2>,2> trac_disp_mats;
    std::array<std::array<std::vector<double>,2>,2> disp_disp_mats;
    for (int k = 0; k < 2; k++) {
        for (int j = 0; j < 2; j++) {
            Problem<2> p_trac_disp = {trac_mesh, disp_mesh, 
                                      ek.traction_mat[k][j], {}};
            trac_disp_mats[k][j] = interact_matrix(p_trac_disp, qs);

            Problem<2> p_disp_disp = {disp_mesh, disp_mesh, 
                                      ek.displacement_mat[k][j], {}};
            disp_disp_mats[k][j] = interact_matrix(p_disp_disp, qs);
            // TODO: For situations like this, it would be useful to slightly
            // encapsulate the std::vector storage mechanism
            for (std::size_t el = 0; el < disp_disp_mats[k][j].size(); el++) {
                disp_disp_mats[k][j][el] = -disp_disp_mats[k][j][el];
            }
        }
    }

    std::vector<double> all_rhs(2 * (n_reduced_trac_dofs + n_reduced_disp_dofs));

    auto trac_rhs_begin = all_rhs.begin();
    auto disp_rhs_begin = all_rhs.begin() + 2 * n_reduced_trac_dofs;
    for (int d = 0; d < 2; d++) {
        std::copy(trac_rhs[d].begin(), trac_rhs[d].end(),
                  trac_rhs_begin + d * n_reduced_trac_dofs);
        std::copy(disp_rhs[d].begin(), disp_rhs[d].end(),
                  disp_rhs_begin + d * n_reduced_disp_dofs);
    }

    // Solve the linear system.
    // return 0;
    int count = 0;
    auto reduced_soln = solve_system(all_rhs, 1e-5,
        [&] (std::vector<double>& x, std::vector<double>& y) {
            std::cout << "iteration " << count << std::endl;
            count++;
            std::array<std::vector<double>,2> trac_temp;
            std::array<std::vector<double>,2> disp_temp;
            auto trac_begin = x.begin();
            auto disp_begin = x.begin() + 2 * n_reduced_trac_dofs;
            for (int i = 0; i < 2; i++) {
                auto trac_temp_reduced = std::vector<double>(n_reduced_trac_dofs);
                auto disp_temp_reduced = std::vector<double>(n_reduced_disp_dofs);
                std::copy(trac_begin + i * n_reduced_trac_dofs,
                          trac_begin + (i + 1) * n_reduced_trac_dofs,
                          trac_temp_reduced.begin());
                std::copy(disp_begin + i * n_reduced_disp_dofs,
                          disp_begin + (i + 1) * n_reduced_disp_dofs,
                          disp_temp_reduced.begin());
                trac_temp[i] = constraints.get_all(trac_temp_reduced, n_trac_dofs);
                disp_temp[i] = constraints.get_all(disp_temp_reduced, n_disp_dofs);
            }

            for (int k = 0; k < 2; k++) {
                std::vector<double> trac_out_full(n_trac_dofs, 0.0);
                std::vector<double> disp_out_full(n_disp_dofs, 0.0);
                for (int j = 0; j < 2; j++) {
                    for (int mi = 0; mi < n_trac_dofs; mi++) {
                        for (int ni = 0; ni < n_trac_dofs; ni++) {
                            trac_out_full[mi] += 
                                trac_trac_mats[k][j][mi * n_trac_dofs + ni] 
                                * trac_temp[j][ni];
                        }
                        for (int ni = 0; ni < n_disp_dofs; ni++) {
                            trac_out_full[mi] += 
                                disp_trac_mats[k][j][mi * n_disp_dofs + ni] 
                                * disp_temp[j][ni];
                        }
                    }

                    for (int mi = 0; mi < n_disp_dofs; mi++) {
                        for (int ni = 0; ni < n_trac_dofs; ni++) {
                            disp_out_full[mi] += 
                                trac_disp_mats[k][j][mi * n_trac_dofs + ni] 
                                * trac_temp[j][ni];
                        }
                        for (int ni = 0; ni < n_disp_dofs; ni++) {
                            disp_out_full[mi] += 
                                disp_disp_mats[k][j][mi * n_disp_dofs + ni] 
                                * disp_temp[j][ni];
                        }
                    }
                }

                auto trac_out_begin = y.begin();
                auto disp_out_begin = y.begin() + 2 * n_reduced_trac_dofs;
                auto trac_out_reduced = constraints.get_reduced(trac_out_full);
                auto disp_out_reduced = constraints.get_reduced(disp_out_full);
                std::copy(trac_out_reduced.begin(), trac_out_reduced.end(),
                          trac_out_begin + k * n_reduced_trac_dofs);
                std::copy(disp_out_reduced.begin(), disp_out_reduced.end(),
                          disp_out_begin + k * n_reduced_disp_dofs);
            }
        });

    // TODO: Calculate the tractions on the slip_mesh.

    auto trac_begin = reduced_soln.begin();
    auto disp_begin = reduced_soln.begin() + 2 * n_reduced_trac_dofs;
    std::array<std::vector<double>,2> disp_soln;
    std::array<std::vector<double>,2> trac_soln;
    for (int i = 0; i < 2; i++) {
        auto reduced_disp = std::vector<double>(n_reduced_trac_dofs);
        auto reduced_trac = std::vector<double>(n_reduced_trac_dofs);
        std::copy(trac_begin + i * n_reduced_trac_dofs,
                  trac_begin + (i + 1) * n_reduced_trac_dofs,
                  reduced_trac.begin());
        std::copy(disp_begin + i * n_reduced_disp_dofs,
                  disp_begin + (i + 1) * n_reduced_disp_dofs,
                  reduced_disp.begin());
        trac_soln[i] = constraints.get_all(reduced_trac, n_trac_dofs);
        disp_soln[i] = constraints.get_all(reduced_disp, n_disp_dofs);
    }

    //TODO: More output?
    std::string file_root = remove_extension(filename);
    if (trac_mesh.facets.size() > 0) {
        hdf_out_surface<2>(file_root + ".trac_out",
                           trac_mesh, {trac_soln[0], trac_soln[1]});
    }

    if (disp_mesh.facets.size() > 0) {
        hdf_out_surface<2>(file_root + ".disp_out",
                           disp_mesh, {disp_soln[0], disp_soln[1]});
    }
}
