#include "elastic.h"
#include "3bem/kernels.h"
#include "3bem/quadrature.h"
#include "3bem/bem.h"
#include "3bem/petsc_interface.h"
#include "3bem/util.h"

using namespace tbem;

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage is 'elastic_process filename'" << std::endl;
        return 1;
    }
    
    auto filename = argv[1];
    auto doc = parse_json(load_file(filename));
    auto elements = collect_elements(doc);
    auto elastic_prob = build_problem(doc, elements);

    auto traction_mesh = elastic_prob.traction_mesh;
    auto slip_mesh = elastic_prob.slip_mesh;

    int n_trac_dofs = 2 * traction_mesh.facets.size();
    std::cout << n_trac_dofs << std::endl;
    int n_slip_dofs = 2 * slip_mesh.facets.size();

    // Gather the imposed boundary conditions
    auto displacement_continuity = 
        ConstraintMatrix::from_constraints(mesh_continuity<2>(traction_mesh));

    //Average constraint.
    // Constraint

    // Remove continuity constraints at the intersection of the fault and the 
    // surface mesh.
    // TODO: Here, the problem is more complex than in the past because I 
    //       should allow the vertices not to match up.
    auto constraints = apply_discontinuities<2>(
        traction_mesh, slip_mesh, displacement_continuity
    );//.add_constraints({average_constraint});

    
    // Setup the kernels that are necessary.
    //TODO: Make these a parameter in the input file/ElasticProblem
    double shear_modulus = 30e9;
    double poisson_ratio = 0.25;
    ElasticKernels<2> ek(shear_modulus, poisson_ratio);

    // Setup the quadrature
    // TODO: Parameters for the file/ElasticProblem!
    QuadStrategy<2> qs(2);

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
    for (std::size_t i = 0; i < traction_mesh.facets.size(); i++) {
        for (int c = 0; c < 2; c++) {
            trac_bcs[0][i * 2 + c] = elastic_prob.traction_bcs.facets[i].vertices[c][0];
            trac_bcs[1][i * 2 + c] = elastic_prob.traction_bcs.facets[i].vertices[c][1];
        }
    }
    //END UNNECESSARY
    
    // Build the RHS for the traction_mesh DOFs 
    std::array<std::vector<double>,2> all_trac_rhs = {
        std::vector<double>(n_trac_dofs, 0.0),   
        std::vector<double>(n_trac_dofs, 0.0)
    };

    for (int k = 0; k < 2; k++) {
        Problem<2> p_mass = {traction_mesh, traction_mesh,
                             one<2>, trac_bcs[k]};
        auto rhs_mass = mass_term(p_mass, qs);
        for (unsigned int i = 0; i < rhs_mass.size(); i++) {
            all_trac_rhs[k][i] += rhs_mass[i];
        }

        for (int j = 0; j < 2; j++) {
            Problem<2> p_slip_trac = {slip_mesh, traction_mesh,
                                      ek.hypersingular_mat[k][j], slip_bcs[j]};
            auto int0 = direct_interact(p_slip_trac, qs);
            Problem<2> p_trac_trac = {traction_mesh, traction_mesh,
                                      ek.adjoint_traction_mat[k][j], trac_bcs[j]};
            auto int1 = direct_interact(p_trac_trac, qs);
            for (unsigned int i = 0; i < int0.size(); i++) {
                all_trac_rhs[k][i] += int0[i] + int1[i];
                std::cout << int1[i] << std::endl;
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
    // return 0;
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

    //TODO: Strip the ".in" from the filename
    hdf_out_surface<2>(std::string(filename) + ".out", traction_mesh, {soln[0], soln[1]});
}
