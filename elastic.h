#ifndef __QWENMBBBZBLAHEWJPQ_ELASTIC_H
#define __QWENMBBBZBLAHEWJPQ_ELASTIC_H

#include "rapidjson/document.h"
#include "constraint.h"

template <int dim>
struct ElasticProblem {
    // The mesh on which the displacement is known
    // Mesh<dim> displacement_mesh;

    // // The mesh on which the traction is known
    // Mesh<dim> traction_mesh;

    // // The mesh on which the slip is known
    // Mesh<dim> slip_mesh;

    // std::vector<Vec<double,dim>> displacement_bcs;
    // std::vector<Vec<double,dim>> traction_bcs;
    // std::vector<Vec<double,dim>> slip_bcs;
};

enum BC {DISPLACEMENT, TRACTION, SLIP, CRACK};

std::string load_file(std::string filename);
rapidjson::Document parse_json(std::string json);
// ElasticProblem<2> load_file(std::string filename);
ConstraintMatrix build_constraints(const ElasticProblem<2>& elastic_prob);

// int main(int argc, char* argv[]) {
//     if (argc < 2) {
//         std::cout << "Usage is 'elastic_process filename'" << std::endl;
//         return 1;
//     }
// 
//     auto elastic_prob = load_file(argv[1]);
// 
//     // Gather the imposed boundary conditions
//     auto constraints = build_constraints(elastic_prob);
// 
//     // Setup the kernels that are necessary.
//     //TODO: Make these a parameter in the input file/ElasticProblem
//     double shear_modulus = 30e9;
//     double poisson_ratio = 0.25;
//     //TODO: template on dimension for the elastic kernels.
//     //TODO: finish the plane strain elastic kernels
//     ElasticKernels<2> ek(shear_modulus, poisson_ratio);
// 
//     // Setup the quadrature
//     // TODO: Parameters for the file/ElasticProblem!
//     QuadStrategy<2> qs(2);
    
    // Build the RHS for the traction_mesh DOFs 

    // Build the RHS for the displacment_mesh DOFs
    
    // Build the LHS matrix for the traction_mesh DOFs

    // Build the LHS matrix for the displacement_mesh DOFs.

    // Solve the linear system.

    // Calculate the tractions on the slip_mesh.
// }

#endif
