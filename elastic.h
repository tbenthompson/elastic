#ifndef __QWENMBBBZBLAHEWJPQ_ELASTIC_H
#define __QWENMBBBZBLAHEWJPQ_ELASTIC_H

#include "rapidjson/document.h"
#include "3bem/constraint.h"
#include "3bem/mesh.h"

struct Parameters {
    int obs_quad_order;
    int src_far_quad_order;
    int n_singular_steps;
    double far_threshold;
    double near_tol;

    double poisson_ratio;
    double shear_modulus;
};

const Parameters default_params{2, 2, 6, 3.0, 1e-2, 0.25, 30e9};

template <int dim>
struct ElasticProblem {

    // The mesh on which the displacement is known
    tbem::Mesh<dim> displacement_mesh;

    // The mesh on which the traction is known
    tbem::Mesh<dim> traction_mesh;

    // The mesh on which the slip is known
    tbem::Mesh<dim> slip_mesh;

    // Boundary conditions can also be defined by points on a mesh.
    tbem::Mesh<dim> displacement_bcs;
    tbem::Mesh<dim> traction_bcs;
    tbem::Mesh<dim> slip_bcs;
};

enum BCType {DISPLACEMENT, TRACTION, SLIP, CRACK};

template <int dim>
struct Element {
    tbem::Vec<tbem::Vec<double,dim>,dim> pts; 
    BCType bc_type;
    tbem::Vec<tbem::Vec<double,dim>,dim> bc;
    int n_refines;
};

std::string load_file(std::string filename);
rapidjson::Document parse_json(std::string json);
std::vector<Element<2>> collect_elements(const rapidjson::Document& doc);
Parameters parse_parameters(const rapidjson::Document& doc);
ElasticProblem<2> build_problem(const std::vector<Element<2>>& elements);


#endif
