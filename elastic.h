#ifndef __QWENMBBBZBLAHEWJPQ_ELASTIC_H
#define __QWENMBBBZBLAHEWJPQ_ELASTIC_H

#include "rapidjson/document.h"
#include "3bem/constraint.h"
#include "3bem/mesh.h"

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

// TODO: Is rapidjson parallelizable?
std::string load_file(std::string filename);
rapidjson::Document parse_json(std::string json);
std::vector<Element<2>> collect_elements(const rapidjson::Document& doc);
ElasticProblem<2> build_problem(const rapidjson::Document& doc, 
                                const std::vector<Element<2>>& elements);


#endif
