#ifndef __QWENMBBBZBLAHEWJPQ_ELASTIC_H
#define __QWENMBBBZBLAHEWJPQ_ELASTIC_H

#include "rapidjson/document.h"
#include "constraint.h"
#include "mesh.h"

template <int dim>
struct ElasticProblem {
    // The mesh on which the displacement is known
    Mesh<dim> displacement_mesh;

    // The mesh on which the traction is known
    Mesh<dim> traction_mesh;

    // The mesh on which the slip is known
    Mesh<dim> slip_mesh;

    std::vector<Vec<Vec<double,dim>,dim>> displacement_bcs;
    std::vector<Vec<Vec<double,dim>,dim>> traction_bcs;
    std::vector<Vec<Vec<double,dim>,dim>> slip_bcs;
};

enum BCType {DISPLACEMENT, TRACTION, SLIP, CRACK};

template <int dim>
struct Element {
    Vec<Vec<double,dim>,dim> pts; 
    BCType bc_type;
    Vec<Vec<double,dim>,dim> bc;
};

std::string load_file(std::string filename);
rapidjson::Document parse_json(std::string json);
std::vector<Element<2>> collect_elements(const rapidjson::Document& doc);
ElasticProblem<2> build_problem(const rapidjson::Document& doc, 
                                const std::vector<Element<2>>& elements);

#endif
