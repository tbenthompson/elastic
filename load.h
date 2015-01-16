#ifndef __QWENMBBBZBLAHEWJPQ_LOAD_H
#define __QWENMBBBZBLAHEWJPQ_LOAD_H

#include <map>
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


std::string load_file(const std::string& filename);
rapidjson::Document parse_json(const std::string& json);
Parameters get_parameters(const rapidjson::Document& doc);

template <size_t dim>
struct Element {
    tbem::Vec<tbem::Vec<double,dim>,dim> pts; 
    std::string bc_type;
    tbem::Vec<tbem::Vec<double,dim>,dim> bc;
    int n_refines;
};
std::vector<Element<2>> get_elements(const rapidjson::Document& doc);

template <size_t dim>
using MeshSet = std::map<std::string,tbem::Mesh<dim>>;

typedef std::vector<std::vector<double>> BC;

typedef std::map<std::string,BC> BCSet;

template <size_t dim>
MeshSet<dim> get_meshes(const std::vector<Element<dim>>& elements);

template <size_t dim>
BCSet get_bcs(const std::vector<Element<dim>>& elements);


#endif
