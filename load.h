#ifndef __QWENMBBBZBLAHEWJPQ_LOAD_H
#define __QWENMBBBZBLAHEWJPQ_LOAD_H

#include <map>
#include "rapidjson/document.h"
#include "3bem/constraint.h"
#include "3bem/mesh.h"

std::string remove_extension(const std::string& filename); 
std::string load_file(const std::string& filename);
rapidjson::Document parse_json(const std::string& json);

struct Parameters {
    int obs_quad_order;
    int src_far_quad_order;
    int n_singular_steps;
    double far_threshold;
    double near_tol;

    double poisson_ratio;
    double shear_modulus;
};

Parameters get_parameters(const rapidjson::Document& doc);

template <size_t dim>
struct Element {
    tbem::Vec<tbem::Vec<double,dim>,dim> pts; 
    std::string bc_type;
    tbem::Vec<tbem::Vec<double,dim>,dim> bc;
    int n_refines;
};

template <size_t dim>
std::vector<Element<dim>> get_elements(const rapidjson::Document& doc);

const std::string mesh_types[3] = {
    "traction", "displacement", "slip"
};
template <size_t dim>
using MeshSet = std::map<std::string, tbem::Mesh<dim>>;

typedef std::vector<std::vector<double>> BC;

struct FieldDescriptor 
{
    std::string where;
    std::string what;
    bool operator<(const FieldDescriptor& fd) const {
        if (where < fd.where) {
            return true;
        }
        if (where > fd.where) {
            return false;
        }
        if (what < fd.what) {
            return true;
        }
        return false;
    }
};

typedef std::map<FieldDescriptor, BC> BCSet;

template <size_t dim>
MeshSet<dim> get_meshes(const std::vector<Element<dim>>& elements);

template <size_t dim>
BCSet get_bcs(const std::vector<Element<dim>>& elements);


#endif
