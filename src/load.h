#ifndef __QWENMBBBZBLAHEWJPQ_LOAD_H
#define __QWENMBBBZBLAHEWJPQ_LOAD_H

#include <map>
#include "rapidjson/document.h"
#include "3bem/constraint.h"
#include "3bem/mesh.h"

std::string load_file(const std::string& filename);
rapidjson::Document parse_json(const std::string& json);

struct Parameters {
    int obs_quad_order;
    int src_far_quad_order;
    int n_singular_steps;
    double far_threshold;
    double near_tol;

    double solver_tol;

    double poisson_ratio;
    double shear_modulus;
};
const Parameters default_params{2, 2, 8, 3.0, 1e-4, 1e-6, 0.25, 30e9};

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

template <size_t dim>
using MeshMap = std::map<std::string, tbem::Mesh<dim>>;

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

typedef std::vector<std::vector<double>> BlockFunction;
typedef std::map<FieldDescriptor,BlockFunction> BCMap;

template <size_t dim>
MeshMap<dim> get_meshes(const std::vector<Element<dim>>& elements);

template <size_t dim>
BCMap get_bcs(const std::vector<Element<dim>>& elements);

template <size_t dim>
std::vector<tbem::Vec<double,dim>> get_pts(const rapidjson::Document& doc);
#endif
