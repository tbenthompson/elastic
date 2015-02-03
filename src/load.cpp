#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include "load.h"
#include "spec.h"
#include "3bem/vec.h"
#include "3bem/vertex_iterator.h"

using namespace tbem;

std::string load_file(const std::string& filename) {
    std::ifstream file;
    file.open(filename);
    if (!file.is_open()) {
        throw std::invalid_argument("Filename: " + filename + " could not be" + 
                                    " opened (does it exist?).");
    }

    // Read the file into a string
    std::stringstream buffer;
    buffer << file.rdbuf();

    file.close();

    return buffer.str();
}

rapidjson::Document parse_json(const std::string& json) {
    rapidjson::Document doc;
    bool parse_error = doc.Parse(json.c_str()).HasParseError();
    if (parse_error) {
        throw std::invalid_argument("There are json errors.");
    }
    return doc;
}

template <size_t dim>
Vec<Vec<double,dim>,dim> parse_tensor(const rapidjson::Value& e_json, 
    std::string field_name, std::string except_text) 
{
    if (!e_json.HasMember(field_name.c_str())) {
        throw std::invalid_argument(except_text);
    }

    const auto& V = e_json[field_name.c_str()];
    if (!V.IsArray() || V.Capacity() != dim) {
        throw std::invalid_argument(except_text);
    }

    Vec<Vec<double,dim>,dim> out;
    for (int d0 = 0; d0 < dim; d0++) {
        if (!V[d0].IsArray() || V[d0].Capacity() != dim) {
            throw std::invalid_argument(except_text);
        }
        for (int d1 = 0; d1 < dim; d1++) {
            if(!V[d0][d1].IsDouble()) {
                throw std::invalid_argument(except_text);
            }
            out[d0][d1] = V[d0][d1].GetDouble();
        }    
    }
    
    return out;
}

//TODO: Should this be configurable?
#define GETPARAM(TYPE, NAME) {\
    if (doc.HasMember(#NAME) && doc[#NAME].Is##TYPE()) {\
        out.NAME = doc[#NAME].Get##TYPE();\
    }}
Parameters get_parameters(const rapidjson::Document& doc) {
    Parameters out = default_params;
    GETPARAM(Int, obs_quad_order);
    GETPARAM(Int, src_far_quad_order);
    GETPARAM(Int, n_singular_steps);
    GETPARAM(Double, far_threshold);
    GETPARAM(Double, near_tol);
    GETPARAM(Double, poisson_ratio);
    GETPARAM(Double, shear_modulus);
    return out;
}

template <size_t dim>
std::vector<Element<dim>> get_elements(const rapidjson::Document& doc) {
    const auto& element_list = doc["elements"];

    std::vector<Element<dim>> out;

    std::string except_text = "An element object must have a (dim x dim) array of float 'pts' and 'bc' array, a string 'bc_type' field, and an integer 'refine' field.";

    for (std::size_t i = 0; i < element_list.Size(); i++) {
        auto& e_json = element_list[i];

        auto corners = parse_tensor<dim>(e_json, "pts", except_text);

        if (!e_json.HasMember("bc_type") || !e_json["bc_type"].IsString()) {
            throw std::invalid_argument(except_text);
        }
        std::string bc_type = e_json["bc_type"].GetString();

        auto bc = parse_tensor<dim>(e_json, "bc", except_text);

        if (!e_json.HasMember("refine") || !e_json["refine"].IsInt()) {
            throw std::invalid_argument(except_text);
        }
        int n_refines = e_json["refine"].GetInt();

        out.push_back({corners, bc_type, bc, n_refines});
    }
    return out;
}

template 
std::vector<Element<2>> get_elements(const rapidjson::Document& doc);
template 
std::vector<Element<3>> get_elements(const rapidjson::Document& doc);

template <size_t dim>
MeshMap<dim> get_meshes(const std::vector<Element<dim>>& elements) {
    std::unordered_map<std::string,std::vector<Mesh<dim>>> facet_sets;
    for (auto e: elements) {
        auto facet = Mesh<dim>{{e.pts}};
        auto refined_facet = facet.refine_repeatedly(e.n_refines);
        facet_sets[e.bc_type].push_back(refined_facet);
    }

    std::vector<std::pair<std::string,Mesh<dim>>> meshes;
    for (auto it = facet_sets.begin(); it != facet_sets.end(); ++it) {
        auto union_mesh = Mesh<dim>::create_union(it->second);
        meshes.push_back(std::make_pair(it->first, union_mesh));
    }

    MeshMap<dim> out(meshes.begin(), meshes.end());

    for (const auto& type: get_mesh_types()) {
        if (out.find(type) == out.end()) {
            out.insert(std::make_pair(type, Mesh<dim>{{}}));
        };
    }

    return out;
}

template 
MeshMap<2> get_meshes(const std::vector<Element<2>>& elements);
template 
MeshMap<3> get_meshes(const std::vector<Element<3>>& elements);

template <size_t dim>
BCMap get_bcs(const std::vector<Element<dim>>& elements) {
    std::unordered_map<std::string,std::vector<Mesh<dim>>> bc_sets;
    for (auto e: elements) {
        auto bc = Mesh<dim>{{e.bc}};
        // BCs need to be refined to match up with the refined mesh.
        auto refined_bc = bc.refine_repeatedly(e.n_refines);
        bc_sets[e.bc_type].push_back(refined_bc);
    }

    BCMap bc_map;
    for (auto it = bc_sets.begin(); it != bc_sets.end(); ++it) {
        FieldDescriptor key{it->first, it->first};
        auto union_mesh = Mesh<dim>::create_union(it->second);
        bc_map[key].resize(dim);
        for (size_t d = 0; d < dim; d++) {
            bc_map[key][d].resize(union_mesh.n_dofs());
            for (auto it = union_mesh.begin(); it != union_mesh.end(); ++it) {
                const auto& vertex_val = *it;
                bc_map[key][d][it.absolute_index()] = vertex_val[d];
            }
        }
    }

    for (const auto& type: get_mesh_types()) {
        FieldDescriptor key{type, type};
        if (bc_map.find(key) == bc_map.end()) {
            bc_map.insert(std::make_pair(key, BlockFunction(dim, std::vector<double>(0))));
        };
    }

    return bc_map;
}

template 
BCMap get_bcs<2>(const std::vector<Element<2>>& elements);
template 
BCMap get_bcs<3>(const std::vector<Element<3>>& elements);

template <size_t dim>
std::vector<Vec<double,dim>> get_pts(const rapidjson::Document& doc) {
    auto n_pts = doc.Size();
    std::vector<Vec<double,dim>> out(n_pts);

    std::string except_text = "Points must be in the form [Double, Double]";

    for (std::size_t i = 0; i < n_pts; i++) {
        auto& pt_json = doc[i];
        if (!pt_json.IsArray() || pt_json.Capacity() != dim) {
            throw std::invalid_argument(except_text);
        }
        for (int d = 0; d < dim; d++) {
            if(!pt_json[d].IsDouble()) {
                throw std::invalid_argument(except_text);
            }
            out[i][d] = pt_json[d].GetDouble();
        }
    }
    return out;
}
template 
std::vector<Vec<double,2>> get_pts(const rapidjson::Document& doc);

template 
std::vector<Vec<double,3>> get_pts(const rapidjson::Document& doc);
