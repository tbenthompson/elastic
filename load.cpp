#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include "load.h"
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

Vec2<Vec2<double>> parse_tensor(const rapidjson::Value& e_json, 
                                std::string field_name,
                                std::string except_text) {
    if (!e_json.HasMember(field_name.c_str())) {
        throw std::invalid_argument(except_text);
    }

    const auto& V = e_json[field_name.c_str()];
    if (!V.IsArray() || V.Capacity() != 2) {
        throw std::invalid_argument(except_text);
    }
    for (int d0 = 0; d0 < 2; d0++) {
        if (!V[d0].IsArray() || V[d0].Capacity() != 2) {
            throw std::invalid_argument(except_text);
        }
        for (int d1 = 0; d1 < 2; d1++) {
            if(!V[d0][d1].IsDouble()) {
                throw std::invalid_argument(except_text);
            }
        }    
    }

    return {{
        {{V[0][0].GetDouble(), V[0][1].GetDouble()}},
        {{V[1][0].GetDouble(), V[1][1].GetDouble()}}
    }};
}

const Parameters default_params{2, 2, 6, 3.0, 1e-2, 0.25, 30e9};
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

std::vector<Element<2>> get_elements(const rapidjson::Document& doc) {
    const auto& element_list = doc["elements"];

    std::vector<Element<2>> out;

    std::string except_text = "An element object must have a (dim x dim) array of float";
    except_text += "'pts' and 'bc' array, a string 'bc_type' field, and an integer ";
    except_text += "'refine' field.";

    for (std::size_t i = 0; i < element_list.Size(); i++) {
        auto& e_json = element_list[i];

        auto corners = parse_tensor(e_json, "pts", except_text);

        if (!e_json.HasMember("bc_type") || !e_json["bc_type"].IsString()) {
            throw std::invalid_argument(except_text);
        }
        std::string bc_type = e_json["bc_type"].GetString();

        auto bc = parse_tensor(e_json, "bc", except_text);

        if (!e_json.HasMember("refine") || !e_json["refine"].IsInt()) {
            throw std::invalid_argument(except_text);
        }

        int n_refines = e_json["refine"].GetInt();
        Element<2> e{corners, bc_type, bc, n_refines};
        out.push_back(e);
    }
    return out;
}

template <size_t dim>
MeshSet<dim> get_meshes(const std::vector<Element<dim>>& elements) {
    std::map<std::string,std::vector<Mesh<dim>>> facet_sets;
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

    MeshSet<dim> out(meshes.begin(), meshes.end());
    return out;
}

template 
MeshSet<2> get_meshes(const std::vector<Element<2>>& elements);
template 
MeshSet<3> get_meshes(const std::vector<Element<3>>& elements);

template <size_t dim>
BCSet get_bcs(const std::vector<Element<dim>>& elements) {
    BCSet bc_sets;
    for (auto e: elements) {
        auto bc = Mesh<dim>{{e.bc}};
        // BCs need to be refined to match up with the refined mesh.
        auto refined_bc = bc.refine_repeatedly(e.n_refines);
        bc_sets[e.bc_type].resize(dim);
        for (size_t d = 0; d < dim; d++) {
            bc_sets[e.bc_type][d].resize(refined_bc.n_dofs());
            for (auto it = refined_bc.begin(); it != refined_bc.end(); ++it) {
                const auto& vertex_val = *it;
                bc_sets[e.bc_type][d][it.absolute_index()] = vertex_val[d];
            }
        }
    }

    return bc_sets;
}

template 
BCSet get_bcs<2>(const std::vector<Element<2>>& elements);
template 
BCSet get_bcs<3>(const std::vector<Element<3>>& elements);
