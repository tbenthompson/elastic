#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include "elastic.h"
#include "3bem/vec.h"

using namespace tbem;
std::string load_file(std::string filename) {
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

rapidjson::Document parse_json(std::string json) {
    rapidjson::Document doc;
    bool parse_error = doc.Parse(json.c_str()).HasParseError();
    if (parse_error) {
        throw std::invalid_argument("There are json errors.");
    }
    return doc;
}

BCType parse_bc_type(std::string bc_type_str) {
    BCType bc_type = DISPLACEMENT;
    if (bc_type_str == "displacement") {
        bc_type = DISPLACEMENT;
    } else if (bc_type_str == "traction") {
        bc_type = TRACTION;
    } else if (bc_type_str == "slip") {
        bc_type = SLIP;
    } else if (bc_type_str == "crack") {
        bc_type = CRACK;
    } else {
        throw std::invalid_argument("bc_type must be one of \
                                     (displacement, traction, slip, crack)");
    }
    return bc_type;
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
        {V[0][0].GetDouble(), V[0][1].GetDouble()},
        {V[1][0].GetDouble(), V[1][1].GetDouble()}
    }};
}

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
struct Element {
    tbem::Vec<tbem::Vec<double,dim>,dim> pts; 
    BCType bc_type;
    tbem::Vec<tbem::Vec<double,dim>,dim> bc;
    int n_refines;
};

std::vector<Element<2>> get_elements(const rapidjson::Document& doc) {
    const auto& element_list = doc["elements"];

    std::vector<Element<2>> out;

    std::string except_text = "An element object must have a (dim x dim) \
                'pts' and 'bc' array, a string 'bc_type' field, and an integer \
                'refine' field.";

    for (std::size_t i = 0; i < element_list.Size(); i++) {
        auto& e_json = element_list[i];

        auto corners = parse_tensor(e_json, "pts", except_text);

        if (!e_json.HasMember("bc_type") || !e_json["bc_type"].IsString()) {
            throw std::invalid_argument(except_text);
        }
        std::string bc_type_str = e_json["bc_type"].GetString();
        BCType bc_type = parse_bc_type(bc_type_str);

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

MeshesAndBCs<2> get_meshes_bcs(const rapidjson::Document& doc) {
    auto elements = get_elements(doc);

    //One facet list per BC type.
    std::vector<Mesh<2>> facet_lists[4]; 
    std::vector<Mesh<2>> bc_lists[4];

    for (auto e: elements) {
        auto facet = Mesh<2>{{Facet<2>{e.pts}}};
        auto refined_facet = facet.refine_repeatedly(e.n_refines);
        facet_lists[e.bc_type].push_back(refined_facet);

        auto bc = Mesh<2>{{Facet<2>{e.bc}}};
        auto refined_bc = bc.refine_repeatedly(e.n_refines);
        bc_lists[e.bc_type].push_back(refined_bc);
    }

    return {
        Mesh<2>::form_union(facet_lists[DISPLACEMENT]),
        Mesh<2>::form_union(facet_lists[TRACTION]),
        Mesh<2>::form_union(facet_lists[SLIP]),
        Mesh<2>::form_union(bc_lists[DISPLACEMENT]),
        Mesh<2>::form_union(bc_lists[TRACTION]),
        Mesh<2>::form_union(bc_lists[SLIP])
    };
}
