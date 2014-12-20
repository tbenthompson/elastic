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

std::vector<Element<2>> collect_elements(const rapidjson::Document& doc) {
    const auto& element_list = doc["elements"];

    std::vector<Element<2>> out;

    for (std::size_t i = 0; i < element_list.Size(); i++) {
        const auto& pts = element_list[i]["pts"];
        //TODO: Better error handling for reading in file.
        //TODO: How to check if entry exists in rapidjson?
        const Vec2<Vec2<double>> corners = {{
            {pts[0][0].GetDouble(), pts[0][1].GetDouble()},
            {pts[1][0].GetDouble(), pts[1][1].GetDouble()}
        }};

        std::string bc_type_str = element_list[i]["bc_type"].GetString();
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

        const auto& bc_element = element_list[i]["bc"];
        Vec2<Vec2<double>> bc = {{
            {bc_element[0][0].GetDouble(), bc_element[0][1].GetDouble()},
            {bc_element[1][0].GetDouble(), bc_element[1][1].GetDouble()}
        }};
        Element<2> e{corners, bc_type, bc};
        out.push_back(e);
    }
    return out;
}

ElasticProblem<2> build_problem(const rapidjson::Document& doc, 
                                const std::vector<Element<2>>& elements) {
    //One facet list per BC type.
    std::vector<Facet<2>> facet_lists[4]; 
    std::vector<Vec<Vec<double,2>,2>> bc_lists[4];
    for (auto e: elements) {
        facet_lists[e.bc_type].push_back(Facet<2>{e.pts});
        bc_lists[e.bc_type].push_back(e.bc);
    }

    return {
        Mesh<2>{facet_lists[DISPLACEMENT], false, nullptr},
        Mesh<2>{facet_lists[TRACTION], false, nullptr},
        Mesh<2>{facet_lists[SLIP], false, nullptr},
        bc_lists[DISPLACEMENT],
        bc_lists[TRACTION],
        bc_lists[SLIP]
    };
}
