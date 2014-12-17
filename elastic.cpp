// #include "bem.h"
// #include "mesh_gen.h"
// #include "quadrature.h"
// #include "kernels.h"
// #include "petsc_interface.h"
// #include "util.h"
#include "vec.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include "elastic.h"

// A handy debugging function.
template <typename T>
void print_members(const T& obj) {
    static const char* kTypeNames[] = { 
        "Null", "False", "True", "Object",
        "Array", "String", "Number" 
    };
    for (rapidjson::Value::ConstMemberIterator itr = obj.MemberBegin();
         itr != obj.MemberEnd();
         ++itr) 
    {
        printf("Type of member %s is %s\n",
               itr->name.GetString(),
               kTypeNames[itr->value.GetType()]);
    }
} 

std::string load_file(std::string filename) {
    std::ifstream file;
    file.open(filename);
    if (!file.is_open()) {
        /* std::cout << "Filename: " << filename << " does not exist." << std::endl; */
        // TODO: Better way to exit from deeper in the stack.
        throw std::invalid_argument("Filename: " + filename + " could not be" + 
                                    " opened (does it exist?).");
    }

    // Read the file into a string
    std::stringstream buffer;
    buffer << file.rdbuf();

    return buffer.str();
}

rapidjson::Document parse_json(std::string json) {
    rapidjson::Document doc;
    bool parse_error = doc.Parse(json.c_str()).HasParseError();
    if (parse_error) {
        /* std::cout << "Filename: " << filename << " has JSON errors." << std::endl; */
        // TODO: Better way to exit from deeper in the stack.
        throw std::invalid_argument("There are json errors.");
    }
    return doc;
}
// 
// 
//     //TODO: Element input format:
//     //{pt1: {x1,y1,z1}, pt2: {x2,y2,z2}, bc_type, bc_val: {bcx, bcy, bcz}}
//     const auto& element_list = doc["elements"];
//     for (std::size_t i = 0; i < element_list.Size(); i++) {
//         const auto& pts = element_list[i]["pts"];
//         //TODO: Better error handling for reading in file.
//         const Vec2<Vec2<double>> corners = {{
//             {pts[0][0].GetDouble(), pts[0][1].GetDouble()},
//             {pts[1][0].GetDouble(), pts[1][1].GetDouble()}
//         }};
// 
//         std::string bc_type_str = element_list[i]["bc_type"].GetString();
//         BC bc_type = DISPLACEMENT;
//         if (bc_type_str == "displacement") {
//             bc_type = DISPLACEMENT;
//         } else if (bc_type_str == "traction") {
//             bc_type = TRACTION;
//         } else if (bc_type_str == "slip") {
//             bc_type = SLIP;
//         } else if (bc_type_str == "crack") {
//             bc_type = CRACK;
//         } else {
//             //TODO: Error handling here.
//         }
// 
//         const auto& bc_element = element_list[i]["bc"];
//         Vec2<double> bc = {bc_element[0].GetDouble(), bc_element[1].GetDouble()};
// 
//         std::cout << corners << std::endl;
//         std::cout << bc_type << std::endl;
//         std::cout << bc << std::endl;
// 
//         // std::cout << element_list[i]["pts"][0][0].GetDouble() << std::endl;
//         // assert(element_list[i].IsArray());
//     }
// 
//     std::cout << "Finished reading file." << std::endl;
// 
//     file.close();
// 
//     ElasticProblem<2> ep;
//     return ep;
// }

ConstraintMatrix build_constraints(const ElasticProblem<2>& elastic_prob) {
//     // Impose continuity constraints on the output displacement field.
//     // The displacement field is unknown on the traction_mesh.
//     auto displacement_continuity = ConstraintMatrix::from_constraints(
//         mesh_continuity(elastic_prob.traction_mesh)
//         );
// 
//     // Remove continuity constraints at the intersection of the fault and the 
//     // surface mesh.
//     // TODO: Here, the problem is more complex than in the past because I 
//     //       should allow the vertices not to match up.
//     auto constraints = apply_discontinuities(
//         elastic_prob.traction_mesh,
//         elastic_prob.slip_mesh,
//         displacement_continuity
//         );
// 
//     return constraints;
}
