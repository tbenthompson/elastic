#ifndef __QWENMBBBZBLAHEWJPQ_ELASTIC_H
#define __QWENMBBBZBLAHEWJPQ_ELASTIC_H

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

const Parameters default_params{2, 2, 6, 3.0, 1e-2, 0.25, 30e9};

template <size_t dim>
struct MeshesAndBCs {
    tbem::Mesh<dim> displacement_mesh;
    tbem::Mesh<dim> traction_mesh;
    tbem::Mesh<dim> slip_mesh;

    tbem::Mesh<dim> displacement_bcs;
    tbem::Mesh<dim> traction_bcs;
    tbem::Mesh<dim> slip_bcs;
};

enum BCType {DISPLACEMENT, TRACTION, SLIP, CRACK};

std::string load_file(std::string filename);
rapidjson::Document parse_json(std::string json);
Parameters get_parameters(const rapidjson::Document& doc);
MeshesAndBCs<2> get_meshes_bcs(const rapidjson::Document& doc);


#endif
