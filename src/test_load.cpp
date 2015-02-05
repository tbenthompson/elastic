#include "UnitTest++.h"
#include <iostream>
#include "load.h"
#include "spec.h"
#include "3bem/util.h"
#include "3bem/vec_ops.h"

using namespace tbem;

TEST(SimpleLoadGood) {
    auto doc = parse_json(load_file("test_data/good.in"));
}

TEST(SimpleLoadBad) {
    CHECK_THROW(parse_json(load_file("test_data/bad.in")), std::invalid_argument);
    CHECK_THROW(parse_json(load_file("doesnotexist")), std::invalid_argument);
}

TEST(RapidJSONBigFile) {
    TIC
    auto file = load_file("test_data/reallybig.in");
    TOC("Loading file");
    TIC2
    auto doc = parse_json(file);
    TOC("Parsing file");
    TIC2
    auto meshes = get_meshes(get_elements<2>(doc));
    int n_elements = meshes["traction"].facets.size();
    TOC("Meshing " + std::to_string(n_elements) + " elements");
}


TEST(GetMeshes) {
    auto doc = parse_json(load_file("test_data/one.in"));
    auto meshes = get_meshes(get_elements<2>(doc));
    CHECK_EQUAL(meshes["traction"].facets.size(), 1);
    CHECK_EQUAL(meshes["traction"].facets[0][0],
                (Vec2<double>{-20.0, 0.0}));
    CHECK_EQUAL(meshes["slip"].facets.size(), 1);
    CHECK_EQUAL(meshes["slip"].facets[0][0],
                (Vec2<double>{-2.0, -2.0}));
}

TEST(AllMeshesCreated) {
    auto doc = parse_json(load_file("test_data/empty.in"));
    auto meshes = get_meshes(get_elements<2>(doc));
    for (const auto& name: get_mesh_types()) {
        CHECK(meshes.find(name) != meshes.end());
    }
}

TEST(AllBCs) {
    auto doc = parse_json(load_file("test_data/empty.in"));
    auto bcs = get_bcs(get_elements<2>(doc));
    for (const auto& name: get_bc_types()) {
        CHECK(bcs.find(FieldDescriptor{name,name}) != bcs.end());
    }
}

TEST(Refinement) {
    auto doc = parse_json(load_file("test_data/refine.in"));
    auto meshes = get_meshes(get_elements<2>(doc));
    CHECK_EQUAL(meshes["traction"].facets.size(), 8);
    for (int i = 0; i < 8; i++) {
        CHECK_EQUAL(meshes["traction"].facets[i][0],
                    (Vec2<double>{(double)i, -(double)i}));
    }
}

TEST(BCRefinement) {
    auto doc = parse_json(load_file("test_data/refine.in"));
    auto bcs = get_bcs(get_elements<2>(doc));
    for (int i = 0; i < 8; i++) {
        CHECK_EQUAL(bcs[(FieldDescriptor{"traction", "traction"})][0][2 * i], i);
        CHECK_EQUAL(bcs[(FieldDescriptor{"traction", "traction"})][1][2 * i], -i);
    }
}

TEST(MalformedElementException) {
    CHECK_THROW(get_elements<2>(parse_json(load_file("test_data/bad2.in"))),
                std::invalid_argument);
    CHECK_THROW(get_elements<2>(parse_json(load_file("test_data/bad3.in"))),
                std::invalid_argument);
    CHECK_THROW(get_elements<2>(parse_json(load_file("test_data/bad4.in"))),
                std::invalid_argument);
    CHECK_THROW(get_elements<2>(parse_json(load_file("test_data/bad5.in"))),
                std::invalid_argument);
}

TEST(LoadParametersDefault) {
    auto doc = parse_json(load_file("test_data/one.in"));
    auto p = get_parameters(doc);
    // All should be defaults.
    CHECK_EQUAL(p.obs_quad_order, default_params.obs_quad_order);
    CHECK_EQUAL(p.src_far_quad_order, default_params.src_far_quad_order);
    CHECK_EQUAL(p.n_singular_steps, default_params.n_singular_steps);
    CHECK_EQUAL(p.far_threshold, default_params.far_threshold);
    CHECK_EQUAL(p.near_tol, default_params.near_tol);
    CHECK_EQUAL(p.solver_tol, default_params.solver_tol);
    CHECK_EQUAL(p.poisson_ratio, default_params.poisson_ratio);
    CHECK_EQUAL(p.shear_modulus, default_params.shear_modulus);
}

TEST(LoadParametersNotDefault) {
    auto doc = parse_json(load_file("test_data/params.in"));
    auto p = get_parameters(doc);
    CHECK_EQUAL(p.obs_quad_order, 4); CHECK_EQUAL(p.src_far_quad_order, 6);
    CHECK_EQUAL(p.n_singular_steps, 9); CHECK_EQUAL(p.far_threshold, 4.0);
    CHECK_EQUAL(p.near_tol, 1e-4);
    CHECK_EQUAL(p.solver_tol, 1e-6);
    CHECK_EQUAL(p.poisson_ratio, 0.28);
    CHECK_EQUAL(p.shear_modulus, 22e9);
}

TEST(Load3D) {
    auto doc = parse_json(load_file("test_data/3d_test.in"));
    auto meshes = get_meshes(get_elements<3>(doc));
    CHECK_EQUAL(meshes["traction"].n_facets(), 64);
}

TEST(LoadPts) {
    auto doc = parse_json(load_file("test_data/pts.in"));
    auto points = get_pts<2>(doc);
    CHECK_EQUAL(points.size(), 2);
    CHECK_ARRAY_EQUAL(&points[0], 
        (Vec<Vec<double,2>,2>{{{{0.0, 1.0}}, {{1.0, 0.0}}}}), 2);
}

TEST(LoadPtsBad) {
    auto doc = parse_json(load_file("test_data/bad_pts.in"));
    CHECK_THROW(get_pts<2>(doc), std::invalid_argument);
}

int main() {
    return UnitTest::RunAllTests();
}
