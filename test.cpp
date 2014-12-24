#include "UnitTest++.h"
#include <iostream>
#include "elastic.h"
#include "3bem/util.h"

using namespace tbem;

TEST(SimpleLoadGood) {
    auto doc = parse_json(load_file("test_data/good.in"));
}

TEST(SimpleLoadBad) {
    CHECK_THROW(parse_json(load_file("test_data/bad.in")), std::invalid_argument);
    CHECK_THROW(parse_json(load_file("doesnotexist")), std::invalid_argument);
}

TEST(LoadElements) {
    auto doc = parse_json(load_file("test_data/one.in"));
    auto elements = collect_elements(doc);
    CHECK_EQUAL(elements[0].pts[0], (Vec2<double>{-20.0,0.0}));
    CHECK_EQUAL(elements[0].pts[1], (Vec2<double>{20.0,0.0}));
    CHECK_EQUAL(elements[0].bc_type, TRACTION);
    CHECK_EQUAL(elements[0].bc[0], (Vec2<double>{0.0,0.0}));
    CHECK_EQUAL(elements[0].bc[1], (Vec2<double>{0.0,0.0}));
    CHECK_EQUAL(elements[1].pts[0], (Vec2<double>{-2.0,-2.0}));
    CHECK_EQUAL(elements[1].pts[1], (Vec2<double>{0.0,0.0}));
    CHECK_EQUAL(elements[1].bc_type, SLIP);
    CHECK_EQUAL(elements[1].bc[0], (Vec2<double>{1.0,1.0}));
    CHECK_EQUAL(elements[1].bc[1], (Vec2<double>{1.0,1.0}));
}

TEST(RapidJSONBigFile) {
    UNITTEST_TIME_CONSTRAINT(100);
    TIC
    auto file = load_file("test_data/reallybig.in");
    TOC("Loading file");
    TIC2
    auto doc = parse_json(file);
    TOC("Parsing file");
    TIC2
    auto elements = collect_elements(doc);
    TOC("Collecting " + std::to_string(elements.size()) + " elements");
    TIC2
    auto prob = build_problem(elements);
    TOC("Meshing " + std::to_string(elements.size()) + " elements");
}

TEST(BuildProblem) {
    auto doc = parse_json(load_file("test_data/one.in"));
    auto elements = collect_elements(doc);
    auto elast_prob = build_problem(elements);
    CHECK_EQUAL(elast_prob.traction_mesh.facets.size(), 1);
    CHECK_EQUAL(elast_prob.traction_mesh.facets[0].vertices[0],
                (Vec2<double>{-20.0, 0.0}));
    CHECK_EQUAL(elast_prob.slip_mesh.facets.size(), 1);
    CHECK_EQUAL(elast_prob.slip_mesh.facets[0].vertices[0],
                (Vec2<double>{-2.0, -2.0}));
}

TEST(Refinement) {
    auto doc = parse_json(load_file("test_data/refine.in"));
    auto elements = collect_elements(doc);
    auto elast_prob = build_problem(elements);
    CHECK_EQUAL(elast_prob.traction_mesh.facets.size(), 8);
    for (int i = 0; i < 8; i++) {
        CHECK_EQUAL(elast_prob.traction_mesh.facets[i].vertices[0],
                    (Vec2<double>{(double)i, -(double)i}));
        CHECK_EQUAL(elast_prob.traction_bcs.facets[i].vertices[0],
                    (Vec2<double>{(double)i, -(double)i}));
    }
}

TEST(MalformedElementException) {
    CHECK_THROW(collect_elements(parse_json(load_file("test_data/bad2.in"))),
                std::invalid_argument);
    CHECK_THROW(collect_elements(parse_json(load_file("test_data/bad3.in"))),
                std::invalid_argument);
    CHECK_THROW(collect_elements(parse_json(load_file("test_data/bad4.in"))),
                std::invalid_argument);
    CHECK_THROW(collect_elements(parse_json(load_file("test_data/bad5.in"))),
                std::invalid_argument);
}

TEST(LoadParametersDefault) {
    auto doc = parse_json(load_file("test_data/one.in"));
    auto p = parse_parameters(doc);
    // All should be defaults.
    CHECK_EQUAL(p.obs_quad_order, 2); CHECK_EQUAL(p.src_far_quad_order, 2);
    CHECK_EQUAL(p.n_singular_steps, 6); CHECK_EQUAL(p.far_threshold, 3.0);
    CHECK_EQUAL(p.near_tol, 1e-2);
    CHECK_EQUAL(p.poisson_ratio, 0.25);
    CHECK_EQUAL(p.shear_modulus, 30e9);
}
TEST(LoadParametersNotDefault) {
    auto doc = parse_json(load_file("test_data/params.in"));
    auto p = parse_parameters(doc);
    CHECK_EQUAL(p.obs_quad_order, 4); CHECK_EQUAL(p.src_far_quad_order, 6);
    CHECK_EQUAL(p.n_singular_steps, 9); CHECK_EQUAL(p.far_threshold, 4.0);
    CHECK_EQUAL(p.near_tol, 1e-4);
    CHECK_EQUAL(p.poisson_ratio, 0.28);
    CHECK_EQUAL(p.shear_modulus, 22e9);
}

int main() {
    return UnitTest::RunAllTests();
}
