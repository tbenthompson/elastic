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
    TIC
    auto doc = parse_json(load_file("test_data/reallybig.in"));
    auto elements = collect_elements(doc);
    TOC("Parsing " + std::to_string(elements.size()) + " elements.");
}

TEST(BuildProblem) {
    auto doc = parse_json(load_file("test_data/one.in"));
    auto elements = collect_elements(doc);
    auto elast_prob = build_problem(doc, elements);
    CHECK_EQUAL(elast_prob.traction_mesh.facets.size(), 1);
    CHECK_EQUAL(elast_prob.traction_mesh.facets[0].vertices[0],
                (Vec2<double>{-20.0, 0.0}));
    CHECK_EQUAL(elast_prob.slip_mesh.facets.size(), 1);
    CHECK_EQUAL(elast_prob.slip_mesh.facets[0].vertices[0],
                (Vec2<double>{-2.0, -2.0}));
}


int main() {
    return UnitTest::RunAllTests();
}
