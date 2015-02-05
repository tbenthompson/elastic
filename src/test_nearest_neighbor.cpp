#include "UnitTest++.h"
#include "3bem/mesh.h"
#include "3bem/vec_ops.h"
#include "nearest_neighbor.h"

using namespace tbem;

TEST(NearestPointOneFacetEndpoint) {
    Facet<2> f{{{1,1},{2,1}}};
    Mesh<2> m{{f}};
    auto result = nearest_pt({0, 0}, m);
    CHECK_EQUAL(result.facets[0], m.facets[0]);
    CHECK_EQUAL(result.pt, m.facets[0][0]);
    CHECK_EQUAL(result.distance, std::sqrt(2));
}

TEST(NearestPointOneFacetNotEndpoint) {
    Facet<2> f{{{-1,-1},{1,-1}}};
    Mesh<2> m{{f}};
    auto result = nearest_pt({0, 0}, m);
    CHECK_EQUAL(result.facets[0], m.facets[0]);
    CHECK_EQUAL(result.pt, (Vec<double,2>{0,-1}));
    CHECK_EQUAL(result.distance, 1.0);
}

TEST(NearestPointTwoFacetsNotEndpoint) {
    Facet<2> f{{{0.5,-1},{1,-1}}};
    Facet<2> f2{{{-1,-1},{0.5,-1}}};
    Mesh<2> m{{f,f2}};
    auto result = nearest_pt({0, 0}, m);
    CHECK_EQUAL(result.facets[0], m.facets[1]);
    CHECK_EQUAL(result.pt, (Vec<double,2>{0,-1}));
    CHECK_EQUAL(result.distance, 1.0);
}

TEST(NearestPointIntersection) {
    Facet<2> f{{{0,0},{1,0}}};
    Facet<2> f2{{{0,1},{0,0}}};
    Mesh<2> m{{f,f2}};
    auto result = nearest_pt({0, 0}, m);
    CHECK_EQUAL(result.facets.size(), 2);
    CHECK(result.facets[0] == f || result.facets[1] == f);
    CHECK(result.facets[0] == f2 || result.facets[1] == f2);
    CHECK_EQUAL(result.pt, (Vec<double,2>{0,0}));
    CHECK_EQUAL(result.distance, 0.0);
}

TEST(DecideRichardsonDirFar) {
    Facet<2> f{{{0,-1},{1,-1}}};
    Vec<double,2> p{0,1};
    NearestPoint<2> np{{f}, p, 0.0};
    auto dir = decide_richardson_dir(p, np);
    CHECK_EQUAL(dir, (Vec<double,2>{0, 1}));
}

TEST(DecideRichardsonDirTouching) {
    Facet<2> f{{{0,0},{0,1}}};
    Vec<double,2> p{0,0};
    NearestPoint<2> np{{f}, p, 0.0};
    auto dir = decide_richardson_dir(p, np);
    CHECK_EQUAL(dir, (Vec<double,2>{-1, 0}));
}

TEST(DecideRichardsonDirOppositeNormal) {
    Facet<2> f{{{0,-1},{1,-1}}};
    Vec<double,2> p{0,-2};
    NearestPoint<2> np{{f}, p, 0.0};
    auto dir = decide_richardson_dir(p, np);
    CHECK_EQUAL(dir, (Vec<double,2>{0, -1}));
}

TEST(DecideRichardsonDirIntersection) {
    Facet<2> f{{{0,0},{1,0}}};
    Facet<2> f2{{{0,1},{0,0}}};
    Vec<double,2> p{0,0};
    NearestPoint<2> np{{f,f2}, p, 0.0};
    auto dir = decide_richardson_dir(p, np);
    CHECK_EQUAL(dir, (Vec<double,2>{0.5, 0.5}));
}
