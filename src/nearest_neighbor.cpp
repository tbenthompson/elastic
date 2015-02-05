#include <limits>
#include "3bem/mesh.h"
#include "3bem/vec_ops.h"
#include "3bem/vertex_iterator.h"
#include "nearest_neighbor.h"

using namespace tbem;

Vec<double,3> closest_pt_seg(const Vec<double,3>& pt, const Facet<3>& seg) {

    return seg[0];
}

Vec<double,2> closest_pt_seg(const Vec<double,2>& pt, const Facet<2>& seg)
{
     auto v = seg[1] - seg[0];
     auto w = pt - seg[0];

     double c1 = dot_product(w, v);
     if (c1 <= 0) {
          return seg[0];
     }

     double c2 = dot_product(v, v);
     if (c2 <= c1) {
          return seg[1];
     }

     double b = c1 / c2;
     return seg[0] + b * v;
}

template <size_t dim>
NearestPoint<dim> nearest_pt(const Vec<double,dim>& pt, const Mesh<dim>& mesh) 
{
    std::vector<Facet<dim>> closest_facets;
    Vec<double,dim> closest_pt;
    auto min_dist2 = std::numeric_limits<double>::max();
    for (size_t facet_idx = 0; facet_idx < mesh.n_facets(); facet_idx++) {
        auto mesh_pt = closest_pt_seg(pt, mesh.facets[facet_idx]);
        auto dist2_to_mesh = dist2(pt, mesh_pt);
        if (dist2_to_mesh < min_dist2) {
            closest_facets.clear();
            closest_facets.push_back(mesh.facets[facet_idx]);
            min_dist2 = dist2_to_mesh;
            closest_pt = mesh_pt;
        } else if (dist2_to_mesh == min_dist2) {
            closest_facets.push_back(mesh.facets[facet_idx]);
        }
    }
    return {closest_facets, closest_pt, std::sqrt(min_dist2)};
}

template 
NearestPoint<2> nearest_pt(const Vec<double,2>& pt, const Mesh<2>& mesh);
template 
NearestPoint<3> nearest_pt(const Vec<double,3>& pt, const Mesh<3>& mesh);

template <size_t dim>
Vec<double,dim> decide_richardson_dir(const Vec<double,dim>& pt,
    const NearestPoint<dim>& mesh_pt) 
{
    auto facet_normal = zeros<Vec<double,dim>>::make();
    for (size_t i = 0; i < mesh_pt.facets.size(); i++) {
        facet_normal += unscaled_normal(mesh_pt.facets[i]); 
    }
    facet_normal /= (double)mesh_pt.facets.size();

    auto which_side = which_side_point(mesh_pt.facets[0], pt);
    if (which_side == INTERSECT) {
        return facet_normal;
    } else if (which_side == FRONT) {
        return facet_normal;
    } else {
        return -facet_normal;
    }
}

template Vec<double,2> 
decide_richardson_dir(const Vec<double,2>& pt, const NearestPoint<2>& mesh_pt); 
template Vec<double,3>
decide_richardson_dir(const Vec<double,3>& pt, const NearestPoint<3>& mesh_pt);
