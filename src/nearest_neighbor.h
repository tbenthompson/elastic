#ifndef __QWELJLADJJN_NEAREST_NEIGHBOR_H
#define __QWELJLADJJN_NEAREST_NEIGHBOR_H

#include "3bem/vec.h"

namespace tbem {
    template <size_t dim>
    struct Mesh;
}

template <size_t dim>
struct NearestPoint {
    const std::vector<tbem::Facet<dim>> facets;
    const tbem::Vec<double,dim> pt;
    const double distance;
};

template <size_t dim>
NearestPoint<dim>
nearest_pt(const tbem::Vec<double, dim>& pt, const tbem::Mesh<dim>& mesh);


template <size_t dim> tbem::Vec<double,dim>
decide_richardson_dir(const tbem::Vec<double,dim>& pt, const NearestPoint<dim>& mesh_pt);

#endif
