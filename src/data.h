#ifndef __reakWHLA123Jjkla_BEM_BUILDER_H
#define __reakWHLA123Jjkla_BEM_BUILDER_H

#include "load.h"
#include "kernels.h"
#include "3bem/quadrature.h"

template <size_t dim>
struct IntegralEquationSpec;

template <size_t dim>
struct BEM {
    const Parameters params;
    const MeshMap<dim> meshes;
    const BCMap bcs;
    const KernelMap<dim> kernels;
    const tbem::QuadStrategy<dim> quad_strategy;
    const std::vector<IntegralEquationSpec<dim>> eqtn_specs;
};

template <size_t dim>
BEM<dim> parse_into_bem(const std::string& filename);

#endif
