#ifndef __RRRTYKHJASK_SYSTEM_H
#define __RRRTYKHJASK_SYSTEM_H
#include "3bem/3bem.h"

template <size_t dim>
struct BEM {
    const MeshMap<dim> meshes;
    const BCMap bcs;
    const KernelMap<dim> kernels;
    const tbem::QuadStrategy<dim> quad_strategy;
    const std::vector<IntegralEquationSpec> eqtn_specs;
    const std::vector<tbem::ConstraintMatrix> constraints;
};


#endif
