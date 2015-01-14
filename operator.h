#ifndef __AQQQHHQHQHHEKLWLKEJR_OPERATOR_H
#define __AQQQHHQHQHHEKLWLKEJR_OPERATOR_H
#include "function.h"
#include "3bem/quadrature.h"

template <size_t dim>
struct Operator 
{
    const tbem::Mesh<dim>& src_mesh;
    const tbem::Mesh<dim>& obs_mesh;
    const std::vector<tbem::Vec<tbem::Vec<double,dim>,dim>> matrix; 

    //this function should check that obs_mesh == fnc.mesh
    Function<dim> operator()(const Function<dim>& fnc);

    template <typename KT>
    static Operator<dim> form_operator(
        const tbem::Mesh<dim>& src_mesh, const tbem::Mesh<dim>& obs_mesh,
        const KT& K, const tbem::QuadStrategy<dim>& quad_strategy);
};

template <size_t dim>
struct BlockOperator 
{
    const std::vector<Operator<dim>> operators;
    const size_t n_rows;
    const size_t n_cols;

    BlockFunction<dim> operator()(const BlockFunction<dim>& fnc);
};


#endif
