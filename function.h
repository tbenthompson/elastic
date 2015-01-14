#ifndef __AQweQWEJL_FUNCTION_H
#define __AQweQWEJL_FUNCTION_H

#include "3bem/mesh.h"
#include "3bem/constraint.h"

// functions are made from data or from operators applied to other functions
// Replace with typedef of MeshField?
template <size_t dim>
struct Function
{
    const tbem::Mesh<dim>& mesh;
    const std::vector<tbem::Vec<double,dim>> values;

    // operations like these should assert mesh == other.mesh
    Function<dim> operator+(const Function<dim>& other);
    Function<dim> operator-(const Function<dim>& other);
};


template <size_t dim>
struct ReducedFunction {
    const size_t expanded_n_dofs;   
    const tbem::ConstraintMatrix constraint_matrix; 
    const std::vector<tbem::Vec<double,dim>> values;

    // ReducedFunction<dim> scale_rows(double factor);
};

template <size_t dim>
struct BlockFunction 
{
    const std::vector<Function<dim>> functions;
};

template <size_t dim>
struct ReducedBlockFunction
{
    const std::vector<ReducedFunction<dim>> functions;

    // BlockFunction<dim> scale_rows(double factor) const;
    std::vector<double> stack_functions() const;

    static ReducedBlockFunction<dim> unstack_functions(const std::vector<double>& rows);
};

#endif
