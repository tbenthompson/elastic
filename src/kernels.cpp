#include "kernels.h"
#include "3bem/elastic_kernels.h"

template <size_t dim>
KernelMap<dim> get_elastic_kernels(double shear_modulus, double poisson_ratio) {
    KernelMap<dim> out;
    out.insert(std::make_pair(
        "displacement",
        KernelPtr<dim>(
            new tbem::ElasticDisplacement<dim>(shear_modulus, poisson_ratio)
        ))
    );
    out.insert(std::make_pair(
        "traction",
        KernelPtr<dim>(
            new tbem::ElasticTraction<dim>(shear_modulus, poisson_ratio)
        ))
    );
    out.insert(std::make_pair(
        "adjoint_traction",
        KernelPtr<dim>(
            new tbem::ElasticAdjointTraction<dim>(shear_modulus, poisson_ratio)
        ))
    );
    out.insert(std::make_pair(
        "hypersingular",
        KernelPtr<dim>(
            new tbem::ElasticHypersingular<dim>(shear_modulus, poisson_ratio)
        ))
    );
    return out;
}

template 
KernelMap<2> get_elastic_kernels(double shear_modulus, double poisson_ratio);
template 
KernelMap<3> get_elastic_kernels(double shear_modulus, double poisson_ratio);

