#ifndef __aaaaaaaaaaaaa_KERNELS_H
#define __aaaaaaaaaaaaa_KERNELS_H

#include <memory>
#include <map>
#include "3bem/kernel.h"

template <size_t dim>
using KernelPtr = 
    std::shared_ptr<tbem::Kernel<dim,
        tbem::Vec<double,dim>,
        tbem::Vec<double,dim>,
        tbem::Vec<tbem::Vec<double,dim>,dim>
    >>;

template <size_t dim>
using KernelMap = std::map<std::string, KernelPtr<dim>>;

template <size_t dim>
KernelMap<dim> get_elastic_kernels(double shear_modulus, double poisson_ratio);

#endif
