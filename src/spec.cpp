#include "spec.h"
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
    return std::move(out);
}

template 
KernelMap<2> get_elastic_kernels(double shear_modulus, double poisson_ratio);
template 
KernelMap<3> get_elastic_kernels(double shear_modulus, double poisson_ratio);

IntegralEquationSpec get_displacement_BIE() {
    //CORRECT -1
    IntegralSpec uut{"displacement", "displacement", "traction", "displacement", -1};
    //CORRECT -1
    IntegralSpec uuu{"displacement", "displacement", "displacement", "traction", -1};
    //CORRECT -1
    IntegralSpec utt{"displacement", "traction", "traction", "displacement", -1};
    //UNKNOWN
    IntegralSpec utu{"displacement", "traction", "displacement", "traction", 1};

    //UNKNOWN
    IntegralSpec ust{"displacement", "slip", "traction", "slip", -1};
    return {
    //CORRECT 1
        {"displacement", "displacement", 1},
        {uut, utt, ust, uuu, utu}
    };
}

IntegralEquationSpec get_traction_BIE() {
    //CORRECT -1
    IntegralSpec tuh{"traction", "displacement", "hypersingular", "displacement", -1};
    //UNKNOWN
    IntegralSpec tua{"traction", "displacement", "adjoint_traction", "traction", -1};
    //CORRECT -1
    IntegralSpec tth{"traction", "traction", "hypersingular", "displacement", -1};
    //CORRECT -1
    IntegralSpec tta{"traction", "traction", "adjoint_traction", "traction", -1};
    
    //UNKNOWN
    IntegralSpec tsh{"traction", "slip", "hypersingular", "slip", 1};
    return {
    //CORRECT
        {"traction", "traction", 1},
        {tuh, tth, tsh, tua, tta}
    };
}

std::vector<std::string> get_mesh_types() {
    return {
        "traction", "displacement", "slip"
    };
}


std::vector<std::string> get_bc_types() {
    return get_mesh_types();
}
