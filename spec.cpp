#include "spec.h"
#include "3bem/elastic_kernels.h"

template <size_t dim>
KernelSet<dim> get_elastic_kernels(double shear_modulus, double poisson_ratio) {
    return {
        {
            "displacement",
            tbem::ElasticDisplacement<dim>(shear_modulus, poisson_ratio)
        },
        {
            "traction",
            tbem::ElasticTraction<dim>(shear_modulus, poisson_ratio)
        },
        {
            "adjoint_traction",
            tbem::ElasticAdjointTraction<dim>(shear_modulus, poisson_ratio)
        },
        {
            "hypersingular",
            tbem::ElasticHypersingular<dim>(shear_modulus, poisson_ratio)
        }
    };
}

template 
KernelSet<2> get_elastic_kernels(double shear_modulus, double poisson_ratio);
template 
KernelSet<3> get_elastic_kernels(double shear_modulus, double poisson_ratio);

IntegralEquationSpec get_displacement_BIE() {
    OperatorSpec uut{"displacement", "displacement", "traction", "displacement", 1};
    OperatorSpec utt{"displacement", "traction", "traction", "displacement", 1};
    OperatorSpec ust{"displacement", "slip", "traction", "slip", 1};
    OperatorSpec uuu{"displacement", "displacement", "displacement", "traction", -1};
    OperatorSpec utu{"displacement", "traction", "displacement", "traction", -1};
    return {
        {"displacement", "displacement", 1},
        {uut, utt, ust, uuu, utu}
    };
}

IntegralEquationSpec get_traction_BIE() {
    OperatorSpec tuh{"traction", "displacement", "hypersingular", "displacement", 1};
    OperatorSpec tth{"traction", "traction", "hypersingular", "displacement", 1};
    OperatorSpec tsh{"traction", "slip", "hypersingular", "displacement", 1};
    OperatorSpec tua{"traction", "displacement", "adjoint_traction", "traction", -1};
    OperatorSpec tta{"traction", "traction", "adjoint_traction", "traction", -1};
    return {
        {"traction", "traction", 1},
        {tuh, tth, tsh, tua, tta}
    };
}
