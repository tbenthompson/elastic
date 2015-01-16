#ifndef __YYYCYYCXYZYZYYZI_SPEC_H
#define __YYYCYYCXYZYZYYZI_SPEC_H

#include <vector>
#include <string>
#include <map>
#include "3bem/kernel.h"

template <size_t dim>
using KernelSet = std::map<
    std::string,
    const tbem::Kernel<dim,
        tbem::Vec<double,dim>,
        tbem::Vec<double,dim>,
        tbem::Vec<tbem::Vec<double,dim>,dim>
    >&
>;

struct OperatorSpec 
{
    const std::string obs_mesh;
    const std::string src_mesh;
    const std::string kernel;
    const std::string function;
    const double multiplier;
};

// Operators for the displacement BIE
OperatorSpec uut("displacement", "displacement", "traction", "displacement", 1);
OperatorSpec utt("displacement", "traction", "traction", "displacement", 1);
OperatorSpec ust("displacement", "slip", "traction", "slip", 1);
OperatorSpec uuu("displacement", "displacement", "displacement", "traction", -1);
OperatorSpec utu("displacement", "traction", "displacement", "traction", -1);

// Operators for the traction BIE
OperatorSpec tuh("traction", "displacement", "hypersingular", "displacement", 1);
OperatorSpec tth("traction", "traction", "hypersingular", "displacement", 1);
OperatorSpec tsh("traction", "slip", "hypersingular", "displacement", 1);
OperatorSpec tua("traction", "displacement", "adjoint_traction", "traction", -1);
OperatorSpec tta("traction", "traction", "adjoint_traction", "traction", -1);

struct MassSpec {
    const std::string obs_mesh;
    const std::string function;
    const double multiplier;
};

struct IntegralEquation {
    MassSpec mass;
    std::vector<OperatorSpec> terms;
};

IntegralEquation displacement_BIE{
    {"displacement", "displacement", 1},
    {uut, utt, ust, uuu, utu}
};

IntegralEquation traction_BIE{
    {"traction", "traction", 1},
    {tuh, tth, tsh, tua, tta}
};

#endif
