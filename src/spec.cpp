#include "spec.h"

IntegralEquationSpec get_displacement_BIE(const std::string& obs_mesh) {
    IntegralSpec uut{obs_mesh, "displacement", "traction", "displacement", -1};
    IntegralSpec uuu{obs_mesh, "displacement", "displacement", "traction", -1};
    IntegralSpec utt{obs_mesh, "traction", "traction", "displacement", -1};
    IntegralSpec utu{obs_mesh, "traction", "displacement", "traction", -1};
    IntegralSpec ust{obs_mesh, "slip", "traction", "slip", -1};

    return {
        {obs_mesh, "displacement", 1},
        {uut, utt, ust, uuu, utu}
    };
}

IntegralEquationSpec get_traction_BIE(const std::string& obs_mesh) {
    IntegralSpec tuh{obs_mesh, "displacement", "hypersingular", "displacement", -1};
    IntegralSpec tua{obs_mesh, "displacement", "adjoint_traction", "traction", -1};
    IntegralSpec tth{obs_mesh, "traction", "hypersingular", "displacement", -1};
    IntegralSpec tta{obs_mesh, "traction", "adjoint_traction", "traction", -1};
    IntegralSpec tsh{obs_mesh, "slip", "hypersingular", "slip", -1};

    return {
        {obs_mesh, "traction", -1},
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
