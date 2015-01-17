#include "spec.h"
#include "load.h"
#include "compute.h"

using namespace tbem;

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage is 'elastic_process filename'" << std::endl;
        return 1;
    }
    
    auto filename = argv[1];
    auto bem_input = parse_into_bem<2>(filename);

    auto disp_BIE_ops = bem_input.compute_integral_equation(get_displacement_BIE());
    auto trac_BIE_ops = bem_input.compute_integral_equation(get_traction_BIE());
    std::cout << "Initial ops" << std::endl;
    std::cout << disp_BIE_ops.terms.size() << std::endl;
    std::cout << trac_BIE_ops.terms.size() << std::endl;

    auto disp_system = separate(disp_BIE_ops, bem_input.bcs);
    auto trac_system = separate(trac_BIE_ops, bem_input.bcs);
    std::cout << "after rhs eval ops" << std::endl;
    std::cout << disp_system.lhs.size() << std::endl;
    std::cout << trac_system.lhs.size() << std::endl;

    // //prep:
    // scale rows
    auto disp_scaled_system = scale_rows(disp_system);
    auto trac_scaled_system = scale_rows(trac_system);

    //reduce using constraints
    //stack rows
    const auto& disp_constraints = bem_input.displacement_constraints;
    auto stacked_rhs = concatenate({
        disp_constraints.get_reduced(disp_scaled_system.rhs[0]),
        disp_constraints.get_reduced(disp_scaled_system.rhs[1]),
        trac_scaled_system.rhs[0],
        trac_scaled_system.rhs[1]
    });


    //solve:
    int count = 0;
    auto reduced_soln = solve_system(stacked_rhs.data, 1e-5,
        [&] (std::vector<double>& x, std::vector<double>& y) {
            std::cout << "iteration " << count << std::endl;
            count++;

            //unstack_rows
            auto unstacked_rhs = expand(stacked_rhs, x);

            //expand using constraints
            int n_disp_dofs = disp_scaled_system.rhs[0].size();
            Function unknown_disp = {
                disp_constraints.get_all(unstacked_rhs[0], n_disp_dofs),
                disp_constraints.get_all(unstacked_rhs[1], n_disp_dofs)
            };
            Function unknown_trac = {
                unstacked_rhs[2],
                unstacked_rhs[3],
            };

            //apply operators
            //TODO: better name than BCMap, FunctionMap?
            BCMap unknowns; 
            unknowns[FieldDescriptor{"traction", "displacement"}] = unknown_disp;
            unknowns[FieldDescriptor{"displacement", "traction"}] = unknown_trac;

            std::cout << "Inside solution." << std::endl;
            auto disp_eval = separate(ComputedIntegralEquation{disp_scaled_system.lhs}, unknowns);
            std::cout << disp_eval.lhs.size() << std::endl;
            // assert(disp_eval.lhs.size() == 0);
            auto trac_eval = separate(ComputedIntegralEquation{trac_scaled_system.lhs}, unknowns);
            std::cout << trac_eval.lhs.size() << std::endl;
            assert(trac_eval.lhs.size() == 0);
            //reduce using constraints
            //stack rows
        }
    );
    
    //handle soln:
    //unstack rows
    //expand using constraints

    //output

    //interior eval
}
