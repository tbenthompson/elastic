#include "data.h"
#include "spec.h"
#include "filenames.h"
#include "nearest_neighbor.h"
#include "reload_soln.h"
#include "3bem/3bem.h"

    // find the nearest boundary point and make that the path to follow away
    // for nearfield integration --> nearest neighbor problems are a real theme!

using namespace tbem;

template <size_t dim>
std::vector<std::vector<double>>
compute_interior(const std::vector<tbem::ObsPt<dim>>& pts,
                 const BEM<dim>& bem,
                 const IntegralEquationSpec<dim>& eqtn_spec,
                 const BCMap& bcs) {
    std::vector<std::vector<double>> results(dim, std::vector<double>(pts.size(), 0.0));
#pragma omp parallel for
    for (size_t i = 0; i < pts.size(); i++) {
        for (const auto& term: eqtn_spec.terms) {
            const auto& src_mesh = bem.meshes.at(term.src_mesh);
            const auto& kernel = bem.kernels.at(term.kernel);
            const auto& field = bcs.at(FieldDescriptor{term.src_mesh, term.function});

            auto problem = make_problem(src_mesh, tbem::Mesh<dim>{{}}, *kernel);
            auto op = mesh_to_point_operator(problem, bem.quad_strategy, pts[i]);
            auto evaluated = apply_operator(op, field);
            for (size_t d = 0; d < dim; d++) {
                results[d][i] += evaluated[d][0];
            }
        }
    }
    return results;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cout << "Usage is './interior input.in pts.in" 
                  << std::endl;
        return 1;
    }

    auto input_filename = argv[1];
    auto bem_input = parse_into_bem<2>(input_filename);

    BlockFunction soln_disp(2);
    if (bem_input.meshes.at("traction").n_facets() > 0) {
        soln_disp = load_surface(disp_out_filename(input_filename));
    }
    BlockFunction soln_trac(2);
    if (bem_input.meshes.at("displacement").n_facets() > 0) {
        soln_trac = load_surface(trac_out_filename(input_filename));
    }
    BlockFunction soln_slip(2);
    if (bem_input.meshes.at("crack_traction").n_facets() > 0) {
        soln_slip = load_surface(slip_out_filename(input_filename));
    }

    BCMap fields = bem_input.bcs; 
    fields[FieldDescriptor{"displacement", "traction"}] = soln_trac;
    fields[FieldDescriptor{"traction", "displacement"}] = soln_disp;
    fields[FieldDescriptor{"crack_traction", "slip"}] = soln_slip;

    auto pts_filename = argv[2];
    auto pts_parsed = parse_json(load_file(pts_filename));
    auto pts = get_pts<2>(pts_parsed);

    auto whole_mesh = Mesh<2>::create_union({
        bem_input.meshes.at("traction"), 
        bem_input.meshes.at("displacement"), 
        bem_input.meshes.at("slip"),
        bem_input.meshes.at("crack_traction")
    });

    std::vector<ObsPt<2>> obs_ptsx;
    std::vector<ObsPt<2>> obs_ptsy;
    for (size_t i = 0; i < pts.size(); i++) {
        //find the nearest neighbor edge, len_scale = distance to the edge
        //while direction = direction away from the edge
        auto mesh_pt = nearest_pt(pts[i], whole_mesh);
        auto dir = decide_richardson_dir(pts[i], mesh_pt);
        const double length_factor = 5.0;
        auto length_scale = hypot(dir) / length_factor;
        obs_ptsx.push_back({length_scale, pts[i], {1, 0}, normalized(dir)});
        obs_ptsy.push_back({length_scale, pts[i], {0, 1}, normalized(dir)});
    }

    auto interior_disp = compute_interior(
        obs_ptsx, bem_input, get_displacement_BIE<2>("displacement"), fields
    );

    auto interior_sxx_sxy = compute_interior(
        obs_ptsx, bem_input, get_traction_BIE<2>("traction"), fields
    );

    auto interior_sxy_syy = compute_interior(
        obs_ptsy, bem_input, get_traction_BIE<2>("traction"), fields
    );

    BlockFunction interior_stress{
        interior_sxx_sxy[0], interior_sxx_sxy[1], interior_sxy_syy[1]
    };

    HDFOutputter interior_disp_file(interior_disp_out_filename(input_filename));
    out_volume(interior_disp_file, pts, interior_disp);

    HDFOutputter interior_stress_file(interior_stress_out_filename(input_filename));
    out_volume(interior_stress_file, pts, interior_stress);
}
