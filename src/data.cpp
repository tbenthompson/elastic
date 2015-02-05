#include "data.h"
#include "spec.h"

template <size_t dim>
BEM<dim> parse_into_bem(const std::string& filename)
{
    auto doc = parse_json(load_file(filename));
    auto params = get_parameters(doc);
    auto elements = get_elements<dim>(doc);
    auto meshes = get_meshes(elements);
    auto bcs = get_bcs(elements);

    return {
        params,
        meshes, 
        bcs,
        get_elastic_kernels<dim>(params.shear_modulus, params.poisson_ratio),
        tbem::QuadStrategy<dim>(params.obs_quad_order, params.src_far_quad_order,
            params.n_singular_steps, params.far_threshold, params.near_tol),
        {
            get_displacement_BIE<dim>("displacement"),
            get_traction_BIE<dim>("traction")
        }
    };
}

template
BEM<2> parse_into_bem(const std::string& filename);
template
BEM<3> parse_into_bem(const std::string& filename);
