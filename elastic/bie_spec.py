from tbempy import *
import tbempy.TwoD
import tbempy.ThreeD

def get_tbem(dim):
    if dim == 2:
        tbem = tbempy.TwoD
    else:
        tbem = tbempy.ThreeD
    return tbem

def default_params():
    #TODO: Load this from a file? Or move it to spec?
    return dict(
        obs_order = 3,
        singular_steps = 8,
        far_threshold = 3.0,
        near_tol = 1e-4,
        solver_tol = 1e-6,
        poisson_ratio = 0.25,
        shear_modulus = 30e9
    )

def bc_types():
    return [
        'traction',
        'displacement',
        'slip',
        'crack_traction',
        'free_slip_traction'
    ]

def mesh_types():
    return bc_types()

def get_all_BIEs():
    return [
        get_displacement_BIE('displacement'),
        get_traction_BIE('traction'),
        get_crack_traction_BIE('crack_traction')
    ]

def get_displacement_BIE(obs_mesh_name):
    return dict(
        obs_mesh = obs_mesh_name,
        unknown_field = 'traction',
        mass_term = dict(
            obs_mesh = obs_mesh_name,
            function = 'displacement',
            multiplier = 1.0
        ),
        terms = displacement_BIE_terms(obs_mesh_name),
        constraint_builder = form_traction_constraints,
        scaling = lambda p: p['shear_modulus']
    )

def form_traction_constraints(tbem, dof_map, meshes):
    return []

def get_traction_BIE(obs_mesh_name):
    return dict(
        obs_mesh = obs_mesh_name,
        unknown_field = 'displacement',
        mass_term = dict(
            obs_mesh = obs_mesh_name,
            function = 'traction',
            multiplier = 1.0
        ),
        terms = traction_BIE_terms(obs_mesh_name),
        constraint_builder = form_displacement_constraints,
        scaling = lambda p: 1.0
    )

def form_displacement_constraints(tbem, dof_map, meshes):
    continuity = tbem.mesh_continuity(meshes['traction'].begin())
    combined_crack_mesh = tbem.Mesh.create_union([
        meshes['crack_traction'],
        meshes['slip']
    ])
    cut_continuity = tbem.cut_at_intersection(
        continuity,
        meshes['traction'].begin(),
        combined_crack_mesh.begin()
    )
    one_component = tbem.convert_to_constraints(cut_continuity)
    all_components = []
    for d in range(tbem.dim):
        start_dof = dof_map[('traction', 'displacement')][d]
        all_components.extend(tbem.shift_constraints(one_component, start_dof))
    return all_components

def get_crack_traction_BIE(obs_mesh_name):
    return dict(
        obs_mesh = obs_mesh_name,
        unknown_field = 'slip',
        mass_term = dict(
            obs_mesh = obs_mesh_name,
            function = 'crack_traction',
            multiplier = 1.0
        ),
        terms = traction_BIE_terms(obs_mesh_name),
        constraint_builder = form_slip_constraints,
        scaling = lambda p: 1.0
    )

def form_slip_constraints(tbem, dof_map, meshes):
    continuity = tbem.mesh_continuity(meshes['crack_traction'].begin())
    one_component = tbem.convert_to_constraints(continuity)
    all_components = []
    for d in range(tbem.dim):
        start_dof = dof_map[('crack_traction', 'slip')][d]
        all_components.extend(tbem.shift_constraints(one_component, start_dof))
    return all_components

def displacement_BIE_terms(obs_mesh_name):
    return [
        dict(
            obs_mesh = obs_mesh_name,
            src_mesh = 'displacement',
            kernel = 'traction',
            function = 'displacement',
            multiplier = -1.0
        ),
        dict(
            obs_mesh = obs_mesh_name,
            src_mesh = 'displacement',
            kernel = 'displacement',
            function = 'traction',
            multiplier = -1.0
        ),
        dict(
            obs_mesh = obs_mesh_name,
            src_mesh = 'traction',
            kernel = 'traction',
            function = 'displacement',
            multiplier = -1.0
        ),
        dict(
            obs_mesh = obs_mesh_name,
            src_mesh = 'traction',
            kernel = 'displacement',
            function = 'traction',
            multiplier = -1.0
        ),
        dict(
            obs_mesh = obs_mesh_name,
            src_mesh = 'slip',
            kernel = 'traction',
            function = 'slip',
            multiplier = -1
        ),
        dict(
            obs_mesh = obs_mesh_name,
            src_mesh = 'crack_traction',
            kernel = 'traction',
            function = 'slip',
            multiplier = -1
        ),
        # dict(
        #     obs_mesh = obs_mesh_name,
        #     src_mesh = 'free_slip_traction',
        #     kernel = 'traction',
        #     function = 'free_slip',
        #     multiplier = -1
        # )
    ]


def traction_BIE_terms(obs_mesh_name):
    return [
        dict(
            obs_mesh = obs_mesh_name,
            src_mesh = 'displacement',
            kernel = 'hypersingular',
            function = 'displacement',
            multiplier = 1.0
        ),
        dict(
            obs_mesh = obs_mesh_name,
            src_mesh = 'displacement',
            kernel = 'adjoint_traction',
            function = 'traction',
            multiplier = 1.0
        ),
        dict(
            obs_mesh = obs_mesh_name,
            src_mesh = 'traction',
            kernel = 'hypersingular',
            function = 'displacement',
            multiplier = 1.0
        ),
        dict(
            obs_mesh = obs_mesh_name,
            src_mesh = 'traction',
            kernel = 'adjoint_traction',
            function = 'traction',
            multiplier = 1.0
        ),
        dict(
            obs_mesh = obs_mesh_name,
            src_mesh = 'slip',
            kernel = 'hypersingular',
            function = 'slip',
            multiplier = 1
        ),
        dict(
            obs_mesh = obs_mesh_name,
            src_mesh = 'crack_traction',
            kernel = 'hypersingular',
            function = 'slip',
            multiplier = 1
        ),
        # dict(
        #     obs_mesh = obs_mesh_name,
        #     src_mesh = 'free_slip_traction',
        #     kernel = 'hypersingular',
        #     function = 'free_slip',
        #     multiplier = 1
        # )
    ]

