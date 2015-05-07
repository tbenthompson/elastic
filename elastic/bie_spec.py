from tbempy import *
import tbempy.TwoD
import tbempy.ThreeD

def get_tbem(dim):
    if dim == 2:
        tbem = tbempy.TwoD
    else:
        tbem = tbempy.ThreeD
    return tbem

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
        get_traction_BIE('traction')
    ]

def get_displacement_BIE(obs_mesh_name):
    return dict(
        obs_mesh = obs_mesh_name,
        unknown_field = 'traction',
        mass_term = dict(
            obs_mesh = obs_mesh_name,
            function = 'displacement',
            multiplier = 1
        ),
        terms = displacement_BIE_terms(obs_mesh_name),
        constraint_builder = form_traction_constraints
    )

def form_traction_constraints(tbem, component_map, dof_map, meshes, d):
    return []

def get_traction_BIE(obs_mesh_name):
    # TODO: Make a namedtuple for the BIE
    return dict(
        obs_mesh = obs_mesh_name,
        unknown_field = 'displacement',
        mass_term = dict(
            obs_mesh = obs_mesh_name,
            function = 'traction',
            multiplier = 1
        ),
        terms = displacement_BIE_terms(obs_mesh_name),
        constraint_builder = form_displacement_constraints
    )

def form_displacement_constraints(tbem, component_map, dof_map, meshes, d):
    # TODO: I'm not sure if the meshes['traction'] is correct when this
    # function is called for a the 'crack_traction' or 'free_slip_traction'
    # meshes
    continuity = tbem.mesh_continuity(meshes['traction'].begin())
    constraints = tbem.convert_to_constraints(continuity)
    start_dof = dof_map[component_map['traction'] + d]
    return tbem.shift_constraints(constraints, start_dof)

# TODO: Make a namedtuple for the terms
def displacement_BIE_terms(obs_mesh_name):
    return [
        dict(
            obs_mesh = obs_mesh_name,
            src_mesh = 'displacement',
            kernel = 'traction',
            function = 'displacement',
            multiplier = -1
        ),
        dict(
            obs_mesh = obs_mesh_name,
            src_mesh = 'displacement',
            kernel = 'displacement',
            function = 'traction',
            multiplier = -1
        ),
        dict(
            obs_mesh = obs_mesh_name,
            src_mesh = 'traction',
            kernel = 'traction',
            function = 'displacement',
            multiplier = -1
        ),
        dict(
            obs_mesh = obs_mesh_name,
            src_mesh = 'traction',
            kernel = 'displacement',
            function = 'traction',
            multiplier = -1
        ),
        # dict(
        #     obs_mesh = obs_mesh_name,
        #     src_mesh = 'slip',
        #     kernel = 'traction',
        #     function = 'slip',
        #     multiplier = -1
        # ),
        # dict(
        #     obs_mesh = obs_mesh_name,
        #     src_mesh = 'crack_traction',
        #     kernel = 'traction',
        #     function = 'slip',
        #     multiplier = -1
        # ),
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
            multiplier = 1
        ),
        dict(
            obs_mesh = obs_mesh_name,
            src_mesh = 'displacement',
            kernel = 'adjoint_traction',
            function = 'traction',
            multiplier = 1
        ),
        dict(
            obs_mesh = obs_mesh_name,
            src_mesh = 'traction',
            kernel = 'hypersingular',
            function = 'displacement',
            multiplier = 1
        ),
        dict(
            obs_mesh = obs_mesh_name,
            src_mesh = 'traction',
            kernel = 'adjoint_traction',
            function = 'traction',
            multiplier = 1
        ),
        # dict(
        #     obs_mesh = obs_mesh_name,
        #     src_mesh = 'slip',
        #     kernel = 'hypersingular',
        #     function = 'slip',
        #     multiplier = 1
        # ),
        # dict(
        #     obs_mesh = obs_mesh_name,
        #     src_mesh = 'crack_traction',
        #     kernel = 'hypersingular',
        #     function = 'slip',
        #     multiplier = 1
        # ),
        # dict(
        #     obs_mesh = obs_mesh_name,
        #     src_mesh = 'free_slip_traction',
        #     kernel = 'hypersingular',
        #     function = 'free_slip',
        #     multiplier = 1
        # )
    ]

