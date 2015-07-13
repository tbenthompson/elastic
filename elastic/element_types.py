from constraints import Constraint, BCConstraint, Term

def displacement(pts, bc):
    dim = len(pts)
    return dict(
        pts = pts,
        type = 'continuous',
        constraints = [
            BCConstraint('displacement', d1, [bc[d2][d1] for d2 in range(dim)])
            for d1 in range(dim)
        ]
    )

def traction(pts, bc):
    dim = len(pts)
    return dict(
        pts = pts,
        type = 'continuous',
        constraints = [
            BCConstraint('traction', d1, [bc[d2][d1] for d2 in range(dim)])
            for d1 in range(dim)
        ]
    )

def mixed(pts, types, bc):
    dim = len(pts)
    return dict(
        pts = pts,
        type = 'continuous',
        constraints = [
            BCConstraint(types[d1], d1, [bc[d2][d1] for d2 in range(dim)])
            for d1 in range(dim)
        ]
    )

def free_slip(pts, normal_displacement, shear_stress):
    cs = [
        BCConstraint('crack_traction', 'tangential0', shear_stress[0]),
        BCConstraint('slip', 'normal', normal_displacement)
    ]
    dim = len(pts)
    if dim == 3:
        cs.append(BCConstraint('crack_traction', 'tangential1', shear_stress[1]))
    return dict(pts = pts, type = 'discontinuous', constraints = cs)

def slip(pts, bc):
    dim = len(pts)
    return dict(
        pts = pts,
        type = 'discontinuous',
        constraints = [
            BCConstraint('slip', d1, [bc[d2][d1] for d2 in range(dim)])
            for d1 in range(dim)
        ]
    )

def crack(pts, bc):
    dim = len(pts)
    return dict(
        pts = pts,
        type = 'discontinuous',
        constraints = [
            BCConstraint('crack_traction', d1, [bc[d2][d1] for d2 in range(dim)])
            for d1 in range(dim)
        ]
    )

def static_friction(pts, normal_displacement, friction_coefficient):
    dim = len(pts)
    if dim == 3:
        return 'unimplemented!'
    cs = [
        Constraint([
            Term('crack_traction', 'normal', friction_coefficient),
            Term('crack_traction', 'tangential0', -1.0),
        ], [0, 0]),
        BCConstraint('slip', 'normal', normal_displacement)
    ]
    return dict(pts = pts, type = 'discontinuous', constraints = cs)
