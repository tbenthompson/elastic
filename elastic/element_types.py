from constraints import Constraint, BCConstraint

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

# def simple_friction(pts, normal_displacement, friction_coefficient):
#     dim = len(pts)
#     cs = [
#         BCConstraint('traction', 'tangential0', shear_stress[0]),
#         BCConstraint('displacement', 'normal', normal_displacement)
#     ]
#     dim = len(pts)
#     if dim == 3:
#         cs.append(BCConstraint('traction', 'tangential1', shear_stress[1]))
#     return dict(
#         pts = pts,
#         type = 'discontinuous',
#         constraints = [
#             BCConstraint('displacement', 'normal'
#         ]
#     )
