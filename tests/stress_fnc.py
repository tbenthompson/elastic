import numpy as np
import sympy
from elastic import traction

def form_stress_lambdas(stress_fnc, body_force):
    x, y = sympy.symbols('x, y')
    f = stress_fnc(x, y)
    fx = sympy.diff(f, x)
    fy = sympy.diff(f, y)
    fxx = sympy.diff(fx, x)
    fxy = sympy.diff(fx, y)
    fyy = sympy.diff(fy, y)
    sxx = fyy - body_force[0] * x
    syy = fxx - body_force[1] * y
    sxy = -fxy
    fsxx = sympy.lambdify((x, y), sxx)
    fsxy = sympy.lambdify((x, y), sxy)
    fsyy = sympy.lambdify((x, y), syy)
    return [fsxx, fsxy, fsyy]

def calc_stress_builder(stress_lambdas):
    def calc_stress(pt, empty):
        x = pt[0]
        y = pt[1]
        sxx = stress_lambdas[0](x, y)
        sxy = stress_lambdas[1](x, y)
        syy = stress_lambdas[2](x, y)
        return [sxx, sxy, sxy, syy]
    return calc_stress

def calc_traction_builder(calc_stress):
    def calc_traction(pt, normal):
        S = calc_stress(pt, [])
        tx = S[0] * normal[0] + S[1] * normal[1]
        ty = S[2] * normal[0] + S[3] * normal[1]
        return tx, ty
    return calc_traction

def traction_bc_builder(traction_fnc):
    def trac_bc(pts):
        nx = -(pts[1, 1] - pts[0, 1])
        ny = pts[1, 0] - pts[0, 0]
        nmag = np.sqrt(nx ** 2 + ny ** 2)
        nx /= nmag
        ny /= nmag
        normal = [nx, ny]
        return traction(pts, [
            traction_fnc(pts[0, :], normal),
            traction_fnc(pts[1, :], normal)
        ])
    return trac_bc
