from elastic import *
import numpy as np
from stress_fnc import traction_bc_builder, form_stress_evaluators,\
        calc_traction_builder

import logging
logging.basicConfig(filename = 'log.txt', filemode = 'w', level = logging.DEBUG)

alpha = np.radians(3.5) # surface slope
beta = np.radians(10.0) # angle of basal thrust
theta = alpha + beta # taper

x0 = 100 * 1000.0 # width of wedge (m)

def basal_depth(x):
    return x * np.tan(theta)
d0 = basal_depth(x0) # depth at backstop (m)

rho = 2500 #density
g = 9.8 #strength of gravity
gravity_vector = [
    -rho * g * np.sin(alpha),
    rho * g * np.cos(alpha)
]

mu_b = 0.7 # basal friction

k1 = 1.0
k2 = 1.0
k3 = 1.0
k4 = 1.0

def stress_fnc(x, y):
    return rho * g / 2.0 * (k1 * x ** 2 * y + k2 * x * y ** 2) + \
           rho * g / 6.0 * (k3 * x ** 3 + k4 * y ** 3)

def main():
    calc_stress = form_stress_evaluators(stress_fnc, [0, 0])
    trac_calc = calc_traction_builder(calc_stress)
    trac_bc = traction_bc_builder(trac_calc)

    refine = 6
    pts = [[0, 0], [100, -100], [100, 0], [0, 0]]
    es = []
    for i in range(len(pts) - 1):
        es.extend(line([pts[i], pts[i + 1]], refine, trac_bc))
    params = dict(
        dense = True,
        shear_modulus = 1.0,
        poisson_ratio = 0.3,
        singular_steps = 8,
        obs_order = 5,
        sinh_order = 10
    )
    result = execute(2, es, params)
    result.params['obs_order'] = 5

    check_error(result, 'continuous', 'traction', calc_traction, 1e-2)

    n_interior = 80
    x_hat, y_hat = np.meshgrid(
        np.linspace(0, 1, n_interior), np.linspace(0, 1, n_interior)
    )

    x_scale = pts[1][0]
    y_scale = pts[1][1]
    x = x_hat * x_scale
    y = y_hat * (x / x_scale) * y_scale
    pts = np.array([x.flatten(), y.flatten()]).T
    normals_x = np.array([np.ones(pts.shape[0]), np.zeros(pts.shape[0])]).T
    normals_y = np.array([np.zeros(pts.shape[0]), np.ones(pts.shape[0])]).T
    tx_interior = result.interior_traction(pts, normals_x)
    ty_interior = result.interior_traction(pts, normals_y)
    stress_interior = [tx_interior[0], tx_interior[1], ty_interior[0], ty_interior[1]]
    check_interior_error(pts, normals_x, stress_interior, calc_stress, 1e-3)

if __name__ == '__main__':
    main()






























#
# def sxx(x, y, k3, k4, k8):
#     return k3 * x + k4 * y + k8 + rho_e * g * x * np.sin(alpha)
#
# def sxy(x, y, k3, k4, k8):
#     return -k3 * y
#
# def syy(x, y, k3, k4, k8):
#     return -rho_e * g * y * np.cos(alpha)
#
# def internal_traction(x, y, nx, ny, k3, k4, k8):
#     tx = sxx(x, y, k3, k4, k8) * nx + sxy(x, y, k3, k4, k8) * ny
#     ty = sxy(x, y, k3, k4, k8) * nx + syy(x, y, k3, k4, k8) * ny
#     return [tx, ty]
#
# def normal_traction(x, y, nx, ny, k3, k4, k8):
#     tx, ty = internal_traction(x, y, nx, ny, k3, k4, k8)
#     return tx * nx + ty * ny
#
# def shear_traction(x, y, nx, ny, k3, k4, k8):
#     tx, ty = internal_traction(x, y, nx, ny, k3, k4, k8)
#     return -tx * ny + ty * nx
#
# def calc_k3(k4, k8):
#     check_pt = x0
#     k3 = sympy.symbols('k3')
#     tau_b = shear_traction(
#         check_pt, basal_depth(check_pt), basal_nx, basal_ny, k3, k4, k8
#     )
#     sigma_b = normal_traction(
#         check_pt, basal_depth(check_pt), basal_nx, basal_ny, k3, k4, k8
#     )
#     eqtn = tau_b + mu_b * (1 - lambda_b) * sigma_b
#     soln = sympy.solve(eqtn, k3)
#     return float(soln[0])
#
#
# k8 = -100e6
# k4 = -1.5 * rho_e * g
# k3 = calc_k3(k4, k8)
#
#
# refine = 5
# boundary_pts = [
#     [0, 0],
#     [x0, 0],
#     [x0, d0]
# ]
#
# def trac_bc_fnc(pts):
#     nx = -(pts[1,1] - pts[0,1])
#     ny = pts[1,0] - pts[0,0]
#     #TODO: NEED TO NORMALIZE HERE OR ALL WILL FAIL IN EPIC FASHION!
#     return traction(pts,
#         [
#             internal_traction(pts[0, 0], pts[0, 1], nx, ny, k3, k4, k8),
#             internal_traction(pts[1, 0], pts[1, 1], nx, ny, k3, k4, k8)
#         ]
#     )
#
# es = line([boundary_pts[2], boundary_pts[0]], refine, trac_bc_fnc)
# es += line([boundary_pts[0], boundary_pts[1]], refine, trac_bc_fnc)
# es += line([boundary_pts[1], boundary_pts[2]], refine, trac_bc_fnc)
#
# params = dict(
#     dense = True,
#     obs_order = 6,
#     sinh_order = 22,
#     singular_steps = 8,
#     gravity = True,
#     gravity_vector = gravity_vector,
# )
# result = execute(2, es, params)
#
# def plot_field(mesh_name, field_name, component, x_axis):
#     f = result.meshes[mesh_name].facets
#     xs = f[:,:,x_axis].reshape(f.shape[0] * f.shape[1])
#     s = result.soln[(mesh_name, field_name)][component]
#     plt.plot(xs, s, '-')
#     plt.show()
#
# from shapely.geometry import Point
# from shapely.geometry.polygon import Polygon
# import time
# s = time.time()
# nx = 100
# ny = 100
# x_vec = np.linspace(0, x0, nx)
# y_vec = np.linspace(0, d0, ny)
# x_mat, y_mat = np.meshgrid(x_vec, y_vec)
# interior_pts = np.array([x_mat.flatten(), y_mat.flatten()]).T
# normals_x = np.array([np.ones(interior_pts.shape[0]), np.zeros(interior_pts.shape[0])]).T
# normals_y = np.array([np.zeros(interior_pts.shape[0]), np.ones(interior_pts.shape[0])]).T
# # ux, uy = result.interior_displacement(interior_pts)
# sxx_est, sxy_est = result.interior_traction(interior_pts, normals_x)
# sxy_est, syy_est = result.interior_traction(interior_pts, normals_y)
# print('interior took ' + str(time.time() - s) + ' secs')
#
# def interior_plot(plt_me, title):
#     polygonBuffer = Polygon(boundary_pts)
#     def not_in_polygon(pts):
#         delete_these = []
#         for i in range(0, interior_pts.shape[0]):
#             candidatePoint = Point(interior_pts[i, 0], interior_pts[i, 1])
#             isIn = polygonBuffer.contains(candidatePoint)
#             if isIn == False:
#                 delete_these.append(i)
#         return delete_these
#
#     # ignore points not in polygon of interest
#     delete_indices = not_in_polygon(interior_pts)
#     plt_me[delete_indices] = np.nan
#     plt_me_mat = plt_me.reshape(x_mat.shape)
#
#
#     fig = plt.figure(facecolor = 'white', dpi = 100)
#     plt.contourf(x_mat, y_mat, plt_me_mat)
#     cbar = plt.colorbar()
#     plt.xlim([0, x0])
#     plt.ylim([d0, 0])
#     plt.title(title)
#
# sxx_exact = sxx(interior_pts[:, 0], interior_pts[:, 1], k3, k4, k8)
# sxy_exact = sxy(interior_pts[:, 0], interior_pts[:, 1], k3, k4, k8)
# syy_exact = syy(interior_pts[:, 0], interior_pts[:, 1], k3, k4, k8)
#
# interior_plot(sxx_exact, 'sxx exact')
# interior_plot(sxy_exact, 'sxy exact')
# interior_plot(syy_exact, 'syy exact')
# interior_plot(sxx_est, 'sxx est')
# interior_plot(sxy_est, 'sxy est')
# interior_plot(syy_est, 'syy est')
# plt.show()
# interior_plot(np.log10(np.abs(sxy)))
#
# # interior_plot(np.abs(sxy) > 50e6)
# # interior_plot(ux)
# # interior_plot(uy)
# # interior_plot(np.log10(np.abs(sxx)))
# # interior_plot(np.log10(np.abs(syy)))
# import ipdb; ipdb.set_trace()
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# #
# # a11 = x0 * (np.sin(2 * theta) - np.tan(theta) * np.sin(theta) ** 2 +
# #     mu_b * (1 - lambda_b) * np.sin(theta) ** 2)
# # a12 = x0 * np.sin(theta) ** 2 * (1 + np.tan(theta) * mu_b * (1 - lambda_b))
# # b1 = rho_e * g * x0 * np.sin(theta) * (np.sin(alpha + theta) -
# #     mu_b * (1 - lambda_s) * (np.cos(alpha + theta))) + \
# #     k8 * np.sin(theta) * (np.cos(theta) + mu_b * (1 - lambda_s) * np.sin(theta))
# #
# # # k3 = (-k4 * a12 + b1) / a11
