import tbempy.TwoD as tbem
from elastic import Element, execute
from elastic.mesh_gen import circle
import matplotlib.pyplot as plt
import numpy as np

def test_gravity():
    r = 1.0
    es = []
    es.extend(circle([0, 0], r, 6, "displacement", lambda x: np.zeros_like(x), False))
    es.extend(circle([0, 0], r, 6, "gravity", lambda x: np.zeros_like(x), False))

    params = dict(
        shear_modulus = 30e9,
        poisson_ratio = 0.25,
        gravity = True,
        gravity_vector = [0.0, -9.8 * 2700],
        dense = True
    )
    result = execute(2, es, params)

    pts = []
    x = y = np.linspace(-r, r, 100)
    for i in range(x.shape[0]):
        for j in range(y.shape[0]):
            if x[i] ** 2 + y[j] ** 2 < r ** 2:
                pts.append([x[i], y[j]])
    pts = np.array(pts)
    normalsx = np.array([[1, 0]] * pts.shape[0])
    normalsy = np.array([[0, 1]] * pts.shape[0])

    ux, uy = result.interior_displacement(pts)
    sxx, sxy = result.interior_traction(pts, normalsx)
    sxy, syy = result.interior_traction(pts, normalsy)
    assert(np.allclose(np.min(uy), -1.10e-7, atol = 1e-9))
    assert(np.allclose(np.max(uy), 0.0, atol = 1e-9))
    assert(np.allclose(np.min(syy), -19500, atol = 50.0))
    assert(np.allclose(np.max(syy), 19500, atol = 50.0))
    assert(np.allclose(np.min(sxy), -6500, atol = 50.0))
    assert(np.allclose(np.max(sxy), 6500, atol = 50.0))
    assert(np.allclose(np.min(sxx), -6500, atol = 50.0))
    assert(np.allclose(np.max(sxx), 6500, atol = 50.0))

if __name__ == '__main__':
    test_gravity()
