from elastic import execute, Element
import numpy as np

def freeslip_run(es, constrained_dir):
    problem = execute(2, es, dict())
    s = problem.soln[('free_slip_traction', 'free_slip')]
    sx = s[0]
    sy = s[1]
    slip_has_happened = np.all(np.array(sx) > 0.0) or np.all(np.array(sy) > 0.0)
    assert(slip_has_happened)
    assert(np.all(constrained_dir(sx, sy) == 0.0))

def test_vertical_free_slip():
    freeslip_run([
        Element([[0, -1], [0, 1]], [[0, -1e10], [0, -1e10]], "free_slip_traction", 4)
        ], lambda sx, sy: sx)

def test_diag_free_slip():
    freeslip_run([
        Element([[-1, -1], [1, 1]], [[-1e10, 0], [-1e10, 0]], "free_slip_traction", 4)
        ], lambda sx, sy: sx - sy)
