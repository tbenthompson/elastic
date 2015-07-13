from elastic import execute, line, free_slip
import numpy as np

def freeslip_run(es, constrained_dir):
    problem = execute(2, es, dict())
    s = problem.soln[('discontinuous', 'slip')]
    sx = s[0]
    sy = s[1]
    slip_has_happened = np.all(np.array(sx) > 0.0) or np.all(np.array(sy) > 0.0)
    assert(slip_has_happened)
    assert(np.all(constrained_dir(sx, sy) == 0.0))

def test_vertical_free_slip():
    freeslip_run(
        line([[0, -1], [0, 1]], 4,
            lambda pts: free_slip(pts, [0, 0], [[-1e10, -1e10]])
        ),
        lambda sx, sy: sx
    )

def test_diag_free_slip():
    freeslip_run(
        line([[-1, -1], [1, 1]], 4,
            lambda pts: free_slip(pts, [0, 0], [[-1e10, -1e10]])
        ),
        lambda sx, sy: sx - sy
    )
