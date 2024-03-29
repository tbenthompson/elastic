from elastic import execute, circle, displacement
from errors import check_error
import numpy as np

def check_zero_traction(disp_bc):
    es = circle([0, 0], 1.0, 3,
        lambda pts: displacement(pts, [disp_bc(pts[0, :]), disp_bc(pts[1, :])]),
    False)
    solved = execute(2, es, dict(dense = True))
    for d in range(2):
        t_component = solved.soln[('continuous', 'traction')][d]
        mu = solved.params['shear_modulus']
        error = np.sqrt(np.sum((t_component / mu) ** 2))
        assert(error < 0.006)

def test_rigid_body_translation():
    # The results from this test diverge with increasing refinement, rapidly
    # failing. This is interesting. The divergence can be counteracted by
    # improving the integration tolerances. It will be interesting to see
    # if this divergence disappears with the upcoming improvements to
    # the robustness of the integration system. I'm suspicious that the
    # instability relates to the inconsistent placement of richardson points.
    def disp_bc(pt):
        return 1.0 * np.ones_like(pt[0]), 1.0 * np.ones_like(pt[0])
    check_zero_traction(disp_bc)

def test_rigid_body_rotation():
    def disp_bc(pt):
        return -pt[1], pt[0]
    check_zero_traction(disp_bc)
