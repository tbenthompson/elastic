from elastic.mesh_gen import circle
from elastic.solver import execute
from errors import check_error
import numpy as np

def test_rigid_body():
    # The results from this test diverge with increasing refinement, rapidly
    # failing. This is interesting. The divergence can be counteracted by
    # improving the integration tolerances. It will be interesting to see
    # if this divergence disappears with the upcoming improvements to
    # the robustness of the integration system. I'm suspicious that the
    # instability relates to the inconsistent placement of richardson points.
    def disp_bc(pt):
        return 1.0 * np.ones_like(pt[0]), 1.0 * np.ones_like(pt[0])
    es = circle([0, 0], 1.0, 4, 'displacement', disp_bc, False)
    solved = execute(2, es, dict(dense = True))
    for d in range(2):
        t_component = solved.soln[('displacement', 'traction')][d]
        mu = solved.input.params['shear_modulus']
        error = np.sqrt(np.sum((t_component / mu) ** 2))
        assert(error < 0.01)