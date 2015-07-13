from elastic import line, slip, traction, execute, Result
import numpy as np
import os

import logging
logging.basicConfig(filename = 'log.txt', filemode = 'w', level = logging.DEBUG)

def test_save_load():
    es = line([[20, 0], [-20, 0]], 6, lambda pts: traction(pts, [[0, 0], [0, 0]])) +\
        line([[-1, -1], [0, 0]], 6, lambda pts: slip(pts, [[1, 1], [1, 1]]))
    result = execute(2, es, dict(
        dense = True
    ))
    filename = os.path.join('tests', 'saved_data')
    result.save(filename)

    result2 = Result.load(filename)

    ux1, uy1 = result.interior_displacement(np.array([[-1, -2]]))
    ux2, uy2 = result2.interior_displacement(np.array([[-1, -2]]))
    assert(ux1 == ux2)
    assert(uy1 == uy2)
