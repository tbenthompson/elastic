from elastic import Element, execute, Result
import numpy as np
import os

def test_save_load():
    es = [
        Element([[20, 0], [-20, 0]], [[0, 0], [0, 0]], "traction", 6),
        Element([[-1, -1], [1, 1]], [[1, 1], [1, 1]], "slip", 6)
    ]
    result = execute(2, es, dict(dense = True))
    filename = os.path.join('tests', 'saved_data')
    result.save(filename)

    result2 = Result.load(filename)

    ux1, uy1 = result.interior_displacement(np.array([[-1, -2]]))
    ux2, uy2 = result2.interior_displacement(np.array([[-1, -2]]))
    assert(ux1 == ux2)
    assert(uy1 == uy2)
