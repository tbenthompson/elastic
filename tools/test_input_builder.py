from input_builder import *

def ones_bc(x,y):
    return np.ones_like(x), np.ones_like(y)

def test_line():
    a = [0, 0]
    b = [1, 0]
    es = line([a, b], 2, "traction", ones_bc)
    assert(len(es) == 4)
    np.testing.assert_almost_equal(es[0].pts[0], a)
    for i in range(3):
        np.testing.assert_almost_equal(es[i].pts[1], es[i + 1].pts[0])
    np.testing.assert_almost_equal(es[3].pts[1], b)

def test_circle():
    es = circle([0, 0], 1.0, 2, "traction", ones_bc, False)
    assert(len(es) == 4)
    for i in range(3):
        np.testing.assert_almost_equal(es[i].pts[1], es[i + 1].pts[0])
    np.testing.assert_almost_equal(es[3].pts[1], es[0].pts[0])
