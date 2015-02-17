from tools.input_builder import *

def freeslip_run(es, constrained_dir):
    fileroot = 'free_slip'
    input_filepath = test_data_dir + fileroot + '.in'
    pts_filepath = test_data_dir + fileroot + '.pts'
    slip_filepath = test_data_dir + fileroot + '.free_slip_out'


    bem_template(input_filepath, es = es)
    run(input_filepath, stdout_dest = subprocess.PIPE)
    points_grid([-2, 2, 20], [-2, 2, 20], pts_filepath)
    interior_run(input_filepath, pts_filepath)
    f = h5py.File(slip_filepath)
    sx = f['values0'][:, 0]
    sy = f['values1'][:, 0]
    slip_has_happened = np.all(np.array(sx) > 0.0) or np.all(np.array(sy) > 0.0)
    assert(slip_has_happened)
    assert(np.all(constrained_dir(sx, sy) == 0.0))

def test_vertical_free_slip():
    freeslip_run([
        Element([[0, -1], [0, 1]], "free_slip_traction", [[0, 1e10], [0, 1e10]], 5)
        ], lambda sx, sy: sx)

def test_diag_free_slip():
    freeslip_run([
        Element([[-1, -1], [1, 1]], "free_slip_traction", [[1e10, 0], [1e10, 0]], 5)
        ], lambda sx, sy: sx - sy)
