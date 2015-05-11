from elastic.input_builder import Element
from elastic.solver import Controller
import multiprocessing
import subprocess
import time

def memory_usage_ps(pid):
    cmd = ['ps', 'v', '-p', str(pid)]
    P = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    out = P.communicate()[0].split(b'\n')
    vsz_index = out[0].split().index(b'RSS')
    mem_kb = float(out[1].split()[vsz_index])
    mem_mb = mem_kb / 1024.0
    return mem_mb

def max_mem_requirements_estimate(refine, dim):
    elements = 2 ** refine
    dofs = 2 * elements
    bytes_per_entry = 8
    n_blocks = dim ** 2
    # The mass matrix and the primary matrix and space for copying one matrix
    n_matrices = 3
    leeway = 1.1
    matrix_mem = leeway * n_matrices * n_blocks * bytes_per_entry * (dofs ** 2)
    matrix_MB = matrix_mem / (2.0 ** 20)
    return matrix_MB

def runner(refine):
    dim = 2
    fileroot = 'memory_stress'
    es = [
        Element([[0, -1], [0, 1]], [[0, 1e10], [0, 1e10]], "displacement", refine)
    ]
    problem = Controller(dim, es, dict(solver_tol = 1e-2))

def measure_memory(refine):
    p = multiprocessing.Process(target = runner, args = [refine])
    p.start()
    most_mem_used = 0
    while p.is_alive():
        cur_mem = memory_usage_ps(p.pid)
        most_mem_used = max(cur_mem, most_mem_used)
    p.join()
    return most_mem_used

def test_memory():
    level1 = 7
    level2 = 8
    mem1 = measure_memory(level1)
    mem2 = measure_memory(level2)
    n_vertices_diff = 2 ** (level2 + 1) - 2 ** (level1 + 1)
    mem_diff = mem2 - mem1
    mem_per_vert = mem_diff * 1e6 / n_vertices_diff
    max_vertices_meade02 = 800e9 / mem_per_vert
    assert(max_vertices_meade02 > 10e6)

if __name__ == '__main__':
    test_memory()
