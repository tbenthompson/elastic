from tools.input_builder import *
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

def max_mem_requirements_estimate(refine):
    elements = 2 ** refine
    dofs = 2 * elements
    bytes_per_entry = 8
    matrix_mem = bytes_per_entry * dofs ** 2
    matrix_MB = matrix_mem / (2.0 ** 20)
    leeway = 3
    return leeway * matrix_MB

def test_memory():
    fileroot = 'memory_stress'
    input_filepath = test_data_dir + fileroot + '.in'
    refine = 8
    es = [
        Element([[0, -1], [0, 1]], "free_slip_traction", [[0, 1e10], [0, 1e10]], refine)
    ]
    bem_template(input_filepath, es = es)
    cmd = run_command(input_filepath, 2)
    P = subprocess.Popen(cmd, stdout = subprocess.PIPE)
    max_mem = 0.0
    while P.poll() is None:
        cur_mem = memory_usage_ps(P.pid)
        max_mem = max(cur_mem, max_mem)

    max_allowed = max_mem_requirements_estimate(refine)
    assert(max_mem < max_allowed)
