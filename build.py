from tools.fabricate import *
import subprocess
import sys
import os

def prepend_dir(files, dirname = 'src'):
    return [os.path.join(dirname, f) for f in files]

shared_sources = [
    'load', 'spec', 'kernels', 'compute', 'data',
    'filenames', 'reload_soln', 'nearest_neighbor'
]
test_sources = prepend_dir(['test_filenames', 'test_load', 'test_nearest_neighbor'] + shared_sources)
run_sources = prepend_dir(['main'] + shared_sources)
interior_sources = prepend_dir(['interior'] + shared_sources)

run_exec_name = 'solve'
interior_exec_name = 'interior'

tbem_loc = '../stablelib'

def get_cpp_flags(tbem_loc):
    cpp_flags = '-Wall -std=c++11 -O3 -DDEBUG'.split()
    cpp_flags.extend([
        '-I' + tbem_loc,
        '-I../lib/unittest-cpp/UnitTest++',
        '-I../lib/rapidjson/include',
        '-fopenmp'
    ])
    return cpp_flags

cpp_flags = get_cpp_flags(tbem_loc)

link_flags = [
    '-Wl,-rpath=' + tbem_loc + '/build',
    '-L' + tbem_loc + '/build/',
    '-l3bem',
    '-lhdf5',
    '-fopenmp'
]

test_link_flags = link_flags + [
    '-L../lib/unittest-cpp/builds',
    '-lUnitTest++'
]

run_link_flags = link_flags

compiler = 'mpic++'

def build():
    compile()
    link()

def compile():
    compile_list(test_sources)
    after()
    compile_list(run_sources)
    after()
    compile_list(interior_sources)

def compile_list(srces):
    for s in srces:
        run(compiler, '-c', s + '.cpp', '-o', s + '.o', cpp_flags)

def obj_files(srces):
    return [s + '.o' for s in srces]

def link_exec(srces, name, flags):
    objs = obj_files(srces)
    run(compiler, '-o', name, objs, flags)

def link():
    after()
    link_exec(test_sources, 'test', test_link_flags)
    link_exec(run_sources, run_exec_name, run_link_flags)
    link_exec(interior_sources, interior_exec_name, run_link_flags)

def tests():
    assert(sys.argv[1] == 'tests')
    subprocess.call('./test', shell = True)
    args = ' '.join(sys.argv[2:])
    cmd = ' py.test -s --tb=short\
        tools/test_input_builder.py\
        tools/pressured_cylinder.py\
        tools/beam_bend.py\
        ' + args
    print cmd
    subprocess.call(cmd, shell = True)

def clean():
    autoclean()

if __name__ == "__main__":
    main(parallel_ok = True, jobs = 12)
