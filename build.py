from tools.fabricate import *
import subprocess
import sys
import shutil
import os
import glob

def prepend_dir(files, dirname = 'src'):
    return [os.path.join(dirname, f) for f in files]

tbem_loc = '../stablelib'
shared_link_flags = [
    '-Wl,-rpath=' + tbem_loc + '/build',
    '-L' + tbem_loc + '/build/',
    '-l3bem',
    '-lhdf5',
    '-fopenmp'
]
compiler = 'mpic++'


shared_sources = [
    'load', 'spec', 'kernels', 'compute', 'data',
    'filenames', 'reload_soln', 'nearest_neighbor'
]

test_files = [
    'test_filenames',
    'test_load',
    'test_nearest_neighbor',
    'test_compute'
]
test_link_flags = [
    '-L../lib/unittest-cpp/builds',
    '-lUnitTest++'
]

executables = []

test = dict()
test['sources'] = prepend_dir(test_files + shared_sources)
test['exec_name'] = 'test'
test['link_flags'] = shared_link_flags + test_link_flags
executables.append(test)

for d in ['2','3']:
    solver = dict()
    solver['sources'] = prepend_dir(['solve' + d + 'd'] + shared_sources)
    solver['exec_name'] = 'solve' + d + 'd'
    solver['link_flags'] = shared_link_flags
    executables.append(solver)

interior = dict()
interior['sources'] = prepend_dir(['interior'] + shared_sources)
interior['exec_name'] = 'interior'
interior['link_flags'] = shared_link_flags
executables.append(interior)

def get_cpp_flags(tbem_loc):
    cpp_flags = '-Wall -std=c++11 -fopenmp \
        -DDEBUG=1 -g -Ofast -funroll-loops'.split()

    cpp_flags.extend([
        '-I' + tbem_loc,
        '-I../lib/unittest-cpp/UnitTest++',
        '-I../lib/rapidjson/include',
    ])
    return cpp_flags


def build():
    compile()
    link()

def compile():
    for e in executables:
        compile_list(e['sources'])
        after()

def compile_list(srces):
    for s in srces:
        run(compiler, '-c', s + '.cpp', '-o', s + '.o', get_cpp_flags(tbem_loc))

def obj_files(srces):
    return [s + '.o' for s in srces]

def link_exec(srces, name, flags):
    objs = obj_files(srces)
    run(compiler, '-o', name, objs, flags)

def link():
    after()
    for e in executables:
        link_exec(e['sources'], e['exec_name'], e['link_flags'])

def tests():
    assert(sys.argv[1] == 'tests')
    subprocess.call(['./test'])
    if os.path.exists('tools/__pycache__'):
        shutil.rmtree('tools/__pycache__')
    cmd = ['py.test']
    cmd.extend(glob.glob('tools/test_*.py'))
    cmd.extend(glob.glob('acctests/*.py'))
    cmd.extend(sys.argv[2:])
    subprocess.call(cmd)

def clean():
    autoclean()

def rebuild():
    clean()
    build()

if __name__ == "__main__":
    main(parallel_ok = True, jobs = 12)
