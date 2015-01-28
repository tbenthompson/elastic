from tools.fabricate import *
import subprocess
import sys

test_sources = ['test_compute', 'test_function', 'test_load', 'load', 'spec', 'compute']
test_sources = ['src/' + f for f in test_sources]
run_sources = ['main', 'load', 'spec', 'compute']
run_sources = ['src/' + f for f in run_sources]

run_name = 'run'

tbem_loc = '../stablelib'
cpp_flags = '-Wall -std=c++11 -Og -DDEBUG'.split()
cpp_flags.extend([
    '-I' + tbem_loc,
    '-I../lib/unittest-cpp/UnitTest++',
    '-I../lib/rapidjson/include',
    '-fopenmp'
])

link_flags = [
    '-Wl,-rpath=' + tbem_loc + '/build',
    '-L' + tbem_loc + '/build/',
    '-l3bem',
    '-fopenmp'
]

test_link_flags = link_flags + [
    '-L../lib/unittest-cpp/builds',
    '-lUnitTest++'
]

run_link_flags = link_flags

def build():
    compile()
    link()

def compile():
    for source in (test_sources + run_sources):
        run('g++', '-c', source + '.cpp', '-o', source + '.o', cpp_flags)

def link():
    after()
    test_objs = [s + '.o' for s in test_sources]
    run('g++', '-o', 'test', test_objs, test_link_flags)

    run_objs = [s + '.o' for s in run_sources]
    run('g++', '-o', run_name, run_objs, run_link_flags)

def tests():
    assert(sys.argv[1] == 'tests')
    args = ' '.join(sys.argv[2:])
    cmd = ' py.test -s \
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
