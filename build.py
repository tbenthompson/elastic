from fabricate import *

test_sources = ['test', 'elastic']
run_sources = ['elastic2d', 'elastic']

tbem_loc = '../3bem_stable'

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
        run('g++', '-c', source+'.cpp', cpp_flags)

def link():
    after()
    test_objs = [s+'.o' for s in test_sources]
    run('g++', '-o', 'test', test_objs, test_link_flags)

    run_objs = [s+'.o' for s in run_sources]
    run('g++', '-o', 'elastic2d', run_objs, run_link_flags)

def clean():
    autoclean()

main(parallel_ok = True, jobs = 12)
