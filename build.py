from fabricate import *

test_sources = ['test', 'elastic']
run_sources = ['elastic2d', 'elastic']

cpp_flags = '-Wall -std=c++11 -O3 -DDEBUG'.split()
cpp_flags.extend([
    '-I../3bem_stable/src',
    '-I../lib/unittest-cpp/src',
    '-I../lib/rapidjson/include',
    '-fopenmp'
])

link_flags = [
    '-Wl,-rpath=../3bem_stable/build',
    '-L../3bem_stable/build/',
    '-l3bem',
    '-fopenmp'
]

test_link_flags = link_flags + [
    '-L../lib/unittest-cpp',
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
    test_objs = [s+'.o' for s in test_sources]
    run('g++', '-o', 'test', test_objs, test_link_flags)

    run_objs = [s+'.o' for s in run_sources]
    run('g++', '-o', 'elastic2d', run_objs, run_link_flags)

def clean():
    autoclean()

main()
