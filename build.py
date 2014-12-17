from fabricate import *

sources = ['test', 'elastic']

cpp_flags = '-Wall -std=c++11 -O3 -DDEBUG'.split()
cpp_flags.extend([
    '-I../3bem_stable/src',
    '-I../lib/unittest-cpp/src',
    '-I../lib/rapidjson/include',
])
link_flags = [
    '-L../3bem_stable/build/',
    '-l3bem',
    '-L../lib/unittest-cpp',
    '-lUnitTest++'
]

def build():
    compile()
    link()

def compile():
    for source in sources:
        run('g++', '-c', source+'.cpp', cpp_flags)

def link():
    objs = [s+'.o' for s in sources]
    run('g++', '-o', 'test', objs, link_flags)

def clean():
    autoclean()

main()
