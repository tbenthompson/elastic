from fabricate import *

sources = ['test']

cpp_flags = '-Wall -std=c++11 -O3 -DDEBUG'.split()
cpp_flags.append('-I../3bem_stable/src')
link_flags = [
    '-L../3bem_stable/build/',
    '-l3bem'
]

def build():
    compile()
    link()

def compile():
    for source in sources:
        run('g++', '-c', source+'.cpp', cpp_flags)

def link():
    objs = [s+'.o' for s in sources]
    run('g++', '-o', 'test', objs)

def clean():
    autoclean()

main()
