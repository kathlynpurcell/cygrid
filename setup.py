#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup
from setuptools.extension import Extension
from Cython.Distutils import build_ext
import numpy
import platform
import os


EX_COMP_ARGS = []
if 'mac' in platform.system().lower():
    EX_COMP_ARGS += ['-stdlib=libc++', '-mmacosx-version-min=10.7']
    os.environ["CC"] = "g++-6"
    # os.environ["CPP"] = "cpp-6"
    # os.environ["CXX"] = "g++-6"
    # os.environ["LD"] = "gcc-6"

GRID_EXT = Extension(
    'cygrid.cygrid',
    ['cygrid/cygrid.pyx'],
    extra_compile_args=['-fopenmp', '-O3', '-std=c++11'] + EX_COMP_ARGS,
    extra_link_args=['-fopenmp'],
    language='c++',
    include_dirs=[
        numpy.get_include(),
    ]
)

HELPER_EXT = Extension(
    'cygrid.helpers',
    ['cygrid/helpers.pyx'],
    extra_compile_args=['-fopenmp', '-O3', '-std=c++11'] + EX_COMP_ARGS,
    extra_link_args=['-fopenmp'],
    language='c++',
    include_dirs=[
        numpy.get_include(),
    ]
)

HPX_EXT = Extension(
    'cygrid.healpix',
    ['cygrid/healpix.pyx'],
    extra_compile_args=['-fopenmp', '-O3', '-std=c++11'] + EX_COMP_ARGS,
    extra_link_args=['-fopenmp'],
    language='c++',
    include_dirs=[
        numpy.get_include(),
    ]
)

HPHASH_EXT = Extension(
    'cygrid.hphashtab',
    ['cygrid/hphashtab.pyx'],
    extra_compile_args=['-fopenmp', '-O3', '-std=c++11'] + EX_COMP_ARGS,
    extra_link_args=['-fopenmp'],
    language='c++',
    include_dirs=[
        numpy.get_include(),
    ]
)

KERNEL_EXT = Extension(
    'cygrid.kernels',
    ['cygrid/kernels.pyx'],
    extra_compile_args=['-fopenmp', '-O3', '-std=c++11'] + EX_COMP_ARGS,
    extra_link_args=['-fopenmp'],
    language='c++',
    include_dirs=[
        numpy.get_include(),
    ]
)

setup(
    name='cygrid',
    version='0.9.5',
    author='Benjamin Winkel, Lars Flöer, Daniel Lenz',
    author_email='bwinkel@mpifr.de',
    description=(
        'Cygrid is a cython-powered convolution-based gridding module '
        'for astronomy'
        ),
    install_requires=[
        'setuptools',
        'cython>=0.20.2',
        'numpy>=1.8',
        'astropy>=1.0',
    ],
    packages=['cygrid'],
    cmdclass={'build_ext': build_ext},
    ext_modules=[
        KERNEL_EXT,
        HELPER_EXT,
        HPX_EXT,
        HPHASH_EXT,
        GRID_EXT,
    ],
    url='https://github.com/bwinkel/cygrid/',
    download_url='https://github.com/bwinkel/cygrid/tarball/0.9.4',
    keywords=['astronomy', 'gridding', 'fits/wcs']
)
