#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup
from setuptools.extension import Extension
from Cython.Distutils import build_ext
import numpy
import platform
import os


VERSION = '0.9.8'


COMP_ARGS = {
    'extra_compile_args': ['-fopenmp', '-O3', '-std=c++11'],
    'extra_link_args': ['-fopenmp'],
    'language': 'c++',
    'include_dirs': [numpy.get_include()],
    }

# need to handle compilation on Windows and MAC-OS machines:
if platform.system().lower() == 'windows':
    COMP_ARGS['extra_compile_args'] = ['/openmp']
    del COMP_ARGS['extra_link_args']

if 'darwin' in platform.system().lower():
    COMP_ARGS['extra_compile_args'].append('-mmacosx-version-min=10.7')
    os.environ["CC"] = "g++-6"
    os.environ["CPP"] = "cpp-6"
    os.environ["CXX"] = "g++-6"
    os.environ["LD"] = "gcc-6"

# Cython extensions:
GRID_EXT = Extension(
    'cygrid.cygrid',
    ['cygrid/cygrid.pyx'],
    **COMP_ARGS
    )

HELPER_EXT = Extension(
    'cygrid.helpers',
    ['cygrid/helpers.pyx'],
    **COMP_ARGS
    )

HPX_EXT = Extension(
    'cygrid.healpix',
    ['cygrid/healpix.pyx'],
    **COMP_ARGS
    )

HPHASH_EXT = Extension(
    'cygrid.hphashtab',
    ['cygrid/hphashtab.pyx'],
    **COMP_ARGS
    )

KERNEL_EXT = Extension(
    'cygrid.kernels',
    ['cygrid/kernels.pyx'],
    **COMP_ARGS
    )

setup(
    name='cygrid',
    version=VERSION,
    author='Benjamin Winkel, Lars Flöer, Daniel Lenz',
    author_email='bwinkel@mpifr.de',
    description=(
        'Cygrid is a cython-powered convolution-based gridding module '
        'for astronomy'
        ),
    install_requires=[
        'setuptools',
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
    download_url='https://github.com/bwinkel/cygrid/tarball/{}'.format(VERSION),
    keywords=['astronomy', 'gridding', 'fits/wcs'],
    classifiers=[
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: MacOS',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python',
        'Programming Language :: Cython',
        'Topic :: Scientific/Engineering :: Astronomy',
        ],
)
