#!/usr/bin/env python
# -*- coding: utf-8 -*-

from distutils.core import setup
from distutils.extension import Extension
from distutils.version import StrictVersion, LooseVersion
from Cython.Distutils import build_ext
import warnings
import sys
import numpy
import platform

EX_COMP_ARGS = []
if 'mac' in platform.system().lower():
    EX_COMP_ARGS += ['-stdlib=libc++', '-mmacosx-version-min=10.7']


# Dependency checking
dependencies = [['cython', '0.20.2'], ['numpy', '1.8'], ['astropy', '1.0']]

for (pkg, minversion) in dependencies:
    try:
        m = __import__(pkg)
        if minversion is not None:
            if StrictVersion(m.__version__) < StrictVersion(minversion):
                if LooseVersion(m.__version__) < LooseVersion(minversion):
                    raise ValueError
                warnings.warn(
                    'Version', m.__version__,
                    'of package', pkg,
                    'might not be sufficient'
                    )
    except ImportError:
        print 'Package', pkg, 'not present.'
        sys.exit(1)
    except ValueError:
        print 'Package', pkg, 'has version', m.__version__
        print 'Version', minversion, 'required.'
        sys.exit(1)


GRID_EXT = Extension(
    "cygrid.cygrid",
    ["cygrid/cygrid.pyx"],
    extra_compile_args=['-fopenmp', '-O3', '-std=c++11'] + EX_COMP_ARGS,
    extra_link_args=['-fopenmp'],
    language='c++',
    include_dirs=[
        numpy.get_include(),
    ]
)

HELPER_EXT = Extension(
    "cygrid.helpers",
    ["cygrid/helpers.pyx"],
    extra_compile_args=['-fopenmp', '-O3', '-std=c++11'] + EX_COMP_ARGS,
    extra_link_args=['-fopenmp'],
    language='c++',
    include_dirs=[
        numpy.get_include(),
    ]
)

HPX_EXT = Extension(
    "cygrid.healpix",
    ["cygrid/healpix.pyx"],
    extra_compile_args=['-fopenmp', '-O3', '-std=c++11'] + EX_COMP_ARGS,
    extra_link_args=['-fopenmp'],
    language='c++',
    include_dirs=[
        numpy.get_include(),
    ]
)

HPRB_EXT = Extension(
    "cygrid.hprainbow",
    ["cygrid/hprainbow.pyx"],
    extra_compile_args=['-fopenmp', '-O3', '-std=c++11'] + EX_COMP_ARGS,
    extra_link_args=['-fopenmp'],
    language='c++',
    include_dirs=[
        numpy.get_include(),
    ]
)

KERNEL_EXT = Extension(
    "cygrid.kernels",
    ["cygrid/kernels.pyx"],
    extra_compile_args=['-fopenmp', '-O3', '-std=c++11'] + EX_COMP_ARGS,
    extra_link_args=['-fopenmp'],
    language='c++',
    include_dirs=[
        numpy.get_include(),
    ]
)

setup(
    name='cygrid',
    version="0.9",
    author="Benjamin Winkel, Lars Flöer, Daniel Lenz",
    author_email="bwinkel@mpifr.de",
    description=(
        "Cygrid is a cython-powered convolution-based gridding module "
        "for astronomy"
        ),
    packages=['cygrid'],
    cmdclass={'build_ext': build_ext},
    ext_modules=[
        KERNEL_EXT,
        HELPER_EXT,
        HPX_EXT,
        HPRB_EXT,
        GRID_EXT,
    ],
    url = 'https://github.com/bwinkel/cygrid/',
    download_url = 'https://github.com/bwinkel/cygrid/tarball/0.9',
    keywords=['astronomy', 'gridding', 'fits/wcs']
)
