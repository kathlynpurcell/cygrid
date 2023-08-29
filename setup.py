#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

# Note: This file needs to be Python 2 / <3.6 compatible, so that the nice
# "This package only supports Python 3.x+" error prints without syntax errors etc.
from setuptools import setup
import glob
import os
import sys
try:
    from configparser import ConfigParser
except ImportError:
    from ConfigParser import ConfigParser
# adapted from https://github.com/astropy/astropy


import sys


TEST_HELP = """
Note: running tests is no longer done using 'python setup.py test'. Instead
you will need to run:
    tox -e test
If you don't already have tox installed, you can install it with:
    pip install tox
If you only want to run part of the test suite, you can also use pytest
directly with::
    python -m pip install --no-build-isolation --no-deps -v -v -v -e .
    python -m pytest -rsx --ignore-glob="*/setup_package.py" cygrid
    # or individual tests:
    python -m pytest -rsx --ignore-glob="*/setup_package.py" cygrid -k <test_func_name/module_name/etc.>
    # doc tests
    python -m pytest -rsx --doctest-plus --doctest-glob="*.rst" --doctest-ignore-import-errors docs

For more information, see:
  https://docs.astropy.org/en/latest/development/testguide.html#running-tests
"""

if 'test' in sys.argv:
    print(TEST_HELP)
    sys.exit(1)


# Import ah_bootstrap after the python version validation

import ah_bootstrap


# A dirty hack to get around some early import/configurations ambiguities
if sys.version_info[0] >= 3:
    import builtins
else:
    import __builtin__ as builtins
builtins._ASTROPY_SETUP_ = True

from astropy_helpers.setup_helpers import (register_commands, get_debug_option,
                                           get_package_info)
from astropy_helpers.git_helpers import get_git_devstr
from astropy_helpers.version_helpers import generate_version_py


# order of priority for long_description:
#   (1) set in setup.cfg,
#   (2) load LONG_DESCRIPTION.rst,
#   (3) load README.rst,
#   (4) package docstring
readme_glob = 'README*'
_cfg_long_description = metadata.get('long_description', '')
if _cfg_long_description:
    LONG_DESCRIPTION = _cfg_long_description

elif os.path.exists('LONG_DESCRIPTION.rst'):
    with open('LONG_DESCRIPTION.rst') as f:
        LONG_DESCRIPTION = f.read()

elif len(glob.glob(readme_glob)) > 0:
    with open(glob.glob(readme_glob)[0]) as f:
        LONG_DESCRIPTION = f.read()

else:
    # Get the long description from the package's docstring
    __import__(PACKAGENAME)
    package = sys.modules[PACKAGENAME]
    LONG_DESCRIPTION = package.__doc__

# Store the package name in a built-in variable so it's easy
# to get from other parts of the setup infrastructure
builtins._ASTROPY_PACKAGE_NAME_ = PACKAGENAME


DOCS_HELP = """
Note: building the documentation is no longer done using
'python setup.py build_docs'. Instead you will need to run:
    tox -e build_docs
If you don't already have tox installed, you can install it with:
    pip install tox
You can also build the documentation with Sphinx directly using::
    python -m pip install --no-build-isolation --no-deps -v -v -v -e .
    cd docs
    # make clean  # to rebuild everything
    make html
    # alternatively (in project dir):
    sphinx-build docs docs/_build/html -b html
    sphinx-build docs docs/_build/html -b html -W  # fail on warnings

For more information, see:
  https://docs.astropy.org/en/latest/install.html#builddocs
"""

if 'build_docs' in sys.argv or 'build_sphinx' in sys.argv:
    print(DOCS_HELP)
    sys.exit(1)


# Only import these if the above checks are okay
# to avoid masking the real problem with import error.
from setuptools import setup  # noqa
from extension_helpers import get_extensions  # noqa

# import numpy as np
# np.import_array()

setup(ext_modules=get_extensions())
