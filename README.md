# Introduction #

- *Version*: 0.9
- *Authors*: Benjamin Winkel, Lars Flöer, Daniel Lenz

[![PyPI version](https://img.shields.io/pypi/v/cygrid.svg)](https://pypi.python.org/pypi/cygrid)
[![PyPI downloads](https://img.shields.io/pypi/dm/cygrid.svg)](https://pypi.python.org/pypi/cygrid)
[![Build Status](https://travis-ci.org/bwinkel/cygrid.svg?branch=master)](https://travis-ci.org/bwinkel/cygrid)
[![Publication](http://img.shields.io/badge/arXiv-1604.06667-blue.svg)](http://arxiv.org/abs/1604.06667)
[![License](https://img.shields.io/badge/license-GPL-blue.svg)](https://www.github.com/bwinkel/cygrid/blob/master/COPYING)

# Purpose#

`cygrid` allows to resample a number of spectra (or data points) to a regular grid - a data cube - using any valid astronomical FITS/WCS projection (see http://docs.astropy.org/en/stable/wcs/).

The method is a based on serialized convolution with finite gridding kernels. Currently, only Gaussian (radial-symmetric or elliptical) kernels are provided (which has the drawback of slight degradation of the effective resolution). The algorithm has very small memory footprint, allows easy parallelization, and is very fast.

A detailed description of the algorithm is given in [Winkel, Lenz & Flöer (2016)](http://adsabs.harvard.edu/abs/2016A%26A...591A..12W), which we kindly ask to be used as reference if you found `cygrid` useful for your research.

# Features

* Supports any WCS projection system as target.
* Conserves flux.
* Low memory footprint.
* Scales very well on multi-processor/core platforms.

# Usage #

### Installation ###

The easiest way to install cygrid is via `pip`:

```
pip install cygrid
```

The installation is also possible from source. Download the tar.gz-file, extract (or clone from GitHub) and simply execute

```
python setup.py install
```

### Dependencies ###

We kept the dependencies as minimal as possible. The following packages are
required:
* `numpy 1.10` or later
* `cython 0.23.4` or later
* `astropy 1.0` or later
(Older versions of these libraries may work, but we didn't test this!)

If you want to run the notebooks yourself, you will also need the Jupyter server, matplotlib and wcsaxes packages. To run the tests, you'll need HealPy.

Note, for compiling the C-extension, openmp is used for parallelization and some C++11 language features are necessary. If you use gcc, for example, you need at least version 4.7 otherwise the setup-script will fail. (If you have absolutely no possibility to upgrade gcc, older version may work if you replace `-std=c++11` with `-std=c++0x` in `setup.py`. Thanks to bs538 for pointing this out.)


### Examples ###

Check out the [`ipython notebooks`](http://nbviewer.jupyter.org/github/bwinkel/cygrid/blob/master/notebooks/index.ipynb) in the repository for some examples of how to use `cygrid`. Note that you only view them on the nbviewer service, and will have to clone the repository to run them on your machine.

### Who do I talk to? ###

If you encounter any problems or have questions, do not hesitate to raise an
issue or make a pull request. Moreover, you can contact the devs directly:

* <bwinkel@mpifr.de>
* <mail@daniellenz.org>
