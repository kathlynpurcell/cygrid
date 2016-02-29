# Introduction #

- *Version*: 0.9
- *Authors*: Benjamin Winkel, Lars Flöer, Daniel Lenz

# Purpose#

`cygrid` allows to resample a number of spectra (or data points) to a regular grid - a data cube - using any valid astronomical FITs/WCS projection (TODO link).

The method is a based on serialized convolution with finite gridding kernels. Currently, only Gaussian (radial-symmetric or elliptical) kernels are provided (which has the drawback of slight degradation of the effective resolution). The algorithm has very small memory footprint, allows easy parallelization, and is very fast.

# Features

* Supports any WCS projection system as target.
* Conserves flux.
* Low memory footprint.
* Scales very well on multi-processor/core platforms.

# Usage #

### Installation ###

Right now, the installation should be done from source. Download the repository
and simply execute

```
python setup.py install
```

In the future, we will make sure you can use pip to install the package.

### Dependencies ###

We kept the dependencies as minimal as possible. The following packages are
required:
* `numpy 1.10` or later
* `cython 0.23.4` or later
* `astropy 1.0` or later

(Older versions of these libraries may work, but we didn't test this!)

### Examples ###

Check out the `ipython notebooks` in the repository for some examples of how to
use `cygrid`. In addition to the `cygrid` dependencies, you will need `healpy`
to run these.

### Who do I talk to? ###

If you encounter any problems or have questions, do not hesitate to raise an
issue or make a pull request. Moreover, you can contact the devs directly:

* <bwinkel@mpifr.de>
* <mail@lfloeer.de>
* <dlenz.bonn@gmail.com>
