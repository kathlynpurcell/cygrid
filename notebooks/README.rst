*************************
Cygrid tutorial notebooks
*************************

This directory contains various `Jupyter <http://jupyter.org>`_ notebooks that
demonstrate how one can use the `cygrid` Python package. GitHub will render
the notebooks, but note that some of the tutorials can be best viewed with
the `nbviewer <https://nbviewer.jupyter.org/>`_ service (just copy-paste the
URL to the notebook into their search field.)

How to use cygrid
=================

Minimal example
---------------

In the `minimal example <01_minimal_example.ipynb>`_
notebook, we show the basic use of cygrid. Moreover, we illustrate some
applications of the WCS package from astropy.

Re-gridding vs. reprojection
----------------------------

Sometimes, it is desired to change the underlying coordinate system or
projection of a map. In such cases, the raw data is already present on a
regular grid. `cygrid` can be used to do such transformation in a very fast
and reliable manner, at the cost of a slight degradation of the angular
resolution of the data. Therefore, the `reproject package
<https://github.com/astrofrog/reproject>`_ might be a good alternative.

In this `notebook <02_regrid_from_healpix.ipynb>`_, we demonstrate how to transform a FITS-image to another projection, which is very useful if you want to combine data sets
of different instruments, e.g., HI and far-infrared data. The two approaches,
re-gridding and reprojection are compared.

Re-gridding (from HEALPix data)
-------------------------------

This is basically the same use case, as above, but shows a use case where
`cygrid` excels over re-project: if you have very large input maps and
want the results fast or want to down-sample the angular resolution while
strictly conserving the integrated flux density.

The `notebook <03_regrid_from_healpix.ipynb>`_
shows, how the healpy library can be used to create a random Gaussian field
from a given power spectrum. Then, `cygrid` is used to re-grid a portion
of the all-sky map to a a smaller FITS image.


Sight-line gridding
-------------------

Certainly a less common use case is to employ `cygrid` to grid a data set to
a list of (irregular) coordinates - which makes only sense for true data
gridding (and not re-projecting data) or if one needs to change the angular
resolution (down-sampling). One use case would be to create a HEALPix map
from a large raw data set. The `cygrid` sight-line gridder can do that, which
is demonstrated in `this notebook <04_sightline_gridding.ipynb>`_.

.. note::

    This notebook needs a complete overhaul. (Regridding from HPX to sight
    lines is not a good example, because interpolation would be better
    suited for this. Alternatively, we could high-light the down-sampling
    aspect.)

Additional material
===================

The convolutional-gridding technique
------------------------------------
`Cygrid` is based on the convolution of the raw data samples with a kernel
(a weighting function with finite support, e.g., a Gaussian). To visualize
what's going on under the hood, we made this `little video
<A01_convolutional_gridding.ipynb>`_.


An on-the-fly map exposure calculator based on cygrid
-----------------------------------------------------
As the title suggests, `cygrid` is like a swiss-army knife for real-life
issues of an astronomer,e.g., if your favorite observatory doesn't provide
you with a proper `exposure calculator <B01_OTF-map_exposure_calculator.ipynb>`_ ;-).
