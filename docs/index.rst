Documentation
=============

This is the documentation for cygrid.
Cygrid allows to resample a number of spectra (or data points) to a regular
grid - a data cube - using any valid astronomical FITS/WCS projection. The
method is a based on serialized convolution with finite gridding kernels.
Currently, only Gaussian (radial-symmetric or elliptical) kernels are provided
(which has the drawback of slight degradation of the effective resolution). The
algorithm has very small memory footprint, allows easy parallelization, and is
very fast.

.. toctree::
  :maxdepth: 2

  cygrid/index.rst

.. note:: The layout of this directory is simply a suggestion.  To follow
          traditional practice, do *not* edit this page, but instead place
          all documentation for the package inside ``cygrid/``.
          You can follow this practice or choose your own layout.
