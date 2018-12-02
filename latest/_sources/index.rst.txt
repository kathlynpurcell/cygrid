
:tocdepth: 3

#####################
cygrid Documentation
#####################

This is the documentation for cygrid.

Cygrid allows to resample a number of spectra (or data points) to a regular
grid - a data cube - using any valid astronomical FITS/WCS projection. The
method is a based on serialized convolution with finite gridding kernels.
Currently, only Gaussian (radial-symmetric or elliptical) kernels are provided
(which has the drawback of slight degradation of the effective resolution). The
algorithm has very small memory footprint, allows easy parallelization, and is
very fast.

***************
Getting Started
***************

.. toctree::
   :maxdepth: 1

   install
   importing_cygrid
   quick_tour

******************
User Documentation
******************

.. toctree::
   :maxdepth: 1

   user_manual
   Tutorials (Jupyter notebooks) <http://nbviewer.jupyter.org/github/bwinkel/cygrid/blob/master/notebooks/>


***************
Project details
***************

.. toctree::
   :maxdepth: 1

   license

***************
Acknowledgments
***************

This code makes use of the excellent work provided by the
`Astropy <http://www.astropy.org/>`__ community. pycraf uses the Astropy package and also the
`Astropy Package Template <https://github.com/astropy/package-template>`__
for the packaging.
