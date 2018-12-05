.. user-manual:

******************
Cygrid user manual
******************

.. currentmodule:: cygrid



Introduction
============


`Cygrid` allows to resample a number of spectra (or data points) to a regular
grid - a data cube - using any valid astronomical FITS/WCS projection. The
method is a based on serialized convolution with finite gridding kernels.
Currently, only Gaussian (radial-symmetric or elliptical) kernels are provided
(which has the drawback of slight degradation of the effective resolution).
The algorithm has very small memory footprint, allows easy parallelization,
and is very fast.

Cygrid is already used in several "production" systems, for example it was
utilized for two major 21-cm HI surveys, EBHIS and HI4PI. Nevertheless,
we cannot guarantee that it's completely bug-free. We kindly invite you to
use the library and we are grateful for feedback. Note, that work on the
documentation is still ongoing.

The `~cygrid` package is available for Linux, Windows, and MacOS operating
systems. To improve computation speed, the `OpenMP technology
<https://www.openmp.org/>`_ is used. Therefore, a suitable C++ compiler must
be installed on your system, if you build from source. For convenience, we
also provide packages for the `Anaconda Python distribution
<https://www.anaconda.com/>`_.



Data gridding
-------------

Gridding data is a ubiquitous task in many scientific applications, for
example to make maps from measured raw data sets in astronomy or geographic
information systems (GIS). Therefore, it is not surprising that also the
popular `~scipy` library provides a function, `~scipy.interpolate.griddata`,
which can be used for this purpose:


.. from https://docs.scipy.org/doc/scipy-0.18.1/reference/generated/scipy.interpolate.griddata.html

.. plot::
    :include-source:

    # Example adapted from "scipy.interpolate.griddata" docs

    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.interpolate import griddata


    def func(x, y):

        return (
            x * (1 - x) *
            np.cos(4 * np.pi * x) *
            np.sin(4 * np.pi * y ** 2) ** 2
            )

    grid_x, grid_y = np.mgrid[0:1:100j, 0:1:200j]

    # generate some data
    points = np.random.rand(1000, 2)
    values = func(points[:,0], points[:,1])

    gridded = griddata(
        points,
        values,
        (grid_x, grid_y),
        method='cubic'
        )

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    ax1.imshow(
        func(grid_x, grid_y).T,
        extent=(0,1,0,1), origin='lower'
        )
    ax1.plot(points[:,0], points[:,1], 'k.', ms=1)
    ax1.set_title('Original')
    ax2.imshow(gridded.T, extent=(0,1,0,1), origin='lower')
    ax2.set_title('Cubic interpolation')
    plt.show()


The points in the left panel show the positions, where the underlying
function was sampled. As can be seen in the right panel,
`~scipy.interpolate.griddata` does a good job in estimating the function
values (on a regular grid) from these irregularly sampled points.

But what happens, if we add noise to the sampled data? Obviously, the gridded
data will somehow be affected by the decreased signal-to-noise ratio and may
not be a good description of the original (noise-free) function anymore. One
possible counter could be to sample (aka observe) the function at more
positions, because this should increase the signal-to-noise ratio in the
gridded map. Unfortunately, the `~scipy.interpolate.griddata` algorithm
doesn't work this way::

    # generate some data
    sigma = 0.01  # or 0.1
    nsize = 1000  # or 100000
    points = np.random.rand(nsize, 2)
    noise = np.random.normal(0, sigma, nsize)

    values = func(points[:,0], points[:,1]) + noise

    grid_and_plot(...)

.. plot::

    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.interpolate import griddata


    def func(x, y):
        return x*(1-x)*np.cos(4*np.pi*x) * np.sin(4*np.pi*y**2)**2

    grid_x, grid_y = np.mgrid[0:1:100j, 0:1:200j]

    fig, axes = plt.subplots(2, 2, figsize=(12, 12))
    for ax, (sigma, nsize) in zip(axes.flatten(), [
            (0.01, 1000), (0.01, 100000), (0.1, 1000), (0.1, 100000)
            ]):

        points = np.random.rand(nsize, 2)
        noise = np.random.normal(0, sigma, nsize)

        values = func(points[:,0], points[:,1]) + noise
        gridded = griddata(points, values, (grid_x, grid_y), method='cubic')

        ax.imshow(gridded.T, extent=(0,1,0,1), origin='lower', vmin=-0.3, vmax=0.3)
        ax.set_title('nsize: {:d}, sigma: {:.2f}'.format(nsize, sigma))

    plt.show()

With a significantly small number of samples or too much noise (e.g., lower
left panel), the gridded data is a very bad representation of the underlying
function. But even with a large input sample size, the
`~scipy.interpolate.griddata` algorithm produces strange outliers in the
gridded map (bottom right panel), but at least one can still see traces of the
original function.

Convolution based gridding
--------------------------

One gridding algorithm, which does a better job at this, is so-called
convolutional gridding. The idea behind this method is relatively simple. The
raw data samples (which are located at irregular coordinates) are convolved
with a kernel function, and then the result is computed at the position of the
pixel centers of the desired regular grid cells. If a Gaussian is used as
kernel, one can safe a lot of computing time, because each raw data sample
will only influence the output grid cells in within a certain sphere around
the raw-sample position.


Mathematically, this approach can be described by the following formula:

.. math::

    R_{i,j}[s]=\frac{1}{W_{i,j}}\sum_n R_n[s](\alpha_n,\delta_n)w(\alpha_{i,j},\delta_{i,j};\alpha_n,\delta_n)\,.

where:

.. math::

    W_{i,j}\equiv\sum_n w(\alpha_{i,j},\delta_{i,j};\alpha_n,\delta_n)\,,

is called the weight map.

Here, :math:`R_n [s]` and :math:`R_{i,j}[s]` are two different representations
of the true signal :math:`s`. The index :math:`n` runs over the list of all
input samples, with respective coordinates :math:`(\alpha_n,\delta_n)`, while
the regular output grid can be parametrized via pixel coordinates
:math:`(i,j)` with associated world coordinates :math:`(\alpha_{i,j},
\delta_{i,j})`. The value of the weighting kernel :math:`w` depends only on
the input and output coordinates. In most cases a radially symmetric kernel is
applied, such that:

.. math::

    w(\alpha_{i,j},\delta_{i,j};\alpha_n,\delta_n)=w\left(\mathcal{D}(\alpha_{i,j},\delta_{i,j};\alpha_n,\delta_n)\right)

with the distance :math:`\mathcal{D}` between a pixel center world coordinate,
:math:`(\alpha_{i,j}, \delta_{i,j})`, and the :math:`n`-th input coordinate,
:math:`(\alpha_n,\delta_n)`.

Since once usually doesn't want to keep all raw data samples in memory, we use
two helper maps. While we iterate over the samples (sum over :math:`n`), we
compute the products :math:`R_n[s](\alpha_n,\delta_n)w(\alpha_{i,j},\delta_{i,
j};\alpha_n,\delta_n)` for all pixels :math:`(i, j)` within a sphere around
the :math:`(\alpha_n,\delta_n)` position, and add them to an empty array
`A[i,j]`, which can be identified with :math:`R_{i,j}[s]\cdot W_{i,j}`.
Simultaneously, we do the same for the weighting factors and add these to
another array `B[i,j]` aka the weight map, :math:`W_{i,j}`. When all raw data
samples have been processed, the two resulting arrays are divided (`A[i,j] /
B[i,j]`), which gives :math:`R_{i,j}`.

.. note::

    In fact, `~cygrid` follows a slightly different approach internally,
    caching some intermediary values, to allow multi-threaded processing. For
    details, please see the `Cygrid paper`_.


.. note::

    In one of our `Jupyter tutorial notebooks <https://github.com/bwinkel/cygrid/blob/master/notebooks/A01_convolutional_gridding.ipynb>`_
    we made an animation that demonstrates how the algorithm works
    step-by-step.

The cygrid algorithm
====================

The great advantage of the approach described above, is, that it can be
straightforwardly applied to spherical geometry, which is necessary for map
making in geographic information systems and astronomy. The sampled
coordinates are often the longitude and latitude on a sphere, while one wants
to display the map in rectangular coordinates. Therefore, a certain map
projections has to be specified, as well.

However, there are libraries to easily convert between the world coordinates
(longitude and latitude) and pixel coordinates (:math:`i` and :math:`j`), such
as `~astropy.wcs`. Then one only needs to use the true-angular distance
function instead instead of the simple Cartesian distance and the above
formulas can be applied.

Of course, in practice there is a little bit more to it. One has to deal with
edge effects and the double-sum over all :math:`n` and :math:`(i, j)` is a
deal breaker for non-trivial map sizes. One can use a kernel function with
finite support, though, and restrict the summation to those input samples
:math:`n` that contribute significantly to :math:`(i, j)` (in spherical
coordinates, it is not an easy task to identify pixels within a certain
sphere!). On top of that, cygrid uses some clever cashing techniques and can
profit from multi-core CPUs to further increase the computational speed of the
gridding process. For details, we refer to the `Cygrid paper`_.


Comparison with scipy.griddata
-------------------------------

After reading all of the above, you're now probably curious, how `~cygrid`
compares to `~scipy.interpolate.griddata`. Here is the result:

.. plot::

    import numpy as np
    import matplotlib.pyplot as plt
    import cygrid


    def func(x, y):
        return x*(1-x)*np.cos(4*np.pi*x) * np.sin(4*np.pi*y**2)**2

    grid_x, grid_y = np.meshgrid(
        np.linspace(0, 1, 100), np.linspace(0, 1, 100)
        )

    fig, axes = plt.subplots(2, 2, figsize=(12, 12))
    for ax, (sigma, nsize) in zip(axes.flatten(), [
            (0.01, 1000), (0.01, 100000), (0.1, 1000), (0.1, 100000)
            ]):

        points = np.random.rand(2, nsize)
        noise = np.random.normal(0, sigma, nsize)

        values = func(points[0], points[1]) + noise

        gridder = cygrid.SlGrid(grid_x.flatten(), grid_y.flatten(), 1)
        kernelsize_fwhm = 0.05
        kernelsize_sigma = kernelsize_fwhm / 2.355
        sphere_radius = 3. * kernelsize_sigma

        gridder.set_kernel(
            'gauss1d',
            (kernelsize_sigma,),
            sphere_radius,
            kernelsize_sigma / 2.
            )
        gridder.grid(points[0], points[1], values[:, None])
        gridded = gridder.get_datacube().squeeze().reshape(grid_x.shape)

        ax.imshow(gridded, extent=(0,1,0,1), origin='lower', vmin=-0.3, vmax=0.3)
        ax.set_title('nsize: {:d}, sigma: {:.2f}'.format(nsize, sigma))

    plt.show()

We think, it is fair to say that `~cygrid` produces much better result
for these complicated situations.

.. note::

    Here we have made the approximation, that the input coordinates (ranging
    from 0 to 1 are in angles on the sphere (in degrees), and not the linear
    Cartesian coordinates. However, for small angles the distortion is
    negligible.)


.. _angular-resolution-label:

Angular resolution
------------------


Because `~cygrid` is based on convolution with a (typically) Gaussian kernel,
the effective angular resolution of the data is degraded:

.. math::

    \sigma_\mathrm{out}^2 = \sigma_\mathrm{in}^2 + \sigma_k^2

where :math:`\sigma_\mathrm{out}` is the standard deviation of the resulting
effective beam width (given it's Gaussian), :math:`\sigma_\mathrm{in}` is the
width of the input beam, and :math:`\sigma_{k}` is the width of the kernel.

.. note::

    One can convert between the standard-deviation width and the FWHM via
    :math:`\sigma = \mathrm{FWHM} / \sqrt{8\ln2}`.

A good practical compromise is often to let the kernel be about half the
size of the input resolution, because then the resulting resolution is
only degraded by about 10%.

There are several additional considerations, which play a role. The kernel
must not be too narrow, especially if the input data are not sufficiently
densely sampled. Each raw data sample can only "reach" output grid pixels
that are well within an angular distance of :math:`<3\sigma_k`. Furthermore,
the output pixel grid must be fine enough to warrant full sampling of the
output signal: :math:`\Delta p\lesssim\mathrm{FWHM}_\mathrm{out}/2\sqrt{2}=\sigma_\mathrm{out}/\sqrt{\ln2}`.
We'd even recommend to chose the output grid such that
:math:`\Delta p\lesssim\mathrm{FWHM}_k/2\sqrt{2}=\sigma_k/\sqrt{\ln2}`.
This is because for performance reasons `~cygrid` does compute the kernel
function values only for the pixel centers, and doesn't perform numerical
integration over the pixel area. Of course, this leads to a certain numerical
inaccuracy, but if the output pixel grid supports full sampling of the kernel,
the errors are on an acceptable level. For further details we refer to
the `Cygrid paper`_.

Under some circumstances, the decrease in angular resolution is even
desirable. If one wants to compare two data sets with different angular
resolution, it is possible to smooth the higher-resolution data to
the lower-resolution map by choosing :math:`\sigma_k^2=\sigma_\mathrm{low}^2 - \sigma_\mathrm{high}^2`.

.. _kernel-parameters-label:

Kernel parameters
---------------------------------

Apart from choosing a proper kernel size (see
:ref:`angular-resolution-label`), one can use three different weighting
functions, `gauss1d`, `gauss2d`, and `tapered_sinc`. Despite its name,
`gauss1d` refers to a two-dimensional radial-symmetric Gaussian - in
contrast to `gauss2d`, which is for an elliptical 2D Gaussian. More
information, e.g., which parameters are required for each of the three
functions, can be found in the `~cygrid.WcsGrid.set_kernel` method
documentation.

There are two other necessary parameters, to be provided to the
`~cygrid.WcsGrid.set_kernel` method: the kernel sphere radius (or support
radius) and the resolution of the internally used `HEALPix
<https://healpix.jpl.nasa.gov/>`_ grid. The sphere radius should be
:math:`3\ldots5\sigma_k` depending on the desired accuracy. As the area of a
sphere goes with radius squared, using :math:`5\sigma_k` will take roughly
three times longer. The internally HPX resolution defines only some details
of the used caches. A good value for this is :math:`\sigma_k/2`.



How to use cygrid? Simple tasks
======================================

TODO: this was shifted to quick tour

Simple example::

    from astropy.io import fits
    import cygrid

    # read-in data
    lon, lat, signal = get_data()

    # define target FITS/WCS header
    header = create_fits_header()

    # prepare gridder
    kernelsize_sigma = 0.2

    kernel_type = 'gauss1d'
    kernel_params = (kernelsize_sigma, )  # must be a tuple
    kernel_support = 3 * kernelsize_sigma
    hpx_maxres = kernelsize_sigma / 2

    mygridder = cygrid.WcsGrid(header)
    mygridder.set_kernel(
        kernel_type,
        kernel_params,
        kernel_support,
        hpx_maxres,
        )

    # do the actual gridding
    mygridder.grid(lon, lat, signal)

    # query result and store to disk
    data_cube = mygridder.get_datacube()
    fits.writeto('example.fits', header=header, data=data_cube)


How does cygrid compare to scipy.griddata or reproject?
--------------------------------------------------------


How to use cygrid? Advanced tasks
=========================================

.. _serialization-label:

Decrease memory footprint
-------------------------


Sight-line gridding
-------------------




Benchmarking
============
see paper


See Also
========

- `Cygrid paper <http://adsabs.harvard.edu/abs/2016A%26A...591A..12W>`_: B. Winkel, D. Lenz & L. Flöer: *Cygrid: A fast Cython-powered convolution-based gridding module for Python*, Astronomy & Astrophysics, Volume 591, A 12, 2016.
- `HEALPix paper <http://adsabs.harvard.edu/abs/2005ApJ...622..759G>`_: K. M. Górski, E. Hivon, A. J. Banday, B. D. Wandelt, F. K. Hansen, M. Reinecke, M. Bartelmann: *HEALPix: A Framework for High-Resolution Discretization and Fast Analysis of Data Distributed on the Sphere*, The Astrophysical Journal, Volume 622, Issue 2, 2005.
- `WCSlib paper I <http://adsabs.harvard.edu/abs/2002A%26A...395.1061G>`_: E. W. Greisen & M. R. Calabretta: *Representations of world coordinates in FITS*, Astronomy & Astrophysics, Volume 395, p.1061, 2002.
- `WCSlib paper II <http://adsabs.harvard.edu/abs/2002A%26A...395.1061G>`_: M. R. Calabretta & E. W. Greisen: *Representations of celestial coordinates in FITS*, Astronomy & Astrophysics, Volume 395, p.1077, 2002.
- `Astropy World Coodinate System package <http://docs.astropy.org/en/stable/
  wcs/index.html>`_, which is used extensively in cygrid.


Reference/API
=============

.. automodapi:: cygrid
    :inherited-members:
    :no-inheritance-diagram:
    :no-main-docstr:
..    :no-heading:

