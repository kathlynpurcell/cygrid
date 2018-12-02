.. user-manual:

******************
Cygrid user manual
******************

.. currentmodule:: cygrid



Introduction
============

The `~cygrid` package (see `Reference/API`_).



Data gridding
-------------

from https://docs.scipy.org/doc/scipy-0.18.1/reference/generated/scipy.interpolate.griddata.html

.. plot::
    :include-source:

    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.interpolate import griddata


    def func(x, y):
        return x*(1-x)*np.cos(4*np.pi*x) * np.sin(4*np.pi*y**2)**2

    grid_x, grid_y = np.mgrid[0:1:100j, 0:1:200j]

    # generate some data
    points = np.random.rand(1000, 2)
    values = func(points[:,0], points[:,1])

    gridded = griddata(points, values, (grid_x, grid_y), method='cubic')

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    ax1.imshow(func(grid_x, grid_y).T, extent=(0,1,0,1), origin='lower')
    ax1.plot(points[:,0], points[:,1], 'k.', ms=1)
    ax1.set_title('Original')
    ax2.imshow(gridded.T, extent=(0,1,0,1), origin='lower')
    ax2.set_title('Cubic interpolation')
    plt.show()


But what happens, if we add noise to the sampled data? Obviously, the gridded
data will not be a good description of the original function anymore. One
possible counter is to sample (aka observer) the function at more positions,
because this should decrease the noise in the gridded map (the
signal-to-noise ratio gets better). Unfortunately, the `~scipy.griddata` algorithm works in a way, which doesn't work in this way::

    # generate some data
    sigma = 0.01
    nsize = 1000
    points = np.random.rand(nsize, 2)
    noise = np.random.normal(0, sigma, nsize)

    values = func(points[:,0], points[:,1]) + noise

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

While one can see traces of the original function in the bottom right panel,
it is evident that there are strange outliers in the gridded map.

Convolution based gridding
--------------------------

One gridding algorithm, which does a better job at this, is so-called
convolutional gridding. Mathematically, it can be described by the following
formula:

.. math::

    R_{i,j}[s]=\frac{1}{W_{i,j}}\sum_n R_n[s](\alpha_n,\delta_n)w(\alpha_{i,j},\delta_{i,j};\alpha_n,\delta_n)\,.

where:

.. math::

    W_{i,j}\equiv\sum_n w(\alpha_{i,j},\delta_{i,j};\alpha_n,\delta_n)\,,

is called the weight map.

Here, :math:`R_n [s]` and :math:`R_{i,j}[s]` are two different representations
of the true signal :math:`s`. The index :math:`n` runs over the list of all
input samples, with respective coordinates :math:`(\alpha_n,\delta_n)`, while
the regular output grid can be parametrized via pixel coordinates :math:`(i,j)`
with associated world coordinates :math:`(\alpha_{i,j}, \delta_{i,j})`. The
value of the weighting kernel :math:`w` depends only on the input and output
coordinates. In most cases a radially symmetric kernel is applied, such that:

.. math::

    w(\alpha_{i,j},\delta_{i,j};\alpha_n,\delta_n)=w\left(\mathcal{D}(\alpha_{i,j},\delta_{i,j};\alpha_n,\delta_n)\right)

with the distance :math:`\mathcal{D}` between a pixel center world coordinate,
:math:`(\alpha_{i,j}, \delta_{i,j})`, and the :math:`n`-th input coordinate,
:math:`(\alpha_n,\delta_n)`.

TODO: add a simple python gridder, which visualizes this (maybe with animated
gif?)

The cygrid algorithm
--------------------
The great advantage of the approach described above, is, that it can be
straightforwardly applied to spherical geometry, which is necessary for map
making in geographic information systems and astronomy. The sampled coordinates
are often the longitude and latitude on a sphere, while one wants to display
the map in rectangular coordinates. Therefore, a certain map projections has
to be specified, as well.

However, there are libraries to easily convert between the world coordinates
(longitude and latitude) and pixel coordinates (:math:`i` and :math:`j`), such
as WCS (TODO REF). Then one only needs to use the true-angular distance
function instead of simple cartesian distance and the above formula will work.

Of course, in practice there is a little bit more to it. One has to deal with
edge effects and the double-sum over all :math:`n` and :math:`(i, j)` is a
dealbreaker for non-trivial map sizes. One can use a kernel function with
finite support, though, and restrict the summation to those input samples
:math:`n` that contribute significantly to :math:`(i, j)` (in spherical
coordinates, this is not an easy task!). On top of that, cygrid uses some
clever cashing techniques and can profit from multi-core CPUs to further
increase the computational speed of the gridding process. For details, we
refer to the Paper TODO REF.

.. _kernel-parameters-label:

Kernel parameters
---------------------------------



Comparison with scipy.griddata
-------------------------------

You're probably curious, how this compares to `~scipy.griddata`. Here is
the result:

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

.. note::

    Here we have made the approximation, that the input coordinates (ranging
    from 0 to 1 are in angles on the sphere (in degrees), and not the
    cartesian description. However, for small angles the distortion is
    negligible.)

Angular resolution
------------------

Drawback: angular resolution is degraded.

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

- `Astropy Units and Quantities package <http://docs.astropy.org/en/stable/
  units/index.html>`_, which is used extensively in cygrid.

Reference/API
=============

.. automodapi:: cygrid
    :inherited-members:
    :no-inheritance-diagram:
    :no-main-docstr:
..    :no-heading:

