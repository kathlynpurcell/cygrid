#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np


__all__ = ['produce_mock_data']


def produce_mock_data(
        mapcenter, mapsize, beamsize_fwhm, num_samples, num_sources
        ):
    '''
    Produce mock raw data (including coordinate pairs) for testing.

    The raw data mimics a continuum observation (i.e., one flux density
    value per position) containing only Gaussian noise and some point
    sources.

    Parameters
    ----------
    mapcenter : tuple of floats
        Center coordinates, `(lon, lat)` of the generated (rectangular)
        map [deg].
    mapsize : tuple of floats
        Size of map, `(width, height)` of the generated map [deg].
    beamsize_fwhm : float
        FWHM beam size [deg].
        The artificial point sources, which are sampled, will be convolved
        with this.
    num_samples : int
        Number of samples to produce. Note that this must be large enough
        to ensure full sampling of the map (otherwise artifacts will
        appear, which are caused by so-called aliasing).
    num_sources : int
        Number of artificial point sources to be placed into the map.

    Returns
    -------
    lons : `~numpy.array` of float64 [1D]
        Longitude coordinates for each sample.
    lats : `~numpy.array` of float64 [1D]
        Latitude coordinates for each sample.
    signal : `~numpy.array` of float64 [1D]
        The artifical raw data signal for each position.

    Notes
    -----
    As in real astronomical measurements, the point sources are convolved
    with the instrument's response function (PSF), or telescope beam.
    The longitude map size is scaled with :math:`\\cos b`, where :math:`b`
    is the latitude. This is to ensure that the true-angular size roughly
    matches the given map width.
    '''

    lon_scale = np.cos(np.radians(mapcenter[1]))
    map_l, map_r = (
        mapcenter[0] - 1.1 * mapsize[0] / 2. / lon_scale,
        mapcenter[0] + 1.1 * mapsize[0] / 2. / lon_scale
        )
    map_b, map_t = (
        mapcenter[1] - 1.1 * mapsize[1] / 2.,
        mapcenter[1] + 1.1 * mapsize[1] / 2.,
        )

    # coordinates are drawn from a uniform distribution
    lons = np.random.uniform(map_l, map_r, num_samples).astype(np.float64)
    lats = np.random.uniform(map_b, map_t, num_samples).astype(np.float64)

    # add Gaussian noise
    signal = np.random.normal(0., 1., len(lons))

    beamsize_sigma = beamsize_fwhm / np.sqrt(8 * np.log(2))

    # put in artifical point source, with random amplitudes
    # we'll assume a Gaussian-shaped PSF

    def gauss2d(x, y, x0, y0, A, s):
        return A * np.exp(-((x - x0) ** 2 + (y - y0) ** 2) / 2. / s ** 2)

    for _ in range(num_sources):

        sou_x = np.random.uniform(map_l, map_r, 1).astype(np.float64)
        sou_y = np.random.uniform(map_b, map_t, 1).astype(np.float64)
        A = np.random.uniform(0, 10, 1)
        signal += gauss2d(lons, lats, sou_x, sou_y, A, beamsize_sigma)

    return lons, lats, signal
