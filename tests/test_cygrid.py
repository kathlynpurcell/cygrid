#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import pytest
import numpy as np
# from numpy.testing import assert_equal, assert_allclose
# from astropy.utils.misc import NumpyRNGContext
from cygrid import WcsGrid, SlGrid


class TestWcsGrid:

    def setup(self):

        mapcenter = (131., 50.)
        mapsize = (1., 1.)  # degrees
        self.beamsize_fwhm = 0.05  # degrees

        avg_num_pixels_x = 5 * mapsize[0] / self.beamsize_fwhm  # per dimension
        avg_num_pixels_y = 5 * mapsize[1] / self.beamsize_fwhm  # per dimension

        self.xcoords = np.random.uniform(
            mapcenter[0] - mapsize[0]/2.,
            mapcenter[0] + mapsize[0]/2.,
            int(avg_num_pixels_x * avg_num_pixels_y),
            ).astype(np.float64)
        self.ycoords = np.random.uniform(
            mapcenter[1] - mapsize[1]/2.,
            mapcenter[1] + mapsize[1]/2.,
            int(avg_num_pixels_x * avg_num_pixels_y),
            ).astype(np.float64)
        self.signal = np.random.normal(0., .01, len(self.xcoords))

        pixsize = self.beamsize_fwhm / 3.
        dnaxis1 = int(mapsize[0] / pixsize + 0.5)
        dnaxis2 = int(mapsize[1] / pixsize + 0.5)
        projection = 'SIN'

        self.header = {
            'NAXIS': 3,
            'NAXIS1': dnaxis1,
            'NAXIS2': dnaxis2,
            'NAXIS3': 1,  # need dummy spectral axis
            'CTYPE1': 'RA---{}'.format(projection),
            'CTYPE2': 'DEC--{}'.format(projection),
            'CUNIT1': 'deg',
            'CUNIT2': 'deg',
            'CDELT1': -pixsize,
            'CDELT2': pixsize,
            'CRPIX1': dnaxis1 / 2.,
            'CRPIX2': dnaxis2 / 2.,
            'CRVAL1': mapcenter[0],
            'CRVAL2': mapcenter[1],
            }

    def teardown(self):

        pass

    def test_gridding(self):

        kernelsize_fwhm = self.beamsize_fwhm / 2
        kernelsize_sigma = kernelsize_fwhm / 2.35
        support_radius = 3. * kernelsize_sigma
        hpx_min_res = kernelsize_sigma / 2.

        mygridder = WcsGrid(self.header)
        mygridder.set_kernel(
            'gauss1d',
            (0.5 / kernelsize_sigma ** 2,),
            support_radius,
            hpx_min_res,
            )

        mygridder.grid(self.xcoords, self.ycoords, self.signal[:, np.newaxis])


class TestSlGrid:

    def setup(self):

        mapcenter = (131., 50.)
        mapsize = (1., 1.)  # degrees
        self.beamsize_fwhm = 0.05  # degrees

        avg_num_pixels_x = 5 * mapsize[0] / self.beamsize_fwhm  # per dimension
        avg_num_pixels_y = 5 * mapsize[1] / self.beamsize_fwhm  # per dimension

        self.input_x = np.random.uniform(
            mapcenter[0] - mapsize[0]/2.,
            mapcenter[0] + mapsize[0]/2.,
            int(avg_num_pixels_x * avg_num_pixels_y),
            ).astype(np.float64)
        self.input_y = np.random.uniform(
            mapcenter[1] - mapsize[1]/2.,
            mapcenter[1] + mapsize[1]/2.,
            int(avg_num_pixels_x * avg_num_pixels_y),
            ).astype(np.float64)
        self.signal = np.random.normal(0., .01, len(self.input_x))

        self.target_x = np.random.uniform(
            mapcenter[0] - mapsize[0]/2.,
            mapcenter[0] + mapsize[0]/2.,
            1000,
            ).astype(np.float64)
        self.target_y = np.random.uniform(
            mapcenter[0] - mapsize[0]/2.,
            mapcenter[0] + mapsize[0]/2.,
            1000,
            ).astype(np.float64)

    def teardown(self):

        pass

    def test_gridding(self):

        kernelsize_fwhm = self.beamsize_fwhm / 2
        kernelsize_sigma = kernelsize_fwhm / 2.35
        support_radius = 3. * kernelsize_sigma
        hpx_min_res = kernelsize_sigma / 2.

        mygridder = SlGrid(self.target_x, self.target_y, 1)
        mygridder.set_kernel(
            'gauss1d',
            (0.5 / kernelsize_sigma ** 2,),
            support_radius,
            hpx_min_res,
            )

        mygridder.grid(self.input_x, self.input_y, self.signal[:, np.newaxis])
