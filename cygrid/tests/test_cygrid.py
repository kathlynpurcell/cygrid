#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import pytest
import numpy as np
from numpy.testing import assert_allclose
from astropy.utils.misc import NumpyRNGContext
from astropy.utils.data import get_pkg_data_filename
import cygrid


class TestWcsGrid:

    def setup(self):

        mapcenter = (131., 50.)
        mapsize = (1., 1.)  # degrees
        self.beamsize_fwhm = 0.05  # degrees

        avg_num_pixels_x = 5 * mapsize[0] / self.beamsize_fwhm
        avg_num_pixels_y = 5 * mapsize[1] / self.beamsize_fwhm

        scale = np.cos(np.radians(mapcenter[1]))

        with NumpyRNGContext(0):
            self.xcoords = np.random.uniform(
                mapcenter[0] - mapsize[0] / 2. / scale,
                mapcenter[0] + mapsize[0] / 2. / scale,
                int(avg_num_pixels_x * avg_num_pixels_y),
                ).astype(np.float64)
            self.ycoords = np.random.uniform(
                mapcenter[1] - mapsize[1] / 2.,
                mapcenter[1] + mapsize[1] / 2.,
                int(avg_num_pixels_x * avg_num_pixels_y),
                ).astype(np.float64)
            self.signal = np.random.normal(0., .01, len(self.xcoords))

        self.signal2 = np.column_stack([self.signal, self.signal ** 2])
        pixsize = self.beamsize_fwhm / 3.
        dnaxis1 = int(mapsize[0] / pixsize + 0.5)
        dnaxis2 = int(mapsize[1] / pixsize + 0.5)
        projection = 'SIN'

        self.header = {
            'NAXIS': 2,
            'NAXIS1': dnaxis1,
            'NAXIS2': dnaxis2,
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

        kernelsize_fwhm = self.beamsize_fwhm / 2
        kernelsize_sigma = kernelsize_fwhm / 2.35
        support_radius = 3. * kernelsize_sigma
        hpx_min_res = kernelsize_sigma / 2.

        self.kernel_args = (
            'gauss1d',
            (0.5 / kernelsize_sigma ** 2,),
            support_radius,
            hpx_min_res,
            )

        self.test_maps = np.load(get_pkg_data_filename(
            'data/cygrid_test_maps.npy'
            ))

    def teardown(self):

        pass

    def _do_gridding2d(self, signal, header, **kwargs):

        mygridder = cygrid.WcsGrid(header, **kwargs)
        mygridder.set_kernel(*self.kernel_args)

        mygridder.grid(self.xcoords, self.ycoords, signal[:, np.newaxis])

        assert_allclose(
            mygridder.get_datacube()[0],
            self.test_maps[0],
            atol=1.e-6
            )

        # test automatic patching of input signal array size (to 2D)
        mygridder.grid(self.xcoords, self.ycoords, signal)

        assert_allclose(
            mygridder.get_datacube()[0],
            self.test_maps[0],
            atol=1.e-6
            )

    def _do_gridding3d(self, signal, header, **kwargs):

        mygridder = cygrid.WcsGrid(header, **kwargs)
        mygridder.set_kernel(*self.kernel_args)

        mygridder.grid(self.xcoords, self.ycoords, signal)

        if kwargs.get('do_store', False):
            np.save('/tmp/cygrid_test_maps.npy', mygridder.get_datacube())

        assert_allclose(mygridder.get_datacube(), self.test_maps, atol=1.e-6)

    def test_gridding2d_naxis3_implicit(self):
        '''
        Implicitly set 3rd (redundant) axis in header
        '''

        header = self.header.copy()
        header['NAXIS'] = 3
        header['NAXIS3'] = 1
        self._do_gridding2d(self.signal, header)

    def test_gridding2d_naxis3_explicit(self):
        '''
        Explicitly set 3rd via kwargs
        '''

        self._do_gridding2d(self.signal, self.header, naxis3=1)

    def test_gridding2d_naxis3_auto(self):
        '''
        Automatic handling of naxis3 (aka set to One, if not present)
        '''

        self._do_gridding2d(self.signal, self.header)

    def test_gridding3d_naxis3_implicit(self):
        '''
        Implicitly set 3rd (redundant) axis in header
        '''

        header = self.header.copy()
        header['NAXIS'] = 3
        header['NAXIS3'] = 2
        self._do_gridding3d(self.signal2, header)

    def test_gridding3d_naxis3_explicit(self):
        '''
        Explicitly set 3rd via kwargs
        '''

        self._do_gridding3d(
            self.signal2, self.header, naxis3=2, do_store=False
            )

    def test_shape_error(self):
        '''
        Implicitly set 3rd (redundant) axis in header
        '''

        header = self.header.copy()
        header['NAXIS'] = 3
        header['NAXIS3'] = 2
        with pytest.raises(cygrid.ShapeError):
            self._do_gridding2d(self.signal, header)

        header['NAXIS3'] = 1
        with pytest.raises(cygrid.ShapeError):
            self._do_gridding3d(self.signal2, header)

        header = self.header.copy()
        with pytest.raises(cygrid.ShapeError):
            self._do_gridding2d(self.signal, header, naxis3=2)

        with pytest.raises(cygrid.ShapeError):
            self._do_gridding3d(self.signal2, header, naxis3=1)

    def test_c_contiguous(self):
        '''
        Cygrid should autocast to C-contiguous if necessary
        '''

        signal2_f_cont = np.require(self.signal2, self.signal2.dtype, 'F')
        self._do_gridding3d(signal2_f_cont, self.header, naxis3=2)


class TestSlGrid:

    def setup(self):

        mapcenter = (131., 50.)
        mapsize = (1., 1.)  # degrees
        self.beamsize_fwhm = 0.05  # degrees

        avg_num_pixels_x = 5 * mapsize[0] / self.beamsize_fwhm
        avg_num_pixels_y = 5 * mapsize[1] / self.beamsize_fwhm

        scale = np.cos(np.radians(mapcenter[1]))

        with NumpyRNGContext(0):
            self.xcoords = np.random.uniform(
                mapcenter[0] - mapsize[0] / 2. / scale,
                mapcenter[0] + mapsize[0] / 2. / scale,
                int(avg_num_pixels_x * avg_num_pixels_y),
                ).astype(np.float64)
            self.ycoords = np.random.uniform(
                mapcenter[1] - mapsize[1] / 2.,
                mapcenter[1] + mapsize[1] / 2.,
                int(avg_num_pixels_x * avg_num_pixels_y),
                ).astype(np.float64)
            self.signal = np.random.normal(0., .01, len(self.xcoords))

            self.target_x = np.random.uniform(
                mapcenter[0] - mapsize[0] / 2. / scale,
                mapcenter[0] + mapsize[0] / 2. / scale,
                1000,
                ).astype(np.float64)
            self.target_y = np.random.uniform(
                mapcenter[1] - mapsize[1] / 2.,
                mapcenter[1] + mapsize[1] / 2.,
                1000,
                ).astype(np.float64)

        kernelsize_fwhm = self.beamsize_fwhm / 2
        kernelsize_sigma = kernelsize_fwhm / 2.35
        support_radius = 3. * kernelsize_sigma
        hpx_min_res = kernelsize_sigma / 2.

        self.kernel_args = (
            'gauss1d',
            (0.5 / kernelsize_sigma ** 2,),
            support_radius,
            hpx_min_res,
            )

        self.test_sls = np.load(get_pkg_data_filename(
            'data/cygrid_test_sightlines.npy'
            ))

    def teardown(self):

        pass

    def _do_gridding(self, signal, naxis3, **kwargs):

        mygridder = cygrid.SlGrid(
            self.target_x,
            self.target_y,
            naxis3
            )
        mygridder.set_kernel(*self.kernel_args)
        mygridder.grid(self.xcoords, self.ycoords, signal)

        if kwargs.get('do_store', False):
            np.save(
                '/tmp/cygrid_test_sightlines.npy',
                mygridder.get_datacube()
                )

        if naxis3 == 1:
            assert_allclose(
                mygridder.get_datacube()[0],
                self.test_sls[0],
                atol=1.e-6
                )
        else:
            assert_allclose(
                mygridder.get_datacube(),
                self.test_sls,
                atol=1.e-6
                )

    def test_gridding_naxis3_1(self):

        self._do_gridding(self.signal, 1)

    def test_gridding_naxis3_2(self):

        signal2 = np.column_stack([self.signal, self.signal ** 2])
        self._do_gridding(signal2, 2, do_store=False)

    def test_shape_error(self):

        signal2 = np.column_stack([self.signal, self.signal ** 2])

        with pytest.raises(cygrid.ShapeError):
            self._do_gridding(self.signal, 2)

        with pytest.raises(cygrid.ShapeError):
            self._do_gridding(signal2, 1)
