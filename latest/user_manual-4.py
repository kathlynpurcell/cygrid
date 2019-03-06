import numpy as np
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.utils.misc import NumpyRNGContext
import cygrid

mapcenter = 60., 30.  # all in degrees
mapsize = 5., 5.
beamsize_fwhm = 0.1
num_samples = 10 ** 6
num_sources = 20

with NumpyRNGContext(1):
    lons, lats, signal = cygrid.produce_mock_data(
        mapcenter, mapsize, beamsize_fwhm, num_samples, num_sources
        )

pixsize = beamsize_fwhm / 3.
dnaxis1 = int(mapsize[0] / pixsize)
dnaxis2 = int(mapsize[1] / pixsize)

target_header = {
    'NAXIS': 2,
    'NAXIS1': dnaxis1,
    'NAXIS2': dnaxis2,
    'CTYPE1': 'RA---SIN',
    'CTYPE2': 'DEC--SIN',
    'CUNIT1': 'deg',
    'CUNIT2': 'deg',
    'CDELT1': -pixsize,
    'CDELT2': pixsize,
    'CRPIX1': dnaxis1 / 2.,
    'CRPIX2': dnaxis2 / 2.,
    'CRVAL1': mapcenter[0],
    'CRVAL2': mapcenter[1],
    }

gridder = cygrid.WcsGrid(target_header)

kernelsize_fwhm = 2.5 / 60.  # degrees
kernelsize_sigma = kernelsize_fwhm / np.log(8 * np.sqrt(2))
sphere_radius = 3. * kernelsize_sigma

gridder.set_kernel(
    'gauss1d',
    (kernelsize_sigma,),
    sphere_radius,
    kernelsize_sigma / 2.
    )
gridder.grid(lons, lats, signal)
gridded_map = gridder.get_datacube()

target_wcs = gridder.get_wcs()

fig = plt.figure(figsize=(7, 7))
ax = fig.add_subplot(111, projection=target_wcs.celestial)
im = ax.imshow(
    gridded_map, vmin=-0.5, vmax=8.0,
    origin='lower', interpolation='nearest'
    )
lon, lat = ax.coords
lon.set_axislabel('R.A. [deg]')
lat.set_axislabel('Dec [deg]')
plt.show()