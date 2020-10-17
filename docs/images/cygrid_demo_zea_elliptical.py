#!/usr/bin/python
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import numpy as np
from kapteyn import maputils
import matplotlib.pyplot as plt
import cygrid
from astropy.io import fits as pf
from astropy import wcs

plt.rc('text', usetex=True)


def true_angular_distance(l1, b1, l2, b2):

    sin_diff_lon = np.sin(np.radians(l2 - l1))
    cos_diff_lon = np.cos(np.radians(l2 - l1))
    sin_lat1 = np.sin(np.radians(b1))
    sin_lat2 = np.sin(np.radians(b2))
    cos_lat1 = np.cos(np.radians(b1))
    cos_lat2 = np.cos(np.radians(b2))

    num1 = cos_lat2 * sin_diff_lon
    num2 = cos_lat1 * sin_lat2 - sin_lat1 * cos_lat2 * cos_diff_lon
    denominator = sin_lat1 * sin_lat2 + cos_lat1 * cos_lat2 * cos_diff_lon

    return np.degrees(np.arctan2(
        np.sqrt(num1 ** 2 + num2 ** 2), denominator
        ))


def astro_signal(xcoords, ycoords, beamsize_fwhm):
    '''
    Create an example signal (noise+baseline+sources)
    '''

    beamsize_sigma = beamsize_fwhm / 2.35

    def gauss2d(x, y, x0, y0, A, s):
        r = true_angular_distance(x, y, x0, y0)
        return A * np.exp(-(r**2) / 2. / s**2)

    s_xcoords = np.arange(0., 359., 45.)
    s_ycoords = np.ones_like(s_xcoords) * 88.5
    A = 10.
    s = beamsize_sigma

    signal = np.random.normal(0., 1., xcoords.size)
    # signal += np.sin(np.radians(xcoords) / 360. * 2. * np.pi * 10 + 0.3) * \
        # np.sin(np.radians(ycoords) / 360. * 2. * np.pi * 2 + 0.1)

    for s_x, s_y in zip(s_xcoords, s_ycoords):
        signal += gauss2d(xcoords, ycoords, s_x, s_y, A, s)

    return signal

def do_plot(header, data, figname, showplot=False):
    # now plot
    f = maputils.FITSimage(externalheader=header, externaldata=data)

    cmap = str('YlGnBu')
    cmap = str('Greys')

    X = np.linspace(0., 300., 6)
    Y = np.linspace(87., 89., 3)

    plt.close()
    fig = plt.figure(figsize=(6, 6))
    frame = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    mplim = f.Annotatedimage(frame, cmap=cmap)
    mplim.Image(interpolation='nearest')
    grat = mplim.Graticule(
        axnum=(1, 2), wylim=(87, 90), wxlim=(-180, 180), startx=X, starty=Y
        )
    grat.setp_gratline(color='0.5', lw=0.5)
    grat.Insidelabels(
        wcsaxis=0, world=X+2.5, constval=87.75, fmt='%g^\circ',
        color='k', fontsize=10, fun=lambda x: x - 2.5
        # deltapx=1, deltapy=1,
        )
    grat.Insidelabels(
        wcsaxis=1, world=Y+0.06, constval=200, fmt='Dms',
        color='k', fontsize=10,
        )
    for wcsaxis in range(4):
        grat.setp_plotaxis(wcsaxis, visible=False)
        grat.setp_ticklabel(wcsaxis, visible=False)
        grat.setp_tickmark(wcsaxis, visible=False)

    mplim.plot()
    plt.savefig(figname, dpi=75)
    if showplot:
        plt.show()


def main():
    mapcenter = (0., 90.)
    mapsize = (5., 5.)  # degrees
    beamsize_fwhm = 0.1  # degrees
    pixsize = beamsize_fwhm / 4.
    dnaxis1 = int(mapsize[0] / pixsize)
    dnaxis2 = int(mapsize[1] / pixsize)

    # define target grid (via fits header according to WCS convention)
    projection = 'ZEA'
    header = {
        'NAXIS': 3,
        'NAXIS1': dnaxis1,
        'NAXIS2': dnaxis2,
        'NAXIS3': 1,  # need dummy spectral axis
        'CTYPE1': 'GLON-{}'.format(projection),
        'CTYPE2': 'GLAT-{}'.format(projection),
        'CUNIT1': 'deg',
        'CUNIT2': 'deg',
        'CDELT1': -pixsize,
        'CDELT2': pixsize,
        'CRPIX1': dnaxis1 / 2.,
        'CRPIX2': dnaxis2 / 2.,
        'CRVAL1': mapcenter[0],
        'CRVAL2': mapcenter[1],
        }

    my_wcs = wcs.WCS(header, naxis=[1, 2])

    xpix, ypix = np.meshgrid(
        np.linspace(0, dnaxis1, 6 * dnaxis1 + 3),
        np.linspace(0, dnaxis2, 6 * dnaxis2 + 3)
        )

    xcoords, ycoords = my_wcs.all_pix2world(xpix, ypix, 0)
    xcoords, ycoords = xcoords.flatten(), ycoords.flatten()

    signal = astro_signal(xcoords, ycoords, beamsize_fwhm)
    signal = signal[:, np.newaxis]  # need dummy spectral axis

    header = pf.Header(header.items())
    mygridder = cygrid.WcsGrid(header, dbg_messages=True)
    kernelsize_fwhm = beamsize_fwhm / 2.
    kernelsize_sigma = kernelsize_fwhm / 2.35

    # mygridder.set_kernel(
    #     'gauss1d',
    #     kernelsize_sigma,
    #     kernelsize_sigma * 3,
    #     kernelsize_sigma / 2.
    #     )
    mygridder.set_kernel(
        'gauss2d',
        (kernelsize_sigma * 10, kernelsize_sigma, np.radians(45.)),
        kernelsize_sigma * 10 * 3,  # major axis is important here
        kernelsize_sigma / 2.  # minor axis is important here
        )
    mygridder.grid(xcoords, ycoords, signal)
    data = mygridder.get_datacube()

    pf.writeto(
        'cygrid_demo_zea_elliptical.fits',
        header=header, data=data, clobber=True
        )
    do_plot(header, data, 'cygrid_demo_zea_elliptical.pdf', showplot=False)

if __name__ == '__main__':
    main()


