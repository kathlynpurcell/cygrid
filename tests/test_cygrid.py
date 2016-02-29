import numpy as np
import pyfits as pf
import cygrid, cygrid2
import matplotlib.pyplot as plt

def gauss2D(x, y, x0, y0, w, A):
    return A * np.exp(-((x-x0)**2 + (y-y0)**2) / 2. / w**2)

N = 50000
in_center_x, in_center_y = 180., 30.
in_width_x, in_width_y = 10., 10.
in_lons = np.random.normal(in_center_x, in_width_x, N)
in_lats = np.random.normal(in_center_y, in_width_y, N)
in_data = np.random.normal(0., 1., N).astype(np.float32)[:, np.newaxis]
in_weights = np.ones_like(in_data)

in_data += gauss2D(in_lons, in_lats, in_center_x, in_center_y, 2., 10.)[:, np.newaxis]

naxis = 200
hdr = pf.Header()
hdr.update({
    'NAXIS': 3,
    'NAXIS1': naxis,
    'NAXIS2': naxis,
    'NAXIS3': 1,
    'CRVAL1': in_center_x,
    'CRVAL2': in_center_y,
    'CRVAL3': 0.,
    'CRPIX1': (naxis+1)/2.,
    'CRPIX2': (naxis+1)/2.,
    'CRPIX3': 1.,
    'CDELT1': -in_width_x / naxis,
    'CDELT2': in_width_y / naxis,
    'CDELT3': 1.,
    'CTYPE1': 'RA---SIN',
    'CTYPE2': 'DEC--SIN',
    'CTYPE3': 'VELO-LSR',
    })




kernelsizefwhm = 30. / 60.
kernelsigma = np.float32(kernelsizefwhm / np.sqrt(8.*np.log(2.)))
kernelradius = np.float32(5. * kernelsigma)
kernelsigmasq = np.float32(kernelsigma**2)


pygridder = cygrid.pygrid(hdr)
# %timeit -n 3 -r 3 pygridder.grid(in_lons, in_lats, in_data, in_weights, kernelradius, kernelsigmasq)
# N = 50000, in_width_x, in_width_y = 10., 10., naxis = 200: 3 loops, best of 3: 58.4 s per loop
pygridder.grid(in_lons, in_lats, in_data, in_weights, kernelradius, kernelsigmasq)
# number of target pixels used 39601
# number of input-output pixel combinations 2000000000
# number of good input-output pixel combinations 11972441
# average number of input pixels per output pixel 50503.7751572
# average good number of input pixels per output pixel 302.326734173
datacube = pygridder.getDatacube()

pygridder2 = cygrid2.pygrid(hdr)
# %timeit -n 2 -r 2 pygridder2.grid(in_lons, in_lats, in_data, in_weights, kernelradius, kernelsigmasq)
# %timeit -n 3 -r 3 pygridder2.grid(in_lons, in_lats, in_data, in_weights, kernelradius, kernelsigmasq)
# N = 50000, in_width_x, in_width_y = 10., 10., naxis = 200: 3 loops, best of 3: 6.3 s per loop
pygridder2.grid(in_lons, in_lats, in_data, in_weights, kernelradius, kernelsigmasq)
# number of target pixels used 40000
# number of input-output pixel combinations 11972909
# number of good input-output pixel combinations 11815510
# average number of input pixels per output pixel 299.322725
# average good number of input pixels per output pixel 295.38775
datacube2 = pygridder2.getDatacube()

print np.std(datacube[0]-datacube2[0])
print np.max(datacube[0]-datacube2[0])
print np.min(datacube[0]-datacube2[0])

plt.close()
plt.imshow(datacube[0], interpolation='nearest')
plt.show()

plt.close()
plt.imshow(datacube2[0], interpolation='nearest')
plt.show()

plt.close()
plt.imshow(datacube[0]-datacube2[0], interpolation='nearest')
plt.show()











import numpy as np
import pyfits as pf
import cygrid, cygrid2
import matplotlib.pyplot as plt

def gauss2D(x, y, x0, y0, w, A):
    return A * np.exp(-((x-x0)**2 + (y-y0)**2) / 2. / w**2)


N = 200000
in_center_x, in_center_y = 180., 30.
in_width_x, in_width_y = 1., 1.
in_lons = np.random.normal(in_center_x, in_width_x, N)
in_lats = np.random.normal(in_center_y, in_width_y, N)
in_data = np.random.normal(0., 1., N).astype(np.float32)[:, np.newaxis]
in_weights = np.ones_like(in_data)

in_data += gauss2D(in_lons, in_lats, in_center_x, in_center_y, 0.1, 10.)[:, np.newaxis]

in_lons = np.array([in_center_x+0.5])
in_lats = np.array([in_center_y])
in_data = np.array([1.]).astype(np.float32)[:, np.newaxis]
in_weights = np.ones_like(in_data)

naxis = 200
hdr = pf.Header()
hdr.update({
    'NAXIS': 3,
    'NAXIS1': naxis,
    'NAXIS2': naxis,
    'NAXIS3': 1,
    'CRVAL1': in_center_x,
    'CRVAL2': in_center_y,
    'CRVAL3': 0.,
    'CRPIX1': (naxis+1)/2.,
    'CRPIX2': (naxis+1)/2.,
    'CRPIX3': 1.,
    'CDELT1': -in_width_x / naxis,
    'CDELT2': in_width_y / naxis,
    'CDELT3': 1.,
    'CTYPE1': 'RA---SIN',
    'CTYPE2': 'DEC--SIN',
    'CTYPE3': 'VELO-LSR',
    })




kernelsizefwhm = hdr['CDELT2'] * 2.
kernelsigma = np.float32(kernelsizefwhm / np.sqrt(8.*np.log(2.)))
kernelradius = np.float32(5. * kernelsigma)
kernelsigmasq = np.float32(kernelsigma**2)


pygridder = cygrid.pygrid(hdr)
# %timeit -n 3 -r 3 pygridder.grid(in_lons, in_lats, in_data, in_weights, kernelradius, kernelsigmasq)
# N = 50000, in_width_x, in_width_y = 10., 10., naxis = 200: 3 loops, best of 3: 58.4 s per loop
pygridder.grid(in_lons, in_lats, in_data, in_weights, kernelradius, kernelsigmasq)
#datacube = pygridder.getDatacube()


pygridder2 = cygrid2.pygrid(hdr)
# %timeit -n 3 -r 3 pygridder2.grid(in_lons, in_lats, in_data, in_weights, kernelradius, kernelsigmasq)
pygridder2.grid(in_lons, in_lats, in_data, in_weights, kernelradius, kernelsigmasq)
#datacube2 = pygridder2.getDatacube()


print np.std(datacube[0]-datacube2[0])
print np.nanmax(datacube[0]-datacube2[0])
print np.nanmin(datacube[0]-datacube2[0])

plt.close()
plt.imshow(datacube[0]-datacube2[0], interpolation='nearest')
plt.show()
