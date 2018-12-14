import numpy as np
from astropy.coordinates import SkyCoord


ebhis_glon_pix, ebhis_glat_pix = np.meshgrid(
    np.arange(ebhis_header['NAXIS1']),
    np.arange(ebhis_header['NAXIS2'])
    )
ebhis_lon_world, ebhis_lat_world = ebhis_wcs.all_pix2world(
    ebhis_glon_pix, ebhis_glat_pix, 0
    )
ebhis_coords_gal = SkyCoord(
    ebhis_lon_world, ebhis_lat_world, frame='galactic', unit='deg'
    )
ebhis_coords_eq = ebhis_coords_gal.icrs