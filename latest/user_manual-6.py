from astropy.wcs import WCS
import matplotlib.pyplot as plt


planck_wcs, ebhis_wcs = WCS(planck_header), WCS(ebhis_header)

fig = plt.figure(figsize=(12, 6))
ax1 = fig.add_subplot(1,2,1, projection=planck_wcs)
ax1.imshow(planck_data, origin='lower')
ax1.set_title('Planck 857 GHz')
ax1.coords['ra'].set_axislabel('Right Ascension')
ax1.coords['dec'].set_axislabel('Declination')

ax2 = fig.add_subplot(1,2,2, projection=ebhis_wcs)
ax2.imshow(ebhis_data, origin='lower')
ax2.set_title('EBHIS 21-cm HI line')
ax2.coords['glon'].set_axislabel('Galactic Longitude')
ax2.coords['glat'].set_axislabel('Galactic Latitude')
ax2.coords['glat'].set_axislabel_position('r')
ax2.coords['glat'].set_ticklabel_position('r')