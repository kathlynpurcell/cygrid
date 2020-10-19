fig.clear()
ax1 = fig.add_subplot(1,2,1, projection=planck_wcs)
ax1.imshow(planck_data, origin='lower')
ax1.set_title('Planck 857 GHz')

ax2 = fig.add_subplot(1,2,2, projection=planck_wcs)
ax2.imshow(ebhis_data_regridded, origin='lower')
ax2.set_title('EBHIS 21-cm HI line')

for ax in [ax1, ax2]:
    ax.coords['ra'].set_axislabel('Right Ascension')
    ax.coords['dec'].set_axislabel('Declination')

ax2.coords['ra'].set_axislabel_position('r')
ax2.coords['dec'].set_ticklabel_position('r')