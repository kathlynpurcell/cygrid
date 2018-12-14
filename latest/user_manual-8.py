import cygrid


# need to add dummy axis again
ebhis_target_header = planck_header.copy()
ebhis_target_header['NAXIS'] = 3
ebhis_target_header['NAXIS3'] = 1

gridder = cygrid.WcsGrid(ebhis_target_header)

# EBHIS resolution is 10.8', but image is strongly oversampled, so we use
# a kernel of 0.012 (twice the pixel size)
kernelsize_fwhm = 0.012
kernelsize_sigma = kernelsize_fwhm / np.log(8 * np.sqrt(2))
sphere_radius = 4. * kernelsize_sigma

gridder.set_kernel(
    'gauss1d',
    (kernelsize_sigma,),
    sphere_radius,
    kernelsize_sigma / 2.
    )

gridder.grid(
    ebhis_coords_eq.ra.value.flatten(),
    ebhis_coords_eq.dec.value.flatten(),
    ebhis_data.reshape(-1, 1)
    )
ebhis_data_regridded = gridder.get_datacube().squeeze()