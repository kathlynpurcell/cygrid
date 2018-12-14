from reproject import reproject_interp


ebhis_data_reprojected, footprint = reproject_interp(
    (ebhis_data, ebhis_header), planck_header
    )