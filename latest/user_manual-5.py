from astropy import units as u
from astropy.utils.data import conf
from astroquery.skyview import SkyView


# Loading data from SkyView can take a while, hence the longer timeout
conf.remote_timeout = 60.

kwargs = dict(
    radius=3 * u.deg,
    pixels='500',
    scaling='Linear',
    )
paths = SkyView().get_images(
    position='00 42 44.330 +41 16 07.50',
    coordinates='J2000', survey=['Planck 857'], **kwargs
    )  # doctest: +ELLIPSIS
planck_header, planck_data = paths[0][0].header, paths[0][0].data

paths = SkyView().get_images(
    position='121.174322 -21.573311',
    coordinates='Galactic', survey=['EBHIS'], **kwargs
    )  # doctest: +ELLIPSIS
ebhis_header, ebhis_data = paths[0][0].header, paths[0][0].data