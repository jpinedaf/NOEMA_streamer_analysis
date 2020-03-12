import urllib.request
from spectral_cube import SpectralCube
import astropy.units as u
from config import file_12CO_blue, file_12CO_red
import os

# import file from Dataverse
url = 'https://dataverse.harvard.edu/api/access/datafile/:persistentId?persistentId=doi:10.7910/DVN/3ZYWRC/KUVRKD'
tmp_file = 'test.fits'
urllib.request.urlretrieve(url, tmp_file)
# Load and create integrated intensity maps
cube = SpectralCube.read(tmp_file).with_spectral_unit(u.km/u.s)
cube_blue = cube.spectral_slab(-1.5 * u.km / u.s, 5 * u.km / u.s)
cube_red = cube.spectral_slab(9 * u.km / u.s, 16 * u.km / u.s)
TdV_red = cube_red.moment(order=0)
TdV_blue = cube_blue.moment(order=0)
# write out FITS files
TdV_red.write(file_12CO_red, overwrite=True)
TdV_blue.write(file_12CO_blue, overwrite=True)
os.remove(tmp_file)