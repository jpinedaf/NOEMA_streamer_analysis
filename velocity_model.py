import numpy as np
import astropy.units as u
# from astropy.constants import G
import velocity_tools.stream_lines as SL
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from astropy.wcs import WCS
from astropy.io import fits
from config import *

# Generate a srteamline
theta0_1 = 80 * u.deg

r0 = 1e4 * u.au
Mstar = 3.2 * u.Msun * 20
# Omega1 = (np.sqrt(G * Mstar * 300*u.au) / r0**2).decompose() #  6.16883729e-14/u.s
# Omega1 = 1e-13/u.s
# Omega1 = 1e-12/u.s

# theta0 = 80 * u.deg
r0 = 1e4 * u.au
Mstar = 3.2 * u.Msun * 20
Omega1 = 1e-12/u.s
# Inclination angle is not well constrained
inc = 45*u.deg
PA_ang = 130*u.deg

# Create Per-emb-2 reference coordinate system
from astropy.coordinates import SkyCoord, FK5
Per2_c = SkyCoord(ra_Per2, dec_Per2, frame='fk5')
Per2_ref = Per2_c.skyoffset_frame()

plt.ion()
xrange = np.array([0, 10])
# configure plotting windows
plt.close()
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

hdu = fits.open(HC3N_TdV_10_9)[0]
wcs = WCS(hdu.header)
fig2 = plt.figure()
ax2 = fig2.add_subplot(111, projection=wcs)
ax2.imshow(hdu.data, vmin=0, vmax=160.e-3, origin='lower', cmap='Greys')
ax2.set_autoscale_on(False)
ax2.plot(ra_Per2, dec_Per2, transform=ax2.get_transform('fk5'), marker='*', color='red') 

fig3 = plt.figure()
ax3 = fig3.add_subplot(111)

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
# 
# theta0 = 80 * u.deg
theta0_list = np.array([100, 120.])*u.deg
phi0_list = np.linspace( 100., 180, num=3)*u.deg
inc=0.*u.deg
for theta0 in theta0_list:
    for phi0 in phi0_list:
        x1, y1, z1 = SL.xyz_stream(M=Mstar, theta0=theta0, phi0=phi0,
                    Omega=-Omega1, r0=r0, inc=inc, PA=PA_ang)
        ax.plot(x1, y1, z1, marker='o', markersize=3)
        # Stream line into arcsec
        dra_stream = -x1.value / distance
        ddec_stream = z1.value / distance
        fil = SkyCoord(dra_stream * u.arcsec, ddec_stream * u.arcsec, 
                         frame=Per2_ref).transform_to(FK5)
        ax2.plot(fil.ra, fil.dec, transform=ax2.get_transform('fk5'))
    x1, y1, z1 = SL.xyz_stream(M=Mstar, theta0=theta0, phi0=phi0,
                    Omega=-Omega1, r0=r0, inc=0*u.deg, PA=0*u.deg)
    ax3.plot(x1, z1, marker='o', markersize=3)