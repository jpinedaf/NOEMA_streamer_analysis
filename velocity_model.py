import numpy as np
import astropy.units as u
import velocity_tools.stream_lines as SL
import matplotlib.pyplot as plt
# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
from astropy.wcs import WCS
from astropy.io import fits
from config import *

# Main parameters to generate a streamline
# r0 = 1.2e4*u.au
Mstar = 3.2*u.Msun
# Omega1 = 1e-13/u.s
# Inclination angle is not well constrained
#
inc = 43*u.deg
PA_ang = 130*u.deg

# Create Per-emb-2 reference coordinate system
from astropy.coordinates import SkyCoord, FK5
Per2_c = SkyCoord(ra_Per2, dec_Per2, frame='fk5')
Per2_ref = Per2_c.skyoffset_frame()
#
# First we setup the plotting windows: 3d plot and one with wcs coordinates
#
plt.ion()
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
# get data and plot background image for streamer
hdu = fits.open(HC3N_TdV_10_9)[0]
wcs = WCS(hdu.header)
fig2 = plt.figure()
ax2 = fig2.add_subplot(111, projection=wcs)
ax2.imshow(hdu.data, vmin=0, vmax=160.e-3, origin='lower', cmap='Greys')
ax2.set_autoscale_on(False)
ax2.plot(ra_Per2, dec_Per2, transform=ax2.get_transform('fk5'), marker='*',
         color='red')
fig3 = plt.figure()
ax3 = fig3.add_subplot(111)
# axes labels
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax2.set_xlabel('Right Ascension (J2000)')
ax2.set_ylabel('Declination (J2000)')
ax3.set_xlabel('Projected distance (au)')
ax3.set_ylabel(r"V$_{lsr}$ (km s$^{-1}$)")
#
# Define the different initial theta and phi to generate the stream lines
#
# theta0_list = np.array([100, 120.])*u.deg
theta0_list = [115.]*u.deg
theta0 = 115.*u.deg
r0_list = [1.6e4]*u.au
r0 = 1.6e4*u.au
# phi0_list = np.linspace(160., 200, num=3)*u.deg
phi0_list = [190.]*u.deg
# phi0 = 200.*u.deg
v_r0_list = np.linspace(0.025, 0.05, num=2)*u.km/u.s
# omega_list = np.logspace(-13, -12, num=2, base=10.)/u.s
omega_list = [1e-13]/u.s
# theta0 = 120.*u.deg
# phi0 = 185.*u.deg
lw_i = 3
ls_list = ['-', '--', ':']
ls_i = '-'
for v_r0 in v_r0_list:
    for omega0 in omega_list:
        # for r0, theta0, phi0, ls_i in zip(r0_list, theta0_list, phi0_list, ls_list):
        for phi0 in phi0_list:
            (x1, y1, z1), (vx1, vy1, vz1) = SL.xyz_stream(
                        mass=Mstar, r0=r0, theta0=theta0, phi0=phi0,
                        omega=omega0, v_r0=v_r0, inc=inc, pa=PA_ang)
            ax.plot(x1, y1, z1, marker='o', markersize=3)
            # Stream line into arcsec
            dra_stream = -x1.value / distance
            ddec_stream = z1.value / distance
            fil = SkyCoord(dra_stream*u.arcsec, ddec_stream*u.arcsec,
                             frame=Per2_ref).transform_to(FK5)
            ax2.plot(fil.ra, fil.dec, transform=ax2.get_transform('fk5'),
                     ls=ls_i, lw=lw_i,
                     label=r"$V_r=${0} $\Omega=${1} $\theta_0=${2} $\phi_0=${3}".format(v_r0, omega0, theta0, phi0))
            ax3.scatter(np.sqrt(x1**2 + z1**2), 6.9*u.km/u.s + vy1, marker='o')
    lw_i -= 1
# plot legend at the end
ax2.legend()