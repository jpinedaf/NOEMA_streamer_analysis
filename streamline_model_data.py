import numpy as np
import astropy.units as u
import velocity_tools.stream_lines as SL
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.io import fits
from config import *

# Main parameters to generate a streamline
# r0 = 1.2e4*u.au
Mstar = 3.2*u.Msun
# Omega1 = 1e-13/u.s
# Inclination angle is not well constrained
#
inc = -43*u.deg
PA_ang = -(180. - 130)*u.deg
PA_ang = 130*u.deg

# Create Per-emb-2 reference coordinate system
from astropy.coordinates import SkyCoord, FK5
Per2_c = SkyCoord(ra_Per2, dec_Per2, frame='fk5')
Per2_ref = Per2_c.skyoffset_frame()
#
# First we setup the plotting windows: 3d plot and one with wcs coordinates
#
plt.ion()
# get data and plot background image for streamer
hdu = fits.open(HC3N_TdV_10_9)[0]
wcs = WCS(hdu.header)
fig2 = plt.figure()
ax2 = fig2.add_subplot(111, projection=wcs)
ax2.imshow(hdu.data, vmin=0, vmax=160.e-3, origin='lower', cmap='Greys')
ax2.set_autoscale_on(False)
ax2.plot(ra_Per2, dec_Per2, transform=ax2.get_transform('fk5'), marker='*',
         color='red')
ax2.plot(ra_poly, dec_poly, transform=ax2.get_transform('fk5'), color='red')
#
fig3 = plt.figure()
ax3 = fig3.add_subplot(111)
# axes labels
ax2.set_xlabel('Right Ascension (J2000)')
ax2.set_ylabel('Declination (J2000)')
ax3.set_xlabel('Projected distance (au)')
ax3.set_ylabel(r"V$_{lsr}$ (km s$^{-1}$)")

r_proj, v_los = per_emb_2_get_vc_r()
ax3.scatter(r_proj, v_los, marker='*', color='k')
#
# Define the different initial theta and phi to generate the stream lines
#
theta0 = 130.*u.deg
r0 = 0.9e4*u.au
phi0 = 365.*u.deg
v_r0 = 0*u.km/u.s
omega0 = 4e-13/u.s
v_lsr = 6.98*u.km/u.s
#
(x1, y1, z1), (vx1, vy1, vz1) = SL.xyz_stream(
            mass=Mstar, r0=r0, theta0=theta0, phi0=phi0,
            omega=omega0, v_r0=v_r0, inc=inc, pa=PA_ang)
d_sky_au = np.sqrt(x1**2 + z1**2)
gd_2000 = (d_sky_au > 2.4e3*u.au)
# Stream line into arcsec
dra_stream = -x1.value / distance
ddec_stream = z1.value / distance
fil = SkyCoord(dra_stream*u.arcsec, ddec_stream*u.arcsec,
                 frame=Per2_ref).transform_to(FK5)
ax2.plot(fil.ra, fil.dec, transform=ax2.get_transform('fk5'),
         ls='-', lw=2, label=my_label)
ax3.scatter(d_sky_au[gd_2000], v_lsr + vy1[gd_2000], marker='o')
ax3.set_ylim(6.9, 7.6)