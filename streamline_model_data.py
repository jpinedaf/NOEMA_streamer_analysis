import numpy as np
import astropy.units as u
import velocity_tools.stream_lines as SL
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.io import fits
from config import *

# Main parameters to generate a streamline
Mstar = 3.2*u.Msun
# Inclination angle is not well constrained
#
inc = -43*u.deg
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
# hdu = fits.open(HC3N_TdV_10_9)[0]
# wcs = WCS(hdu.header)
# fig2 = plt.figure(figsize=(4, 4))
# ax2 = fig2.add_subplot(111, projection=wcs)
# ax2.imshow(hdu.data, vmin=0, vmax=160.e-3, origin='lower', cmap='Greys')
# ax2.set_autoscale_on(False)
# ax2.plot(ra_Per2, dec_Per2, transform=ax2.get_transform('fk5'), marker='*',
#          color='red')
# ax2.plot(ra_poly, dec_poly, transform=ax2.get_transform('fk5'), color='black')
# hdu = fits.open(HC3N_TdV_10_9)[0]
# wcs = WCS(hdu.header)
import aplpy
fig2 = aplpy.FITSFigure(HC3N_TdV_10_9, figsize=(4, 4))
# ax2 = fig2.add_subplot(111, projection=wcs)
fig2.show_grayscale(vmin=0, vmax=160.e-3, invert=True)#, cmap='Greys')
# ax2.set_autoscale_on(False)
# ax2.plot(ra_Per2, dec_Per2, transform=ax2.get_transform('fk5'), marker='*',
#          color='red')
# ax2.plot(ra_poly, dec_poly, transform=ax2.get_transform('fk5'), color='black')
# fig2.show_polygons(np.array([ra_poly, dec_poly]))
fig2.show_regions(region_file)
fig2.show_markers(ra_Per2, dec_Per2, marker='*')
fig2.add_beam()

#
fig3 = plt.figure(figsize=(4, 4))
ax3 = fig3.add_subplot(111)
# axes labels
# ax2.set_xlabel('Right Ascension (J2000)')
# ax2.set_ylabel('Declination (J2000)')
ax3.set_xlabel('Projected distance (au)')
ax3.set_ylabel(r"V$_{lsr}$ (km s$^{-1}$)")
#
r_proj, v_los = per_emb_2_get_vc_r()
from scipy import stats
xmin=0; xmax=12000; ymin=6.3; ymax=7.7
xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
positions = np.vstack([xx.ravel(), yy.ravel()])
#
gd_vlos = np.isfinite(r_proj*v_los)
values = np.vstack([r_proj[gd_vlos].value, v_los[gd_vlos].value])
kernel = stats.gaussian_kde(values)
zz = np.reshape(kernel(positions).T, xx.shape)
zz /= zz.max()
ax3.contourf(xx, yy, zz, cmap='Greys', levels=np.arange(0.1, 1.2, 0.1), vmin=0., vmax=1.1)
#
# Define the different initial theta and phi to generate the stream lines
#
theta0 = 130.*u.deg
r0 = 0.9e4*u.au
phi0 = 365.*u.deg
v_r0 = 0*u.km/u.s
omega0 = 4e-13/u.s
# v_lsr = 6.98*u.km/u.s
v_lsr = 7.05*u.km/u.s
#
(x1, y1, z1), (vx1, vy1, vz1) = SL.xyz_stream(
            mass=Mstar, r0=r0, theta0=theta0, phi0=phi0,
            omega=omega0, v_r0=v_r0, inc=inc, pa=PA_ang, rmin=5.5e3*u.au)
d_sky_au = np.sqrt(x1**2 + z1**2)
# gd_2000 = (d_sky_au > 2.4e3*u.au)
# Stream line into arcsec
dra_stream = -x1.value / distance
ddec_stream = z1.value / distance
fil = SkyCoord(dra_stream*u.arcsec, ddec_stream*u.arcsec,
                 frame=Per2_ref).transform_to(FK5)
fig2.show_polygons([fil.ra.value*u.deg, fil.dec.value*u.deg], ls='-', lw=2, color='red')
ax3.plot(d_sky_au, v_lsr + vy1, color='red')
ax3.xaxis.set_ticks(np.arange(3e3, 12e3, 2e3))
ax3.yaxis.set_ticks(np.arange(6.9, 7.6, 0.2))
ax3.set_ylim(6.9, 7.6)
ax3.set_xlim(2.0e3, 11e3)
# save files
fig2.savefig('figures/Per-emb-2_HC3N_10-9_TdV_streamline.pdf', bbox_inches='tight')
fig3.savefig('figures/Per-emb-2_HC3N_10-9_Vlsr_streamline.pdf', bbox_inches='tight')
