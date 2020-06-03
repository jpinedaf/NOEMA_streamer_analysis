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
import aplpy
fig2 = aplpy.FITSFigure(HC3N_TdV_10_9, figsize=(4, 4))
fig2.show_grayscale(vmin=0, vmax=160.e-3, invert=True)
# setup and colorbar
setup_plot_noema(fig2, label_col='black', star_col='yellow')

fig2.show_regions(region_file)
#
fig3 = plt.figure(figsize=(5, 4))
ax3 = fig3.add_subplot(111)
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
v_lsr = 7.05*u.km/u.s
#
(x1, y1, z1), (vx1, vy1, vz1) = SL.xyz_stream(
            mass=Mstar, r0=r0, theta0=theta0, phi0=phi0,
            omega=omega0, v_r0=v_r0, inc=inc, pa=PA_ang, rmin=5.5e3*u.au)
d_sky_au = np.sqrt(x1**2 + z1**2)
# Stream line into arcsec
dra_stream = -x1.value / distance
ddec_stream = z1.value / distance
fil = SkyCoord(dra_stream*u.arcsec, ddec_stream*u.arcsec,
               frame=Per2_ref).transform_to(FK5)
#
fig2.show_markers(fil.ra.value*u.deg, fil.dec.value*u.deg,
                  marker='o', color='red', s=3)
fig2.add_label(0.75, 0.9, r"HC$_3$N ($10-9$)", color='black', relative=True,
               size=14, weight=60)
freq_HC5N_10_9 = fits.getheader(HC3N_TdV_10_9)['RESTFREQ'] * u.Hz
fig2.show_circles(ra_Per2, dec_Per2, pb_noema(freq_HC5N_10_9).to(u.deg)*0.5,
                  ls=':', color='black')
ax3.plot(d_sky_au, v_lsr + vy1, color='red')
ax3.xaxis.set_ticks(np.arange(3e3, 12e3, 2e3))
ax3.yaxis.set_ticks(np.arange(6.9, 7.6, 0.2))
ax3.set_ylim(6.9, 7.55)
ax3.set_xlim(2.0e3, 9e3)
ax3.text(2400, 7.5, r"HC$_3$N ($10-9$)")
ax3.text(2400, 7.45, r"Streamline model", color='red')
# save files
fig2.add_colorbar(axis_label_text=r'Integrated Intensity (Jy beam$^{-1}$ km s$^{-1}$)',
                  ticks=[0, 0.05, 0.1, 0.15])
fig2.colorbar.hide()

import pickle
stream_model = {'ra': fil.ra.value*u.deg, 'dec': fil.dec.value*u.deg,
                'd_sky_au': d_sky_au, 'vlsr': v_lsr + vy1}
with open(stream_pickle, 'wb') as f:
    pickle.dump(stream_model, f)

KDE_vel_rad = {'radius': xx, 'v_lsr': yy, 'dens': zz}
with open(vlsr_rad_kde_pickle, 'wb') as f:
    pickle.dump(KDE_vel_rad, f)

