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
ax2.plot(ra_poly, dec_poly, transform=ax2.get_transform('fk5'), color='red')
#
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

r_proj, v_los = per_emb_2_get_vc_r()
ax3.scatter(r_proj, v_los, marker='*', color='k')
#
# Define the different initial theta and phi to generate the stream lines
#
# theta0_list = np.array([100, 120.])*u.deg
theta0_list = [130., 132.]*u.deg
theta0 = 130.*u.deg
#r0 = 1.0e4*u.au
r0 = 0.9e4*u.au
# phi0_list = [190.]*u.deg
phi0_list = [365., 370.]*u.deg
phi0 = 365.*u.deg
# phi0 = 200.*u.deg
#v_r0_list = np.linspace(0.0, 0.5, num=3)*u.km/u.s
v_r0_list = [0.0]*u.km/u.s
v_r0 = 0*u.km/u.s
omega_list = [4e-13, 1e-14]/u.s
omega0 = 4e-13/u.s
#lw_i = 2
#ls_i = '-'
v_lsr = 7.05*u.km/u.s
#
def stream_label(v_r=None, omega=None, theta=None, phi=None):
    my_label = ""
    if v_r is not None:
        my_label = r"{0} $V_r=${1}".format(my_label, v_r)
    if omega is not None:
        my_label = r"{0} $\Omega=${1}".format(my_label, omega)
    if theta is not None:
        my_label = r"{0} $\theta_0=${1}".format(my_label, theta)
    if phi is not None:
        my_label = r"{0} $\phi_0=${1}".format(my_label, phi)
    return my_label

#for v_r0 in v_r0_list:
for omega0 in omega_list:
        # for r0, theta0, phi0, ls_i in zip(r0_list, theta0_list, phi0_list, ls_list):
        # for phi0 in phi0_list:
 #       for theta0, phi0 in zip(theta0_list, phi0_list):
    my_label = stream_label(v_r=v_r0, omega=omega0,
                            theta=theta0, phi=phi0)
    (x1, y1, z1), (vx1, vy1, vz1) = SL.xyz_stream(
                mass=Mstar, r0=r0, theta0=theta0, phi0=phi0,
                omega=omega0, v_r0=v_r0, inc=inc, pa=PA_ang)
    d_sky_au = np.sqrt(x1**2 + z1**2)
    gd_2000 = (d_sky_au > 2.4e3*u.au)
    ax.plot(x1, y1, z1, marker='o', markersize=3)
    # Stream line into arcsec
    dra_stream = -x1.value / distance
    ddec_stream = z1.value / distance
    fil = SkyCoord(dra_stream*u.arcsec, ddec_stream*u.arcsec,
                   frame=Per2_ref).transform_to(FK5)
    ax2.plot(fil.ra, fil.dec, transform=ax2.get_transform('fk5'),
             ls='-', lw=2, label=my_label)
    ax3.scatter(d_sky_au[gd_2000], v_lsr + vy1[gd_2000], marker='o')

# Add axes to see rotations
x_b = np.array([1, 0, 0])*8e3/distance
y_b = np.array([0, 0, 0])*8e3/distance  # This is the axes origin
z_b = np.array([0, 0, 1])*8e3/distance
nx_b, ny_b, nz_b = SL.rotate_xyz(x_b, y_b, z_b, inc=inc, pa=PA_ang)

# original axes
my_axis = SkyCoord(-x_b*u.arcsec, z_b*u.arcsec,
                   frame=Per2_ref).transform_to(FK5)
ax2.plot(my_axis.ra, my_axis.dec, transform=ax2.get_transform('fk5'),
         color='k')
# new axes
my_axis_new = SkyCoord(-nx_b*u.arcsec, nz_b*u.arcsec,
                       frame=Per2_ref).transform_to(FK5)
ax2.plot(my_axis_new.ra[0:2], my_axis_new.dec[0:2],
         transform=ax2.get_transform('fk5'),
         color='red')
if ny_b[-1] > 0:
    new_ax_color = 'red'
else:
    new_ax_color = 'blue'
ax2.plot(my_axis_new.ra[1:], my_axis_new.dec[1:],
         transform=ax2.get_transform('fk5'),
         color=new_ax_color)
# plot legend at the end
ax2.legend()