import numpy as np
import astropy.units as u
# from astropy.constants import G
import velocity_tools as VT
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
Omega1 = 1e-12/u.s
# Rc = VT.stream_lines.R_cent(M=Mstar, Omega=Omega1, r0=r0)
# #
# r = np.arange(r0.value, Rc.value, step=-10) * u.au
# angles = np.linspace(0, 2 * np.pi, 100)
# theta1 = VT.stream_lines.stream_line(r, M=Mstar, theta0=theta0_1, 
#     Omega=Omega1, r0=r0)
# dphi1 = VT.stream_lines.dphi(theta1, theta0=theta0_1)
# phi1 = 15 * u.deg
# dphi_f = 10 * u.deg
# # Convert from sphereical into cartesian coordinates
# z1 = r * np.cos(theta1)
# y1 = r * np.sin(theta1) * np.sin(phi1 + dphi1)
# x1 = r * np.cos(theta1) * np.cos(phi1 + dphi1)

# z1_p = z1
# y1_p = r * np.sin(theta1) * np.sin(phi1 + dphi1 + dphi_f)
# x1_p = r * np.cos(theta1) * np.cos(phi1 + dphi1 + dphi_f)

# z1_m = z1
# y1_m = r * np.sin(theta1) * np.sin(phi1 + dphi1 - dphi_f)
# x1_m = r * np.cos(theta1) * np.cos(phi1 + dphi1 - dphi_f)
# #
# inc = 30*u.deg
# PA_ang = -50*u.deg
# inc = 0*u.deg
# PA_ang = 50*u.deg
# #
# x1_new, y1_new, z1_new = VT.stream_lines.rotate_xyz(x1, y1, z1, inc=inc, PA=PA_ang)
# x1_p_new, y1_p_new, z1_p_new = VT.stream_lines.rotate_xyz(x1_p, y1_p, z1_p, inc=inc, PA=PA_ang)
# x1_m_new, y1_m_new, z1_m_new = VT.stream_lines.rotate_xyz(x1_m, y1_m, z1_m, inc=inc, PA=PA_ang)

theta0 = 80 * u.deg
r0 = 1e4 * u.au
Mstar = 3.2 * u.Msun * 20
Omega1 = -1e-12/u.s
inc = 30*u.deg
PA_ang = 50*u.deg
x1, y1, z1 = VT.stream_lines.xyz_stream(M=Mstar, theta0=theta0, phi0=115*u.deg,
    Omega=Omega1, r0=r0, inc=inc, PA=PA_ang)
x2, y2, z2 = VT.stream_lines.xyz_stream(M=Mstar, theta0=theta0, phi0=145*u.deg,
    Omega=Omega1, r0=r0, inc=inc, PA=PA_ang)
x3, y3, z3 = VT.stream_lines.xyz_stream(M=Mstar, theta0=theta0, phi0=175*u.deg,
    Omega=Omega1, r0=r0, inc=inc, PA=PA_ang)

plt.ion()
xrange = np.array([0, 10])
# plt.plot(xrange, xrange)
plt.close()
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot(x1, y1, z1, color='r', marker='o', markersize=3)
ax.plot(x2, y2, z2, color='b', marker='o', markersize=3)
ax.plot(x3, y3, z3, color='k', marker='o', markersize=3)
# # ax.plot(x1_m, y1_m, z1_m, color='r', marker='o', markersize=3)
# # ax.plot(x1_p, y1_p, z1_p, color='r', marker='o', markersize=3)
# # #
# ax.plot(x1_new, y1_new, z1_new, color='b', marker='o', markersize=3)
# ax.plot(x1_p_new, y1_p_new, z1_p_new, color='b', marker='o', markersize=3)
# ax.plot(x1_m_new, y1_m_new, z1_m_new, color='b', marker='o', markersize=3)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

# Stream line into arcsec
dra_stream1 = -x1.value / distance
ddec_stream1 = z1.value / distance

dra_stream2 = -x2.value / distance
ddec_stream2 = z2.value / distance

dra_stream3 = -x3.value / distance
ddec_stream3 = z3.value / distance

# Create Per-emb-2 reference coordinate system
from astropy.coordinates import SkyCoord, FK5
Per2_c = SkyCoord(ra_Per2, dec_Per2, frame='fk5')
Per2_ref = Per2_c.skyoffset_frame()
#
fil_1 = SkyCoord(dra_stream1 * u.arcsec, ddec_stream1 * u.arcsec, frame=Per2_ref).transform_to(FK5)
fil_2 = SkyCoord(dra_stream2 * u.arcsec, ddec_stream2 * u.arcsec, frame=Per2_ref).transform_to(FK5)
fil_3 = SkyCoord(dra_stream3 * u.arcsec, ddec_stream3 * u.arcsec, frame=Per2_ref).transform_to(FK5)
# other_fk5  = other.transform_to(FK5)
# 
hdu = fits.open(HC3N_TdV_10_9)[0]
wcs = WCS(hdu.header)

fig2 = plt.figure()
ax2 = fig2.add_subplot(111, projection=wcs)
ax2.imshow(hdu.data, vmin=0, vmax=160.e-3, origin='lower', cmap='Greys')
ax2.set_autoscale_on(False)
ax2.plot(fil_1.ra, fil_1.dec, transform=ax2.get_transform('fk5'), color='r') 
ax2.plot(fil_2.ra, fil_2.dec, transform=ax2.get_transform('fk5'), color='b') 
ax2.plot(fil_3.ra, fil_3.dec, transform=ax2.get_transform('fk5'), color='orange') 
ax2.plot(ra_Per2, dec_Per2, transform=ax2.get_transform('fk5'), marker='*', color='red') 
# ax.plot(x1, y1, z1, color='r', marker='o', markersize=3)
# ax.plot(x1_m, y1_m, z1_m, color='r', marker='o', markersize=3)
# ax.plot(x1_p, y1_p, z1_p, color='r', marker='o', markersize=3)
# #
# ax.plot(x1_new, y1_new, z1_new, color='b', marker='o', markersize=3)
# ax.plot(x1_p_new, y1_p_new, z1_p_new, color='b', marker='o', markersize=3)
# ax.plot(x1_m_new, y1_m_new, z1_m_new, color='b', marker='o', markersize=3)
# ax.set_xlabel('X')
# ax2.set_xlabel('RA (J2000)')
# ax2.set_ylabel('Dec (J2000)')
