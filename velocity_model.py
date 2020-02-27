import numpy as np
import astropy.units as u
from astropy.constants import G
import velocity_tools as VT
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Generate a srteamline
theta0_1 = 70*u.deg

r0 = 1e4*u.au
Mstar = 3.2*u.Msun
# Omega1 = (np.sqrt(G * Mstar * 300*u.au) / r0**2).decompose() #  6.16883729e-14/u.s
Omega1 = 1e-13/u.s
Rc = VT.stream_lines.R_cent(M=Mstar, Omega=Omega1, r0=r0)
#
r = np.arange(r0.value, Rc.value, step=-10) * u.au
angles = np.linspace(0, 2 * np.pi, 100)
theta1 = VT.stream_lines.stream_line(r, M=Mstar, theta0=theta0_1, 
    Omega=Omega1, r0=r0)
dphi1 = VT.stream_lines.dphi(theta1, theta0=theta0_1)
phi1 = 15 * u.deg
dphi_f = 10 * u.deg
# Convert from sphereical into cartesian coordinates
z1 = r * np.cos(theta1)
y1 = r * np.sin(theta1) * np.sin(phi1 + dphi1)
x1 = r * np.cos(theta1) * np.cos(phi1 + dphi1)

z1_p = z1
y1_p = r * np.sin(theta1) * np.sin(phi1 + dphi1 + dphi_f)
x1_p = r * np.cos(theta1) * np.cos(phi1 + dphi1 + dphi_f)

z1_m = z1
y1_m = r * np.sin(theta1) * np.sin(phi1 + dphi1 - dphi_f)
x1_m = r * np.cos(theta1) * np.cos(phi1 + dphi1 - dphi_f)
#
inc = 30*u.deg
PA_ang = -50*u.deg
inc = 0*u.deg
PA_ang = 50*u.deg
#
x1_new, y1_new, z1_new = VT.stream_lines.rotate_xyz(x1, y1, z1, inc=inc, PA=PA_ang)
x1_p_new, y1_p_new, z1_p_new = VT.stream_lines.rotate_xyz(x1_p, y1_p, z1_p, inc=inc, PA=PA_ang)
x1_m_new, y1_m_new, z1_m_new = VT.stream_lines.rotate_xyz(x1_m, y1_m, z1_m, inc=inc, PA=PA_ang)
plt.ion()
xrange = np.array([0, 10])
# plt.plot(xrange, xrange)
plt.close()
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot(x1, y1, z1, color='r', marker='o', markersize=3)
ax.plot(x1_m, y1_m, z1_m, color='r', marker='o', markersize=3)
ax.plot(x1_p, y1_p, z1_p, color='r', marker='o', markersize=3)
#
ax.plot(x1_new, y1_new, z1_new, color='b', marker='o', markersize=3)
ax.plot(x1_p_new, y1_p_new, z1_p_new, color='b', marker='o', markersize=3)
ax.plot(x1_m_new, y1_m_new, z1_m_new, color='b', marker='o', markersize=3)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
