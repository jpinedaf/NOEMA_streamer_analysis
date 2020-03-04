import matplotlib.pyplot as plt
import velocity_tools.stream_lines as SL
import astropy.units as u
import numpy as np
plt.ion()

inc = -43*u.deg
PA_ang = 130*u.deg

x_b = np.array([0, 0])
y_b = np.array([0, 0])
z_b = np.array([0, 1])
nx_b, ny_b, nz_b = SL.rotate_xyz(x_b, y_b, z_b, inc=inc, pa=PA_ang)

x_r = np.array([0, 0])
y_r = np.array([0, 0])
z_r = np.array([0, -1])
nx_r, ny_r, nz_r = SL.rotate_xyz(x_r, y_r, z_r, inc=inc, pa=PA_ang)

if ny_r[1] > 0:
    plt.plot(nx_r, nz_r, color='r')
else:
    plt.plot(nx_r, nz_r, color='b')

if ny_b[1] > 0:
    plt.plot(nx_b, nz_b, color='r', ls='--')
else:
    plt.plot(nx_b, nz_b, color='b', ls='--')
plt.axis('equal')