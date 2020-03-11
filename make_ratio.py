from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import aplpy
from config import *

def get_pb_corr(radius, pb_fwhm):
    pb_sigma = pb_fwhm / 2.355
    return np.exp(-radius * radius/2/pb_sigma**2)

# data8, hd8 = fits.getdata(HC3N_TdV_8_7, header=True)
# data10, hd10 = fits.getdata(HC3N_TdV_10_9_match, header=True)
tdv8, N8 = col_dens_hc3n(HC3N_TdV_8_7, nup=8, get_tdv=True)
N10 = col_dens_hc3n(HC3N_TdV_10_9_match, nup=10)

Ratio = N10 / N8
mask = (tdv8 > 0.45*u.K)
Ratio[~mask] = np.nan

hd8 = fits.getheader(HC3N_TdV_8_7)
hd10 = fits.getheader(HC3N_TdV_10_9_match)
hdu_R = fits.PrimaryHDU(Ratio, hd8)

Pop_Up = 1.963E-02

import velocity_tools.coordinate_offsets as c_offset
results = c_offset.generate_offsets(hd8, ra_Per2, dec_Per2, pa_angle=0*u.deg, inclination=0*u.deg)

PB_10 = get_pb_corr(results.r, pb_noema(fits.getheader(HC3N_TdV_10_9)['RESTFREQ'] * u.Hz))
PB_8 = get_pb_corr(results.r, pb_noema(fits.getheader(HC3N_TdV_8_7)['RESTFREQ'] * u.Hz))

Ratio_pb_corr = Ratio / PB_10 * PB_8
hdu_R_PB_corr = fits.PrimaryHDU(Ratio_pb_corr, hd8)

plt.ion()

fig = aplpy.FITSFigure(hdu_R, figsize=(4,4))
fig.show_colorscale(vmin=0, vmax=1, cmap='inferno')
setup_plot_noema(fig, label_col='black', star_col='white')
fig.add_colorbar(axis_label_text='Ratio of integrated intensity maps')
fig.show_regions(region_file_white)
fig.add_label(0.75, 0.9, r'HC$_3$N (10-9)/(8-7)', relative=True, color='black')
#
fig.show_contour(fits.PrimaryHDU(PB_10, hd8), levels=[0.5], colors='black', linewidths=0.7)
fig.show_contour(fits.PrimaryHDU(PB_8, hd8), levels=[0.5], colors='yellow', linewidths=0.7)


fig2 = aplpy.FITSFigure(hdu_R_PB_corr, figsize=(4,4))
fig2.show_colorscale(vmin=0, vmax=1, cmap='inferno')
setup_plot_noema(fig2, label_col='black', star_col='white')
fig2.add_colorbar(axis_label_text='Ratio of integrated intensity maps')
fig2.show_regions(region_file_white)
fig2.add_label(0.75, 0.9, r'HC$_3$N (10-9)/(8-7)', relative=True, color='black')
#
fig2.show_contour(fits.PrimaryHDU(PB_10, hd8), levels=[0.5], colors='black', linewidths=0.7)
fig2.show_contour(fits.PrimaryHDU(PB_8, hd8), levels=[0.5], colors='yellow', linewidths=0.7)

from astropy.wcs import WCS
from regions import read_ds9

regions = read_ds9(region_file)
wcs_ratio = WCS(hdu_R_PB_corr)

mask_ratio = (regions[0].to_pixel(wcs_ratio)).to_mask()
N10_cutout = mask_ratio.cutout(N10)
PB10_cutout = mask_ratio.cutout(PB_10)
#
ratio_cutout = mask_ratio.cutout(Ratio)
ratio_pb_corr_cutout = mask_ratio.cutout(Ratio_pb_corr)
#
my_mask = (mask_ratio.data == 1)
# <N_10> and <N_10>_PB_corr
N10_cutout_mean = np.mean(N10_cutout[my_mask])
N10_pb_cutout_mean = np.mean(N10_cutout[my_mask] / PB10_cutout[my_mask])
print("Mean col={0:3.2}x10^13, Mean col PB corr={1:3.2}x10^13".format(N10_cutout_mean*1e-13/Pop_Up,
       N10_pb_cutout_mean*1e-13/Pop_Up))
print("Difference of {0:2.0} %".format((N10_cutout_mean - N10_pb_cutout_mean) * 100 / N10_cutout_mean))
#
ratio_cutout_mean = np.mean(ratio_cutout[my_mask])
ratio_pb_corr_cutout_mean = np.mean(ratio_pb_corr_cutout[my_mask])
#
A_pix = hd8['CDELT1'] * hd8['CDELT1'] * u.deg**2
A_beam = np.pi * hd8['BMAJ'] * hd8['BMIN'] / 4 * u.deg**2
N_pix = (np.sum(my_mask) * A_pix / A_beam).decompose()
#
ratio_cutout_std = np.std(ratio_cutout[my_mask])
ratio_pb_corr_cutout_std = np.std(ratio_pb_corr_cutout[my_mask])
ratio_cutout_mean_err = ratio_cutout_std/np.sqrt(N_pix)
ratio_pb_corr_cutout_mean_err = ratio_pb_corr_cutout_std/np.sqrt(N_pix)

print(np.around(ratio_cutout_mean, decimals=2), np.around(ratio_cutout_mean_err, decimals=2))
print(np.around(ratio_pb_corr_cutout_mean, decimals=2), np.around(ratio_pb_corr_cutout_mean_err, decimals=2))


N_tot_mean = N10_cutout_mean / Pop_Up / u.cm**2
N_tot_mean_pb = N10_pb_cutout_mean / Pop_Up / u.cm**2
Area_pix = (A_pix.to(u.rad**2) * (distance*u.pc)**2).to(u.pc**2, equivalencies=u.dimensionless_angles())

N_tot_cutout_sum = (np.sum(N10_cutout[my_mask]) * Area_pix / Pop_Up).decompose()
N_tot_pb_cutout_sum = (np.sum(N10_cutout[my_mask] / PB10_cutout[my_mask]) * Area_pix / Pop_Up).decompose()
print('{0:.2e}'.format(N_tot_cutout_sum))
X_HC3N_min = 3.0e-10
X_HC3N_max = 2.8e-9

from astropy.constants import m_p, G
mu_mass = 2.8
M_total_min = (N_tot_cutout_sum * mu_mass * m_p).to(u.Msun) / X_HC3N_max
M_total_max = (N_tot_cutout_sum * mu_mass * m_p).to(u.Msun) / X_HC3N_min
print(np.around(M_total_min, decimals=3), np.around(M_total_max, decimals=3))

M_total_pb_min = (N_tot_pb_cutout_sum * mu_mass * m_p).to(u.Msun) / X_HC3N_max
M_total_pb_max = (N_tot_pb_cutout_sum * mu_mass * m_p).to(u.Msun) / X_HC3N_min
print(np.around(M_total_pb_min, decimals=3), np.around(M_total_pb_max, decimals=3))
