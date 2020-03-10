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
ratio_cutout = mask_ratio.cutout(Ratio)
ratio_pb_corr_cutout = mask_ratio.cutout(Ratio_pb_corr)
#
my_mask = (mask_ratio.data == 1)
# ratio_streamer = ratio_cutout[gd]
# ratio_pb_corr_streamer = ratio_pb_corr_cutout[gd]

ratio_cutout_mean = np.mean(ratio_cutout[my_mask == 1])
ratio_pb_corr_cutout_mean = np.mean(ratio_pb_corr_cutout[my_mask == 1])

A_beam_pix = np.pi * hd8['BMAJ'] * hd8['BMIN'] / (4 * hd8['CDELT1'] * hd8['CDELT1'])
N_pix = np.sum(my_mask)/A_beam_pix

ratio_cutout_std = np.std(ratio_cutout[my_mask ==1.])
ratio_pb_corr_cutout_std = np.std(ratio_pb_corr_cutout[my_mask ==1.])
ratio_cutout_mean_err = ratio_cutout_std/np.sqrt(N_pix)
ratio_pb_corr_cutout_mean_err = ratio_pb_corr_cutout_std/np.sqrt(N_pix)

print(np.around(ratio_cutout_mean, decimals=2), np.around(ratio_cutout_mean_err, decimals=2))
print(np.around(ratio_pb_corr_cutout_mean, decimals=2), np.around(ratio_pb_corr_cutout_mean_err, decimals=2))
