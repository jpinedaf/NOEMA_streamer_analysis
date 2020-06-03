import astropy.units as u
import numpy as np
from astropy.io import fits
import matplotlib as mpl
from matplotlib import rc
from astropy.constants import c, h, k_B

#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('font',**{'family':'serif'})#,'sans-serif':['Helvetica']})
rc('text', usetex=True)
# set tickmarks inwards
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'

stream_pickle = 'Streamer_model.pickle'
vlsr_rad_kde_pickle = 'Velocity_Radius_KDE.pickle'

# file_HC3N_10_9_TdV_mJy
HC3N_TdV_8_7 = 'data/Per-emb-2-HC3N_8-7_TdV.fits'
HC3N_TdV_10_9 = 'data/Per-emb-2-HC3N_10-9_TdV.fits'
HC3N_TdV_10_9_match = 'data/Per-emb-2-HC3N_10-9_TdV_smooth_regrid.fits'
HC3N_Ratio = 'data/Per-emb-2-HC3N_ratio_10-9_8-7.fits'
HC3N_Vc_10_9 = 'data/Per-emb-2-HC3N_10-9_fit_Vc.fits'
HC3N_dv_10_9 = 'data/Per-emb-2-HC3N_10-9_fit_sigma_v.fits'
ALMA_cont = 'data/Per-emb-2_1.3mm_ALMA_mJy.fits'
# file_CCS_TdV_mJy
CCS_TdV = 'data/Per-emb-2-CCS_8_7-7_6_TdV.fits'
CS_TdV = 'data/Per-emb-2-13CS_2-1_TdV.fits'
N2Hp_TdV = 'data/Per-emb-2-N2Hp-1-0_TdV.fits'
N2Dp_TdV = 'data/Per-emb-2-N2Dp-1-0_TdV.fits'
#
region_file_white = 'Streamer_North.reg'
region_file = 'Streamer_North_v2.reg'
stream_pickle = 'Streamer_model.pickle'
vlsr_rad_kde_pickle = 'Velocity_Radius_KDE.pickle'
ra_Per2 = 15 * (3 + (32 + 17.92/60.) / 60.) * u.deg
dec_Per2 = (30 + (49 + 48.03 / 60.) / 60.) * u.deg
ra_ALMA_zoom = 53.07475*u.deg
dec_ALMA_zoom = 30.8299*u.deg
radius_ALMA_zoom = (2.5*u.arcsec).to(u.deg)
distance = 300.
HC3N_10_9_levels = 8e-3 * np.arange(5, 20, 3)

# CO contour levels to show outflow
file_12CO_blue = 'data/Per-emb-2-CO_2-1-TdV-blue.fits'
file_12CO_red = 'data/Per-emb-2-CO_2-1-TdV-red.fits'
CO_red_levs = 0.2 * np.arange(5, 55, 10)
CO_blue_levs = 0.2 * np.arange(5, 55, 10)
#
# HC3N ratios
R_10_8_mean = 0.48
R_10_8_mean_err = 0.02
# Polygon coordinates added by hand from the region file
# The next line is to close the polygon
poly = np.array([53.079519, 30.83452, 53.078404, 30.836177, 53.076689, 30.83787, 53.075746, 30.837944,
                 53.075017, 30.837208, 53.075067, 30.835867, 53.076296, 30.833681, 53.07735, 30.831875,
                 53.077525, 30.829164, 53.078662, 30.827808, 53.079905, 30.82885, 53.080033, 30.832016])
ra_poly = np.append(np.reshape(poly, [-1, 2])[:, 0], poly[0])*u.deg
dec_poly = np.append(np.reshape(poly, [-1, 2])[:, 1], poly[1])*u.deg


def per_emb_2_get_vc_r():
    """
    Returns the centroid velocity and projected separation in the sky of the
    centroid velocity from Per-emb-2.
    r_proj is in u.au and V_los is in u.km/u.s
    :return: r_proj, V_los
    """
    from regions import read_ds9
    from astropy.wcs import WCS
    from astropy.io import fits
    import velocity_tools.coordinate_offsets as c_offset
    # load region file and WCS structures
    regions = read_ds9(region_file)
    wcs_Vc = WCS(HC3N_Vc_10_9)
    #
    hd_Vc = fits.getheader(HC3N_Vc_10_9)
    results = c_offset.generate_offsets(hd_Vc, ra_Per2, dec_Per2, pa_angle=0*u.deg, inclination=0*u.deg)
    rad_au = (results.r * distance*u.pc).to(u.au, equivalencies=u.dimensionless_angles())
    # Vc_all =
    #
    mask_Vc = (regions[0].to_pixel(wcs_Vc)).to_mask()
    Vc_cutout = mask_Vc.cutout(fits.getdata(HC3N_Vc_10_9))
    rad_cutout = mask_Vc.cutout(rad_au)
    #
    gd = (mask_Vc.data == 1)
    v_los = Vc_cutout[gd]*u.km/u.s
    r_proj = rad_cutout[gd]
    return r_proj, v_los


def pb_noema(freq_obs):
    """
    Primary beam diameter for NOEMA at the observed frequency.
        PB = 64.1 * (72.78382*u.GHz) / freq_obs

    :param freq_obs: is the observed frequency in GHz.
    :return: The primary beam FWHM in arcsec
    """
    return (64.1 * u.arcsec * 72.78382 * u.GHz / freq_obs).decompose()


def pb_sma(freq_obs):
    """
    Primary beam diameter for SMA at the observed frequency.
        PB = 48.0 * (231.0*u.GHz) / freq_obs

    :param freq_obs: is the observed frequency in GHz.
    :return: The primary beam FWHM in arcsec
    """
    return (48.0 * u.arcsec * 231 * u.GHz / freq_obs).decompose()


def setup_plot_noema(fig_i, label_col='black', star_col='red'):
    """
    Setup of NOEMA plots, since they will show all the same format.
    """
    fig_i.set_system_latex(True)
    fig_i.ticks.set_color(label_col)
    fig_i.recenter(53.075, 30.8299, radius=45 * (u.arcsec).to(u.deg))
    fig_i.set_nan_color('0.9')
    fig_i.add_beam(color=label_col)
    distance = 300.  # pc
    ang_size = (5e3 / distance) * u.arcsec
    fig_i.add_scalebar(ang_size, label='5,000 au', color=label_col)
    fig_i.show_markers(ra_Per2.value, dec_Per2.value, marker='*', s=60, layer='star',
                       edgecolor=star_col, facecolor=label_col, zorder=31)
    fig_i.tick_labels.set_xformat('hh:mm:ss')
    fig_i.tick_labels.set_yformat('dd:mm:ss')
    fig_i.ticks.set_length(7)
    fig_i.axis_labels.set_xtext(r'Right Ascension (J2000)')
    fig_i.axis_labels.set_ytext(r'Declination (J2000)')
    return

def convert_into_mili(file_name):
    """
    It converts a file into one rescaled by 1e3.
    This is useful to convert between Jy -> mJy or m/s into km/s
    for plotting purposes (e.g. to use with aplpy).

    Usage:
    fig = aplpy.FITSFigure(convert_into_mili(file_in_Jy), figsize=(4,4))
    fig.show_colorscale(vmin=0, vmax=160, cmap='inferno')
    fig.add_colorbar()
    fig.colorbar.set_axis_label_text(r'Integrated Intensity (mJy beam$^{-1}$ km s$^{-1}$)')

    :param file_name: string with filename to process
    :return: hdu
    """
    data, hd = fits.getdata(file_name, header=True)
    return fits.PrimaryHDU(data=data*1e3, header=hd)

def plot_ratios(ax, dens, ratio, temp, N_col=13, label_i=r'N(HC$_3$N)=10$^{12}$ cm$^{-2}$', file_out='test.pdf'):
    """
    """
    ax.scatter(dens[temp == 10.], ratio[temp == 10.], c='#377eb8', alpha=0.7, label=r'T$_{\rm kin}$=10 K')
    ax.scatter(dens[temp == 12.25], ratio[temp == 12.25], c='#e41a1c', alpha=0.7, label=r'T$_{\rm kin}$=12.25 K')
    dens_min = 500
    dens_max = 3e6
    ax.set_xlim(dens_min, dens_max)
    ax.set_ylim(0, 1)
    ax.set_xscale('log')
    ax.set_xlabel(r'H$_2$ Density (cm$^{-3}$)')
    ax.set_ylabel(r'Ratio of HC$_3$N (10-9)/(8-7)')
    ax.hlines(R_10_8_mean, dens_min, dens_max)
    ax.hlines([R_10_8_mean - R_10_8_mean_err, R_10_8_mean + R_10_8_mean_err], dens_min, dens_max, linestyle="--")
    ax.text(1e3, 0.85, label_i, fontsize=12)
    ax.text(1e3, 0.775, r'$\Delta$v=0.5 km s$^{-1}$', fontsize=12)
    ax.legend(frameon=False, loc=4)
    return


def col_dens_hc3n(file_tdv, nup=8, get_tdv=False):
    """
    converts integrated intensity file into column density of the upper level

    :param file_tdv: FITS file with integrated intensity file
    :param nup: integer (default = 8)
    :return: Column density in units of cm^-2
    """
    from radio_beam import Beam
    if nup == 8:
        nu = 72.78382*u.GHz
        A_ij = 0.294E-04/u.s
    elif nup == 10:
        nu = 90.97902*u.GHz
        A_ij = 0.581e-4/u.s
    else:
        print('This transition is not available yet')
        return np.nan
    data, hd = fits.getdata(file_tdv, header=True)
    beam = Beam.from_fits_header(hd)
    # freq = hd['RESTFREQ'] * u.Hz
    Jy_K = (1*u.Jy).to(u.K, equivalencies=beam.jtok_equiv(nu))
    TdV = data * Jy_K
    if get_tdv:
        return TdV, ((8 * np.pi * k_B * nu ** 2) / (h * c ** 3 * A_ij) * TdV * u.km / u.s).to(1 / u.cm ** 2)
    else:
        return ((8 * np.pi * k_B * nu**2)/(h * c**3 * A_ij) * TdV * u.km / u.s).to(1/u.cm**2)
