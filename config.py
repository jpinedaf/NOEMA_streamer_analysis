import astropy.units as u
import numpy as np
# set tickmarks inwards
import matplotlib as mpl
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'

HC3N_TdV_10_9 = 'data/Per-emb-2-HC3N_10-9_TdV.fits'
HC3N_Vc_10_9 = 'data/Per-emb-2-HC3N_10-9_fit_Vc.fits'
region_file = 'Streamer_North_v2.reg'
ra_Per2 = 15 * (3 + (32 + 17.92/60.) / 60.) * u.deg
dec_Per2 = (30 + (49 + 48.03 / 60.) / 60.) * u.deg
distance = 300.

poly = np.array([53.079519, 30.83452, 53.078404, 30.836177, 53.076689, 30.83787, 53.075746, 30.837944,
                 53.075017, 30.837208, 53.075067, 30.835867, 53.076296, 30.833681, 53.07735, 30.831875,
                 53.077525, 30.829164, 53.078662, 30.827808, 53.079905, 30.82885, 53.080033, 30.832016])
ra_poly = np.append(np.reshape(poly, [-1, 2])[:, 0], poly[0])*u.deg
dec_poly = np.append(np.reshape(poly, [-1, 2])[:, 1], poly[1])*u.deg
# np.array([53.078662, 53.077437, 53.075771, 53.075683, 53.075246, 53.075067, 53.076296, 53.07735,
#                     53.077525, 53.078662, 53.079979, 53.079629, 53.078662])*u.deg
# dec_poly = np.array([30.836617, 30.840306, 30.840231, 30.839853, 30.838425, 30.835867, 30.833681, 30.831875,
#                      30.829164, 30.827808, 30.829392, 30.834283, 30.836617])*u.deg


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


def PB_NOEMA(freq_obs):
    """ Primary beam diameter for NOEMA at the observed frequency.

    freq_obs is the observed frequency in GHz.

    PB = 64.1 * (72.78382*u.GHz) / freq_obs
    """
    return (64.1 * u.arcsec * 72.78382 * u.GHz / freq_obs).decompose()


def PB_SMA(freq_obs):
    """ Primary beam diameter for SMA at the observed frequency.

    freq_obs is the observed frequency in GHz.

    PB = 48.0 * (231.0*u.GHz) / freq_obs
    """
    return (48.0 * u.arcsec * 231 * u.GHz / freq_obs).decompose()


def setup_plot_NOEMA(fig_i, label_col='black', star_col='red'):
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