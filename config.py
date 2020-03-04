import astropy.units as u
import numpy as np

HC3N_TdV_10_9 = 'data/Per-emb-2-HC3N_10-9_TdV.fits'
HC3N_Vc_10_9 = 'data/Per-emb-2-HC3N_10-9_fit_Vc.fits'
region_file = 'Streamer_North.reg'
ra_Per2 = 15 * (3 + (32 + 17.92/60.) / 60.) * u.deg
dec_Per2 = (30 + (49 + 48.03 / 60.) / 60.) * u.deg
distance = 300.

ra_poly = np.array([53.078662, 53.077437, 53.075771, 53.075683, 53.075246, 53.075067, 53.076296, 53.07735,
                    53.077525, 53.078662, 53.079979, 53.079629, 53.078662])*u.deg
dec_poly = np.array([30.836617, 30.840306, 30.840231, 30.839853, 30.838425, 30.835867, 30.833681, 30.831875,
                     30.829164, 30.827808, 30.829392, 30.834283, 30.836617])*u.deg


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
    results = c_offset.generate_offsets(hd_Vc, ra_Per2, dec_Per2, PA_Angle=0*u.deg, inclination=0*u.deg)
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