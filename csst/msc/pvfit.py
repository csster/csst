import numpy as np
import matplotlib.pyplot as plt


def radec_uv(ra, dec, ra0, dec0):
    ra_ref, dec_ref = ra0 / (180. / np.pi), dec0 / (180. / np.pi)
    ra1, dec1 = ra / (180. / np.pi), dec / (180. / np.pi)
    u = np.cos(dec1) * np.sin(ra1 - ra_ref) / (
            np.cos(dec_ref) * np.cos(dec1) * np.cos(ra1 - ra_ref) + np.sin(dec_ref) * np.sin(dec1)) * (180. / np.pi)
    v = (np.cos(dec_ref) * np.sin(dec1) - np.cos(dec1) * np.sin(dec_ref) * np.cos(ra1 - ra_ref)) / (
            np.sin(dec_ref) * np.sin(dec1) + np.cos(dec_ref) * np.cos(dec1) * np.cos(ra1 - ra_ref)) * (180. / np.pi)

    return u, v


def uv_radec(u, v, ra0, dec0):
    ra_ref, dec_ref = ra0 / (180. / np.pi), dec0 / (180. / np.pi)
    u, v = u / 180. * np.pi, v / 180. * np.pi
    ra = ra0 + np.degrees(np.arctan(u / (np.cos(dec_ref) - v * np.sin(dec_ref))))
    dec = np.degrees(
        np.arctan(np.cos(ra / (180. / np.pi) - ra_ref) * (v + np.tan(dec_ref)) / (1 - v * np.tan(dec_ref))))

    return ra, dec


def xy_uv(x, y, header):
    u = (header["CD1_1"] * (x - header["CRPIX1"]) + header["CD1_2"] * (y - header["CRPIX2"]))  # /180.*np.pi
    v = (header["CD2_1"] * (x - header["CRPIX1"]) + header["CD2_2"] * (y - header["CRPIX2"]))  # /180.*np.pi

    return u, v


def uv_xy(u, v, header):
    x = header["CRPIX1"] + (header["CD2_2"] * u - header["CD1_2"] * v) / (
            header["CD2_2"] * header["CD1_1"] - header["CD1_2"] * header["CD2_1"])
    y = header["CRPIX2"] + (header["CD2_1"] * u - header["CD1_1"] * v) / (
            header["CD1_2"] * header["CD2_1"] - header["CD2_2"] * header["CD1_1"])

    return x, y


def rotate(cdelt, corota0):
    corota = corota0 / (180. / np.pi)
    cd11 = cdelt * np.cos(corota)
    cd12 = cdelt * np.sin(corota)
    cd21 = -cdelt * np.sin(corota)
    cd22 = cdelt * np.cos(corota)

    return cd11, cd12, cd21, cd22


def pv_fit(u_xy, v_xy, u_radec, v_radec):
    u_detect_match, v_detect_match = u_xy, v_xy
    r_detect_match = np.sqrt(u_detect_match ** 2 + v_detect_match ** 2)

    u_gaia_match, v_gaia_match = u_radec, v_radec
    # solve the parameters
    list_all, list_u, list_v = [], [], []
    for i in range(len(u_detect_match)):
        list_all.append(
            [1, u_detect_match[i], v_detect_match[i], r_detect_match[i], u_detect_match[i] ** 2,
             u_detect_match[i] * v_detect_match[i], v_detect_match[i] ** 2,
             u_detect_match[i] ** 3, u_detect_match[i] ** 2 * v_detect_match[i],
             u_detect_match[i] * v_detect_match[i] ** 2, v_detect_match[i] ** 3,
             r_detect_match[i] ** 3])  # , u_detect_match[i] ** 4, u_detect_match[i] ** 3 * v_detect_match[i],
        #             u_detect_match[i] ** 2 * v_detect_match[i] ** 2,
        #             u_detect_match[i] * v_detect_match[i] ** 3, v_detect_match[i] ** 4])
        list_u.append(u_gaia_match[i])
        list_v.append(v_gaia_match[i])

    m = np.array(list_all)[:, :]
    v_u = np.array(list_u)
    v_v = np.array(list_v)

    fit_u = np.linalg.lstsq(m, v_u, rcond=None)  # [0]
    fit_v = np.linalg.lstsq(m, v_v, rcond=None)  # [0]
    r_u, residual_u = fit_u[0], fit_u[1]
    r_v, residual_v = fit_v[0], fit_v[1]
    print(r_u)
    print(r_v)

    return r_u, r_v


def uv_cor(u, v, r_u, r_v):
    u_old, v_old = u, v
    r_old = np.sqrt(u_old ** 2 + v_old ** 2)
    u_new = r_u[0] + r_u[1] * u_old + r_u[2] * v_old + r_u[3] * r_old + r_u[
        4] * u_old ** 2 + r_u[5] * u_old * v_old + r_u[6] * v_old ** 2 + r_u[
                7] * u_old ** 3 + r_u[8] * u_old ** 2 * v_old + r_u[
                9] * u_old * v_old ** 2 + r_u[10] * v_old ** 3 + r_u[11] * r_old ** 3  # + r_u[12] * u_old ** 4 + r_u[
    #                13] * u_old ** 3 * v_old + r_u[14] * u_old ** 2 * v_old ** 2 + r_u[15] * u_old * v_old ** 3 + r_u[
    #                16] * v_old ** 4
    v_new = r_v[0] + r_v[1] * u_old + r_v[2] * v_old + r_v[3] * r_old + r_v[
        4] * u_old ** 2 + r_v[5] * u_old * v_old + r_v[6] * v_old ** 2 + r_v[
                7] * u_old ** 3 + r_v[8] * u_old ** 2 * v_old + r_v[
                9] * u_old * v_old ** 2 + r_v[10] * v_old ** 3 + r_v[11] * r_old ** 3  # + r_v[12] * u_old ** 4 + r_v[
    #                13] * u_old ** 3 * v_old + r_v[14] * u_old ** 2 * v_old ** 2 + r_v[15] * u_old * v_old ** 3 + r_v[
    #                16] * v_old ** 4

    return u_new, v_new


def uv_match(x_image, y_image, x_gaia, y_gaia, match_aperture=100.):
    """
    Calculate the relation between x,y and delta_x, delta_y
    :x_gaia: x of gaia stars in pixel from raw wcs, array
    :y_gaia: y of gaia stars in pixel from raw wcs, array
    :x_image: x of stars in pixel from detection catalog, array
    :y_image: y of stars in pixel from detection catalog, array
    :match_aperture: the aperture for star matching, default is 100
    """

    x, y, ra, dec, dis_min, flag_xy, flag_radec = [], [], [], [], [], [], []
    for i in range(len(x_image)):
        distance = (x_gaia - x_image[i]) ** 2 + (y_gaia - y_image[i]) ** 2
        index = np.argmin(np.sqrt(distance))
        if np.min(np.sqrt(distance)) < match_aperture:  # find the closest star within 50 pixels
            x.append(x_image[i])
            y.append(y_image[i])
            ra.append(x_gaia[index])
            dec.append(y_gaia[index])
            dis_min.append(np.sqrt((x_gaia[index] - x_image[i]) ** 2 + (y_gaia[index] - y_image[i]) ** 2))
            flag_xy.append(i)
            flag_radec.append(index)

    x0, y0, ra0, dec0, dis_min0 = np.array(x), np.array(y), np.array(ra), np.array(dec), np.array(dis_min)
    print(len(dis_min0), np.sum(dis_min0) / len(dis_min0), np.std(x0 - ra0), np.std(y0 - dec0))

    return x0, y0, ra0, dec0, flag_xy, flag_radec


def xy_radec(x, y, h, r_u, r_v):
    u_xy, v_xy = xy_uv(x, y, h)
    u_cor, v_cor = uv_cor(u_xy, v_xy, r_u, r_v)
    ra_cor, dec_cor = uv_radec(u_cor, v_cor, h["CRVAL1"], h["CRVAL2"])

    return ra_cor, dec_cor


def radec_xy(ra, dec, h):
    u_radec, v_radec = radec_uv(ra, dec, h["CRVAL1"], h["CRVAL2"])
    x, y = uv_xy(u_radec, v_radec, h)

    return x, y


#######
def pointsource_detection(data, fwhm, threshold, roundlo, roundhi):
    from photutils.detection import DAOStarFinder
    from astropy.stats import sigma_clipped_stats
    from photutils.centroids import centroid_sources, centroid_1dg
    import warnings
    from photutils.aperture import CircularAperture, CircularAnnulus
    from photutils.aperture import aperture_photometry

    mean, median, std = sigma_clipped_stats(data, sigma=3.0)

    daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold * std, ratio=1, exclude_border=True, roundlo=roundlo,
                            roundhi=roundhi)
    sources = daofind(data - median)
    positions = np.transpose((sources['xcentroid'], sources['ycentroid']))

    warnings.filterwarnings("ignore")
    x_new, y_new = centroid_sources(data, sources['xcentroid'], sources['ycentroid'], box_size=11,
                                    centroid_func=centroid_1dg)
    positions_new = list(zip(x_new, y_new))

    aperture = CircularAperture(positions_new, r=6.)
    annulus_aperture = CircularAnnulus(positions_new, r_in=10, r_out=15)
    phot_table = aperture_photometry(data - median, aperture)
    flux = phot_table["aperture_sum"]
    magnitude = -2.5 * np.log10(flux)

    sources_new = list(zip(x_new, y_new, flux, magnitude))

    return sources_new


def cat_match(u_radec, v_radec, u_xy, v_xy, u_xy_base, v_xy_base, match_aperture):
    import torch

    u_radec_match, v_radec_match, u_xy_match, v_xy_match, flag_radec, flag_xy = uv_match(u_radec, v_radec, u_xy,
                                                                                         v_xy,
                                                                                         match_aperture=match_aperture)
    u_xy_new, v_xy_new = torch.index_select(torch.tensor(u_xy_base), 0, torch.tensor(flag_xy)), torch.index_select(
        torch.tensor(v_xy_base), 0, torch.tensor(flag_xy))
    u_radec_new, v_radec_new = torch.index_select(torch.tensor(u_radec), 0,
                                                  torch.tensor(flag_radec)), torch.index_select(torch.tensor(v_radec),
                                                                                                0, torch.tensor(
            flag_radec))
    r_u, r_v = pv_fit(u_xy_new, v_xy_new, u_radec_new, v_radec_new)

    return r_u, r_v


def ref_shrink(ra, dec, head):
    x_gaia, y_gaia = radec_xy(ra, dec, head)
    x_chunk = x_gaia[np.where((abs(x_gaia - 5000) < 6000) & (abs(y_gaia - 5000) < 6000))]
    y_chunk = y_gaia[np.where((abs(x_gaia - 5000) < 6000) & (abs(y_gaia - 5000) < 6000))]
    u_xy, v_xy = xy_uv(x_chunk, y_chunk, head)
    ra, dec = uv_radec(u_xy, v_xy, head["CRVAL1"], head["CRVAL2"])

    return ra, dec

#######
def singlechip_wcsfit(fitsname, ref_cat="gaia.fits"):
    """ Function to calculate the pv parameters of single chip CSST data

    :param fitsname: *_img.fits
    :param ref_cat: reference catalog from gaia
    :param time_interval: time interval between gaia and obs time, default is 2027.67 - 2016.0. If no pm, set 0
    :return: the wcs including PV parameters, and the median and std of delta_ra and delta_dec
    """

    from astropy.io import fits
    from scipy.spatial.distance import cdist
    import torch

    sci = fits.open(fitsname)
    data = sci[1].data
    head_ori = sci[1].header
    exptime = sci[0].header["EXPSTART"]
    time_interval = (exptime - 51543.0)/365.0
    ra_ref = head_ori["CRVAL1"]
    dec_ref = head_ori["CRVAL2"]
    gaia = fits.open(ref_cat)
    gaia_cat = gaia[1].data
    pm_ra = gaia_cat["pmra"]
    pm_dec = gaia_cat["pmdec"]
    ra_gaia = gaia_cat["X_WORLD"] + pm_ra * time_interval / 1000. / 3600. / np.cos(gaia_cat["Y_WORLD"] / 180. * np.pi)
    dec_gaia = gaia_cat["Y_WORLD"] + pm_dec * time_interval / 1000. / 3600.  # correct the pm for ra and dec

    ra_gaia_all = ra_gaia[~np.isnan(ra_gaia)]  # remove the nan value of gaia catalog
    dec_gaia_all = dec_gaia[~np.isnan(dec_gaia)]

    ra_gaia, dec_gaia = ref_shrink(ra_gaia_all, dec_gaia_all ,head_ori)
    print(str(len(ra_gaia)) + " gaia stars have been used")

    # detect point sources by photilus
    sources = pointsource_detection(data, fwhm=2.0, threshold=10.0, roundlo=-0.3, roundhi=0.3)
    x_detect, y_detect, flux_detect, mag_detect = zip(*sources)
    x_detect = np.array(x_detect) + 1  # difference betwewen python and sextractor frame
    y_detect = np.array(y_detect) + 1
    print(str(len(x_detect)) + " objects(stars) have been detected")

    stars = 0
    u_radec, v_radec = radec_uv(ra_gaia, dec_gaia, ra_ref, dec_ref)
    for i in range(-15, 15):
        for j in range(-15, 15):
            u_xy, v_xy = xy_uv(x_detect + i * 10, y_detect + j * 10, head_ori)
            uvp_radec = np.column_stack((u_radec, v_radec))
            uvp_xy = np.column_stack((u_xy, v_xy))
            dist = cdist(uvp_radec, uvp_xy, metric='euclidean')
            star_match = len(dist[np.where(dist < 0.001)])
            if star_match > stars:
                x_shift, y_shift, stars = i * 10, j * 10, star_match

    print(stars, x_shift, y_shift)

    u_xy, v_xy = xy_uv(x_detect + x_shift, y_detect + y_shift, head_ori)
    u_xy_base, v_xy_base = xy_uv(x_detect, y_detect, head_ori)
    u_radec, v_radec = radec_uv(ra_gaia, dec_gaia, ra_ref, dec_ref)

    r_u, r_v = cat_match(u_radec, v_radec, u_xy, v_xy, u_xy_base, v_xy_base, match_aperture=0.001)
    u_xy_new, v_xy_new = uv_cor(u_xy_base, v_xy_base, r_u, r_v)

    r_u, r_v = cat_match(u_radec, v_radec, u_xy_new, v_xy_new, u_xy_base, v_xy_base, match_aperture=0.0005)
    u_xy_new, v_xy_new = uv_cor(u_xy_base, v_xy_base, r_u, r_v)

    r_u, r_v = cat_match(u_radec, v_radec, u_xy_new, v_xy_new, u_xy_base, v_xy_base, match_aperture=0.0003)
    u_xy_new, v_xy_new = uv_cor(u_xy_base, v_xy_base, r_u, r_v)

    u_radec_match, v_radec_match, u_xy_match, v_xy_match, flag_radec, flag_xy = uv_match(u_radec, v_radec,
                                                                                         u_xy_new, v_xy_new,
                                                                                         match_aperture=0.0001)
    u_xy_new, v_xy_new = torch.index_select(torch.tensor(u_xy_base), 0, torch.tensor(flag_xy)), torch.index_select(
        torch.tensor(v_xy_base), 0, torch.tensor(flag_xy))
    u_radec_new, v_radec_new = torch.index_select(torch.tensor(u_radec), 0,
                                                  torch.tensor(flag_radec)), torch.index_select(torch.tensor(v_radec),
                                                                                                0, torch.tensor(
            flag_radec))
    r_u, r_v = pv_fit(u_xy_new, v_xy_new, u_radec_new, v_radec_new)

    ra_old, dec_old = uv_radec(u_radec_new, v_radec_new, ra_ref, dec_ref)
    u_xy_new, v_xy_new = uv_cor(u_xy_new, v_xy_new, r_u, r_v)
    ra_star, dec_star = uv_radec(u_xy_new, v_xy_new, ra_ref, dec_ref)

    delta_ra = np.array(ra_star - ra_old) * 3600. * 1000. * np.cos(np.array(dec_old) / 180. * np.pi)
    delta_dec = np.array(dec_star - dec_old) * 3600. * 1000.
    print(len(delta_ra), np.median(delta_ra), np.std(delta_ra))
    print(np.median(delta_dec), np.std(delta_dec))

    position = head_ori.index("WCSDIM")
    wcs_infor = head_ori[position:position + 13]
    wcs_infor["CTYPE1"] = "RA---TPV"
    wcs_infor["CTYPE2"] = "DEC--TPV"
    wcs_infor.set("RADESYS", "ICRS")

    wcs_infor.set("PV1_0", r_u[0])
    wcs_infor.set("PV1_1", r_u[1])
    wcs_infor.set("PV1_2", r_u[2])
    wcs_infor.set("PV1_3", r_u[3])
    wcs_infor.set("PV1_4", r_u[4])
    wcs_infor.set("PV1_5", r_u[5])
    wcs_infor.set("PV1_6", r_u[6])
    wcs_infor.set("PV1_7", r_u[7])
    wcs_infor.set("PV1_8", r_u[8])
    wcs_infor.set("PV1_9", r_u[9])
    wcs_infor.set("PV1_10", r_u[10])
    wcs_infor.set("PV1_11", r_u[11])

    wcs_infor.set("PV2_0", r_v[0])
    wcs_infor.set("PV2_1", r_v[2])
    wcs_infor.set("PV2_2", r_v[1])
    wcs_infor.set("PV2_3", r_v[3])
    wcs_infor.set("PV2_4", r_v[6])
    wcs_infor.set("PV2_5", r_v[5])
    wcs_infor.set("PV2_6", r_v[4])
    wcs_infor.set("PV2_7", r_v[10])
    wcs_infor.set("PV2_8", r_v[9])
    wcs_infor.set("PV2_9", r_v[8])
    wcs_infor.set("PV2_10", r_v[7])
    wcs_infor.set("PV2_11", r_v[11])
    wcs_infor.set("RAOFF", np.round(np.median(delta_ra), 2), comment="mas in unit")
    wcs_infor.set("DECOFF", np.round(np.median(delta_dec), 2), comment="mas in unit")
    wcs_infor.set("RASTD", np.round(np.std(delta_ra), 2), comment="mas in unit")
    wcs_infor.set("DECSTD", np.round(np.std(delta_dec), 2), comment="mas in unit")

    return wcs_infor


# test
def main():
    wcs = singlechip_wcsfit(fitsname="CSST_MSC_MS_SCI_20270810081950_20270810082220_100000100_13_img.fits",
                            ref_cat="gaia.fits")

    print(wcs)


if __name__ == "__main__":
    main()
