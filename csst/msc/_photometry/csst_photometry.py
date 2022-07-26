#!/usr/local/bin/python3.8
# author Hu Zou
# v1.1.2 first version of csst photometry pipeline based on sextractor
# v1.2.0 change weight and flag names due to the CSST naming
#        add FLAGS_ISO to remove flagged objects in PSFEx
# v1.2.1 change gain to exposure time
#        add options of output images
#        change catalog columns to the defined L2 data format
#        change config directory
# v1.2.2 add pick_psfstars, select PSF stars
#        add photometry types: PSF or MODEL
# v1.2.3 add computing time
# v1.2.4 bug in pick_psfstars with missing max_elp

# from __future__ import (absolute_import, division, print_function, unicode_literals)

import argparse
import os
import sys
import time

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy import table
from astropy import wcs as awcs
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from scipy.interpolate import UnivariateSpline

# import ..magfluxconvert as magf
from .magfluxconvert import asinhpogson, fluxerr2magerr, magerr2fluxerr
from .stats import sigmaclip_limitsig, weighted_mean
# import system
from shutil import which

# prog_dir = sys.path[0]
# CONFIG_PATH = os.path.join(prog_dir, 'config/')

from csst import PACKAGE_PATH
from ..data_manager import CsstMscDataManager

CONFIG_PATH = PACKAGE_PATH + "/msc/phot_config/"

__version__ = "1.2.3"


def fits_center(head):
    nx = head['NAXIS1']
    ny = head['NAXIS2']
    wcs = awcs.WCS(head)
    ra, dec = wcs.all_pix2world(nx / 2.0, ny / 2.0, 1)
    return ra, dec


def valid_coordinates(coordinates, size=None):
    """
    convert tuple, list or np.ndarray of cooridinates to a 2d ndarray
      and check the coordinates within an image size
    INPUT
        coordinates: 2d array for coordinates
        size: size of an image
    OUTPUT:
        coord: reshape coordinates
        indcoord: index of coordinates in a size range
    """
    if isinstance(coordinates, (list, tuple, np.ndarray)):
        coord = np.atleast_2d(coordinates)
        if coord.shape[0] != 2 and coord.shape[1] != 2:
            raise ValueError("coordinates should have at least one axis with 2 elements")
        if coord.shape[1] != 2 and coord.shape[0] == 2:
            coord = coord.transpose()
    else:
        raise TypeError("coordinates should be list or array of (x,y) pixel positions")

    if size is None:
        return coord
    else:
        if len(size) != 2:
            raise ValueError("size should have 2 elements")
        nx, ny = size
        x = coord[:, 0]
        y = coord[:, 1]
        indcoord = np.arange(coord.shape[0])
        good = (x >= 0.5) & (x < nx + 0.5) & (y >= 0.5) & (y < ny + 0.5)
        if np.any(good):
            indcoord = indcoord[good]
        else:
            raise ValueError('coordinates are not in the image range')
        return coord, indcoord


def closest_match(coord1, coord2, min_dist=1.0):
    """
    find closest pairs between two sets of coordinates
    coord1: coordinates to be matched
    coord2: coordinates matched to
    min_dist: separation tolerance

    output: idx1,idx2
        idx1: matched index for coord1
        idx2: matched index for coord2
    """
    coord1 = valid_coordinates(coord1)
    coord2 = valid_coordinates(coord2)

    n1 = len(coord1)
    n2 = len(coord2)
    index1 = []
    index2 = []
    x2 = coord2[:, 0]
    y2 = coord2[:, 1]
    for i in range(n1):
        ix1, iy1 = coord1[i]
        index = np.where((np.abs(x2 - ix1) < min_dist) & (np.abs(y2 - iy1) < min_dist))[0]
        nmatch = len(index)
        if nmatch < 1:
            continue
        if nmatch > 1:
            x2tmp = x2[index];
            y2tmp = y2[index]
            dist = ((x2tmp - ix1) ** 2 + (y2tmp - iy1) ** 2)
            indsort = np.argsort(dist)
            index1.append(i)
            index2.append(index[indsort])
        else:
            index1.append(i)
            index2.append(index)
    return index1, index2


# def get_psf1(fp_img, outdir=None, psf_size=71, degree=3, variability=0.3, fwhm_range=[1.5, 20.0], max_elp=0.3,
#              sampling=0, nthread=1, min_sn=20, detect_thresh=5.0, detect_minarea=3):
#     """
#     get PSF profile for a specified image
#     INPUT
#         fp_img: input image
#         outdir: output directory
#         psf_size: size of the PSF image
#         degree: order of spatially varied PSF
#         variablility: allowed FWHM variability
#         fwhm_range: FWHM range of selected stars
#         max_elp: max (A-B)/(A+B) for selected stars
#         sampling: sampling step in pixels (0 = Auto)
#         nthread: number of threads to get PSF
#         min_sn: minimum S/N for souces
#         detect_thresh: threshold for detections
#         detect_minarea: min pixel area for detections
#     """
#     fitsname, _ = os.path.splitext(fp_img)
#     if outdir is not None:
#         if not os.path.exists(outdir):
#             os.mkdir(outdir)
#         _, rootname = os.path.split(fitsname)
#         fitsname = os.path.join(outdir, rootname)
#     if not system.cmd_exists('sex') or not system.cmd_exists('psfex'):
#         raise OSError('No sex or psfex SHELL command found! Please install')
#
#     sexfile = os.path.join(CONFIG_PATH, 'csst_psfex.sex')
#     covfile = os.path.join(CONFIG_PATH, 'gauss_4.0_7x7.conv')
#     nnwfile = os.path.join(CONFIG_PATH, 'default.nnw')
#     parfile = os.path.join(CONFIG_PATH, 'csst_psfex.param')
#     pexfile = os.path.join(CONFIG_PATH, 'csst.psfex')
#     psffile = fitsname + '_psf.fits'
#
#     # head=fits.getheader(fp_img)
#     # date=head['date-obs']
#     # ccd=head['ccd_no']
#     # flagfile,weightfile=get_flagweight(date,ccd)
#     rootname, _ = os.path.splitext(fp_img)
#     indpos = str.rfind(rootname, '_img')
#     flagfile = rootname[:indpos] + '_flg.fits'
#     weightfile = rootname[:indpos] + '_wht.fits'
#     head = fits.getheader(fp_img, 0)
#     gain = head['exptime']
#
#     command = 'sex ' + fp_img + ' -c ' + sexfile + ' -CATALOG_NAME ' + psffile + ' -FILTER_NAME ' + covfile + \
#               ' -STARNNW_NAME ' + nnwfile + ' -PARAMETERS_NAME ' + parfile + ' -WEIGHT_IMAGE ' + weightfile + \
#               ' -FLAG_IMAGE ' + flagfile + ' -NTHREADS ' + str(nthread) + ' -DETECT_MINAREA ' + str(detect_minarea) + \
#               ' -DETECT_THRESH ' + str(detect_thresh) + ' -GAIN ' + str(gain)
#     print(command)
#     os.system(command)
#     command = 'psfex ' + psffile + ' -c ' + pexfile + ' -PSF_SIZE ' + str(psf_size) + ',' + str(
#         psf_size) + ' -CHECKPLOT_DEV NULL -CHECKPLOT_TYPE NONE -CHECKIMAGE_TYPE NONE -PSFVAR_DEGREES ' + str(
#         degree) + ' -SAMPLE_VARIABILITY ' + str(variability) + ' -SAMPLE_FWHMRANGE ' + str(fwhm_range[0]) + ',' + str(
#         fwhm_range[1]) + ' -PSF_SAMPLING ' + str(sampling) + ' -XML_NAME ' + xmlfile + ' -SAMPLE_MAXELLIP ' + str(
#         max_elp) + ' -NTHREADS ' + str(nthread) + ' -SAMPLE_MINSN ' + str(min_sn)
#     print(command)
#     os.system(command)


def pick_psfstars(psffile, remove_polution=False, min_nstar=10, max_nstar=1500, class_star=0.7, match_dist=3.0,
                  min_separation=20, min_sn=20, max_elp=0.3):
    """
    pick PSF stars (remove PSF stars with neighbour objects, limit maximum number of stars
    INPUT
        psffile, catalog used for selecting PSF stars
        remove_polution: wether to remove objects with pollution from nearby objects
        min_nstar: mininum stars for constructing PSF model
        max_nstar: maximum stars for constructing PSF model
        max_elp: maximum ellpticity
        class_star: class threshhold for PSF stars
        match_dist: matching distance for checking polution in coordinates pixel
        min_separation: minimum separation between objects in pixels
    """

    # read catalog and keep coordinates
    hdulist = fits.open(psffile, mode='update')
    cat = hdulist[2].data

    if len(cat) < min_nstar:
        hdulist.close()
        return False

    if remove_polution:
        print("select isolate stars...")
        coord = np.array((cat['XWIN_IMAGE'], cat['YWIN_IMAGE'])).transpose()
        coord1 = coord.copy()
        index1, index2 = closest_match(coord, coord1, min_dist=match_dist)
        index1 = np.array(index1)
        mask = np.array([len(ind) == 1 for ind in index2])
        index = index1[mask]
        if len(index) >= min_nstar:
            cat = cat[index]

    if len(cat) > max_nstar:
        print('stars in the catalog, needed:', len(cat), max_nstar)
        mask = (cat['flags'] == 0) & (cat['CLASS_STAR'] > class_star) & (cat['ELLIPTICITY'] < max_elp) & (
                    cat['SNR_WIN'] > min_sn)
        if mask.sum() >= min_nstar:
            cat = cat[mask]
        if len(cat) > max_nstar:
            indsort = np.argsort(cat['flux_aper'])
            cat = cat[indsort[-max_nstar:]]

    hdulist[2].data = cat
    hdulist.flush()
    hdulist.close()
    print('isolated stars:', len(index))
    return True


# if save_starpos:
#	f=open(fitsname+'-psfxy.txt','w')
#	for i in range(len(cat)):
#		f.write('%7.2f %7.2f\n' % (cat['XWIN_IMAGE'][i],cat['YWIN_IMAGE'][i]))
#	f.close()

def get_psf(ccd_id, dm, psf_size=101, degree=2, variability=0.3, fwhm_range=[2.0, 20.0], max_elp=0.3,
            sampling=0, min_sn=20.0, detect_thresh=5.0, detect_minarea=5, seeing=0.15, pixel_scale=0.075,
            filter_name='gauss_4.0_7x7.conv', phot_aper=10.0, back_size='400,400', check_plots=False, nthread=1, **kwd):
    # ,save_starpos=False
    """
    get PSF profile for a specified image
    INPUT
        fp_img: input image
        outdir: output directory
        psf_size: size of the PSF image
        degree: order of spatially varied PSF
        variablility: allowed FWHM variability
        fwhm_range: FWHM range of selected stars
        max_elp: max (A-B)/(A+B) for selected stars
        sampling: sampling step in pixels (0 = Auto)
        min_sn: minimum S/N for souces
        detect_thresh: threshold for detections
        detect_minarea: min pixel area for detections
        seeing: PSF fwhm in arcsec
        pixel_scale: pixel scale in arcsec
        filter_name: convolve kernel file
        phot_aper: aperture used for photometry
        back_size: background grid size
        check_plots: output the check plots
        nthread: number of threads to get PSF
        **kwd: keywords for PICK_PSFSTARS
    """

    # if not system.cmd_exists('sex') or not system.cmd_exists('psfex'):
    if which("sex") is None or which("psfex") is None:
        raise OSError('No sex or psfex SHELL command found! Please install')

    # rootname, _ = os.path.splitext(fp_img)
    # # indpos = str.rfind(rootname, '_img')
    # # flagfile = rootname[:indpos] + '_flg.fits'
    # # weightfile = rootname[:indpos] + '_wht.fits'
    # flagfile = fp_img.replace("_img", "_flg")
    # weightfile = fp_img.replace("_img", "_wht")

    fp_img = dm.l1_ccd(ccd_id=ccd_id, post="img_L1.fits")
    fp_wht = dm.l1_ccd(ccd_id=ccd_id, post="wht_L1.fits")
    fp_flg = dm.l1_ccd(ccd_id=ccd_id, post="flg_L1.fits")
    fp_psf = dm.l1_ccd(ccd_id=ccd_id, post="psf.fits")
    fp_cat = dm.l1_ccd(ccd_id=ccd_id, post="cat.fits")

    head = fits.getheader(fp_img, 0)
    gain = head['exptime']
    saturate = 50000.0 / gain

    # fitsname,_=os.path.splitext(fp_img)
    # if outdir is not None:
    #     # _, filename = os.path.split(rootname[:indpos])
    #     # fitsname = os.path.join(outdir, filename)
    #     fitsname = os.path.join(outdir, os.path.basename(fp_img))

    sexfile = os.path.join(CONFIG_PATH, 'csst_psfex.sex')
    covfile = os.path.join(CONFIG_PATH, filter_name)
    nnwfile = os.path.join(CONFIG_PATH, 'default.nnw')
    parfile = os.path.join(CONFIG_PATH, 'csst_psfex.param')
    pexfile = os.path.join(CONFIG_PATH, 'csst.psfex')
    # fp_psf = fitsname + '_psf.fits'

    command = 'sex ' + fp_img + ' -c ' + sexfile + ' -CATALOG_NAME ' + fp_psf + ' -PARAMETERS_NAME ' + parfile + \
              ' -DETECT_MINAREA ' + str(detect_minarea) + ' -DETECT_THRESH ' + str(detect_thresh) +\
              ' -FILTER_NAME ' + covfile + ' -WEIGHT_IMAGE ' + fp_wht + ' -FLAG_IMAGE ' + fp_flg + \
              ' -PHOT_APERTURES ' + str(phot_aper) + ' -GAIN ' + str(gain) + ' -PIXEL_SCALE ' + str(pixel_scale) + \
              ' -SEEING_FWHM ' + str(seeing) + ' -SATUR_LEVEL ' + str(saturate) + ' -STARNNW_NAME ' + nnwfile + \
              ' -BACK_SIZE ' + back_size + ' -NTHREADS ' + str(nthread)
    print(command)
    os.system(command)

    pickflag = pick_psfstars(fp_psf, min_sn=20, **kwd)
    if pickflag is False:
        if os.path.isfile(fp_psf): os.remove(fp_psf)
        return False

    if not check_plots:
        command = 'psfex ' + fp_psf + ' -c ' + pexfile + ' -PSF_SIZE ' + str(psf_size) + ',' + str(
            psf_size) + ' -SAMPLE_MINSN ' + str(
            min_sn) + ' -CHECKPLOT_DEV NULL -CHECKPLOT_TYPE NONE -CHECKIMAGE_TYPE NONE -PSFVAR_DEGREES ' + str(
            degree) + ' -SAMPLE_VARIABILITY ' + str(variability) + ' -SAMPLE_FWHMRANGE ' + str(
            fwhm_range[0]) + ',' + str(fwhm_range[1]) + ' -PSF_SAMPLING ' + str(sampling) + ' -SAMPLE_MAXELLIP ' + str(
            max_elp) + ' -NTHREADS ' + str(nthread)
    else:
        checknames = fp_img + '_psffwhm.png,' + fp_img + '_psfelp.png'
        command = 'psfex ' + fp_psf + ' -c ' + pexfile + ' -PSF_SIZE ' + str(psf_size) + ',' + str(
            psf_size) + ' -SAMPLE_MINSN ' + str(min_sn) + ' -PSFVAR_DEGREES ' + str(
            degree) + ' -SAMPLE_VARIABILITY ' + str(variability) + ' -SAMPLE_FWHMRANGE ' + str(
            fwhm_range[0]) + ',' + str(fwhm_range[1]) + ' -PSF_SAMPLING ' + str(sampling) + ' -SAMPLE_MAXELLIP ' + str(
            max_elp) + ' -NTHREADS ' + str(
            nthread) + ' -CHECKPLOT_TYPE FWHM,ELLIPTICITY -CHECKPLOT_DEV PNG -CHECKPLOT_NAME ' + checknames
    print(command)
    os.system(command)
    return True


def photometry(ccd_id, dm, detect_thresh=1.0, analysis_thresh=1.0, clean='Y', nthread=1, head_zpt=None,
               checkfiles=['sky', 'seg', 'res', 'mod'], aper_cor=True, seeing=0.15, pixel_scale=0.075,
               phot_type='model'):
    """
    perform aperture and model photometry
    INPUT:
    fp_img: input image
    outdir: output directory
    detect_thresh: detecting threshold
    analysis_thresh: analysis threshold
    clean: clean the detections
    nthread: number of threads
    head_zpt: flux zeropoint
    aper_cor: wether do aperture correction
	seeing: in arcsec FWHM
	pixel_scale: pixel scale in arcsec
    phot_type: providing the photometry type, if psf, only perform aperture + psf photometry; if model, perform aperture+psf+model photometry

    checkfiles: check types
        sky: BACKGROUND
        seg: SEGMENTATION
        res: -MODEL
        mod: MODEL
    """

    # fitsname, _ = os.path.splitext(fp_img)
    # indpos = str.rfind(fitsname, '_img')
    # fitsname = fitsname[:indpos]
    # if outdir is not None:
    #     if not os.path.exists(outdir):
    #         os.mkdir(outdir)
    #     _, rootname = os.path.split(fitsname)
    #     fitsname = os.path.join(outdir, rootname)
    fp_img = dm.l1_ccd(ccd_id=ccd_id, post="img_L1.fits")
    fp_wht = dm.l1_ccd(ccd_id=ccd_id, post="wht_L1.fits")
    fp_flg = dm.l1_ccd(ccd_id=ccd_id, post="flg_L1.fits")
    fp_psf = dm.l1_ccd(ccd_id=ccd_id, post="psf.fits")
    fp_cat = dm.l1_ccd(ccd_id=ccd_id, post="cat.fits")

    sexfile = os.path.join(CONFIG_PATH, 'csst_phot.sex')
    covfile = os.path.join(CONFIG_PATH, 'default.conv')
    nnwfile = os.path.join(CONFIG_PATH, 'default.nnw')
    if phot_type == "psf":
        parfile = os.path.join(CONFIG_PATH, 'csst_psfphot.params')
    else:
        parfile = os.path.join(CONFIG_PATH, 'csst_modphot.params')
    # fp_psf = fp_img + '_psf.fits'
    head = fits.getheader(fp_img, 0)
    head1 = fits.getheader(fp_img, 1)
    if 'CCDZP' not in head1:
        print("No zeropoint in the header")
        return False
        # date=head['date-obs']
    # ccd=head['ccd_no']

    # get flag and weight images
    # flagfile,weightfile=get_flagweght(date,ccd)

    # rootname, _ = os.path.splitext(fp_img)
    # indpos = str.rfind(rootname, '_img')
    # flagfile = rootname[:indpos] + '_flg.fits'
    # weightfile = rootname[:indpos] + '_wht.fits'
    # fp_flg = fp_img.replace("_img", "_flg")
    # fp_wht = fp_img.replace("_img", "_wht")

    gain = head['exptime']

    types = {'sky': 'BACKGROUND', 'seg': 'SEGMENTATION', 'mod': 'MODELS', 'res': '-MODELS'}
    fp_cat = dm.l1_ccd(ccd_id, post="fluxadu.fits")
    if checkfiles is None:
        checktype = "NONE"
        checknames = "check.fits"
    else:
        ntypes = len(checkfiles)
        keys = types.keys()
        checktype = "'"
        checknames = "'"
        for i in range(ntypes):
            ikey = checkfiles[i]
            checktype += types[ikey] + " "
            checknames += dm.l1_ccd(ccd_id, post=ikey + ".fits") + " "# fp_img + "_" + ikey + ".fits "
        checktype += "'"
        checknames += "'"

    # checktype="'BACKGROUND -MODELS MODELS -PSFS PSFS'"
    # checknames='"'+fitsname+'-bkg.fits '+fitsname+'-submodel.fits '+fitsname+'-model.fits '+fitsname+'-subpsf.fits '+fitsname+'-psf.fits '+'"'
    command = 'sex ' + fp_img + ' -c ' + sexfile + ' -CATALOG_NAME ' + fp_cat + ' -FILTER_NAME ' + covfile + ' -STARNNW_NAME ' + nnwfile + ' -DETECT_THRESH ' + str(
        detect_thresh) + ' -ANALYSIS_THRESH ' + str(
        analysis_thresh) + ' -CLEAN ' + clean + ' -WEIGHT_IMAGE ' + fp_wht + ' -FLAG_IMAGE ' + fp_flg + ' -PIXEL_SCALE ' + str(
        pixel_scale) + ' -SEEING_FWHM ' + str(seeing) + ' -GAIN ' + str(
        gain) + ' -PARAMETERS_NAME ' + parfile + ' -CHECKIMAGE_TYPE ' + checktype + ' -CHECKIMAGE_NAME ' + checknames + ' -PSF_NAME ' + fp_psf + ' -NTHREADS ' + str(
        nthread)
    print(command)
    os.system(command)

    fluxadu = Table.read(fp_cat)
    # calibrate flux and aperture corrections
    fluxcalib = calibrate_fluxadu(fluxadu, head1, head_zpt=head_zpt)
    # plotname=fitsname
    plotname = None
    # aperture corrections
    if fluxcalib is not False:  # and aper_cor is True:
        fluxcor = magnitude_correction(fluxcalib, head1, elp_lim=0.5, sigma=2.5, plot_name=plotname, magerr_lim=0.1,
                                       sig_limit=0.08, aper_size=[3, 4, 5, 6, 8, 10, 13, 16, 20, 25, 30, 40],
                                       aper_ref=5, class_lim=0.5)
        # corfile = fp_img + '_cat.fits'
        corfile = dm.l1_ccd(ccd_id, post="cat.fits")
        fluxcor.write(corfile, format='fits', overwrite=True)
        ind_head1 = head1.index('NAXIS2') + 1
        ind_head = head.index('NEXTEND') + 1
        f = fits.open(corfile, mode='update')
        f[0].header.extend(head.cards[ind_head:], bottom=True)
        f[1].header.extend(head1.cards[ind_head1:], bottom=True)
        f.flush(output_verify='fix+warn')
        f.close()
    # if os.path.isfile(fp_cat):
    #     os.remove(fp_cat)


def calibrate_fluxadu(fluxadu, head, head_zpt=None):
    # F in nanomaggy = f in adu * 10**((c-22.5)/(-2.5)) where c is zeropoint
    fluxcalib = Table(fluxadu.copy())
    if head_zpt is not None:
        c = head_zpt
    else:
        if 'CCDZP' not in head or head['CCDZP'] <= 0:
            return False
        else:
            c = head['CCDZP']

    # calibrate flux and errors
    keys = ['APER', 'AUTO', 'HYBRID', 'PSF', 'MODEL', 'SPHEROID', 'DISK', 'MAX_MODEL', 'EFF_MODEL', 'MEAN_MODEL',
            'MAX_SPHEROID', 'EFF_SPHEROID', 'MEAN_SPHEROID', 'MAX_DISK', 'EFF_DISK', 'MEAN_DISK']
    fluxkeys = ['FLUX_' + ikey for ikey in keys]
    fluxerrkeys = ['FLUXERR_' + ikey for ikey in keys]

    nkeys = len(keys)
    colnames = fluxadu.colnames
    for i in range(nkeys):
        if fluxkeys[i] not in colnames: continue
        print("calibrate " + fluxkeys[i])
        flux = fluxadu[fluxkeys[i]]
        cflux = np.zeros_like(flux)
        cflux = flux * 10.0 ** ((c - 22.5) / (-2.5))
        fluxcalib[fluxkeys[i]] = cflux
        if i < 7:
            flux_error = fluxadu[fluxerrkeys[i]]
            cflux_error = np.zeros_like(flux)
            cflux_error = flux_error * 10.0 ** ((c - 22.5) / (-2.5))
            fluxcalib[fluxerrkeys[i]] = cflux_error
    return fluxcalib


def magnitude_correction(fluxcalib, head, plot_name=None, magerr_lim=0.05, elp_lim=0.3, sigma=3.0, iters=None,
                         sig_limit=0.1, aper_size=[3, 4, 5, 6, 8, 10, 13, 16, 20, 25, 30, 40], aper_ref=6,
                         class_lim=0.5, min_nstar=5):
    """
    correction for aperture and model magnitudes
    """
    low_errlim = 0.005
    fluxcor = fluxcalib.copy()

    elp = fluxcalib['ELLIPTICITY']
    flux = fluxcalib['FLUX_AUTO']
    fluxerr = fluxcalib['FLUXERR_AUTO']
    radius = fluxcalib['FLUX_RADIUS']
    automag, automagerr = fluxerr2magerr(flux, fluxerr)
    tmpmask = (automagerr < low_errlim) & (automagerr > 0)
    automagerr[tmpmask] = low_errlim
    mask_auto = (fluxcalib['FLAGS'] == 0) & (elp < elp_lim) & (automagerr < magerr_lim) & (automagerr > 0) & (
                automag < 90.0) & (radius > 1.0) & (fluxcalib['CLASS_STAR'] > class_lim)

    print('aperture correction for aperture magnitudes ...')
    apersize = np.array(aper_size)
    naper = len(apersize)
    nobj = len(fluxcalib)
    aperstr = ''
    for iaper in apersize:
        aperstr += str(iaper) + ','
    head["apersize"] = (aperstr, 'aperture radii in pixels')
    indref = np.argmin(np.abs(apersize - aper_ref))

    flux = fluxcalib['FLUX_APER']
    fluxerr = fluxcalib['FLUXERR_APER']
    apermag, apermagerr = fluxerr2magerr(flux, fluxerr)
    apmag8 = apermag[:, indref]
    apmag8err = apermagerr[:, indref]
    tmpmask = (apmag8err > 0) & (apmag8err < low_errlim)
    apmag8err[tmpmask] = low_errlim
    mask = mask_auto & (apmag8 > 0.0) & (apmag8 < 90.0) & (apmag8err < magerr_lim) & (apmag8err > 0)
    for i in range(naper):
        mask = mask & (apermag[:, i] < 90.0) & (apermag[:, i] > 0.0) & (apermagerr[:, i] < magerr_lim) & (
                    apermagerr[:, i] > 0)
    aperflag = True
    apercor = np.zeros(naper)
    apercor_std = np.zeros(naper)
    nstar_aper = 0
    if mask.sum() < min_nstar:
        aperflag = False
        print('not enough stars to do aperture magnitude correction')
    else:
        print('isolated stars: ', mask.sum())
        magdiff = -np.transpose(apermag[mask, :].transpose() - apmag8[mask])
        diff_masked = sigmaclip_limitsig(magdiff, sigma=sigma, maxiters=iters, axis=0)
        mask1 = mask
        mask = np.logical_not(np.any(diff_masked.mask, axis=1))
        nstar_aper = mask.sum()
        diff_masked = diff_masked[mask]
        for i in range(naper):
            weighterr = np.sqrt(apmag8err[mask1][mask] ** 2 + apermagerr[:, i][mask1][mask] ** 2)
            cor, _, corerr = weighted_mean(diff_masked[:, i], weighterr, weight_square=False)
            corerr /= np.sqrt(nstar_aper)
            apercor[i] = cor
            apercor_std[i] = corerr
        # apercor=np.median(diff_masked,axis=0)
        # apercor_std=np.std(diff_masked,axis=0)/np.sqrt(nstar_aper)
        print('aperture correction using stars:', nstar_aper)
        print(np.array([apercor, apercor_std]).transpose())
        if plot_name is not None:
            plt.figure(figsize=(9, 6))
            xsize = apersize.reshape(1, naper).repeat(len(magdiff), axis=0)
            plt.plot(xsize, magdiff, 'k.', markersize=1.0)
            xsize = apersize.reshape(1, naper).repeat(len(diff_masked), axis=0)
            plt.plot(xsize, diff_masked, 'b.', markersize=1.0)
            plt.plot(apersize, apercor, 'r*')
            plt.plot(apersize, apercor, 'r-')
            plt.xlabel('radius (pixels)')
            plt.ylabel('Mag_REF - Mag')
            plt.xlim([0, 43])
            plt.ylim([min(apercor) - 2, max(apercor) + 2])
            plt.title('aperture correction for aperture magnitudes')
            plt.savefig(plot_name + '-APERcor.png', format='png')
            plt.close()
    head['ns_aper'] = (nstar_aper, 'number of stars used in aperture correction')
    for i in np.arange(naper):
        cf = apercor[i]
        # ce=apercor_std[i]
        ce = 0.0
        head['apcor' + str(i)] = (cf, 'mag correction for aperture #{}'.format(i))
        head['aperr' + str(i)] = (ce, 'mag correction error for aperture #{}'.format(i))
        ff = fluxcalib['FLUX_APER'][:, i]
        fe = fluxcalib['FLUXERR_APER'][:, i]
        ff1 = ff * 10.0 ** (-cf / 2.5)
        fe1 = np.sqrt(10.0 ** (-2 * cf / 2.5) * fe ** 2 + ff ** 2 * 10.0 ** (-2 * cf / 2.5) * np.log(
            10.0) ** 2 / 2.5 ** 2 * ce ** 2)
        fluxcor['FLUX_APER'][:, i] = ff1
        fluxcor['FLUXERR_APER'][:, i] = fe1
        mag, magerr = fluxerr2magerr(ff1, fe1)
        apermag[:, i] = mag
        apermagerr[:, i] = magerr
    mag_col = Table.Column(apermag, 'MAG_APER')
    magerr_col = Table.Column(apermagerr, 'MAGERR_APER')
    fluxcor.add_column(magerr_col, fluxcor.index_column('FLUXERR_APER') + 1)
    fluxcor.add_column(mag_col, fluxcor.index_column('FLUXERR_APER') + 1)

    # automatic magnitude correction
    # mag correction for kron mag
    print('correct for kron magnitudes ...')
    # mask=mask_auto #automag<90.0
    # if aperflag:
    if True:
        mask = (fluxcor['FLAGS'] == 0)  # & (elp < elp_lim)
        k = fluxcor['KRON_RADIUS']
        A = fluxcor['A_IMAGE']
        E = fluxcor['ELLIPTICITY']
        B = A * (1.0 - E)
        # B=fluxcor['B_IMAGE']
        kronr = k * np.sqrt(A * B)
        ipt = UnivariateSpline(apersize, apercor, k=3, s=0)
        cf = ipt(kronr)
        ff = fluxcor['FLUX_AUTO']
        fe = fluxcor['FLUXERR_AUTO']
        ff1 = ff * 10.0 ** (-cf / 2.5)
        fluxcor['FLUX_AUTO'] = ff1
        mag, magerr = fluxerr2magerr(ff1, fe)
        mag_col = Table.Column(mag, 'MAG_AUTO')
        magerr_col = Table.Column(magerr, 'MAGERR_AUTO')
        fluxcor.add_column(magerr_col, fluxcor.index_column('FLUXERR_AUTO') + 1)
        fluxcor.add_column(mag_col, fluxcor.index_column('FLUXERR_AUTO') + 1)
        print('total kron sources:', mask.sum())
        if plot_name is not None:
            plt.figure(figsize=(9, 6))
            plt.plot(kronr[mask], apmag8[mask] - automag[mask], 'kx')
            plt.plot(apersize, apercor, 'g--')
            rr = np.linspace(0, 43, 100)
            plt.plot(rr, ipt(rr), 'r-')
            plt.xlabel('radius (pixels)')
            plt.ylabel('Mag_REF - Mag')
            plt.xlim([0, 43])
            plt.ylim([min(apercor) - 2, max(apercor) + 2])
            plt.title('aperture correction for automatic magnitudes')
            plt.savefig(plot_name + '-AUTOcor.png', format='png')
            plt.close()

    # other flux/mag correction
    print('other flux/mag correction: HYBRID PSF MODEL SPHEROID DISK ...')
    keys = ['HYBRID', 'PSF', 'MODEL', 'SPHEROID', 'DISK', 'MAX_MODEL', 'EFF_MODEL', 'MEAN_MODEL', 'MAX_SPHEROID',
            'EFF_SPHEROID', 'MEAN_SPHEROID', 'MAX_DISK', 'EFF_DISK', 'MEAN_DISK']
    fluxkeys = ['FLUX_' + ikey for ikey in keys]
    fluxerrkeys = ['FLUXERR_' + ikey for ikey in keys]
    magkeys = keys.copy()
    for i, ikey in enumerate(keys):
        if i < 5:
            magkeys[i] = 'MAG_' + ikey
        else:
            magkeys[i] = 'MU_' + ikey
    magerrkeys = ['MAGERR_' + ikey for ikey in keys]
    corkeys = ['HYBCOR', 'PSFCOR', 'MODCOR']
    corerrkeys = ['HYBERR', 'PSFERR', 'MODERR']
    nkeys = len(keys)
    colnames = fluxcalib.colnames
    for i in range(nkeys):
        if fluxkeys[i] not in colnames: continue
        print('correct ' + fluxkeys[i])
        flux = fluxcalib[fluxkeys[i]]
        if i < 5:
            fluxerr = fluxcalib[fluxerrkeys[i]]
        else:
            fluxerr = np.zeros_like(flux)
        if i < 3:
            kmag, kmagerr = fluxerr2magerr(flux, fluxerr)
            tmpmask = (kmagerr < low_errlim) & (kmagerr > 0)
            kmagerr[tmpmask] = low_errlim
            mask = mask_auto & (kmag < 90.0) & (apmag8 > 0) & (apmag8 < 90.0) & (kmag > 0) & (kmagerr > 0) & (
                        kmagerr < magerr_lim) & (apmag8err > 0) & (apmag8err < magerr_lim) & (
                               np.abs(kmag - apmag8) < 0.5)
            corflag = True
            if mask.sum() < min_nstar:
                print('not enough stars to correct ' + magkeys[i])
                cor = 0.0;
                corerr = 0.0;
                nstar_cor = 0
                corflag = False
            else:
                print('isolated stars for ' + magkeys[i] + ':', mask.sum())
                magdiff = apmag8[mask] - kmag[mask]
                diff_masked = sigmaclip_limitsig(magdiff, sigma=sigma, maxiters=iters, sig_limit=sig_limit)
                mask1 = np.logical_not(diff_masked.mask)
                nstar_cor = mask1.sum()
                diff_masked = magdiff[mask1]
                weighterr = np.sqrt(kmagerr[mask][mask1] ** 2 + apmag8err[mask][mask1] ** 2)
                cor, _, corerr = weighted_mean(diff_masked, weighterr, weight_square=False)
                corerr /= np.sqrt(nstar_cor)
                print('correction using stars:', nstar_cor)
                print([cor, corerr])
            head['ns_' + keys[i]] = (nstar_cor, 'number of stars used in ' + keys[i] + ' correction')
            head[corkeys[i]] = (cor, 'mag correction for ' + keys[i])
            head[corerrkeys[i]] = (corerr, 'mag correction error')
        cf = cor
        ce = corerr
        ff = flux
        fe = fluxerr
        ff1 = ff * 10.0 ** (-cf / 2.5)
        fe1 = np.sqrt(10.0 ** (-2 * cf / 2.5) * fe ** 2 + ff ** 2 * 10.0 ** (-2 * cf / 2.5) * np.log(
            10.0) ** 2 / 2.5 ** 2 * ce ** 2)
        fluxcor[fluxkeys[i]] = ff1
        if i < 5:
            fluxcor[fluxerrkeys[i]] = fe1
        mag, magerr = fluxerr2magerr(ff1, fe1)
        mag_col = Table.Column(mag, magkeys[i])
        if i < 5:
            magerr_col = Table.Column(magerr, magerrkeys[i])
            fluxcor.add_column(magerr_col, fluxcor.index_column(fluxerrkeys[i]) + 1)
            fluxcor.add_column(mag_col, fluxcor.index_column(fluxerrkeys[i]) + 1)
        else:
            fluxcor.add_column(mag_col, fluxcor.index_column(fluxkeys[i]) + 1)
        if corflag and plot_name is not None and i < 3:
            plt.figure(figsize=(9, 6))
            plt.plot(apmag8, apmag8 - kmag, 'kx')
            plt.plot(apmag8[mask], magdiff, 'bx')
            plt.plot(apmag8[mask][mask1], magdiff[mask1], 'r.')
            plt.plot([14, 25], [cor, cor], 'g')
            plt.xlabel('Mag_REF')
            plt.ylabel('Mag_REF-Mag_' + keys[i])
            plt.xlim([14, 25])
            plt.ylim([cor - 0.5, cor + 0.5])
            plt.title('aperture correction for ' + keys[i] + ' magnitudes')
            plt.savefig(plot_name + '-' + keys[i] + 'cor.png', format='png')
            plt.close()
    # add model magnitude and best magnitude
    # print('calculating the model and best magnitudes ...')
    # chi2psf=fluxcalib['CHI2_PSF']
    # chi2mod=fluxcalib['CHI2_MODEL']
    # nobj=len(fluxcor)
    # best_flux=np.zeros(nobj)
    # best_fluxerr=np.zeros(nobj)
    # best_mag=np.zeros(nobj)
    # best_magerr=np.zeros(nobj)
    # best_type=np.repeat(0,nobj)
    # mask=chi2psf < chi2mod
    # best_mag[mask]=fluxcor['MAG_PSF'][mask]
    # best_magerr[mask]=fluxcor['MAGERR_PSF'][mask]
    # best_flux[mask]=fluxcor['FLUX_PSF'][mask]
    # best_fluxerr[mask]=fluxcor['FLUXERR_PSF'][mask]
    # best_mag[~mask]=fluxcor['MAG_MODEL'][~mask]
    # best_magerr[~mask]=fluxcor['MAGERR_MODEL'][~mask]
    # best_flux[~mask]=fluxcor['FLUX_MODEL'][~mask]
    # best_fluxerr[~mask]=fluxcor['FLUXERR_MODEL'][~mask]
    # best_type[~mask]=1
    # bf_col=Table.Column(best_flux,'FLUX_CHI2MIN')
    # bferr_col=Table.Column(best_fluxerr,'FLUXERR_CHI2MIN')
    # bm_col=Table.Column(best_mag,'MAG_CHI2MIN')
    # bmerr_col=Table.Column(best_magerr,'MAGERR_CHI2MIN')
    # btype_col=Table.Column(best_type,'TYPE_CHI2MIN')
    # fluxcor.add_columns([bf_col,bferr_col,bm_col,bmerr_col,btype_col])
    return fluxcor


def rename_columns(cat, keys, rkeys):
    for i in range(len(keys)):
        cat.rename_column(keys[i], rkeys[i])


def rename_catalog():
    keys = ['', '']


def match_sdss(sexcat, sdsscat, outcat=None):
    sex = Table.read(sexcat, format='fits')
    sdss = Table.read(sdsscat, format='fits')
    rasex = np.asarray(sex['ALPHAWIN_J2000'])
    decsex = np.asarray(sex['DELTAWIN_J2000'])
    rasdss = np.asarray(sdss['ra'])
    decsdss = np.asarray(sdss['dec'])
    sexcoord = SkyCoord(rasex, decsex, frame="icrs", unit='deg')
    sdsscoord = SkyCoord(rasdss, decsdss, frame="icrs", unit='deg')
    idxsdss, idxsex, sep2d, dist3d = sexcoord.search_around_sky(sdsscoord, 1.0 * u.arcsec)
    if len(idxsex) == 0:
        return None
    filters = ['u', 'g', 'r', 'i', 'z']
    types = ['deV', 'exp', 'aper', 'petro', 'model', 'psf', 'cmodel']
    for ifilt in filters:
        if 'devFluxIvar_' + ifilt in sdss.colnames:
            sdss.rename_column('devFluxIvar_' + ifilt, 'deVFluxIvar_' + ifilt)
        for itype in types:
            fkey = itype + 'Flux_' + ifilt
            ferrkey = itype + 'FluxIvar_' + ifilt
            mkey = itype + 'Mag_' + ifilt
            merrkey = itype + 'MagErr_' + ifilt
            if itype == 'aper':
                fkey = itype + 'Flux7_' + ifilt
                ferrkey = itype + 'Flux7Ivar_' + ifilt
                mkey = itype + 'Mag7_' + ifilt
                merrkey = itype + 'Mag7Err_' + ifilt
            flux = sdss[fkey]
            fluxerr = np.sqrt(1.0 / sdss[ferrkey])
            mag, magerr = fluxerr2magerr(flux, fluxerr)
            sdss[fkey] = mag
            sdss[ferrkey] = magerr
            sdss.rename_column(fkey, mkey)
            sdss.rename_column(ferrkey, merrkey)
    keys = sdss.colnames
    rkeys = ['sdss_' + ikey for ikey in keys]
    rename_columns(sdss, keys, rkeys)
    number = Table.Column(np.zeros(len(sdss)).astype('int'), 'NUMBER')
    number[idxsdss] = sex['NUMBER'][idxsex]
    sdss.add_column(number)
    sexsdss = table.join(sex, sdss, join_type='left')
    sexsdss = sexsdss.filled(fill_value=0)
    if outcat: sexsdss.write(outcat, format='fits', overwrite=True)
    return sexsdss


def match_ps1(sexcat, ps1cat, outcat=None):
    sex = Table.read(sexcat, format='fits')
    ps1 = ps1cat
    if isinstance(ps1cat, str):
        ps1 = Table.read(ps1cat, format='fits')
    if 'flags' in ps1.colnames:
        ps1.remove_column('flags')
    rasex = np.asarray(sex['ALPHAWIN_J2000'])
    decsex = np.asarray(sex['DELTAWIN_J2000'])
    raps1 = np.asarray(ps1['RA'])
    decps1 = np.asarray(ps1['DEC'])
    sexcoord = SkyCoord(rasex, decsex, frame="icrs", unit='deg')
    ps1coord = SkyCoord(raps1, decps1, frame="icrs", unit='deg')
    idxps1, idxsex, sep2d, dist3d = sexcoord.search_around_sky(ps1coord, 1.0 * u.arcsec)
    if len(idxsex) == 0:
        return None
    keys = ps1.colnames
    rkeys = ['ps_' + ikey for ikey in keys]
    rename_columns(ps1, keys, rkeys)
    number = Table.Column(np.zeros(len(ps1)).astype('int'), 'NUMBER')
    number[idxps1] = sex['NUMBER'][idxsex]
    ps1.add_column(number)
    sexps1 = table.join(sex, ps1, join_type='left')
    sexps1 = sexps1.filled(fill_value=0)
    if outcat: sexps1.write(outcat, format='fits', overwrite=True)
    return sexps1


def do_phot(ccd_id, dm: CsstMscDataManager, stage=None, ):
    """
    get PSF and do photometry
    ifits: fits image
    outdir: output directory
    stage: psf or phot to run only get_psf or photometry
    """

    fp_psf = dm.l1_ccd(ccd_id, post="psf.fits")
    fp_cat = dm.l1_ccd(ccd_id, post="cat.fits")

    # get psfex
    time0 = time.time()

    time1 = time0
    if stage is None or stage == 'psf':
        # get_psf(fp_img,outdir=outdir,psf_size=71,degree=3,variability=0.3,fwhm_range=[1.5,20.0],max_elp=0.3,sampling=0,min_sn=10.0,detect_thresh=5.0,detect_minarea=5,back_size='400,400',check_plots=False,nthread=0,remove_polution=True,min_nstar=15,max_nstar=1500,class_star=0.7,match_dist=3.0,min_separation=20)
        get_psf(ccd_id, dm, psf_size=71, degree=3, variability=0.3, fwhm_range=[1.5, 20.0], max_elp=0.3,
                sampling=0, min_sn=5.0, detect_thresh=5.0, detect_minarea=5, back_size='400,400', check_plots=False,
                nthread=0, remove_polution=True, min_nstar=9, max_nstar=1500, class_star=0.7, match_dist=3.0,
                min_separation=20)
        time1 = time.time()
        psftime = time1 - time0
    time2 = time1
    if stage is None or stage == 'phot':
        if not os.path.isfile(fp_psf):
            print("PSF file not exist!")
        else:
            photometry(ccd_id, dm, detect_thresh=1.0, analysis_thresh=1.0, clean='Y', head_zpt=None,
                       checkfiles=['sky', 'seg'], nthread=0, phot_type='model')
            time2 = time.time()
            phottime = time2 - time1
            if os.path.isfile(fp_cat):
                f = fits.open(fp_cat, mode='update')
                if 'psftime' in vars().keys():
                    f[1].header['PSFTIME'] = (psftime, 'seconds')
                f[1].header['PHOTTIME'] = (phottime, 'seconds')
                f.flush(output_verify='fix+warn')
                f.close()
    print('Total time for ccd_id={}: {}'.format(ccd_id, time2 - time0))


def phot_main():
    parser = argparse.ArgumentParser(
        description='Photometric pipeline for CSST automatic, aperture, PSF and model photometry.',
        fromfile_prefix_chars='@')
    parser.add_argument('fits', metavar='FITS', type=str, nargs='*',
                        help='Calibrated CSST images to. @filelist to read from file.')
    parser.add_argument('-o', '--outdir', type=str, help='output directory; default same as the input fits image',
                        default=None)
    parser.add_argument('-s', '--stage', type=str, help='stage in the photometry (psf,phot)', default=None)
    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    args = parser.parse_args()

    # get list of fits, psfs, flags or weights
    fitslist = args.fits
    nfits = len(fitslist)
    stage = args.stage
    outdir = args.outdir

    # print help if no parameters provided
    if nfits == 0:
        parser.print_help()
        sys.exit(-1)

    # do photometry
    for i, ifits in enumerate(fitslist):
        start = time.time()
        # check file exist
        if not os.path.isfile(ifits):
            continue
        do_phot(ifits, outdir=outdir, stage=stage)


if __name__ == '__main__':
    phot_main()
