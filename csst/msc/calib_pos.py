import os
import time
from functools import partial
from multiprocessing import Pool
from subprocess import Popen
import joblib

import healpy as hp
import numpy as np
from astropy import table
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS

from .. import PACKAGE_PATH
from ..core.processor import CsstProcessor

CONFIG_PATH = PACKAGE_PATH + "/msc/pos_calib_config/"


class CsstProcMscPositionCalibration(CsstProcessor):

    def join_data(self, img_list, wht_list, flg_list):
        """
        Prepare data for running scamp
        Combine all image data, weight files, flag files to their one frame.

        Parameters
        ----------
        img_list:
            image files to join together, e.g.,MSC_210304093000_0000000_06_img.fits
        wht_list:
            weight files to join together, e.g.,MSC_210304093000_0000000_06_img.fits
        flg_list:
            flag files to join together, e.g.,e.g.,MSC_210304093000_0000000_06_flg.fits

        Returns
        -------
        The joined multi-extension file (not stacked), weight, flag files.
        e.g., MSC_210304093000_0000000_img.fits,MSC_210304093000_0000000_wht.fits, MSC_210304093000_0000000_flg.fits.

        """

        hdul_img = fits.HDUList()
        hdul_wht = fits.HDUList()
        hdul_flg = fits.HDUList()
        for i in range(len(img_list)):
            hdul_img.append(img_list[i][0])
            hdul_img.append(img_list[i][1])
            hdul_wht.append(wht_list[i][0])
            hdul_wht.append(wht_list[i][1])
            hdul_flg.append(flg_list[i][0])
            hdul_flg.append(flg_list[i][1])
        hdul_img.writeto(self.dm.l1_hardcode("combined_img.fits"), overwrite=True)
        hdul_wht.writeto(self.dm.l1_hardcode("combined_wht.fits"), overwrite=True)
        hdul_flg.writeto(self.dm.l1_hardcode("combined_flg.fits"), overwrite=True)

        return

    def run_sextractor(self, ccd_id):
        """
        Run sextractor

        Parameters 
        ----------
        ccd_id:
            file path of image, e.g., MSC_210304093000_0000000_06_img.fits

        Returns
        -------
        The photometric catalog, with position and flux, e.g., MSC_210304093000_0000000_06_img.acat
        """
        # input image
        fp_img = self.dm.l1_sci(ccd_id=ccd_id, suffix="img", ext="fits")
        # output catalog
        fp_cat = self.dm.l1_sci(ccd_id=ccd_id, suffix="img", ext="acat")

        # run sex and output the catalog
        config_sextractor = CONFIG_PATH + "new_csst_realtime.no.weight.sex"
        sex_comd1 = 'sex -c ' + config_sextractor + ' '
        sex_comd2 = fp_img + ' -CATALOG_NAME ' + fp_cat
        sex_comd3 = ' -PARAMETERS_NAME ' + CONFIG_PATH + 'csst_realtime.param' + ' -FILTER_NAME ' \
                    + CONFIG_PATH + 'csst_realtime.conv' + ' -STARNNW_NAME ' + CONFIG_PATH + 'csst_realtime.nnw'
        sex_comd = sex_comd1 + sex_comd2 + sex_comd3
        print(sex_comd)
        p = Popen(sex_comd, shell=True)
        p.wait()

    def combine_catalog(self):
        """
        Combine the sex catalog together

        Output
        ------
        The combined catalog,e.g., MSC_210304093000_0000000.acat.fits

        """

        print("combine catalog ")

        # output catalog
        output_catnm = self.dm.l1_hardcode(hdcd="combined_acat.fits", comment="combine_catalog")
        # fn = path_output + img_list[0][0].header['FILENAME'][0:-7]
        # output_catnm = str(fn + '.acat.fits')

        hdul = fits.HDUList()
        for ccd_id in self.dm.target_ccd_ids:
            this_fp = self.dm.l1_sci(ccd_id=ccd_id, suffix="img", ext="acat")
            print("processing {}".format(this_fp))
            for hdu in fits.open(this_fp):
                hdul.append(hdu)

        print("writing to {}".format(output_catnm))
        hdul.writeto(output_catnm, overwrite=True)

        return

    def run_scamp(self):
        """
        Run scamp

        Returns
        -------
        Image header updated with WCS keywords, MSC_210304093000_0000000.acat.head.
        """
        # image_prefix = img_list[0][0].header['FILENAME'][0:-7]

        fp_comb_cat = self.dm.l1_hardcode(hdcd="combined_acat.fits", comment="run_scamp")
        fp_ref_cat = self.dm.l1_hardcode(hdcd="gaia_ldac.fits", comment="run_scamp")
        fp_merge_cat = self.dm.l1_hardcode(hdcd="merged.cat", comment="run_scamp")
        fp_full_cat = self.dm.l1_hardcode(hdcd="full.cat", comment="run_scamp")

        config_scamp = CONFIG_PATH + "default2.scamp"
        scamp_comd = 'scamp ' + fp_comb_cat + \
                     ' -ASTREFCAT_NAME ' + fp_ref_cat + \
                     ' -MERGEDOUTCAT_NAME ' + fp_merge_cat + \
                     ' -FULLOUTCAT_NAME ' + fp_full_cat + \
                     ' -c ' + config_scamp
        print(scamp_comd)
        p = Popen(scamp_comd, shell=True)
        p.wait()

    @staticmethod
    def convert_hdu_to_ldac(hdu):
        """
        Convert an hdu table to a fits_ldac table (format used by astromatic suite)

        Parameters
        ----------
        hdu:
            `astropy.io.fits.BinTableHDU` or `astropy.io.fits.TableHDU`
            HDUList to convert to fits_ldac HDUList

        Returns
        -------
        tbl1:
            `astropy.io.fits.BinTableHDU`
            Header info for fits table (LDAC_IMHEAD)
        tbl2:
            `astropy.io.fits.BinTableHDU`
            Data table (LDAC_OBJECTS)
        """
        tblhdr = np.array([hdu[1].header.tostring()])
        col1 = fits.Column(name='Field Header Card', array=tblhdr, format='13200A')
        cols = fits.ColDefs([col1])
        tbl1 = fits.BinTableHDU.from_columns(cols)
        tbl1.header['TDIM1'] = '(80, {0})'.format(len(hdu[1].header))
        tbl1.header['EXTNAME'] = 'LDAC_IMHEAD'

        dcol = fits.ColDefs(hdu[1].data)
        tbl2 = fits.BinTableHDU.from_columns(dcol)
        tbl2.header['EXTNAME'] = 'LDAC_OBJECTS'
        return tbl1, tbl2

    def get_refcat(self, search_radius, silent=True):
        """
        Get reference catalog for scamp. The reference cat is GAIA EDR3.

        Parameters
        ----------
        image_prefix:
            a image to get its reference catalog, e.g.,MSC_210304093000_0000000_img.fits.
            Usually the center of the image is the wcs parameters CRVAL1,CRVAL1.
        search_radius:
            circle radius for searching, units: degree. e.g., 2 degree for a 1x1 deg^2 image.
            For large ccd size, use larger radius. csst, r=3 deg.
        path_gaia: directory of the reference catalog.

        Returns
        -------
        outcat:
            filename of the cross matched catalog.
            This catalog is used as a reference catalog for running scamp.
            e.g.,MSC_210304093000_0000000.gaialac.fits

        """
        # image_prefix = img_list[0][0].header['FILENAME'][0:-7]
        # fname = image_prefix + '_img.fits'
        # gaianame = image_prefix + '.gaia.fits'
        # gaialacnm = image_prefix + '.gaialac.fits'
        # outcat = gaianame

        path_gaia = self.dm.dir_pcref

        # combined image
        fname = self.dm.l1_hardcode("combined_img.fits", "get_refcat")
        gaianame = self.dm.l1_hardcode("gaia.fits", "get_refcat")
        gaialacnm = self.dm.l1_hardcode("gaia_ldac.fits", "get_refcat")
        # ref cat name
        # outcat = self.dm.l1_hardcode("ref.cat", "get_refcat")

        # read combined image
        hdu = fits.open(fname)
        header1 = hdu[0].header
        header2 = hdu[1].header
        deltatime = 0.00
        ra = float(header2['CRVAL1'])
        dec = float(header2['CRVAL2'])
        c = SkyCoord(ra, dec, unit=(u.deg, u.deg))
        print('ra, dec ra.deg dec.deg= ', ra, dec, c.ra.deg, c.dec.deg)
        ra = c.ra.deg
        dec = c.dec.deg
        c = SkyCoord(ra, dec, unit=(u.deg, u.deg))
        phi = c.ra.deg / (180. / np.pi)
        theta = (90. - c.dec.deg) / (180 / np.pi)
        vec = hp.ang2vec(theta, phi)
        list1 = hp.query_disc(nside=32, vec=vec, radius=np.radians(search_radius))
        if -1 in list1:
            list1.remove(-1)
        ipring = np.array(list1)
        pix = np.unique(ipring)
        npix = pix.size
        print(ipring, 'ipring', pix, 'pix', npix, 'npix')
        dt = np.dtype(
            [('ra', '>f8'), ('dec', '>f8'), ('ra_error', '>f8'), ('dec_error', '>f8'), ('phot_g_mean_mag', '>f8'),
             ('pmra', '>f8'), ('pmra_error', '>f8'), ('pmdec', '>f8'), ('pmdec_error', '>f8')])
        refcat = table.Table(dtype=dt)
        for i in pix:
            print('i= %5.5d' % i, path_gaia)
            fname = path_gaia + 'healpix-' + '%5.5d' % i + '.fits'
            print('fname=', fname)
            if not silent:
                print('Reading ', fname)
            d = table.Table.read(fname)
            refcat = [refcat, d]
            refcat = table.vstack(refcat, join_type='inner')
        refcat.rename_column('ra', 'X_WORLD')
        refcat.rename_column('dec', 'Y_WORLD')
        print('delta_time between obs_cat and ref_cat:', deltatime)

        mask = (refcat['pmdec'] != refcat['pmdec'])
        refcat['pmdec'][mask] = 0
        mask = (refcat['pmra'] != refcat['pmra'])
        refcat['pmra'][mask] = 0
        refcat['X_WORLD'] = refcat['X_WORLD'] + deltatime * refcat['pmra'] / np.cos(
            refcat['Y_WORLD'] / 180. * np.pi) / 3600.0 / 1000.0
        refcat['Y_WORLD'] = refcat['Y_WORLD'] + deltatime * refcat['pmdec'] / 3600.0 / 1000.0
        refcat['ra_error'] = refcat['ra_error'] / 1000.0 / 3600.0
        refcat['dec_error'] = refcat['dec_error'] / 1000.0 / 3600.0
        refcat.rename_column('ra_error', 'ERRA_WORLD')
        refcat.rename_column('dec_error', 'ERRB_WORLD')
        refcat.rename_column('phot_g_mean_mag', 'MAG')
        refcat.write(gaianame, format='fits', overwrite=True)

        # print('exist')
        hdu = fits.open(gaianame)
        hdu1 = self.convert_hdu_to_ldac(hdu)
        hdup = fits.PrimaryHDU()
        hdu = hdu1[0]
        tbhdu = hdu1[1]
        thdulist = fits.HDUList([hdup, hdu, tbhdu])
        # if os.path.exists(gaialacnm):
        #     os.remove(gaialacnm)
        thdulist.writeto(gaialacnm, overwrite=True)

        print('##################### end #####################')
        return

    @staticmethod
    def rewrite_wcs_head(head, output):
        """
        Rewrite the WCS head from Scamp to the standard fits header

        Parameters
        ----------
        head: scamp head file names

        Returns
        -------
        wcshead: a new head file in fits format

        """
        # wcshead = head + '.fits'
        f = open(head, 'r')
        f1 = open(output, 'w')
        a = ''
        i = 0
        for v in f.readlines():
            sp = ''
            asp = ''
            i += 1
            if len(v) <= 81:
                sp = ' ' * (81 - len(v))
            if 'END' in v:
                asp = ' ' * 80 * (36 - i % 36)
                i = i + (36 - i % 36)
                # print(i)
            a = a + v + sp + asp
        f1.write(a.replace('\n', ''))
        f1.close()
        f.close()
        return

    def check_astrometry(self):
        """
        Check position calibration quality

        Returns
        -------
        ccdraoff:
            ra difference between observed catalog and reference catalog
        ccddecoff:
            dec difference between observed catalog and reference catalog
        ccdraoff_med,ccddecoff_med:
            median of the ra or dec difference
        ccdraoff_rms,ccddecoff_rms:
            rms of the ra or dec difference
        """
        print('############## check the astrometry quality and save files ################')
        r1 = []
        d1 = []

        # image_prefix = (img_list[0][0].header)['FILENAME'][0:-7]
        # fn = path_output + image_prefix

        # combined catalog (head)
        wcshead_original = self.dm.l1_hardcode(hdcd="combined_acat.head", comment="check_astrometry")
        wcshead_rewritten = self.dm.l1_hardcode(hdcd="combined_acat.head.fits", comment="check_astrometry")
        # acat = self.dm.l1_hardcode(hdcd="combined_acat.fits", comment="check_astrometry")
        acat_change = self.dm.l1_hardcode(hdcd="combined_acat.change.fits", comment="check_astrometry")

        self.rewrite_wcs_head(wcshead_original, wcshead_rewritten)
        # acat = fits.open(fn + '.acat.fits')
        # acat_change = str(fn + '.acat.change.fits')
        # cat_suffix = '.acat'

        fp_scamp_coord = self.dm.l1_hardcode(hdcd="scamp_coord.txt", comment="check_astrometry")

        hdul = fits.HDUList()
        # for i in range(0, len(img_list)):
        for i, this_ccd_id in enumerate(self.dm.target_ccd_ids):
            wcshdr = fits.getheader(wcshead_rewritten, i, ignore_missing_simple=True)
            # read headers and change to RA---TPV,DEC--TPV for wcs_transfer package
            wcshdr['CTYPE1'] = 'RA---TPV'
            wcshdr['CTYPE2'] = 'DEC--TPV'
            w = WCS(wcshdr)
            # print(wcshdr)
            # cat_nm = path_output + img_list[i][0].header['FILENAME'] + cat_suffix
            cat_nm = self.dm.l1_sci(ccd_id=this_ccd_id, suffix="img", ext="acat")
            cat_i = fits.open(cat_nm)
            sexcat = cat_i[2].data
            ra_sex = sexcat['ALPHA_J2000']
            dec_sex = sexcat['DELTA_J2000']
            x = sexcat['XWIN_IMAGE']
            y = sexcat['YWIN_IMAGE']
            r, d = w.all_pix2world(x, y, 0)  # convert xwin,ywin to ra,de
            sexcat['ALPHA_J2000'] = r
            sexcat['DELTA_J2000'] = d
            cat_i[2].data = sexcat
            hdul.append(cat_i[0])
            hdul.append(cat_i[1])
            hdul.append(cat_i[2])
            r1 = np.hstack((r1, r))
            d1 = np.hstack((d1, d))
        obsc = SkyCoord(ra=r1 * u.degree, dec=d1 * u.degree)
        tmp_cat = np.zeros((len(obsc), 2))
        tmp_cat[:, 0] = obsc.ra
        tmp_cat[:, 1] = obsc.dec
        np.savetxt(fp_scamp_coord, tmp_cat, fmt="%.10f %.10f", delimiter="\n")
        hdul.writeto(acat_change, overwrite=True)  # update the cat with new ra,dec (from 1st scamp wcs.)

    def write_headers(self):
        """
        Wrtie history to header
        """
        # head_suffix = img_list[0][0].header['FILENAME'][0:-7] + '.acat.head.fits'
        head_suffix = self.dm.l1_hardcode(hdcd="combined_acat.head.fits", comment="write_headers")

        hdul2 = fits.open(head_suffix, ignore_missing_simple=True)
        # for i in range(0, len(img_list)):
        for i, this_ccd_id in enumerate(self.dm.target_ccd_ids):
            # fits_nm = img_list[i][0].header['FILENAME'] + '.head'
            fits_nm = self.dm.l1_sci(this_ccd_id, suffix="img", ext="head")

            hdul1 = fits.open(fits_nm, mode='update', ignore_missing_simple=True)
            hdr = hdul1[0].header
            hdr2 = hdul2[i].header
            hdr.extend(hdr2, unique=True, update=True)
            WCS_S = 0
            WCS_V = '2.0.4'
            WCS_P = 'default.scamp'
            WCS_TOL = time.strftime('%Y-%m-%d %H:%M:%S %p')
            hdr.set('WCS_S', '0', '0=done')
            hdr.set('WCS_V', WCS_V, 'Version of WCS calibration')
            hdr.set('WCS_P', WCS_P, 'Configure file name of WCS')
            hdr.set('WCS_TOL', WCS_TOL, 'Time of last wcs calibration')

        return

    def prepare(self, dm, n_jobs=-1): #, path_gaia, path_output, search_radius=2.0):
        self.dm = dm
        self.n_jobs = n_jobs
        # self.path_gaia = dm.dir_pcref
        # self.path_output = path_output
        # self.search_radius = search_radius

    def run(self, img_list, wht_list, flg_list):

        path_gaia = self.dm.dir_pcref
        path_output = self.dm.dir_l1
        # fn_list = [self.dm.l1_sci(ccd_id=_, suffix="img", ext="fits") for _ in self.dm.target_ccd_ids]

        print('preparing files for position calibration....')
        self.join_data(img_list, wht_list, flg_list)

        print('################## run sextractor ###################')
        joblib.Parallel(n_jobs=self.n_jobs)(joblib.delayed(self.run_sextractor)(_) for _ in self.dm.target_ccd_ids)
        # p = Pool()
        # prod_x = partial(self.run_sextractor, path_output=self.dm.dir_l1)
        # result = p.map(prod_x, self.dm.target_ccd_ids)
        # p.close()
        # p.join()
        print('################## sextractor done ###################')

        print('############### combine sextractor catalog ###############')
        self.combine_catalog()

        print('############### get reference catalog ###############')
        self.get_refcat(search_radius=2.0, silent=True)
        # gaianame gaialacnm
        # Popen('cp ' + refcat + ' ref.cat', shell=True)

        print('############### run scamp ##################')
        self.run_scamp()
        print('################ scamp done #################')

        print('Checking astrometry quality....')
        self.check_astrometry()

        print('################ updating headers.... #############')
        self.write_headers()

        print('#### Position calibration process done ####')

        return

    def cleanup(self):
        # clean up environment
        # image_prefix = img_list[0][0].header['FILENAME'][0:-7]
        # for i in range(0, len(img_list)):
        #     fn = img_list[i][0].header['FILENAME'] + '.acat'
        #     if os.path.isfile(path_output + fn):
        #         os.remove(path_output + fn)
        # if os.path.isfile(path_output + image_prefix + '.gaia.fits'):
        #     os.remove(path_output + image_prefix + '.gaia.fits')
        # if os.path.isfile(path_output + image_prefix + '.gaialac.fits'):
        #     os.remove(path_output + image_prefix + '.gaialac.fits')
        # if os.path.isfile(path_output + 'scamp.xml'):
        #     os.remove(path_output + 'scamp.xml')
        # if os.path.isfile(path_output + 'full_1.cat'):
        #     os.remove(path_output + 'full_1.cat')
        # if os.path.isfile(path_output + 'merged_1.cat'):
        #     os.remove(path_output + 'merged_1.cat')
        # if os.path.isfile(path_output + image_prefix + '_img.fits.back'):
        #     os.remove( path_output + image_prefix + '_img.fits.back')
        # if os.path.isfile(path_output + image_prefix + '_wht.fits'):
        #     os.remove(path_output + image_prefix + '_wht.fits')
        # if os.path.isfile(path_output + image_prefix + '_flg.fits'):
        #     os.remove(path_output + image_prefix + '_flg.fits')
        pass
