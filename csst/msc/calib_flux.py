import os
import time
from subprocess import Popen

import numpy as np
from astropy import table
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.stats import sigma_clip
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales

from .. import PACKAGE_PATH
from ..core.processor import CsstProcessor
from .data_manager import CsstMscDataManager

# Edited on Jun. 17, 2016
# add color term
__author__ = 'ZZM'


class CsstProcFluxCalibration(CsstProcessor):

    def __init__(self, dm : CsstMscDataManager, **kwargs):
        super().__init__(**kwargs)
        self.dm = dm

    def ps1_mags(self, cat, band, obs_region):
        # get ps1 color-corrected magnitude for specific filter
        # ref:https://desi.lbl.gov/trac/wiki/DecamLegacy/Reductions
        color_gi = cat['MEDIAN'][:, 0] - cat['MEDIAN'][:, 2]
        coeff = np.array([0.0, 0.0, 0.0, 0.0])
        if band == 'g':
            ifilt = 0
            magcut = 21.5
            coeff = np.array([0.00225738, -0.02543153, 0.10575001, 0.01077462])
        elif band == 'bokr' or band == 'r':
            ifilt = 1
            magcut = 20.0
            coeff = np.array([-0.00966711, 0.02808938, -0.07650502, 0.0024329])
        elif band == 'z':
            magcut = 20.0
            ifilt = 3
            coeff = np.array([-0.0094723, 0.03040996, -0.08774486, 0.01029133])
            # coeff=coeff=np.array([0.0,0.0,0.0,0.0])
        else:
            print(('no filter:', band))
            return
        colorterm = coeff[0] * color_gi ** 3 + coeff[1] * color_gi ** 2 + coeff[2] * color_gi + coeff[3]
        mags = cat['MEDIAN'][:, ifilt] + colorterm

        # select good phot_mag
        goodid = (cat['NMAG_OK'][:, 0] > 1) & (cat['NMAG_OK'][:, 2] > 1) \
                 & (cat['NMAG_OK'][:, ifilt] > 1) & (cat['STDEV'][:, ifilt] < 0.1) \
                 & (color_gi < 2.7) & (color_gi > 0.4) \
                 & (cat['MEDIAN'][:, ifilt] < magcut) & (cat['MEDIAN'][:, ifilt] > 16) \
                 & (cat['RA'] > obs_region[0]) & (cat['RA'] < obs_region[1]) \
                 & (cat['DEC'] > obs_region[2]) & (cat['DEC'] < obs_region[3])  # 2016.06.03
        print('starPS1.sum()=', goodid.sum())
        if goodid.sum() < 20: goodid = list(range(mags.size))
        mags = mags[goodid]
        magerr = cat['STDEV'][:, ifilt][goodid]
        ra = cat['RA'][goodid]
        dec = cat['DEC'][goodid]
        brmedian = np.median(color_gi[goodid])

        return mags, magerr, ra, dec, brmedian

    def read_ps1cat(self, ra, dec, outcat, path='/line12/Pan-STARRS/chunks-qz-star-v2/', silent=True):
        # Read a piece of the PS1 calibration catalog containing ra,dec
        from healpy import ang2pix
        c = SkyCoord(ra, dec, unit=(u.deg, u.deg))
        if isinstance(ra[0], str):
            c = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
        phi = c.ra.deg / (180. / np.pi)
        theta = (90. - c.dec.deg) / (180 / np.pi)
        ipring = ang2pix(32, theta, phi)
        pix = np.unique(ipring)
        npix = pix.size
        if npix > 8:
            print('too many healpix files:', npix, ' check image wcs!')
            return
        dt = np.dtype([('OBJ_ID', '>i8'), ('RA', '>f8'), ('DEC', '>f8'), ('NMAG_OK', '>i2', (5,)), \
                       ('STDEV', '>f4', (5,)), ('MEAN', '>f4', (5,)), ('MEDIAN', '>f4', (5,)), \
                       ('MEAN_AP', '>f4', (5,)), ('MEDIAN_AP', '>f4', (5,))])
        ps = table.Table(dtype=dt)
        for i in pix:
            fname = path + 'ps1-' + '%5.5d' % i + '.fits'
            if not silent: print(('Reading ', fname))
            if os.path.isfile(fname):
                d = table.Table.read(fname)
                ps = [ps, d]
                ps = table.vstack(ps)
        if outcat: ps.write(outcat, format='fits', overwrite=True)
        return ps

    def read_gaiacat(self, ra, dec, outcat, path='./', silent=True):
        # function: Read a piece of the GAIA calibration catalog containing ra,dec
        from healpy import ang2pix
        c = SkyCoord(ra, dec, unit=(u.deg, u.deg))
        if isinstance(ra[0], str):
            c = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
        phi = c.ra.deg / (180. / np.pi)
        theta = (90. - c.dec.deg) / (180 / np.pi)
        ipring = ang2pix(32, theta, phi)
        pix = np.unique(ipring)
        npix = pix.size
        if npix > 8:
            print('too many healpix files:', npix, ' check image wcs!')
            return
        dt = np.dtype([('source_id', '>i8'), ('ra', '>f8'), ('dec', '>f8'), ('ra_error', '>f4'), ('dec_error', '>f4'),
                       ('phot_g_mean_mag', '>f4'), ('phot_g_mean_flux_over_error', '>f4'), ('phot_g_n_obs', '>i2'),
                       ('phot_bp_mean_mag', '>f4'), ('phot_bp_mean_flux_over_error', '>f4'), ('phot_bp_n_obs', '>i2'),
                       ('phot_rp_mean_mag', '>f4'), ('phot_rp_mean_flux_over_error', '>f4'), ('phot_rp_n_obs', '>i2'),
                       ('astrometric_weight_al', '>f4'), ('astrometric_n_obs_al', '>i2'),
                       ('astrometric_n_good_obs_al', '>i2'), ('astrometric_excess_noise', '>f4'),
                       ('astrometric_excess_noise_sig', '>f4'), ('duplicated_source', '?'), ('ref_epoch', '>f4'),
                       ('parallax', '>f4'), ('parallax_error', '>f4'), ('pmra', '>f4'), ('pmra_error', '>f4'),
                       ('pmdec', '>f4'), ('pmdec_error', '>f4'), ('phot_variable_flag', '?')])
        ps = table.Table(dtype=dt)
        for i in pix:
            fname = path + 'chunk-' + '%5.5d' % i + '.fits'
            if not silent: print(('Reading ', fname))
            if os.path.isfile(fname):
                d = table.Table.read(fname)
                ps = [ps, d]
                ps = table.vstack(ps)
        if outcat: ps.write(outcat, format='fits', overwrite=True)
        return ps

    def gaia_mags(self, cat, band, obs_region):
        # get Gaia color-corrected magnitude for specific filter
        # color transp. parameters
        filterlist = ['N', 'u', 'g', 'r', 'i', 'z', 'y']
        # ffpar=np.array([[-0.12372707,0.87415984,2.68152378,0.46688053],\
        #    [ 0.0715824,0.42018862,1.58311364,0.30184685],\
        #    [-0.02193249,0.22376644,0.51403047,-0.04977904],\
        #    [ 0.0850583,  0.0045054, -0.17416451, 0.02696644],\
        #    [ 0.03638692,  0.0513301, -0.58839035, 0.02328468],\
        #    [-0.14492398, 0.26053385, -0.82220974, -0.0115362 ],\
        #    [-0.27416685, 0.36465315, -0.89341154, -0.02812699]])
        ffpar = np.array([[-1.92165818e+00, 7.88689593e+00, -1.60028226e+00, 1.53590717e+00], \
                          [-2.00544954e+00, 7.17304891e+00, -4.92617803e+00, 3.05285639e+00], \
                          [2.36560800e-02, 1.78153510e-01, 4.17186700e-01, 2.06523985e+00], \
                          [8.04794200e-02, -6.17626400e-02, -1.81644150e-01, 1.89622143e+00], \
                          [2.01861000e-03, 1.22630040e-01, -6.96895960e-01, 9.65536332e-01], \
                          [-4.75608000e-02, 1.86804250e-01, -9.25372450e-01, 2.67927720e+00], \
                          [-7.93956600e-02, 2.28310630e-01, -1.02675552e+00, 1.76769826e+00]])
        coeff = ffpar[filterlist.index(band), :]
        # print('filterlist.index(band):', filterlist.index(band))
        color_br = cat['phot_bp_mean_mag'] - cat['phot_rp_mean_mag']  ##- 0.3161466496 #VEGAMAG -> AB
        colorterm = coeff[0] * color_br ** 3 + coeff[1] * color_br ** 2 + coeff[2] * color_br + coeff[3]
        gaia_G = cat['phot_g_mean_mag']  ## +0.1001113078 #VEGAMAG -> AB
        mags = gaia_G + colorterm
        magcut = 19.0

        # select good phot_mag
        goodid = (cat['phot_g_n_obs'] > 3) & (cat['phot_bp_n_obs'] > 1) & (cat['phot_rp_n_obs'] > 1) \
                 & (color_br < 2.5) & (color_br > 0.5) \
                 & (gaia_G < magcut) & (gaia_G > 16) \
                 & (cat['ra'] > obs_region[0]) & (cat['ra'] < obs_region[1]) \
                 & (cat['dec'] > obs_region[2]) & (cat['dec'] < obs_region[3])
        print('starGAIA.sum()=', goodid.sum())
        if goodid.sum() < 20: goodid = list(range(mags.size))
        mags = mags[goodid]
        magerr = mags * 0.0
        ra = cat['ra'][goodid]
        dec = cat['dec'][goodid]
        gimedian = np.median(color_br[goodid])
        refc = SkyCoord(ra, dec, unit=(u.deg, u.deg), frame='fk5')

        return mags, magerr, gimedian, refc

    def read_inputcat(self, image, outcat='refcat.fits', refdir='', silent=True):
        # imname=os.path.split(image)[-1]
        if refdir == '':
            refdir = os.path.split(image)[0]
        # catname='MSC_'+imname[7:]+'.cat'
        catid = image[image.rfind('_') - 10:image.rfind('_')]
        catname = 'MSC_210525120000_' + catid + '.cat'
        cat = os.path.join(refdir, catname)
        data = table.Table.read(cat, format='ascii')
        data.write(outcat, format='fits', overwrite=True)
        return data

    def read_inputcat_(self, ccd_id):
        cat = self.dm.l0_cat(ccd_id)
        data = table.Table.read(cat, format='ascii')
        return data

    def input_mags(self, cat, usepixcood=False):
        goodid = (cat['mag'] < 20) & (cat['mag'] > 16)
        n = 0
        while (n < 4):
            goodid = (cat['mag'] < 20 + n) & (cat['mag'] > 16)
            n = n + 0.5
            if goodid.sum() > 20:
                n = n + 5
        mags = cat['mag'][goodid]
        magerr = mags * 0.0
        if usepixcood:
            x = cat['xImage'][goodid]
            y = cat['yImage'][goodid]
            z = np.zeros_like(x)
            refc = SkyCoord(x=x, y=y, z=z, unit=u.pix, representation_type='cartesian')
        else:
            ra = cat['ra'][goodid]
            dec = cat['dec'][goodid]
            refc = SkyCoord(ra, dec, unit=(u.deg, u.deg), frame='fk5')
        brmedian = 0.0

        return mags, magerr, brmedian, refc

    def rewrite_sex_cat(self, cat, workdir=''):
        """ split combined_acat.fits """
        data = fits.open(cat)
        # if workdir == '':
        #     workdir = os.path.split(cat)[0]
        index = []
        if (np.size(data) % 3 == 0) & (np.size(data) != 0):
            for i in range(int(np.size(data) / 3)):
                # print(i)
                hdr = data[3 * i + 1].data
                filename = str(hdr)[str(hdr).rfind('FITSFILE') + 11:
                                    str(hdr).rfind(".fits' / File name of the analyse")]
                subfile = os.path.join(workdir, filename + '.acat')
                if not os.path.exists(subfile):
                    hdu = fits.HDUList([data[0], data[3 * i + 1], data[3 * i + 2]])
                    hdu.writeto(subfile, overwrite=True)
                index.append(filename)
        return index

    def rewrite_sex_cat_(self, cat):
        """ split combined_acat.fits into {ccd_id}_img.acat """
        data = fits.open(cat)
        nhdu = len(data)
        # if workdir == '':
        #     workdir = os.path.split(cat)[0]
        print("@rewrite_sex_cat_: {} ({} HDUs)...".format(cat, nhdu))
        # index = []
        # if (np.size(data) % 3 == 0) & (np.size(data) != 0):
        #     for i in range(int(np.size(data) / 3)):
        #         # print(i)
        #         hdr = data[3 * i + 1].data
        #         filename = str(hdr)[str(hdr).rfind('FITSFILE') + 11:
        #                             str(hdr).rfind(".fits' / File name of the analyse")]
        #         subfile = os.path.join(workdir, filename + '.acat')
        #         if not os.path.exists(subfile):
        #             hdu = fits.HDUList([data[0], data[3 * i + 1], data[3 * i + 2]])
        #             hdu.writeto(subfile, overwrite=True)
        #         index.append(filename)
        # return index
        for i, ccd_id in enumerate(self.dm.target_ccd_ids):
            hdu = fits.HDUList([data[0], data[3 * i + 1], data[3 * i + 2]])
            hdu.writeto(self.dm.l1_sci(ccd_id, suffix="img", ext="acat"), overwrite=False)
        return

    def split_wcs_head(self, wcshead, im_index=[], workdir=''):
        if len(im_index) == 0:
            return

        head = fits.open(wcshead, ignore_missing_simple=True)
        if np.size(head) != np.size(im_index):
            return

        for i, filename in enumerate(im_index):
            pathsplit = os.path.split(filename)
            if workdir == '':
                workdir = pathsplit[0]
            headname = os.path.join(workdir, pathsplit[1] + '.whead.fits')
            wheader = head[i].header
            prihdu = fits.PrimaryHDU(header=wheader)
            prihdu.writeto(headname, overwrite=True)

    def split_wcs_head_(self, wcshead):
        """ split combined_acat.head.fits into {ccd_id}_whead.fits """
        head = fits.open(wcshead, ignore_missing_simple=True)
        nhdu = len(head)
        print("@split_wcs_head_: {} ({} HDUs) to ccd_id_img.whead.fits".format(wcshead, nhdu))

        # for i, filename in enumerate(im_index):
        #     pathsplit = os.path.split(filename)
        #     if workdir == '':
        #         workdir = pathsplit[0]
        #     headname = os.path.join(workdir, pathsplit[1] + '.whead.fits')
        #     wheader = head[i].header
        #     prihdu = fits.PrimaryHDU(header=wheader)
        #     prihdu.writeto(headname, overwrite=True)
        for i, ccd_id in enumerate(self.dm.target_ccd_ids):
            prihdu = fits.PrimaryHDU(header=head[i].header)
            prihdu.writeto(self.dm.l1_sci(ccd_id, suffix="img", ext="whead.fits"), overwrite=True)
        return

    def run_sextractor(self, image, cat):
        image0 = image
        header0 = fits.getheader(image, 0)
        header = fits.getheader(image, 1)

        CONFIG_PATH = PACKAGE_PATH + '/msc/flux_calib_config/'
        # sexfile='/home/zhouzm/data/csst/code/cali_default.sex'
        sexfile = CONFIG_PATH + 'cali_default.sex'
        sex_comd3 = ' -PARAMETERS_NAME ' + CONFIG_PATH + 'cali_default.param' + ' -FILTER_NAME ' + CONFIG_PATH + 'cali_default.conv' + ' -STARNNW_NAME ' + CONFIG_PATH + 'cali_default.nnw'
        exptime = header0['EXPTIME']
        gain = header['GAIN1'] * exptime
        # saturate=header['SATURATE']/exptime
        sex_zp = '0.0'
        imwht = image.replace('_img', '_wht')
        # apture='4,9,15,21,27,36,39,44'
        apture = '10'
        sex = 'sex ' + image + ' -c ' + sexfile + ' -CATALOG_NAME ' + cat + \
              ' -PHOT_APERTURES ' + apture + ' -MAG_ZEROPOINT ' + sex_zp + \
              ' -WEIGHT_IMAGE ' + imwht + ' -GAIN ' + str(gain) + \
              ' -GAIN_KEY ' + 'abcdefg' + sex_comd3
        # sex WFST_000001_r_210103150723_9.fits -c cali_default.sex -CATALOG_NAME ./WFST_000001_r_210103150723_9_calibphot.fits -PHOT_APERTURES 36 -MAG_ZEROPOINT 0.0 -WEIGHT_IMAGE WFST_000001_r_210103150723_9.wht.fits -GAIN 2.047026888573 -GAIN_KEY abcdefg
        print(sex)
        p = Popen(sex, shell=True)
        p.wait()

        return cat

    def prepare(self, this_ccd_id, wcsdir, workdir, usewcsresult=False, newcat=False):
        # make calibration files: Obs. cat (SExtractor), PSF model (PSFex)
        # sex_zp=str(25.5+2.5*np.log10(exptime))
        # if not os.path.exists(workdir):
        #     os.mkdir(workdir)
        image = self.dm.l1_sci(this_ccd_id, suffix="img", ext="fits")

        # imag_file_ename = os.path.split(image)[-1]
        imag_file_ename = self.dm.l1_sci(this_ccd_id, suffix="img", ext="fits")
        # imname = imag_file_ename[0:imag_file_ename.find('.fits')]
        # imname0 = imag_file_ename[0:imag_file_ename.find('_img') - 3]
        # wcscat0 = os.path.join(wcsdir, imname0 + '.acat.fits')
        wcscat0 = self.dm.l1_hardcode(hdcd="combined_acat.fits", comment="prepare")
        # wcscat1 = os.path.join(wcsdir, imname + '.acat')
        wcscat1 = self.dm.l1_sci(this_ccd_id, suffix="img", ext="acat")
        # wcshead0 = os.path.join(wcsdir, imname0 + '.acat.head.fits')
        wcshead0 = self.dm.l1_hardcode(hdcd="combined_acat.head.fits") # scamp output head, converted to fits
        # wcshead0txt = wcshead0[:wcshead0.rfind('.fits')]
        wcshead0txt = self.dm.l1_hardcode(hdcd="combined_acat.head") # scamp output head, txt format
        # wcshead1 = os.path.join(wcsdir, imname + '.acat.head.fits')
        wcshead1 = self.dm.l1_sci(this_ccd_id, suffix="img", ext="acat.head.fits")

        # cat = os.path.join(workdir, imname + '.acat')
        cat = self.dm.l1_sci(this_ccd_id, suffix="img", ext="acat") # sex output catalog
        # ref = os.path.join(workdir, imname + '.rcat')
        ref = self.dm.l1_sci(this_ccd_id, suffix="img", ext="rcat") # reference catalog (to be generated)
        # wcshead2 = os.path.join(workdir, imname + '.whead.fits')
        wcshead2 = self.dm.l1_sci(this_ccd_id, suffix="img", ext="whead.fits") # wcs head (to be generated)

        # cali_ref='GAIA' #calibration reference data

        # run SExtractor
        im_index = []
        if (os.path.isfile(cat)) & (newcat):
            os.remove(cat)

        # if splitted sex cat does not exist or combined, no matter splitted or combined,
        # if (not os.path.isfile(wcscat1)) or (os.path.isfile(wcscat0) & usewcsresult):
        #     im_index = self.rewrite_sex_cat(wcscat0, workdir)
        #     if os.path.isfile(wcshead0):
        #         self.split_wcs_head(wcshead0, im_index, workdir)
        #     elif os.path.isfile(wcshead0txt):
        #         wcshead0 = self.rewrite_wcs_head(wcshead0txt)
        #         self.split_wcs_head(wcshead0, im_index, workdir)

        # wcscat1: {ccd_id}_img.acat
        # wcscat0: combined_acat.fits
        # usewcsresult = True --> for test
        # if (not os.path.isfile(wcscat1)) or (os.path.isfile(wcscat0) & usewcsresult):
        #     im_index = self.rewrite_sex_cat_(wcscat0)
        #     if os.path.isfile(wcshead0):
        #         self.split_wcs_head(wcshead0, im_index, workdir)
        #     elif os.path.isfile(wcshead0txt):
        #         wcshead0 = self.rewrite_wcs_head(wcshead0txt)
        #         self.split_wcs_head(wcshead0, im_index, workdir)

        # if sex cat for this image does not exist, run sex (deprecated)
        # if (not os.path.isfile(wcscat1)):
        #     self.run_sextractor(image, wcscat1)

        # image: {ccd_id}_img.fits
        # wcshead1: img.acat.head.fits, scamp wcs head for each image
        # wcshead2: img.whead.fits
        header = self.combine_head_(this_ccd_id, prime=False)

        return header

    # def getebv(image):
    #    ebv=0.0
    #    header=fits.getheader(image, 0)
    #    if image[-3:]=='.fz':
    #        header=fits.getheader(image, 1)
    #    obj=header['OBJECT']
    #    band=header['FILTER']
    #    tilesfile='code/latest/bass-tiles.fits'
    #    tiles=table.Table.read(tilesfile)
    #    tid=(tiles['TID']==str(obj)[:4])
    #    if tid.sum()==1: ebv=tiles['EBV'][tid].data[0]
    #    #header.set('ebv',round(ebv,4),'E(B-V) from SFD1998')
    #    return ebv

    def getmlim(self, fwhm=0.15, avsky=7.0, rdnoise=5.0, zpt=25.8, ebv=0.0, filter='g'):
        # E(B-V): 3.995, 3.214, 2.165, 1.592, 1.211, 1.064 for ugrizY decals
        # https://github.com/dstndstn/tractor/blob/master/tractor/sfd.py#L10
        # if filter=='g':k=3.214; pixsize=0.33
        k = 0.0
        pixsize = 0.075
        snr = 5.0
        # nea=(((fwhm/2.35)**2*4*np.pi)**(1/1.15)+(8.91*(0.45/pixsize)**2)**(1/1.15))**1.15
        # noise=np.abs(avsky)*nea*exptime+7.3**2*nea
        # flim=0.5/exptime*(snr**2+np.sqrt(snr**4+4*snr**2*noise))

        noise = np.sqrt(avsky * np.pi * (fwhm / 2.35) ** 2 + rdnoise ** 2)
        flim = noise * snr
        mlim = -2.5 * np.log10(flim) + zpt - k * ebv
        return mlim

    def getfwhm(self, fwhmsex, ellip, obsme, obsflags):
        nanid = np.isnan(fwhmsex)
        if nanid.sum() > 0:
            fwhmsex[np.isnan(fwhmsex)] = 0.0
        starid1 = (ellip <= 0.15) & (obsme <= 0.15) & (obsme > 0) & (obsflags < 1) & (fwhmsex > 0)
        starid2 = fwhmsex[starid1] < 5.0  # temp. value, 20220314
        starid = starid1.nonzero()[0][starid2]
        if starid2.sum() < 10:
            starid = list(range(fwhmsex.size))

        fwhmsub = fwhmsex[starid]
        fwhmid = sigma_clip(fwhmsub, sigma=2.5)
        f1, f2 = np.histogram(fwhmid.data[~fwhmid.mask])
        fbin = np.where(f1 == f1.max())
        f2l = f2[np.max([fbin[0][0] - 1, 0])]
        f2h = f2[np.min([fbin[0][0] + 1, 10])]
        id = np.where((fwhmid.data[~fwhmid.mask] > f2l) & (fwhmid.data[~fwhmid.mask] < f2h))
        fwhm = np.median(fwhmid.data[~fwhmid.mask][id])
        if np.isnan(fwhm):
            fwhm = np.mean(fwhmsex)
        # print 'fwhm=',fwhm,f1,f2,fbin, f2l,f2h
        # print id,fwhmid.data[~fwhmid.mask][id]
        # print fwhmid, fwhmsex
        return fwhm

    def makedatedir(self, path):
        index = path.rfind('/')
        rootpath = path[0:path.rfind('/', 0, index)]
        walk = os.walk(rootpath)
        datedir = walk.next()[1]
        for i in datedir:
            if not os.path.exists(i):
                os.mkdir(i)
                print(('makedir ', i))
                ##########################

    def rewrite_wcs_head(self, head):
        # rewrite the WCS head from Scamp to the standard fits header.
        wcshead = head + '.fits'
        f = open(head, 'r')
        f1 = open(wcshead, 'w')
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
        return wcshead

    def combine_head(self, image, wcshead1='', wcshead2='', prime=False):
        # combine image head and wcs head keywords
        inst_head_file = image[0:image.find('.fits')] + '.head'
        if os.path.isfile(inst_head_file):
            h0 = fits.getheader(inst_head_file, ignore_missing_simple=True)
        else:
            h0 = fits.getheader(image, 1, ignore_missing_simple=True)

        if os.path.isfile(wcshead1):
            hw = fits.getheader(wcshead1, ignore_missing_simple=True)
        elif os.path.isfile(wcshead2):
            hw = fits.getheader(wcshead2, ignore_missing_simple=True)

        if prime:
            hprime = fits.getheader(image, 0, ignore_missing_simple=True)
            h0.extend(hprime, unique=True, update=True)

        try:
            h0.extend(hw, unique=True, update=True)
        except:
            print(image, ': no wcs header file')

        # update the keywords CTYPE1 & CTYPE2 in header
        if ('CTYPE1' in list(h0.keys())) & ('PV1_0' in list(h0.keys())):
            h0['CTYPE1'] = 'RA---TPV'
            h0['CTYPE2'] = 'DEC--TPV'

        return h0

    def combine_head_(self, ccd_id, prime=False):
        """ combine image head and wcs head keywords """
        # inst_head_file = image[0:image.find('.fits')] + '.head'
        inst_head_file = self.dm.l1_sci(ccd_id, suffix="img", ext="head")
        # if os.path.isfile(inst_head_file):
        h0 = fits.getheader(inst_head_file, ignore_missing_simple=True)
        # else:
        #     h0 = fits.getheader(image, 1, ignore_missing_simple=True)

        # if os.path.isfile(wcshead1):
        #     hw = fits.getheader(wcshead1, ignore_missing_simple=True)
        # elif os.path.isfile(wcshead2):
        wcshead2 = self.dm.l1_sci(ccd_id, suffix="img", ext="whead.fits")
        hw = fits.getheader(wcshead2, ignore_missing_simple=True)

        # if prime:
        #     # if True, include primary header
        #     hprime = fits.getheader(image, 0, ignore_missing_simple=True)
        #     h0.extend(hprime, unique=True, update=True)

        # try:
        h0.extend(hw, unique=True, update=True)
        # except:
        #     print(image, ': no wcs header file')

        # update the keywords CTYPE1 & CTYPE2 in header
        if ('CTYPE1' in list(h0.keys())) & ('PV1_0' in list(h0.keys())):
            h0['CTYPE1'] = 'RA---TPV'
            h0['CTYPE2'] = 'DEC--TPV'
            success_scamp = True
        else:
            success_scamp = False

        h0["PCSCAMP"] = success_scamp
        return h0

    def match_calib(self, obsc, refc, obsm, refm, obsme, refme, obsflags, fwhmsex=np.array([])):
        # matching obs and ref catalogs
        '''
        coeff0=self.match_calib(obsc,refc,obsm,refm)
        obsc,refc: obs and ref coordinates
        obsm,refm: obs and ref magnitudes
        '''
        coeff0 = cstd = cobsm = crefm = -1
        ccdraoff = ccddecoff = -1
        csize = 0
        fwhm = -1
        if obsm.size > 0 and refm.size > 0:
            obsid = (obsme < 0.2)  # 0.2->0.1 #20170518
            obsc = obsc[obsid]
            obsm = obsm[obsid]
            obsme = obsme[obsid]
            obsflags = obsflags[obsid]
            idx, d2d, d3d = obsc.match_to_catalog_sky(refc)
            ref_uid = np.unique(idx)
            obs_uid = np.full_like(ref_uid, -1)
            tmpj = -1
            for i in ref_uid:
                tmpj = tmpj + 1
                iid = (idx == i)
                iiid = (d2d.deg[iid] == d2d.deg[iid].min())
                obs_uid[tmpj] = iid.nonzero()[0][iiid.nonzero()[0]][0]

            uidlim = d2d[obs_uid].arcsecond < 0.5  # set match radius=0.5 arcsec
            if uidlim.sum() > 0:
                # csize=uidlim.sum()
                obs_uidlim = obs_uid[uidlim]
                ref_uidlim = ref_uid[uidlim]

                ccdraoff = np.median(
                    (obsc[obs_uidlim].ra - refc[ref_uidlim].ra).arcsec * np.cos(obsc[obs_uidlim].dec.deg * np.pi / 180))
                ccddecoff = np.median((obsc[obs_uidlim].dec - refc[ref_uidlim].dec).arcsec)
                # if fwhmsex.size >1:
                #    fwhm=np.median(fwhmsex[obs_uidlim])
                # matched magnitude

                # select good photometric value
                goodphot = (obsme[obs_uidlim] < 0.1) & (obsflags[obs_uidlim] < 2)
                csize = goodphot.sum()
                if goodphot.sum() <= 5:
                    goodphot = (obsme[obs_uidlim] < 0.2) & (obsflags[obs_uidlim] < 5)
                if goodphot.sum() >= 1:
                    obsm = obsm[obs_uidlim][goodphot]
                    refm = refm[ref_uidlim][goodphot]
                    obsme = obsme[obs_uidlim][goodphot]
                    refme = refme[ref_uidlim][goodphot]
                    combme = np.sqrt(obsme ** 2 + refme ** 2)
                    clip = sigma_clip(refm - obsm, sigma=2.5)  # sigma=2.5
                    coeff0 = np.median(clip.data[~clip.mask])
                    # wht=1./combme
                    # coeff0=np.average(refm-obsm,weights=wht)
                    cstd = np.std(clip.data[~clip.mask])
                    cobsm = obsm[~clip.mask]
                    crefm = refm[~clip.mask]
                    # if fwhmsex.size >1:
                    #    fwhm=np.median(fwhmsex[obsid][obs_uidlim][goodphot])
                    # print 'fwhm=',fwhm*0.26
        # print(coeff0,cstd,csize,cobsm,crefm,ccdraoff,ccddecoff)
        # print('crefm.min,max=',refm.min(),refm.max())
        return coeff0, cstd, csize, cobsm, crefm, ccdraoff, ccddecoff, fwhm

    def calib(self, this_ccd_id, addhead=False, morehead=True, plot=False, nodel=True, update=False, upcat=True):

        # set directories in case of necessary
        refdir = self.dm.dir_l0
        wcsdir = L1dir = workdir = self.dm.dir_l1

        # read image data
        fp_img = self.dm.l1_sci(this_ccd_id, suffix="img", ext="fits")
        fp_wht = self.dm.l1_sci(this_ccd_id, suffix="wht", ext="fits")
        fp_flg = self.dm.l1_sci(this_ccd_id, suffix="flg", ext="fits")
        imgdata = fits.open(fp_img)
        whtdata = fits.open(fp_wht)
        flgdata = fits.open(fp_flg)

        # imag_file_ename = os.path.split(image)[-1]
        imag_file_ename = os.path.split(self.dm.l1_sci(this_ccd_id, "img", "fits"))[-1]
        imname = imag_file_ename[0:imag_file_ename.rfind('.')]
        # image_head_file = os.path.join(workdir, imname + '.head')
        image_head_file = self.dm.l1_sci(this_ccd_id, suffix="img", ext="head")

        # psname = os.path.join(workdir, imname + '_calib.png')
        psname = self.dm.l1_hardcode(hdcd='flux_calib.png', comment="calib")
        # if plot and os.path.isfile(psname) and (not update):
        #     print((psname, ' is existed, pass'))
        #     return
        ckf = self.dm.l1_hardcode(hdcd='checkwcs.ls', comment="calib")
        ckim = self.dm.l1_hardcode(hdcd='checkim.ls', comment="calib")

        # get astrometry head
        # header=self.combine_head(image)

        cali_ref = 'GAIA'  # calibration reference data
        # Get the photometric catalog,  astrometry head
        cat = newcat = self.dm.l1_sci(this_ccd_id, suffix="img", ext="acat")
        ref = self.dm.l1_sci(this_ccd_id, suffix="img", ext="rcat")
        header = self.prepare(this_ccd_id, wcsdir, workdir, newcat=False)

        # get obs. Information from header
        # k=list(header.keys())
        w = WCS(header)
        band = header['FILTER'][0]  # ['CCDLABEL']
        # !#airmass=header['AIRMASS']
        # exptime=header['EXPTIME']
        exptime = 150.0
        naxis1 = 9216  # header['NAXIS1']
        naxis2 = 9232  # header['NAXIS2']
        # gain=header['GAIN1']
        pixsize = np.mean(proj_plane_pixel_scales(w)) * 3600.0  # pixsize=0.33
        imr0, imd0 = w.all_pix2world(naxis1 / 2., naxis2 / 2., 0)

        # obs catalog: coord & magnitude
        catdata = fits.open(cat)
        sexcat = catdata[2].data
        # sexcat=table.Table.read(cat,2)
        if upcat:
            x = sexcat['XWIN_IMAGE']
            y = sexcat['YWIN_IMAGE']
            r, d = w.all_pix2world(x, y, 0)
            sexcat['ALPHAWIN_J2000'] = r
            sexcat['DELTAWIN_J2000'] = d
            catdata.writeto(newcat, overwrite=True)
        catdata.close()

        select_stars = (sexcat['MAG_APER'] - sexcat['MAG_AUTO'] < 0.05)
        if select_stars.sum() > 0:
            sexcat = sexcat[select_stars]
        r = sexcat['ALPHAWIN_J2000']  # ALPHA_J2000->ALPHAWIN_J2000
        d = sexcat['DELTAWIN_J2000']
        # back=sexcat['BACKGROUND']
        fwhmsex = sexcat['FWHM_IMAGE']
        obsflags = sexcat['FLAGS']
        cstar = sexcat['CLASS_STAR']
        ellip = sexcat['ELLIPTICITY']
        obsm = sexcat['MAG_APER']  # aperture radii=13pix
        obsme = sexcat['MAGERR_APER']
        ccdnstar = obsm.size  # ccdnstar; total number of stars detected on a CCD
        if len(sexcat) <= 0:
            # ckim=image[image.rfind('/')+1:image.rfind('.')][0:5]+'checkim.ls'
            imcheck = open(ckim, "a")
            # print('line501: image needs update. Pass.')
            # print(image, file=imcheck)
            print("CCD_ID={}: sex catalog length = {}, needs check".format(this_ccd_id, len(sexcat)))
            imcheck.close()
            return

        # obs coordinate
        obsc = SkyCoord(r, d, unit=(u.deg, u.deg), frame='fk5')
        obsrmin = r.min()
        obsrmax = r.max()
        obsdmin = d.min()
        obsdmax = d.max()
        obs_region = [obsrmin, obsrmax, obsdmin, obsdmax]
        # measure image seeing:
        seeing = -1
        fwhm = -1
        fwhm = self.getfwhm(fwhmsex, ellip, obsme, obsflags)
        seeing = fwhm * pixsize

        # get reference catalog
        # ps=self.read_gaiacat(r,d,outcat=ref,silent=True)

        # read input catalog, write to img.rcat
        ps = self.read_inputcat_(this_ccd_id)
        ps.write(ref, format='fits', overwrite=True)
        # remove tmp file: cat & ref

        if np.size(ps) < 1:
            wcscheck = open(ckf, "a")
            # print(image, file=wcscheck)
            print("CCD_ID={}: ref catalog length = 0, needs check".format(this_ccd_id, len(ps)))
            wcscheck.close()
            return
        else:
            # refm,refme,gimedian,refc=self.gaia_mags(ps,band,obs_region)
            refm, refme, gimedian, refc = self.input_mags(ps, usepixcood=False)

        # 4: different doors
        im1 = (obsc.ra.deg >= imr0) & (obsc.dec.deg <= imd0)
        im2 = (obsc.ra.deg < imr0) & (obsc.dec.deg <= imd0)
        im3 = (obsc.ra.deg >= imr0) & (obsc.dec.deg > imd0)
        im4 = (obsc.ra.deg < imr0) & (obsc.dec.deg > imd0)
        ref1 = (refc.ra.deg >= imr0) & (refc.dec.deg <= imd0)
        ref2 = (refc.ra.deg < imr0) & (refc.dec.deg <= imd0)
        ref3 = (refc.ra.deg >= imr0) & (refc.dec.deg > imd0)
        ref4 = (refc.ra.deg < imr0) & (refc.dec.deg > imd0)
        # for all chips
        coeff0 = self.match_calib(obsc, refc, obsm, refm, obsme, refme, obsflags, fwhmsex)
        if coeff0[2] <= 0:
            imcheck = open(ckim, "a")
            # print('line560: image needs update. Pass.')
            print("CCD_ID={}: coeff0[2] = {} <=0: ".format(this_ccd_id, coeff0[2]))
            print("this is coeff0:", coeff0)
            imcheck.close()
            # return
            # coeff0, cstd, csize, cobsm, crefm, ccdraoff, ccddecoff, fwhm

        coeff1 = self.match_calib(obsc[im1], refc[ref1], obsm[im1], refm[ref1], obsme[im1], refme[ref1], obsflags[im1])
        coeff2 = self.match_calib(obsc[im2], refc[ref2], obsm[im2], refm[ref2], obsme[im2], refme[ref2], obsflags[im2])
        coeff3 = self.match_calib(obsc[im3], refc[ref3], obsm[im3], refm[ref3], obsme[im3], refme[ref3], obsflags[im3])
        coeff4 = self.match_calib(obsc[im4], refc[ref4], obsm[im4], refm[ref4], obsme[im4], refme[ref4], obsflags[im4])
        coeff = [coeff0[0], coeff1[0], coeff2[0], coeff3[0], coeff4[0]]
        std = [coeff0[1], coeff1[1], coeff2[1], coeff3[1], coeff4[1]]
        match = [coeff0[2], coeff1[2], coeff2[2], coeff3[2], coeff4[2]]
        ccdraoff = coeff0[5]
        ccddecoff = coeff0[6]

        colortpar = '0.0,0.0,0.0,0.0'
        aperture_radii = 10
        aper_r_com = '(pixels) photo-aperture radius'
        colortpar_com = 'par.in [br^3,br^2,br^1,br^0]'

        ccdzp_com = 'zero point for CCD'
        # ccdzpa_com='zero point for CCD ampA'
        # ccdzpb_com='zero point for CCD ampB'
        # ccdzpc_com='zero point for CCD ampC'
        # ccdzpd_com='zero point for CCD ampD'
        ccdphrms_com = 'zpt rms of the matched objects in CCD'
        # ccdphrms_coma='zpt rms of the matched objects in CCD ampA'
        # ccdphrms_comb='zpt rms of the matched objects in CCD ampB'
        # ccdphrms_comc='zpt rms of the matched objects in CCD ampC'
        # ccdphrms_comd='zpt rms of the matched objects in CCD ampD'

        fwhm_com = 'FWHM in pixel'
        seeing_com = 'seeing in arcsec'
        ccdraoff_com = 'median positional offset from GAIA, in arcsec'
        ccddecoff_com = 'median positional offset from GAIA, in arcsec'
        # transparency_com='median transparency.'#
        ccdnstar_com = 'total number of stars detected on a CCD'
        ccdnmatch_com = 'total number of matched stars in 2 arcsec'
        # ccdnmatcha_com='number of matched stars in CCD ampA'
        # ccdnmatchb_com='number of matched stars in CCD ampB'
        # ccdnmatchc_com='number of matched stars in CCD ampC'
        # ccdnmatchd_com='number of matched stars in CCD ampD'
        ccdmdncol_com = 'median (BP-RP)_GAIA of matched stars in CCD'
        cali_ref_com = 'the reference database for calibration'
        vernum = 'FluxCalib_v1.0'
        vernum_com = 'version of calibration code'

        # keys=['cali_ref','ccdzp','ccdzpa','ccdzpb','ccdzpc','ccdzpd','ccdphoff','ccdphrms', 'phrmsA','phrmsB',
        # 'phrmsC','phrmsD','aper_r','fwhm','seeing','raoff','decoff','trans','ccdnstar','nmatch','nmatcha',
        # 'nmatchb','nmatchc','nmatchd', 'mdncol','colt_par','ebv','EXTNAME','cali_v']
        # tmpa=round(coeff0[0],4)

        # set header keywords: flux calibration information##
        header.set('cali_ref', cali_ref, cali_ref_com)
        header.set('COMMENT', '=' * 66, before='cali_ref')
        header.set('COMMENT', 'Flux calibration information', before='cali_ref')
        header.set('COMMENT', '=' * 66, before='cali_ref')
        header.set('ccdzp', round(coeff0[0], 4), ccdzp_com)
        # header.set('ccdzpa',float('%.4f' % coeff1[0]),ccdzpa_com)
        # header.set('ccdzpb',float('%.4f' % coeff2[0]),ccdzpb_com)
        # header.set('ccdzpc',float('%.4f' % coeff3[0]),ccdzpc_com)
        # header.set('ccdzpd',float('%.4f' % coeff4[0]),ccdzpd_com)
        header.set('ccdphrms', round(coeff0[1], 4), ccdphrms_com)
        # header.set('phrmsA',float('%.4f' % coeff1[1]),ccdphrms_coma)
        # header.set('phrmsB',float('%.4f' % coeff2[1]),ccdphrms_comb)
        # header.set('phrmsC',float('%.4f' % coeff3[1]),ccdphrms_comc)
        # header.set('phrmsD',float('%.4f' % coeff4[1]),ccdphrms_comd)

        header.set('aper_r', round(aperture_radii, 4), aper_r_com)
        header.set('fwhm', float('%.4f' % fwhm), fwhm_com)
        # header.set('seeing',round(seeing,4),seeing_com)
        header.set('raoff', round(ccdraoff, 4), ccdraoff_com)
        header.set('decoff', round(ccddecoff, 4), ccddecoff_com)
        # header.set('trans',round(transparency,4),transparency_com)
        header.set('ccdnstar', round(ccdnstar), ccdnstar_com)
        header.set('nmatch', round(coeff0[2]), ccdnmatch_com)
        # header.set('nmatcha',round(coeff1[2]),ccdnmatcha_com)
        # header.set('nmatchb',round(coeff2[2]),ccdnmatchb_com)
        # header.set('nmatchc',round(coeff3[2]),ccdnmatchc_com)
        # header.set('nmatchd',round(coeff4[2]),ccdnmatchd_com)
        header.set('mdncol', round(gimedian, 4), ccdmdncol_com)
        # header.set('colt_par',colortpar,colortpar_com)
        # header.set('ebv',round(ebv,4),'E(B-V) from SFD1998')

        # Calculate and set SKY & magnitude limiting
        if not ('SKYRMS' in list(header.keys())):
            # imdata=fits.getdata(image, 0)
            imdata = imgdata[1].data
            # TODO: this is not the center of CSST CCDs
            skystat = sigma_clipped_stats(imdata[500:1500, 500:1500], sigma=3.)
            ccdskyrms = skystat[2]  # rms/pixel of the sky in CCD
            ccdsky_com = '(e-/s per pixel)'
            ccdskyrms_com = 'rms/pixel of the sky in unit of e-/s'
            ccdsky = skystat[1]
            header.set('sky', float('%.4f' % ccdsky), ccdsky_com)
            header.set('skyrms', float('%.4f' % ccdskyrms), ccdskyrms_com)

        avsky = header['SKYRMS']
        # ebv=getebv(image)
        ebv = 0.0
        ccdzp = coeff0[0]
        mlim = self.getmlim(fwhm=fwhm, avsky=avsky, rdnoise=5.0, zpt=ccdzp, ebv=ebv, filter=band)
        mlim_com = 'magnitude limiting of 5-sigma galaxy detection'
        header.set('mlim', round(mlim, 2), mlim_com)
        # set signals of calibration progress
        opetime = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
        header.set('FLUX_S', 0, 'flux calibration status')
        header.set('FLUX_V', '1.3', vernum_com)
        header.set('FLUX_TOL', opetime, 'flux calibration operation time')

        # if True, write to _img.head
        if morehead:
            # catroot=cat[:cat.rfind('phot.fits')]
            prihdu = fits.PrimaryHDU(header=header)
            prihdu.writeto(image_head_file, overwrite=True)
            # chead = fits.open(image_head_file,mode='update')
            # imh=chead[0].header
            # imh.extend(header,unique=True,update=True)
            # chead.flush()
            # chead.close()

        # crate and open QC1 list file to save the calibrated image names
        qc1list = os.path.join(workdir, 'file_list.tmp')  # the name of QC1 list
        qc1 = open(qc1list, "a")  # open the file

        # if True, add head to L1 image
        if addhead:
            # print(('add head for ',image))
            # imwht = image.replace('_img', '_wht')
            # imflg = image.replace('_img', '_flg')
            for suffix, im in zip(["img", "wht", "flg"], [imgdata, whtdata, flgdata]):
                # newimage = os.path.join(L1dir, os.path.split(img)[-1])
                # newimage = newimage.replace('.fits', '_L1.fits')
                newimage = self.dm.l1_sci(this_ccd_id, suffix=suffix+"_L1")
                # im = fits.open(img,mode='readonly')
                h0 = fits.PrimaryHDU(data=im[0].data, header=im[0].header)
                h1 = fits.ImageHDU(data=im[1].data, header=header) # updates are here!
                hdulist = fits.HDUList([h0, h1])
                hdulist.writeto(newimage, overwrite=True)
                # im.writeto(newimage, overwrite=True)
                # print the calibrated image name to QC1 list file
                print(newimage, file=qc1)
        # close the QC1 list file
        qc1.close()

        # plot
        if plot:
            print('plot calibration chart ...')
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            fig = plt.figure(figsize=(6, 5), dpi=120)
            plt.subplots_adjust(left=0.2, bottom=0.15)
            plt.ylim(coeff0[0] - 5 * coeff0[1] - 0.1, coeff0[0] + 5 * coeff0[1] + 0.1)
            plt.xlim(np.min(coeff0[4]) - 0.2, np.max(coeff0[4]) + 0.2)
            plt.xlabel(r'$\mathrm{MAG_{Gaia}}$')
            plt.ylabel(r'$\mathrm{MAG_{Gaia}-MAG_{inst}}$')
            plt.plot(coeff1[4], -coeff1[3] + coeff1[4], 'r.')
            plt.plot(coeff2[4], -coeff2[3] + coeff2[4], 'y.')
            plt.plot(coeff3[4], -coeff3[3] + coeff3[4], 'b.')
            plt.plot(coeff4[4], -coeff4[3] + coeff4[4], 'g.')
            label0 = 'zpt=%6.3f, rms:%6.3f, matched:%i' % (coeff0[0:3])
            label1 = 'zpta=%6.3f, rms:%6.3f, matched:%i' % (coeff1[0:3])
            label2 = 'zptb=%6.3f, rms:%6.3f, matched:%i' % (coeff2[0:3])
            label3 = 'zptc=%6.3f, rms:%6.3f, matched:%i' % (coeff3[0:3])
            label4 = 'zptd=%6.3f, rms:%6.3f, matched:%i' % (coeff4[0:3])
            xtmin = np.min(coeff0[4]) - 0.2
            xtmax = np.max(coeff0[4]) + 0.2
            l1, = plt.plot([xtmin, xtmax], np.array([0, 0]) + coeff0[0], 'k-', lw=1, label=label0)
            l2, = plt.plot([xtmin, xtmax], np.array([0, 0]) + coeff1[0], 'r-', lw=1, label=label1)
            l3, = plt.plot([xtmin, xtmax], np.array([0, 0]) + coeff2[0], 'y-', lw=1, label=label2)
            l4, = plt.plot([xtmin, xtmax], np.array([0, 0]) + coeff3[0], 'b-', lw=1, label=label3)
            l5, = plt.plot([xtmin, xtmax], np.array([0, 0]) + coeff4[0], 'g-', lw=1, label=label4)
            plt.legend(handles=[l1, l2, l3, l4, l5], frameon=False, loc=2, prop={'size': 10})
            plt.savefig(psname)

        return coeff, std, match

    def run(self, addhead=True, morehead=False, plot=False, nodel=True, update=False, upcat=True):

        # if len(fn_list) == 0:
        #     print('Flux calibration: No input images in img_list!')
        #     return
        # time1=time.time()
        print('############### run flux calibration ###############')
        # split combined_acat.head.fits
        wcshead0 = self.dm.l1_hardcode(hdcd="combined_acat.head.fits", comment="run")
        self.split_wcs_head_(wcshead0)

        # for i in range(len(fn_list)):
        for i, this_ccd_id in enumerate(self.dm.target_ccd_ids):
            # image = workdir + fn_list[i]
            # if len(img_list) == len(fn_list):
            #     imgdata = img_list[i]
            #     whtdata = wht_list[i]
            #     flgdata = flg_list[i]
            # else:
            # print(('flux calibration: '+image))
            # if not os.path.isfile(image):
            #     print(('cannot find the file:' + image))
            # else:
            self.calib(this_ccd_id, addhead=addhead, morehead=morehead, plot=plot, nodel=nodel, update=update)
        # time2=time.time()
        print('\n############### flux calibration done #############\n')

    def cleanup(self, fn_list, workdir, nodel=False):
        # clean up environment
        for image in fn_list:
            cat = os.path.join(workdir, image[:-5] + '.acat')
            ref = os.path.join(workdir, image[:-5] + '.rcat')
            whead = os.path.join(workdir, image[:-5] + '.whead.fits')
            # psname =os.path.join(workdir,image+'_calib.png')
            if not nodel:
                try:
                    os.remove(cat)
                except FileNotFoundError:
                    print()
                try:
                    os.remove(ref)
                except FileNotFoundError:
                    print()
                try:
                    os.remove(whead)
                except FileNotFoundError:
                    print()
