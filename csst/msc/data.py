from collections import OrderedDict
import astropy.io.fits as fits
from astropy.io.fits import HDUList, PrimaryHDU, ImageHDU
from astropy.io.fits.header import Header
from ..core.data import CsstData, INSTRUMENT_LIST
import numpy as np


__all__ = ["CsstMscData", "CsstMscImgData"]


class CsstMscData(CsstData):
    _l1img_types = {'sci': True, 'weight': True, 'flag': True}

    def __init__(self, priHDU, imgHDU, **kwargs):
        super(CsstMscData, self).__init__(priHDU, imgHDU, **kwargs)
        self._l1hdr_global = priHDU.header.copy()
        self._l1data['sci'] = ImageHDU()
        self._l1data['weight'] = ImageHDU()
        self._l1data['flag'] = ImageHDU()

    def set_flat(self, flat):
        """ set flat

        Parameters
        ----------
        flat:
            flat image

        Returns
        -------

        """
        self._auxdata['flat'] = flat

    def set_bias(self, biasimg):
        """ set bias """
        self._auxdata['bias'] = biasimg

    def set_dark(self, darkimg):
        """ set dark """
        self._auxdata['dark'] = darkimg

    def set_badpixel(self, badpixelimg):
        """ set badpixel """
        self._auxdata['badpixel'] = badpixelimg

    def get_flat(self):
        """ get flat """
        return self._auxdata['flat']

    def get_bias(self):
        """ get bias """
        return self._auxdata['bias']

    def get_dark(self):
        """ get dark """
        return self._auxdata['dark']

    def get_badpixel(self):
        """ get badpixel """
        return self._auxdata['badpixel']

    def init_l0data(self):
        """ initialize L0 data """
        pass

    def set_l1keyword(self, key, value, comment=''):
        """ set L1 keyword """
        print('check out whether ' + key + " is a valid key and " + value + " is valid value")
        self._l1hdr_global.set(key, value, comment)

    def set_l1data(self, imgtype, img):
        """ set L1 data """
        try:
            if imgtype == 'sci':
                self._l1data[imgtype].header['EXTNAME'] = 'img'
                self._l1data[imgtype].header['BUNIT'] = 'e/s'
                self._l1data[imgtype].data = img.astype(np.float32) / self._l1hdr_global['exptime']
            elif imgtype == 'weight':
                self._l1data[imgtype].header['EXTNAME'] = 'wht'
                self._l1data[imgtype].data = img.astype(np.float32)
            elif imgtype == 'flag':
                self._l1data[imgtype].header['EXTNAME'] = 'flg'
                self._l1data[imgtype].data = img.astype(np.uint16)
            else:
                raise TypeError('unknow type image')
        except Exception as e:
            print(e)
        print('save image data to l1data')

    def save_l1data(self, imgtype, filename):
        """ save L1 data """
        print('check ' + imgtype + ' is validate')
        try:
            if self._l1img_types[imgtype]:
                super().save_l1data(imgtype, filename)
        except Exception as e:
            print(e)


class CsstMscImgData(CsstMscData):
    def __init__(self, priHDU, imgHDU, **kwargs):
        # print('create CsstMscImgData')
        super(CsstMscImgData, self).__init__(priHDU, imgHDU, **kwargs)

    def __repr__(self):
        return "<CsstMscImgData: {} {}>".format(self.instrument, self.detector)

    @staticmethod
    def read(fp):
        """ read data from fits file

        Parameters
        ----------
        fp:
            the file path of fits file

        Returns
        -------
        CsstMscImgData

        Example
        -------

        >>> fp = "MSC_MS_210527171000_100000279_16_raw.fits"
        >>> from csst.msc import CsstMscImgData
        >>> data = CsstMscImgData.read(fp)
        >>> # print some info
        >>> print("data: ", data)
        >>> print("instrument: ", data.get_l0keyword("pri", "INSTRUME"))
        >>> print("object: ", data.get_l0keyword("pri", "OBJECT"))
        """

        try:
            with fits.open(fp) as hdulist:
                instrument = hdulist[0].header.get('INSTRUME')  # strip or not?
                detector = hdulist[0].header.get('DETECTOR')  # strip or not?
                print("@CsstMscImgData: reading data {} ...".format(fp))
                assert instrument in INSTRUMENT_LIST
                if instrument == 'MSC' and 6 <= int(detector[3:5]) <= 25:
                    # multi-band imaging
                    hdu0 = hdulist[0].copy()
                    hdu1 = hdulist[1].copy()
                    data = CsstMscImgData(hdu0, hdu1, instrument=instrument, detector=detector)
                    return data
        except Exception as e:
            print(e)
