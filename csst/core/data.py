from collections import OrderedDict

import astropy.io.fits as fits
from astropy.io.fits import HDUList, PrimaryHDU

from csst.core.exception import CsstException


__all__ = ["CsstData", "INSTRUMENT_LIST"]

INSTRUMENT_LIST = ["MSC", ]


class CsstData:
    """ General CSST data class """
    _primary_hdu = []
    _l0data = []  # HDUList
    _l1hdr_global = []
    _l1data = OrderedDict()  # dict object
    _l2data = OrderedDict()  #
    _auxdata = OrderedDict()

    def __init__(self, primaryHDU, imgHDU, instrument=None, detector=None):
        # print('create CsstData')
        self._primary_hdu = primaryHDU
        self._l0data = imgHDU
        self.instrument = instrument
        self.detector = detector

    def get_l0data(self, copy=True):
        """ get level 0 data from CsstData class

        Parameters
        ----------
        copy : bool
            if True, return a copy.
        """
        if copy:
            return self._l0data.data.copy()
        else:
            return self._l0data.data

    def get_l0keyword(self, ext="pri", key="INSTRUME"):
        """ get a specific keyword from fits header of level 0 image data

        Parameters
        ----------
        ext: {"pri"| "img"}
            the HDU extension
        key:
            the key
        """
        if ext == 'pri':
            try:
                return self._primary_hdu.header.get(key)
            except Exception as e:
                print(e)
        elif ext == 'img':
            try:
                return self._l0data.header.get(key)
            except Exception as e:
                print(e)
        else:
            raise CsstException

    def set_l1keyword(self, key, value):
        """ set  L1 keyword """
        raise NotImplementedError("Well, not implemented...")

    def set_l1data(self, *args, **kwargs):
        print('save image data to l2data')
        raise NotImplementedError

    def get_auxdata(self, name):
        """ get aux data

        Parameters
        ----------
        """
        print('Parent class returns zero image.')
        # return np.zeros_like(self.get_l0data())
        raise NotImplementedError

    def save_l1data(self, imgtype, filename):
        """ save L1 image and auxilary data to file

        Parameters
        ----------
        imgtype: {}
            image type
        """
        print("save L1 image to a fits file with name " + filename)
        try:
            self._l1hdr_global.set('TYPE', imgtype, 'Type of Level 1 data')
            pri_hdu = PrimaryHDU(header=self._l1hdr_global)
            hdulist = HDUList([pri_hdu, self._l1data[imgtype]])
            hdulist.writeto(filename)
        except Exception as e:
            print(e)

    def read(self, **kwargs):
        """ read data from fits file """
        raise NotImplementedError
