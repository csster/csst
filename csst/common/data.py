from collections import OrderedDict

import astropy.io.fits as fits
from astropy.io.fits import HDUList, PrimaryHDU

from .CsstException import CsstException


INSTRUMENT_LIST = ["MSC", ]


class CsstData:
    """
    general CSST data class
    """
    _primary_hdu = []
    _l0data = []  # HDUList
    _l1hdr_global = []
    _l1data = OrderedDict()  # dict object
    _l2data = OrderedDict()  #
    _auxdata = OrderedDict()

    def __init__(self, primaryHDU, imgHDU, instrument=None, detector=None):
        print('create CsstData')
        self._primary_hdu = primaryHDU
        self._l0data = imgHDU
        self.instrument = instrument
        self.detector = detector

    def get_l0data(self, copy=True):
        """
        obtain level 0 data from CsstData class
        copy: True: if the user want to copy the memory of the data to the new class;
              False: only reference of the data memory is written to the new class
        """
        if copy:
            return self._l0data.data.copy()
        else:
            return self._l0data.data

    def get_l0keyword(self, ext="pri", key="INSTRUME"):
        """
        obtain keywords of the fits header of level 0 image data from the CsstData class
        ext: the index of extension. if it equals to 'pri', looking up keywords from primary session, otherwise from extension sessions
        key: the name of the key
        """
        if ext == 'pri':
            try:
                value = self._primary_hdu.header.get(key)
            except Exception as e:
                print(e)
        elif ext == 'img':
            try:
                value = self._l0data.header.get(key)
            except Exception as e:
                print(e)
        else:
            raise CsstException

    def set_l1keyword(self, key, value):
        print('check out whether ' + key + " is a valid key and " + value + " is valid value")

    def set_l1data(self, img):
        print('save image data to l2data')

    def get_auxdata(self, name):
        print('Parent class returns zero image.')
        # return np.zeros_like(self.get_l0data())
        return

    def save_l1data(self, imgtype, filename):
        """
        asve level 1 image and auxilary data to data file
        imgtype
        """
        print("save L1 image to a fits file with name " + filename)
        try:
            self._l1hdr_global.set('TYPE', imgtype, 'Type of Level 1 data')
            pri_hdu = PrimaryHDU(header=self._l1hdr_global)
            hdulist = HDUList([pri_hdu, self._l1data[imgtype]])
            hdulist.writeto(filename)
        except Exception as e:
            print(e)
