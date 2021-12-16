from collections import OrderedDict

import astropy.io.fits as fits
from astropy.io.fits import HDUList, PrimaryHDU

from csst.msc.mscdata import CsstMscImgData
from .CsstException import CsstException


class CsstDataFactory:
    """
    This class is designed to create CsstData and its inherited classes according to different kinds of input data format.
    """

    def __init__(self):
        pass

    def createData(self, fitsfilename):
        """

        Parameters
        ----------
        fitsfilename:
            the file name of fits files

        Returns
        -------

        """

        try:
            hl = fits.open(fitsfilename)
            instrument = hl[0].header.get('INSTRUME')  # strip or not?
            detector = hl[0].header.get('DETECTOR')    # strip or not?
            print(instrument, detector)
            if instrument == 'MSC' and int(detector[3:5]) >= 6 and int(detector[3:5]) <= 25:
                # multi-band imaging
                data = CsstMscImgData(hl[0], hl[1])

        except Exception as e:
            print(e)
        return data


class CsstData:
    _primary_hdu = []
    _l0data = []  # HDUList
    _l1hdr_global = []
    _l1data = OrderedDict()  # dict object
    _l2data = OrderedDict()  #
    _auxdata = OrderedDict()

    def __init__(self, primaryHDU, imgHDU):
        print('create CsstData')
        self._primary_hdu = primaryHDU
        self._l0data = imgHDU

    def get_l0data(self, *, copy):
        '''
        obtain level 0 data from CsstData class
        copy: True: if the user want to copy the memory of the data to the new class;
              False: only reference of the data memory is written to the new class
        '''
        if copy:
            return self._l0data.data.copy()
        else:
            return self._l0data.data

    def get_l0keyword(self, ext, key):
        '''
        obtain keywords of the fits header of level 0 image data from the CsstData class
        ext: the index of extension. if it equals to 'pri', looking up keywords from primary session, otherwise from extension sessions
        key: the name of the key
        '''
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
        return value

    def set_l1keyword(self, key, value):
        print('check out whether ' + key + " is a valid key and " + value + " is valid value")

    def set_l1data(self, img):
        print('save image data to l2data')

    def get_auxdata(self, name):
        print('Parent class returns zero image.')
        return np.zeros_like(get_l0data())

    def save_l1data(self, imgtype, filename):
        '''
        asve level 1 image and auxilary data to data file 
        imgtype
        '''
        print("save L1 image to a fits file with name " + filename)
        try:
            self._l1hdr_global.set('TYPE', imgtype, 'Type of Level 1 data')
            pri_hdu = PrimaryHDU(header=self._l1hdr_global)
            hdulist = HDUList([pri_hdu, self._l1data[imgtype]])
            hdulist.writeto(filename)
        except Exception as e:
            print(e)
