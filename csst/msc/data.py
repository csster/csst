# from abc import ABC
from collections import OrderedDict
import astropy.io.fits as fits
from astropy.io.fits import HDUList, PrimaryHDU, ImageHDU
from astropy.io.fits.header import Header
from ..core.data import CsstData, INSTRUMENT_LIST
import numpy as np


__all__ = ["CsstMscImgData", ]


class CsstMscImgData(CsstData):
    _l1img_types = {'sci': True, 'weight': True, 'flag': True}

    def __init__(self, hdus=None, file=None):
        """

        Parameters
        ----------
        hdus:
            a list of HDUs
        file:
            open file object
        """
        if hdus is None:
            hdus = []
        super(CsstMscImgData, self).__init__(hdus=hdus, file=file)

        # self._l1hdr_global = self[0].header.copy()
        # self._l1data = dict()
        # self._l1data['sci'] = ImageHDU()
        # self._l1data['weight'] = ImageHDU()
        # self._l1data['flag'] = ImageHDU()

    @property
    def instrument(self):
        return self[0].header["INSTRUME"]

    @property
    def detector(self):
        return self[0].header["DETECTOR"]

    def get_flat(self, fp):
        """ get flat """
        return fits.getdata(fp)

    def get_bias(self, fp):
        """ get bias """
        return fits.getdata(fp)

    def get_dark(self, fp):
        """ get dark """
        return fits.getdata(fp)

    def get_l1data(self):
        """ get L1 raw """
        imgdata = self.get_data(hdu=1)
        exptime = self.get_keyword("EXPTIME", hdu=0)
        # image
        img = self.deepcopy(name="img", data=imgdata.astype(np.float32) / exptime)
        img[1].header['BUNIT'] = 'e/s'
        # weight
        wht = self.deepcopy(name="wht", data=imgdata.astype(np.float32))
        wht[1].header.remove('BUNIT')
        # flag
        flg = self.deepcopy(name="flg", data=imgdata.astype(np.uint16))
        flg[1].header.remove('BUNIT')
        return img, wht, flg

    def __repr__(self):
        return "<CsstMscImgData: {} {}>".format(self.instrument, self.detector)

    # @staticmethod
    # def read(fp):
    #     """ read raw from fits file
    #
    #     Parameters
    #     ----------
    #     fp:
    #         the file path of fits file
    #
    #     Returns
    #     -------
    #     CsstMscImgData
    #
    #     Example
    #     -------
    #
    #     >>> fp = "MSC_MS_210527171000_100000279_16_raw.fits"
    #     >>> from csst.msc import CsstMscImgData
    #     >>> raw = CsstMscImgData.read(fp)
    #     >>> # print some info
    #     >>> print("raw: ", raw)
    #     >>> print("instrument: ", raw.get_l0keyword("pri", "INSTRUME"))
    #     >>> print("object: ", raw.get_l0keyword("pri", "OBJECT"))
    #     """
    #     return CsstMscImgData.fromfile(fp)
