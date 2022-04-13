from collections import OrderedDict

import astropy.io.fits as fits
from astropy.io.fits import HDUList, PrimaryHDU
import numpy as np

from csst.core.exception import CsstException
from astropy.io import fits
from copy import deepcopy

__all__ = ["CsstData", "INSTRUMENT_LIST"]

INSTRUMENT_LIST = ["MSC", ]


class CsstData(fits.HDUList):
    """ General CSST raw class """
    # hdu_pri = fits.PrimaryHDU()
    # hdu_l0 = fits.ImageHDU(raw=None, name="raw")
    # hdu_l1 = fits.ImageHDU(raw=None, name="l1")

    """
    header methods:
    'add_blank', 'add_comment', 'add_history', 'append', 'cards', 'clear', 'comments', 'copy', 'count', 'extend', 
    'fromfile', 'fromkeys', 'fromstring', 'fromtextfile', 'get', 'index', 'insert', 'items', 'keys', 'pop', 'popitem', 
    'remove', 'rename_keyword', 'set', 'setdefault', 'strip', 'tofile', 'tostring', 'totextfile', 'update', 'values'
    """

    def __init__(self, hdus=None, file=None):
        """

        Parameters
        ----------
        hdus:
            a alist of HDUs
        file:
            open file object
        """
        if hdus is None:
            hdus = []
        super(CsstData, self).__init__(hdus=hdus, file=file)

    def get_data(self, copy=True, hdu=1):
        """ get level 0 raw from CsstData class

        Parameters
        ----------
        copy : bool
            if True, return a copy.
        """
        if copy:
            return self[hdu].data.copy()
        else:
            return self[hdu].data

    def get_keyword(self, key="INSTRUME", hdu=0):
        """ get keyword from fits header

        Parameters
        ----------
        key:
            the key
        """
        return self[hdu].header.get(key)

    def set_keyword(self, key, value, hdu=1):
        """ set keyword

        Parameters
        ----------
        key:
            key
        value:
            value
        hdu:
            0 for primary hdu, 1+ for raw hdu

        """
        self[hdu].header[key] = value
        return

    def set_data(self, data, hdu=1):
        """ set image raw """
        self[hdu].data = data
        return

    # def writeto(self, fp, overwrite=False):
    #     """ save L1 image and aux raw to file
    #
    #     Parameters
    #     ----------
    #     fp: str
    #         image type
    #     overwrite : bool
    #         if True, overwrite file
    #     """
    #     self.writeto(fp, overwrite=overwrite)

    def get_auxdata(self):
        """ get aux raw
        In future, this is to automatically get aux raw from database.
        """
        raise NotImplementedError

    @classmethod
    def read(cls, name, ignore_missing_simple=True):
        """ read raw from fits file, should be implemented in child classes """
        return cls.fromfile(name, ignore_missing_simple=ignore_missing_simple)

    def deepcopy(self, name=None, data=None):
        """ generate a deep copy of self """
        cp = self.__class__(deepcopy(self))
        if name is not None:
            cp[1].name = name
        if data is not None:
            cp[1].data = data
        return cp

    @property
    def data(self):
        if len(self) == 1:
            return self[1].data
        return self[1].data

    @property
    def exptime(self):
        return self[0].header["EXPTIME"]
