from astropy.io import fits

from .data.data import INSTRUMENT_LIST
from csst.common.data.mscdata import CsstMscImgData


class CsstDataFactory:
    """
    This class is designed to create CsstData and its inherited classes according to different kinds of input data format.
    """

    def __init__(self):
        pass

    @staticmethod
    def createData(fitsfilename):
        """ create CSST Data instances

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
            detector = hl[0].header.get('DETECTOR')  # strip or not?
            print(instrument, detector)
            assert instrument in INSTRUMENT_LIST
            if instrument == 'MSC' and 6 <= int(detector[3:5]) <= 25:
                # multi-band imaging
                data = CsstMscImgData(hl[0], hl[1], instrument=instrument, detector=detector)
            return data
        except Exception as e:
            print(e)
