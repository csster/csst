from astropy.io import fits


class Header(fits.Header):
    """ a custom Header class

    Conf
    ----
    https://docs.astropy.org/en/stable/io/fits/api/headers.html#astropy.io.fits.Header

    """

    def __init__(self, cards=[], copy=False):
        super(Header, self).__init__(cards=cards, copy=copy)


test_hdr_str = """
SIMPLE  =                    T / conforms to FITS standard
BITPIX  =                    8 / array data type
NAXIS   =                    0 / number of array dimensions
EXTEND  =                    T
"""

if __name__ == "__main__":

    hdr = Header.fromstring(test_hdr_str, sep='\n')
    print(hdr['SIMPLE'], hdr['BITPIX'], len(hdr))
