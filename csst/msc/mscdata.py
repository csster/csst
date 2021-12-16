from collections import OrderedDict
import astropy.io.fits as fits
from astropy.io.fits import HDUList, PrimaryHDU, ImageHDU
from astropy.io.fits.header import Header
from csst.common.data import CsstData

class CsstMscData(CsstData):
    _l1img_types = {'sci':True,'weight':True,'flag':True}
    def __init__(self, primaryHDU, imgHDU):
        print('create CsstMscData')
        super().__init__(primaryHDU, imgHDU)
        self._l1hdr_global = primaryHDU.header.copy()
#         self._l1hdr_global['SIMPLE']  =  'T' #/ conforms to FITS standard
#         self._l1hdr_global['NAXIS']  =  0
        self._l1data['sci'] =ImageHDU()
        self._l1data['weight'] = ImageHDU()
        self._l1data['flag'] = ImageHDU()
        
    def set_flat(self, flatimg):
        self._auxdata['flat'] = flatimg
    
    def set_bias(self, biasimg):
        self._auxdata['bias'] = biasimg
        
    def set_dark(self, darkimg):
        self._auxdata['dark'] = darkimg
    
    def set_badpixel(self, badpixelimg):
        self._auxdata['badpixel'] = badpixelimg
        
    def get_flat(self):
        return self._auxdata['flat']
    
    def get_bias(self):
        return self._auxdata['bias']
    
    def get_dark(self):
        return self._auxdata['dark']
    
    def get_badpixel(self):
        return self._auxdata['badpixel']
    
    def init_l0data(self):
        pass
#         hdr_global = Header(self._primary_hdu.header, copy = True)
#         hdr_l1 = Header(self._l0data.header, copy = True)
#         hdu_pri = PrimaryHDU(hdr = hdr_global)
#         hdu_img = self._l0data.copy()
        
    def set_l1keyword(self, key, value, comment=''):
        print('check out whether '+key+" is a valid key and "+value+" is valid value")
        self._l1hdr_global.set(key,value,comment)
        
    def set_l1data(self, imgtype, img):
        try:
            if self._l1img_types[imgtype]:
                self._l1data[imgtype].data = img.copy()
        except Exception as e:
            print(e)
        print('save image data to l1data')
        
    def save_l1data(self, imgtype, filename):
        print('check '+imgtype+' is validate')
        try:
            if self._l1img_types[imgtype]:
                super().save_l1data(imgtype, filename)
        except Exception as e:
            print(e)

class CsstMscImgData(CsstMscData):
    def __init__(self, primaryHDU, imgHDU):
        print('create CsstMscImgData')
        super().__init__(primaryHDU, imgHDU)
        
