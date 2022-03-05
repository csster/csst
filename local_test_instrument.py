from csst.msc import CsstMscImgData
from csst.msc.instrument import CsstMscInstrumentProc
from astropy.io.fits import getdata


fp = "/data/cali_20211012/L0/150s/MSC_MS_210525120000_100000000_13_raw.fits"
bs = "/data/ref/MSC_CLB_210525190000_100000014_13_combine.fits"
dk = "/data/ref/MSC_CLD_210525192000_100000014_13_combine.fits"
ft = "/data/ref/MSC_CLF_210525191000_100000014_13_combine.fits"
op = "/data/test/"
data = CsstMscImgData.read(fp)
# print(repr(data._l1data.header))
data.set_bias(getdata(bs))
data.set_dark(getdata(dk))
data.set_flat(getdata(ft))

proc = CsstMscInstrumentProc()
proc.run(data)
data.save_l1data('sci', op+'sci.fits')
data.save_l1data('flag', op+'flag.fits')
data.save_l1data('weight', op+'weight.fits')
