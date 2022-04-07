
# def do_one_exposure():

# node name
HOST = "tulip"
# working directory
DIR_WORKING = ""
# gaia catalog directory (for position calibration)
DIR_GAIA_CATALOG = ""


# raw images
fp_list = ["MSC_MS_210527171000_100000279_6_raw.fits",
           "MSC_MS_210527171000_100000279_7_raw.fits",
           "MSC_MS_210527171000_100000279_8_raw.fits",
           "MSC_MS_210527171000_100000279_9_raw.fits",
           "MSC_MS_210527171000_100000279_10_raw.fits",
           "MSC_MS_210527171000_100000279_11_raw.fits",
           "MSC_MS_210527171000_100000279_12_raw.fits",
           "MSC_MS_210527171000_100000279_13_raw.fits",
           "MSC_MS_210527171000_100000279_14_raw.fits",
           "MSC_MS_210527171000_100000279_15_raw.fits",
           "MSC_MS_210527171000_100000279_16_raw.fits",
           "MSC_MS_210527171000_100000279_17_raw.fits",
           "MSC_MS_210527171000_100000279_18_raw.fits",
           "MSC_MS_210527171000_100000279_19_raw.fits",
           "MSC_MS_210527171000_100000279_20_raw.fits",
           "MSC_MS_210527171000_100000279_21_raw.fits",
           "MSC_MS_210527171000_100000279_22_raw.fits",
           "MSC_MS_210527171000_100000279_23_raw.fits",
           "MSC_MS_210527171000_100000279_24_raw.fits",
           "MSC_MS_210527171000_100000279_25_raw.fits",
           ]

from csst.msc.data import CsstMscImgData
from csst.msc.instrument import CsstMscInstrumentProc
from astropy.io import fits

# get aux data
bs = fits.getdata("/data/ref/MSC_CLB_210525190000_100000014_13_combine.fits")
dk = fits.getdata("/data/ref/MSC_CLD_210525192000_100000014_13_combine.fits")
ft = fits.getdata("/data/ref/MSC_CLF_210525191000_100000014_13_combine.fits")

fp_img_list = []
fp_flg_list = []
fp_wht_list = []
data_list = []

for fp in fp_list:
    # read image data
    data = CsstMscImgData.read(fp)

    # set aux data
    data.set_bias(bs)
    data.set_dark(dk)
    data.set_flat(ft)

    # initialize Instrument Processor
    instProc = CsstMscInstrumentProc()
    instProc.prepare()
    instProc.run(data)
    instProc.cleanup()

    # output filepath
    fp_img = fp.replace("raw.fits", "img.fits")
    fp_flg = fp.replace("raw.fits", "flg.fits")
    fp_wht = fp.replace("raw.fits", "wht.fits")

    # save l1 data
    data.save_l1data('sci', fp_img)
    data.save_l1data('flag', fp_flg)
    data.save_l1data('weight', fp_wht)

    # store in a list for position calibration
    fp_img_list.append(fp_img)
    fp_flg_list.append(fp_flg)
    fp_wht_list.append(fp_wht)

    data_list.append(data)


# position calibration
from csst.msc.astrometry import CsstProcMscPositionCalibration
pcProc = CsstProcMscPositionCalibration()
pcProc.prepare(search_radius=2.,)
pcProc.run(data_list)
pcProc.cleanup()



