
# def do_one_exposure():

# node name
HOST = "tulip"
# working directory
DIR_WORK = "/share/HDD7/csstpipeline/msc"
# gaia catalog directory (for position calibration)
DIR_GAIA_CATALOG = ""

DIR_TEST = "/share/Cycle-3-SimuData/multipleBandsImaging/CSST_shearOFF/MSC_0000020/"  #  MSC_MS_210525220000_100000020_06_raw.fits
PATH_BIAS = "/share/HDD7/csstpipeline/ref/MSC_CLB_210525200000_100000016_{:02d}_combine.fits"
PATH_DARK = "/share/HDD7/csstpipeline/ref/MSC_CLD_210525202000_100000016_{:02d}_combine.fits"
PATH_FLAT = "/share/HDD7/csstpipeline/ref/MSC_CLF_210525201000_100000016_{:02d}_combine.fits"

import os
os.chdir(DIR_WORK)

import glob
from csst.msc.data import CsstMscImgData
from csst.msc.instrument import CsstMscInstrumentProc
# from astropy.io import fits

for i_ccd in range(6, 26):
    # i_ccd = 6
    print("processing CCD {}".format(i_ccd))
    fp_raw = glob.glob("{}/MSC_MS_*{:02}_raw.fits".format(DIR_TEST, i_ccd))
    assert len(fp_raw) == 1
    fp_raw = fp_raw[0]

    raw = CsstMscImgData.read(fp_raw)
    bias = raw.get_bias(PATH_BIAS.format(i_ccd))
    dark = raw.get_bias(PATH_DARK.format(i_ccd))
    flat = raw.get_bias(PATH_FLAT.format(i_ccd))

    # initialize Instrument Processor
    instProc = CsstMscInstrumentProc()
    instProc.prepare(n_jobs=2)
    img, wht, flg = instProc.run(raw, bias, dark, flat)
    instProc.cleanup()

    # save img, wht, flg to somewhere
    img.writeto("{}/{}.fits".format(DIR_WORK, img.get_keyword("FILENAME")))
    wht.writeto("{}/{}.fits".format(DIR_WORK, wht.get_keyword("FILENAME")))
    flg.writeto("{}/{}.fits".format(DIR_WORK, flg.get_keyword("FILENAME")))


# TODO:  position calibration
from csst.msc.astrometry import CsstProcMscPositionCalibration
pcProc = CsstProcMscPositionCalibration()
pcProc.prepare(search_radius=2.,)
pcProc.run(data_list)
pcProc.cleanup()



