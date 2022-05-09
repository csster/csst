# def do_one_exposure():

import glob
import os

from csst.msc.data import CsstMscImgData
from csst.msc.calib_pos import CsstProcMscPositionCalibration
from csst.msc.inst_corr import CsstMscInstrumentProc

HOSTNAME = os.uname()[1]
if HOSTNAME == "tulip":
    # on Tulip
    DIR_TEST = "/share/Cycle-3-SimuData/multipleBandsImaging/CSST_shearOFF/MSC_0000020/"  # MSC_MS_210525220000_100000020_06_raw.fits
    PATH_BIAS = "/share/HDD7/csstpipeline/ref/MSC_CLB_210525200000_100000016_{:02d}_combine.fits"
    PATH_DARK = "/share/HDD7/csstpipeline/ref/MSC_CLD_210525202000_100000016_{:02d}_combine.fits"
    PATH_FLAT = "/share/HDD7/csstpipeline/ref/MSC_CLF_210525201000_100000016_{:02d}_combine.fits"
    DIR_WORK = "/share/HDD7/csstpipeline/msc"
    # gaia catalog directory (for position calibration)
    DIR_GAIA_CATALOG = ""

elif HOSTNAME == "Dandelion":
    # on Dandelion
    DIR_TEST = "/home/csstpipeline/L1Pipeline/msc/MSC_0000020"
    PATH_BIAS = "/home/csstpipeline/L1Pipeline/msc/ref/MSC_CLB_210525200000_100000016_{:02d}_combine.fits"
    PATH_DARK = "/home/csstpipeline/L1Pipeline/msc/ref/MSC_CLD_210525202000_100000016_{:02d}_combine.fits"
    PATH_FLAT = "/home/csstpipeline/L1Pipeline/msc/ref/MSC_CLF_210525201000_100000016_{:02d}_combine.fits"
    # working directory
    DIR_WORK = "/home/csstpipeline/L1Pipeline/msc/work/"
    # gaia catalog directory (for position calibration)
    DIR_GAIA_CATALOG = "/home/csstpipeline/L1Pipeline/msc/gaia_dr3/"

else:
    raise ValueError("Invalid HOSTNAME {}!".format(HOSTNAME))

# define CCD ID list
CCD_ID_LIST = [6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 24, 25]

# prohibit multi-threading in backend
os.environ["MKL_NUM_THREADS"] = '1'
os.environ["NUMEXPR_NUM_THREADS"] = '1'
os.environ["OMP_NUM_THREADS"] = '1'

# Step 1. Correct instrumental effect
os.chdir(DIR_WORK)

img_list = []
wht_list = []
flg_list = []
fn_list = []
for i_ccd in CCD_ID_LIST:
    print("processing CCD {}".format(i_ccd))
    fp_raw = glob.glob("{}/MSC_MS_*_{:02}_raw.fits".format(DIR_TEST, i_ccd))
    assert len(fp_raw) == 1
    fp_raw = fp_raw[0]

    # read data with CsstMscImgData.read
    raw = CsstMscImgData.read(fp_raw)
    # in future, get_* functions grab
    bias = raw.get_bias(PATH_BIAS.format(i_ccd))
    dark = raw.get_dark(PATH_DARK.format(i_ccd))
    flat = raw.get_flat(PATH_FLAT.format(i_ccd))

    # initialize Instrument Processor
    instProc = CsstMscInstrumentProc()
    instProc.prepare(n_jobs=2)
    img, wht, flg = instProc.run(raw, bias, dark, flat)
    instProc.cleanup()
    fp_img = img[0].header["FILENAME"] + '.fits'

    # save img, wht, flg to somewhere
    img.writeto("{}/{}.fits".format(DIR_WORK, img.get_keyword("FILENAME")), overwrite=True)
    wht.writeto("{}/{}.fits".format(DIR_WORK, wht.get_keyword("FILENAME")), overwrite=True)
    flg.writeto("{}/{}.fits".format(DIR_WORK, flg.get_keyword("FILENAME")), overwrite=True)
    # save header
    img[1].header.tofile("{}/{}.head".format(DIR_WORK, img.get_keyword("FILENAME").replace(".fits", "")),
                         overwrite=True)

    # append img, wht, flg list
    img_list.append(img)
    wht_list.append(wht)
    flg_list.append(flg)
    fn_list.append(fp_img)


# Step 2. Calibrate Position
pcProc = CsstProcMscPositionCalibration()
pcProc.run(img_list, wht_list, flg_list, fn_list, DIR_GAIA_CATALOG, DIR_WORK, 2.0)
# if img_list:
#     pcProc.run(img_list, wht_list, flg_list, fn_list, DIR_GAIA_CATALOG, DIR_WORK, 2.0)
# else:
#     for i_ccd in CCD_ID_LIST:
#         fp_img = "{}/MSC_MS_*_{:02}_img.fits".format(DIR_WORK, i_ccd)
#         fp_wht = "{}/MSC_MS_*_{:02}_wht.fits".format(DIR_WORK, i_ccd)
#         fp_flg = "{}/MSC_MS_*_{:02}_flg.fits".format(DIR_WORK, i_ccd)
#         img = CsstMscImgData.read(fp_img)
#         wht = CsstMscImgData.read(fp_wht)
#         flg = CsstMscImgData.read(fp_flg)
#         img_list.append(img)
#         wht_list.append(wht)
#         flg_list.append(flg)
#         fn_list.append(fp_img)
#     pcProc.run(img_list, wht_list, flg_list, fn_list, DIR_GAIA_CATALOG, DIR_WORK, 2.0)
pcProc.cleanup(img_list, DIR_WORK)


# Step 3. Calibrate Flux
from csst.msc.calib_flux import CsstProcFluxCalibration
fcProc = CsstProcFluxCalibration()
#fcProc.prepare()
fcProc.run(fn_list,img_list, wht_list, flg_list,wcsdir=DIR_WORK,L1dir=DIR_WORK,workdir=DIR_WORK,refdir=DIR_TEST,addhead=True,morehead=False,plot=False,nodel=False,update=False,upcat=True)
fcProc.cleanup(fn_list,DIR_WORK)
