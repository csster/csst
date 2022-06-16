import glob
import os

from csst.msc.backbone import CCD_ID_LIST, VER_SIMS
from csst.msc.calib_flux import CsstProcFluxCalibration
from csst.msc.calib_pos import CsstProcMscPositionCalibration
from csst.msc.data import CsstMscImgData
from csst.msc.inst_corr import CsstMscInstrumentProc
from csst.msc.phot import CsstMscPhotometryProc

CONFIG_DANDELION = dict(
    # test and working directory
    dir_raw="/home/csstpipeline/L1Pipeline/msc/MSC_0000020",
    dir_work="/home/csstpipeline/L1Pipeline/msc/work/",
    filename_fmt="{}/MSC_MS_*_{:02}_raw.fits",
    # on Dandelion
    path_aux="/home/csstpipeline/L1Pipeline/msc/ref/MSC_{}_*_{:02d}_combine.fits",
    # gaia catalog directory (for position calibration)
    dir_gaia_catalog="/home/csstpipeline/L1Pipeline/msc/gaia_dr3/",
    # version of simulation data
    ver_sim="C3",
    # only 18 cores available in cloud machine from PMO
    n_jobs=18,
    # shut down backend multithreading
    backend_multithreading=False
)


CONFIG_PMO = dict(
    # test and working directory
    dir_raw="/share/simudata/CSSOSDataProductsSims/data/CSSTSimImage_C5/NGP_AstrometryON_shearOFF/MSC_0000100",
    dir_work="/home/user/L1Pipeline/msc/work/",
    filename_fmt='{}/CSST_MSC_MS_SCI_*_{:02d}_L0_1.fits',
    # on PMO
    path_aux="/data/sim_data/MSC_0000100/ref/MSC_{}_*_{:02d}_combine.fits",
    # gaia catalog directory (for position calibration)
    dir_gaia_catalog="/home/user/L1Pipeline/msc/gaia_dr3/",
    # version of simulation data
    ver_sim="C5.1",
    # only 18 cores available in cloud machine from PMO
    n_jobs=18,
    # shut down backend multithreading
    backend_multithreading=False
)


def do_one_exposure(dir_raw="", dir_work="", filename_fmt="", path_aux="", dir_gaia_catalog="", ver_sim="C5.1", ccd_ids=None,
                    n_jobs=18, backend_multithreading=False):

    # currently C3 and C5.1 are tested
    try:
        assert ver_sim in VER_SIMS
    except AssertionError as ae:
        print("Available options for ver_sim: ", VER_SIMS)
        raise ae

    # shut down backend multi-threading
    if not backend_multithreading:
        # prohibit multi-threading in backend
        os.environ["MKL_NUM_THREADS"] = '1'
        os.environ["NUMEXPR_NUM_THREADS"] = '1'
        os.environ["OMP_NUM_THREADS"] = '1'

    if ccd_ids is None:
        ccd_ids = CCD_ID_LIST

    # Step 1. Correct instrumental effect
    os.chdir(dir_work)

    img_list = []
    wht_list = []
    flg_list = []
    fn_list = []
    for i_ccd in ccd_ids:
        print("processing CCD {}".format(i_ccd))
        fp_raw = glob.glob(filename_fmt.format(dir_raw, i_ccd))
        assert len(fp_raw) == 1
        fp_raw = fp_raw[0]
        print(fp_raw)
        # fn_raw = os.path.basename(fp_raw)

        # read data with CsstMscImgData.read
        raw = CsstMscImgData.read(fp_raw)
        # in the future, use get_* functions grab
        bias = raw.get_bias(path_aux.format("CLB", i_ccd))
        dark = raw.get_dark(path_aux.format("CLD", i_ccd))
        flat = raw.get_flat(path_aux.format("CLF", i_ccd))

        # initialize Instrument Processor
        instProc = CsstMscInstrumentProc()
        instProc.prepare(n_jobs=n_jobs)
        img, wht, flg = instProc.run(raw, bias, dark, flat, ver_sim)
        instProc.cleanup()
        fp_img = img[0].header["FILENAME"] + '.fits'

        # save img, wht, flg to somewhere
        img.writeto("{}/{}.fits".format(dir_work, img.get_keyword("FILENAME")), overwrite=True)
        wht.writeto("{}/{}.fits".format(dir_work, wht.get_keyword("FILENAME")), overwrite=True)
        flg.writeto("{}/{}.fits".format(dir_work, flg.get_keyword("FILENAME")), overwrite=True)
        # save header
        img[1].header.tofile(
            "{}/{}.head".format(dir_work, img.get_keyword("FILENAME").replace(".fits", "")), overwrite=True)

        # append img, wht, flg list
        img_list.append(img)
        wht_list.append(wht)
        flg_list.append(flg)
        fn_list.append(fp_img)

    # # parallel step 1
    # def corr_one_img(i_ccd, dir_l0, dir_l1):
    #     fp_raw = glob.glob("{}/MSC_MS_*_{:02}_raw.fits".format(dir_l0, i_ccd))
    #     assert len(fp_raw) == 1
    #     fp_raw = fp_raw[0]
    #
    #     # read data with CsstMscImgData.read
    #     raw = CsstMscImgData.read(fp_raw)
    #
    #     # in future, get_* functions grab
    #     bias = raw.get_bias(PATH_BIAS.format(i_ccd))
    #     dark = raw.get_dark(PATH_DARK.format(i_ccd))
    #     flat = raw.get_flat(PATH_FLAT.format(i_ccd))
    #
    #     # initialize Instrument Processor
    #     instProc = CsstMscInstrumentProc()
    #     instProc.prepare(n_jobs=1, n_threads=1)
    #     img, wht, flg = instProc.run(raw, bias, dark, flat)
    #     instProc.cleanup()
    #     fp_img = img[0].header["FILENAME"] + '.fits'
    #
    #     # save img, wht, flg to somewhere
    #     img.writeto("{}/{}.fits".format(dir_l1, img.get_keyword("FILENAME")), overwrite=True)
    #     wht.writeto("{}/{}.fits".format(dir_l1, wht.get_keyword("FILENAME")), overwrite=True)
    #     flg.writeto("{}/{}.fits".format(dir_l1, flg.get_keyword("FILENAME")), overwrite=True)
    #     # save header
    #     img[1].header.tofile("{}/{}.head".format(dir_l1, img.get_keyword("FILENAME").replace(".fits", "")),
    #                          overwrite=True)
    #     return OrderedDict(img=img, wht=wht, flg=flg, fp_img=fp_img)
    #
    #
    # result = joblib.Parallel(n_jobs=NJOBS, verbose=5)(
    #     joblib.delayed(corr_one_img)(i_ccd, dir_l0, dir_l1) for i_ccd in CCD_ID_LIST)
    # img_list = [_["img"] for _ in result]
    # wht_list = [_["wht"] for _ in result]
    # flg_list = [_["flg"] for _ in result]
    # fn_list = [_["fp_img"] for _ in result]

    # Step 2. Calibrate Position
    pcProc = CsstProcMscPositionCalibration()
    pcProc.run(img_list, wht_list, flg_list, fn_list, dir_gaia_catalog, dir_work, 2.0)
    pcProc.cleanup(img_list, dir_work)
    # if img_list:
    #     pcProc.run(img_list, wht_list, flg_list, fn_list, dir_gaia_catalog, dir_l1, 2.0)
    # else:
    #     for i_ccd in CCD_ID_LIST:
    #         fp_img = "{}/MSC_MS_*_{:02}_img.fits".format(dir_l1, i_ccd)
    #         fp_wht = "{}/MSC_MS_*_{:02}_wht.fits".format(dir_l1, i_ccd)
    #         fp_flg = "{}/MSC_MS_*_{:02}_flg.fits".format(dir_l1, i_ccd)
    #         img = CsstMscImgData.read(fp_img)
    #         wht = CsstMscImgData.read(fp_wht)
    #         flg = CsstMscImgData.read(fp_flg)
    #         img_list.append(img)
    #         wht_list.append(wht)
    #         flg_list.append(flg)
    #         fn_list.append(fp_img)
    #     pcProc.run(img_list, wht_list, flg_list, fn_list, dir_gaia_catalog, dir_l1, 2.0)

    # Step 3. Calibrate Flux
    fcProc = CsstProcFluxCalibration()
    # fcProc.prepare()
    fcProc.run(
        fn_list, img_list, wht_list, flg_list, wcsdir=dir_work, L1dir=dir_work, workdir=dir_work, refdir=dir_raw,
        addhead=True, morehead=False, plot=False, nodel=False, update=False, upcat=True)
    fcProc.cleanup(fn_list, dir_work)

    # Step 4. Photometry
    ptProc = CsstMscPhotometryProc()
    ptProc.prepare()
    ptProc.run(fn_list, out_dir=dir_work, n_jobs=n_jobs)
    ptProc.cleanup()

    return


if __name__ == "__main__":
    # identify where you are
    HOSTNAME = os.uname()[1]
    # you have to run this pipeline in some well-defined servers
    assert HOSTNAME in ["ubuntu", "Dandelion"]
    # get config parameters
    if HOSTNAME == "ubuntu":
        config = CONFIG_PMO
    elif HOSTNAME == "Dandelion":
        config = CONFIG_DANDELION
    else:
        raise ValueError("HOSTNAME {} not known!".format(HOSTNAME))

    # process this exposure
    do_one_exposure(**config)

    for k, v in config.items():
        eval("{}=config[\"{}\"]".format(k, k))

