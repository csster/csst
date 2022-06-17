import glob
import os

from csst.msc.backbone import CCD_ID_LIST, VER_SIMS
from csst.msc.calib_flux import CsstProcFluxCalibration
from csst.msc.calib_pos import CsstProcMscPositionCalibration
from csst.msc.data import CsstMscImgData
from csst.msc.inst_corr import CsstMscInstrumentProc
from csst.msc.phot import CsstMscPhotometryProc
from csst.msc.data_manager import CsstMscDataManager


CONFIG_DANDELION = dict(
    # test and working directory
    dir_l0="/home/csstpipeline/L1Pipeline/msc/MSC_0000020",
    dir_l1="/home/csstpipeline/L1Pipeline/msc/work/",
    # on Dandelion
    path_aux="/home/csstpipeline/L1Pipeline/msc/ref/MSC_{}_*_{:02d}_combine.fits",
    # gaia catalog directory (for position calibration)
    dir_pcref="/home/csstpipeline/L1Pipeline/msc/gaia_dr3/",
    # version of simulation data
    ver_sim="C3",
    # only 18 cores available in cloud machine from PMO
    n_jobs=18,
    # shut down backend multithreading
    backend_multithreading=False
)


CONFIG_PMO = dict(
    # test and working directory
    dir_l0="/share/simudata/CSSOSDataProductsSims/data/CSSTSimImage_C5/NGP_AstrometryON_shearOFF/MSC_0000100",
    dir_l1="/home/user/L1Pipeline/msc/work/",
    # on PMO
    path_aux="/data/sim_data/MSC_0000100/ref/MSC_{}_*_{:02d}_combine.fits",
    # gaia catalog directory (for position calibration)
    dir_pcref="/home/user/L1Pipeline/msc/gaia_dr3/",
    # version of simulation data
    ver_sim="C5.1",
    # only 18 cores available in cloud machine from PMO
    n_jobs=18,
    # shut down backend multithreading
    backend_multithreading=False
)


def do_one_exposure(ver_sim="C5.1", dir_l0="", dir_l1="", dir_pcref="", path_aux="", ccd_ids=None,
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

    # set paths
    dm = CsstMscDataManager(ver_sim=ver_sim, dir_l0=dir_l0, dir_l1=dir_l1, dir_pcref=dir_pcref, path_aux=path_aux)
    # request for ccd_ids
    ccd_ids = dm.get_ccd_ids(ccd_ids)

    # Step 1. Correct instrumental effect
    os.chdir(dir_l1)

    img_list = []
    wht_list = []
    flg_list = []
    fn_list = []
    for this_ccd_id in ccd_ids:
        print("processing CCD {}".format(this_ccd_id))
        fp_raw = dm.l0_sci(ccd_id=this_ccd_id)

        # read data with CsstMscImgData.read
        raw = CsstMscImgData.read(fp_raw)

        # in the future, use get_* functions grab
        bias = dm.get_bias(this_ccd_id)
        dark = dm.get_dark(this_ccd_id)
        flat = dm.get_flat(this_ccd_id)

        # initialize Instrument Processor
        instProc = CsstMscInstrumentProc()
        instProc.prepare(n_jobs=n_jobs)
        img, wht, flg = instProc.run(raw, bias, dark, flat, ver_sim)
        instProc.cleanup()
        fp_img = img[0].header["FILENAME"] + '.fits'

        # save img, wht, flg to somewhere
        img.writeto(dm.l1_sci(ccd_id=this_ccd_id, suffix="img", ext="fits"), overwrite=True)
        wht.writeto(dm.l1_sci(ccd_id=this_ccd_id, suffix="wht", ext="fits"), overwrite=True)
        flg.writeto(dm.l1_sci(ccd_id=this_ccd_id, suffix="flg", ext="fits"), overwrite=True)
        # save header
        img[1].header.tofile(dm.l1_sci(ccd_id=this_ccd_id, suffix="img", ext="head"), overwrite=True)

        # append img, wht, flg list
        img_list.append(img)
        wht_list.append(wht)
        flg_list.append(flg)
        fn_list.append(fp_img)

    # Step 2. Calibrate Position
    pcProc = CsstProcMscPositionCalibration()
    pcProc.run(img_list, wht_list, flg_list, fn_list, dir_pcref, dir_l1, 2.0)
    pcProc.cleanup(img_list, dir_l1)
    """
    pcProc = CsstProcMscPositionCalibration()
    pcProc.prepare(dm)
    pcProc.run(2.0)
    pcProc.cleanup()
    
    get these parameters from dm.l1_sci(*):
    img_list, wht_list, flg_list, fn_list, dir_pcref, dir_l1
    img_list, dir_l1
    """

    # Step 3. Calibrate Flux
    fcProc = CsstProcFluxCalibration()
    # fcProc.prepare()
    fcProc.run(
        fn_list, img_list, wht_list, flg_list, wcsdir=dir_l1, L1dir=dir_l1, workdir=dir_l1, refdir=dir_l0,
        addhead=True, morehead=False, plot=False, nodel=False, update=False, upcat=True)
    fcProc.cleanup(fn_list, dir_l1)

    """
    fcProc = CsstProcFluxCalibration()
    fcProc.prepare(dm)
    fcProc.run(addhead=True, morehead=False, plot=False, nodel=False, update=False, upcat=True)
    fcProc.cleanup()
    
    get these parameters from dm.l1_sci(*):
    fn_list, img_list, wht_list, flg_list, wcsdir=dir_l1, L1dir=dir_l1, workdir=dir_l1, refdir=dir_l0,
    fn_list, dir_l1
    """

    # Step 4. Photometry
    # ptProc = CsstMscPhotometryProc()
    # ptProc.prepare()
    # ptProc.run(fn_list, out_dir=dir_l1, n_jobs=n_jobs)
    # ptProc.cleanup()

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

