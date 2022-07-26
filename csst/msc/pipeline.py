import glob
import joblib
import os

from csst.msc.backbone import CCD_ID_LIST, VER_SIMS
from csst.msc.calib_flux import CsstProcFluxCalibration
from csst.msc.calib_pos import CsstProcMscPositionCalibration
from csst.msc.data import CsstMscImgData
from csst.msc.inst_corr import CsstMscInstrumentProc
from csst.msc.phot import CsstMscPhotometryProc
from csst.msc.data_manager import CsstMscDataManager

# dandelion, csstpipeline, C3
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

# pmo, user, C3
CONFIG_150s = dict(
    # test and working directory
    dir_l0="/data/L1Pipeline/msc/150s",
    dir_l1="{}/L1Pipeline/msc/work_150s/".format(os.getenv("HOME")),
    # on Dandelion
    path_aux="/data/L1Pipeline/msc/ref/MSC_{}_*_{:02d}_combine.fits",
    # gaia catalog directory (for position calibration)
    dir_pcref="/data/L1Pipeline/msc/gaia_dr3/",
    # version of simulation data
    ver_sim="C3",
    # only 18 cores available in cloud machine from PMO
    n_jobs=18,
    # shut down backend multithreading
    backend_multithreading=False,
    ccd_ids=CCD_ID_LIST,
)

# dandelion, cham, C3
CONFIG_CHAM = dict(
    # test and working directory
    dir_l0="/data/L1Pipeline/msc/MSC_0000020",
    dir_l1="{}/L1Pipeline/msc/work/".format(os.getenv("HOME")),
    # on Dandelion
    path_aux="/data/L1Pipeline/msc/ref/MSC_{}_*_{:02d}_combine.fits",
    # gaia catalog directory (for position calibration)
    dir_pcref="/data/L1Pipeline/msc/gaia_dr3/",
    # version of simulation data
    ver_sim="C3",
    # only 18 cores available in cloud machine from PMO
    n_jobs=18,
    # shut down backend multithreading
    backend_multithreading=False
)

# dandelion, anyuser, C5.2
CONFIG_ANYUSER_DDL_C52 = dict(
    # test and working directory
    dir_l0="/nfsdata/share/csst_simulation_data/Cycle-5-SimuData/multipleBandsImaging/NGP_AstrometryON_shearOFF/MSC_0000100",
    dir_l1="{}/L1Pipeline/msc/work/".format(os.getenv("HOME")),
    # on Dandelion
    path_aux="/nfsdata/users/cham/L1Test/ref_C5.2/MSC_{}_*_{:02d}_combine.fits",
    # gaia catalog directory (for position calibration)
    dir_pcref="/data/L1Pipeline/msc/gaia_dr3/",
    # version of simulation data
    ver_sim="C5.2",
    # only 18 cores available in cloud machine from PMO
    n_jobs=18,
    # shut down backend multithreading
    backend_multithreading=False
)

# pmo, user, C5.2
CONFIG_PMO = dict(
    # test and working directory
    # dir_l0="/share/simudata/CSSOSDataProductsSims/data/CSSTSimImage_C5/NGP_AstrometryON_shearOFF/MSC_0000100",
    dir_l0="/data/sim_data/new/MSC_0000100",  # C5.2 new C5 simulation data
    dir_l1="/home/user/L1Pipeline/msc/work_C5.2",
    # on PMO
    path_aux="/data/sim_data/MSC_0000100/ref/MSC_{}_*_{:02d}_combine.fits",
    # gaia catalog directory (for position calibration)
    dir_pcref="/home/user/L1Pipeline/msc/gaia_dr3/",
    # version of simulation data
    ver_sim="C5.2",
    # only 18 cores available in cloud machine from PMO
    n_jobs=18,
    # shut down backend multithreading
    backend_multithreading=False
)


def do_one_exposure(ver_sim="C5.2", dir_l0="", dir_l1="", dir_pcref="", path_aux="", ccd_ids=None,
                    n_jobs=18, backend_multithreading=False, runproc=(1, 0, 0, 0), dumpfile="test.dump"):
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
    # set target ccd_ids
    dm.set_ccd_ids(ccd_ids)

    # Step 1. Correct instrumental effect
    if not os.path.exists(dir_l1):
        print("@pipeline: dir_l1 does not exist, making {}".format(dir_l1))
        os.mkdir(dir_l1)
    print("@pipeline: cd {}".format(dir_l1))
    os.chdir(dir_l1)

    if runproc[0]:
        print("@pipeline: run instrument correction [1/4]")

        img_list = []
        wht_list = []
        flg_list = []
        fn_list = []
        for this_ccd_id in dm.target_ccd_ids:
            print("processing CCD {}".format(this_ccd_id))
            fp_raw = dm.l0_ccd(ccd_id=this_ccd_id)

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
            # fp_img = img[0].header["FILENAME"] + '.fits'

            # save img, wht, flg to somewhere
            img.writeto(dm.l1_ccd(ccd_id=this_ccd_id, post="img.fits"), overwrite=True)
            wht.writeto(dm.l1_ccd(ccd_id=this_ccd_id, post="wht.fits"), overwrite=True)
            flg.writeto(dm.l1_ccd(ccd_id=this_ccd_id, post="flg.fits"), overwrite=True)
            # save header
            img[1].header.tofile(dm.l1_ccd(ccd_id=this_ccd_id, post="img.head"), overwrite=True)

            # append img, wht, flg list
            img_list.append(img)
            wht_list.append(wht)
            flg_list.append(flg)
            # fn_list.append(fp_img)

            joblib.dump((img_list, wht_list, flg_list, dm), dumpfile)
    else:
        print("@pipeline: skip instrument correction [1/4]")
        (img_list, wht_list, flg_list, dm) = joblib.load(dumpfile)

    # Step 2. Calibrate Position
    if runproc[1]:
        print("@pipeline: run position calibration [2/4]")
        pcProc = CsstProcMscPositionCalibration(dm, n_jobs=n_jobs)
        # pcProc.prepare()
        pcProc.run(img_list, wht_list, flg_list)
        # pcProc.cleanup()
    else:
        # this stage saves data to files, no need to dump variables
        print("@pipeline: skip position calibration [2/4]")

    # Step 3. Calibrate Flux
    if runproc[2]:
        print("@pipeline: run flux calibration [3/4]")
        fcProc = CsstProcFluxCalibration(dm)
        # fcProc.prepare()
        fcProc.run(addhead=True, morehead=False, plot=False, nodel=False, update=False, upcat=True)
        # fcProc.cleanup(fn_list, dm.dir_l1)
    else:
        print("@pipeline: skip flux calibration [3/4]")

    # Step 4. Photometry
    if runproc[3]:
        print("@pipeline: run photometry [4/4]")
        # fn_list = [os.path.basename(dm.l1_sci(ccd_id=_, suffix="img_L1", ext="fits")) for _ in dm.target_ccd_ids]
        ptProc = CsstMscPhotometryProc(dm)
        # ptProc.prepare()
        ptProc.run(n_jobs=n_jobs)
        # ptProc.cleanup()
    else:
        print("@pipeline: skip photometry [4/4]")

    # print("@pipeline: dump DM object to dm.dump ...")
    # joblib.dump(dm, "dm.dump")

    return dm


if __name__ == "__main__":
    # identify where you are
    HOSTNAME = os.uname()[1]
    # you have to run this pipeline in some well-defined servers
    assert HOSTNAME in ["ubuntu", "dandelion"]
    # get config parameters
    # if HOSTNAME == "ubuntu":
    #     config = CONFIG_PMO
    # elif HOSTNAME == "Dandelion":
    #     config = CONFIG_DANDELION
    # else:
    #     raise ValueError("HOSTNAME {} not known!".format(HOSTNAME))
    # for k, v in config.items():
    #     eval("{}=config[\"{}\"]".format(k, k))

    # process this exposure
    # do_one_exposure(runproc=(1, 1, 0, 0), **CONFIG_CHAM)
    # do_one_exposure(runproc=(0, 0, 1, 0), **CONFIG_CHAM)

    # # on Dandelion, test another C3 data
    # from csst.msc.pipeline import do_one_exposure, CONFIG_150s
    # print(CONFIG_150s)
    # do_one_exposure(runproc=(1, 1, 1, 1), **CONFIG_150s)
    #
    # # on PMO, test C5.2 data
    # from csst.msc.pipeline import do_one_exposure, CONFIG_PMO
    # print(CONFIG_PMO)
    # do_one_exposure(runproc=(1, 1, 1, 1), **CONFIG_PMO)

    # on DDL, test C5.2 data
    from csst.msc.pipeline import do_one_exposure, CONFIG_ANYUSER_DDL_C52
    print(CONFIG_ANYUSER_DDL_C52)
    do_one_exposure(runproc=(1, 1, 0, 0), **CONFIG_ANYUSER_DDL_C52)

