from pathlib import Path
import joblib
import numpy as np
from ccdproc import cosmicray_lacosmic
from deepCR import deepCR

from ..core.processor import CsstProcessor, CsstProcStatus
from ._photometry.csst_photometry import do_phot
from ..msc import CsstMscImgData
from .. import PACKAGE_PATH


# DEEPCR_MODEL_PATH = PACKAGE_PATH + "/msc/deepcr_model/CSST_2021-12-30_CCD23_epoch20.pth"


class CsstMscPhotometryProc(CsstProcessor):

    def __init__(self):
        super(CsstProcessor, self).__init__()

    def prepare(self, **kwargs):
        # check whether SExtractor and PSFEx exists
        pass

    def cleanup(self):
        pass

    def run(self, fl_in, out_dir, n_jobs=1, verbose=5):
        # assert len(fl_in) == len(fl_out)
        joblib.Parallel(n_jobs=n_jobs, verbose=verbose)(
            joblib.delayed(do_phot)(fl_in[i], outdir=out_dir) for i in range(len(fl_in)))
        return
