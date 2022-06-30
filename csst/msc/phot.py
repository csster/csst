from pathlib import Path
import joblib
import numpy as np
from ccdproc import cosmicray_lacosmic
from deepCR import deepCR

from ..core.processor import CsstProcessor, CsstProcStatus
from ._photometry.csst_photometry import do_phot
from ..msc import CsstMscImgData
from .data_manager import CsstMscDataManager
from .. import PACKAGE_PATH


# DEEPCR_MODEL_PATH = PACKAGE_PATH + "/msc/deepcr_model/CSST_2021-12-30_CCD23_epoch20.pth"


class CsstMscPhotometryProc(CsstProcessor):

    def __init__(self, dm: CsstMscDataManager):
        super(CsstProcessor, self).__init__()
        self.dm = dm

    def prepare(self, **kwargs):
        # check whether SExtractor and PSFEx exists
        pass

    def cleanup(self):
        pass

    def run(self, n_jobs=1, verbose=5):
        joblib.Parallel(n_jobs=n_jobs, verbose=verbose)(
            joblib.delayed(do_phot)(
                fitsfile=self.dm.l1_sci(ccd_id=ccd_id, suffix="img_L1", ext="fits"),
                outdir=self.dm.dir_l1,
                psffile=self.dm.l1_sci(ccd_id=ccd_id, suffix="psf", ext="fits"),
                catfile=self.dm.l1_sci(ccd_id=ccd_id, suffix="cat", ext="fits"),
            ) for ccd_id in self.dm.target_ccd_ids)
        return
