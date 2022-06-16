import glob
import os
import re
import logging

from .backbone import VER_SIMS, CCD_ID_LIST, CCD_FILTER_MAPPING


class CsstMscNamingRules:
    """ this class defines the file format of the input & output of CSST MSC L1 pipeline

    C3:
        MSC_MS_210525220000_100000020_06_raw.fits
        MSC_CRS_210525220000_100000020_06_raw.fits
        MSC_210525120000_0000020_06.cat

    C5.1:
        CSST_MSC_MS_SCI_20270810081950_20270810082220_100000100_06_L0_1.fits
        CSST_MSC_MS_CRS_20270810081950_20270810082220_100000100_06_L0_1.fits
        MSC_10000100_chip_06_filt_y.cat
        MSC_10000100_chip_06_filt_y.log

    """

    def __init__(self, ver_sim="C5.1", dir_l0="", dir_l1="", force_all_ccds=False):
        assert ver_sim in VER_SIMS

        self.dir_l0 = dir_l0
        self.dir_l1 = dir_l1
        self.ver_sim = ver_sim

        print("globbing files ... ", end="")
        fps = self.glob_dir(dir_l0, ver_sim=ver_sim)
        if force_all_ccds:
            assert len(fps) == len(CCD_ID_LIST)
        else:
            assert len(fps) > 0
        print("{} L0 images found!".format(len(fps)))

        if ver_sim == "C3":
            # get info
            # print(re.split(r"[_.]", fps[0]))
            self._instrument, self._survey, \
                self._exp_start, self._exp_id, \
                _ccd_id, self._l0_suffix, _ext = re.split(r"[_.]", fps[0])

            self._exp_start = int(self._exp_start)
            self._exp_id = int(self._exp_id)

            # available CCD IDs
            self.available_ccd_ids = [int(re.split(r"[_.]", fp)[4]) for fp in fps]
            self.available_ccd_ids.sort()

        elif ver_sim == "C5.1":
            # get info
            print(re.split(r"[_.]", fps[0]))
            self._telescope, self._instrument, self._survey, self._imagetype, \
                self._exp_start, self._exp_stop, self._exp_id, \
                _ccd_id, self._l0_suffix, self._version, _ext = re.split(r"[_.]", fps[0])
            self._exp_start = int(self._exp_start)
            self._exp_stop = int(self._exp_stop)
            self._exp_id = int(self._exp_id)
            # available CCD IDs
            self.available_ccd_ids = [int(re.split(r"[_.]", fp)[7]) for fp in fps]
            self.available_ccd_ids.sort()

    @staticmethod
    def glob_dir(dir_raw, ver_sim="C5"):
        if ver_sim == "C3":
            pattern = os.path.join(dir_raw, "MSC_MS_*_raw.fits")
        elif ver_sim == "C5.1":
            pattern = os.path.join(dir_raw, "CSST_MSC_MS_SCI_*.fits")
        fps = glob.glob(pattern)
        fps = [os.path.basename(fp) for fp in fps]
        fps.sort()
        return fps

    def pc_combined_image(self, ccd_id):
        if self.ver_sim == "C3":
            return ""

    @property
    def pc_scamp_coord(self):
        """ SCAMP coord """
        return "scamp_coord.txt"

    def l1_final_output_image(self, ccd_id):
        return

    def l0_cat(self, ccd_id=6):
        if self.ver_sim == "C3":
            fn = "MSC_{}_{:07d}_{:02d}.cat".format(
                self._exp_start, self._exp_id-100000000, ccd_id)
        elif self.ver_sim == "C5.1":
            fn = "MSC_{}_chip_{:02d}_filt_{}.cat".format(
                self._exp_id-90000000, ccd_id, CCD_FILTER_MAPPING[ccd_id])
        return os.path.join(self.dir_l0, fn)

    def l0_log(self, ccd_id=6):
        if self.ver_sim == "C5.1":
            fn = "MSC_{}_chip_{:02d}_filt_{}.log".format(
                self._exp_id - 90000000, ccd_id, CCD_FILTER_MAPPING[ccd_id])
            return os.path.join(self.dir_l0, fn)

    def l0_sci(self, ccd_id=6):
        if self.ver_sim == "C3":
            fn = "{}_{}_{}_{}_{:02d}_raw.fits".format(
                self._instrument, self._survey, self._exp_start, self._exp_id, ccd_id)
        elif self.ver_sim == "C5.1":
            fn = "{}_{}_{}_SCI_{}_{}_{}_{:02d}_L0_1.fits".format(
                self._telescope, self._instrument, self._survey,
                self._exp_start, self._exp_stop, self._exp_id, ccd_id)
        return os.path.join(self.dir_l0, fn)

    def l0_crs(self, ccd_id=6):
        if self.ver_sim == "C3":
            fn = "{}_CRS_{}_{}_{:02d}_raw.fits".format(
                self._instrument, self._exp_start, self._exp_id, ccd_id)
        elif self.ver_sim == "C5.1":
            fn = "{}_{}_{}_CRS_{}_{}_{}_{:02d}_L0_1.fits".format(
                self._telescope, self._instrument, self._survey,
                self._exp_start, self._exp_stop, self._exp_id, ccd_id)
        return os.path.join(self.dir_l0, fn)

    def l1_sci(self, ccd_id=6):
        if self.ver_sim == "C3":
            fn = "{}_{}_{}_{:02d}_L1.fits".format(
                self._instrument, self._survey, self._exp_start, self._exp_id, ccd_id)
        elif self.ver_sim == "C5.1":
            fn = "{}_{}_{}_SCI_{}_{}_{}_{:02d}_L1.fits".format(
                self._telescope, self._instrument, self._survey,
                self._exp_start, self._exp_stop, self._exp_id, ccd_id)
        return os.path.join(self.dir_l1, fn)

    def l0_aux(self, ccd_id=6, suffix="img"):
        if self.ver_sim == "C3":
            fn = "{}_{}_{}_{:02d}_{}.fits".format(
                self._instrument, self._survey,
                self._exp_start, self._exp_id, ccd_id, suffix)
        elif self.ver_sim == "C5.1":
            fn = "{}_{}_{}_SCI_{}_{}_{}_{:02d}_{}.fits".format(
                self._telescope, self._instrument, self._survey,
                self._exp_start, self._exp_stop, self._exp_id, ccd_id, suffix)
        return os.path.join(self.dir_l0, fn)


if __name__ == "__main__":
    from csst.msc.naming import CsstMscNamingRules
    nr = CsstMscNamingRules(ver_sim="C3", dir_l0="/data/L1Pipeline/msc/MSC_0000020", dir_l1="/data/L1Pipeline/msc/work")
    print(nr.l0_sci(ccd_id=6))
    print(os.path.exists(nr.l0_sci(ccd_id=6)))
    print(nr.l0_crs(ccd_id=6))
    print(os.path.exists(nr.l0_sci(ccd_id=8)))
    print(nr.l0_cat(8))
    print(os.path.exists(nr.l0_cat(ccd_id=8)))
    print(nr.l0_log(8))
    print(os.path.exists(nr.l0_log(ccd_id=8)))
    print(nr.l0_aux(8, "img"))
    print(nr.available_ccd_ids)
    print(nr.l1_sci(25))

    from csst.msc.naming import CsstMscNamingRules
    nr = CsstMscNamingRules(ver_sim="C5.1", dir_l0="/data/sim_data/MSC_0000100", dir_l1="/home/user/L1Pipeline/msc/work")
    print(nr.l0_sci(ccd_id=6))
    print(os.path.exists(nr.l0_sci(ccd_id=6)))
    print(nr.l0_crs(ccd_id=6))
    print(os.path.exists(nr.l0_sci(ccd_id=8)))
    print(nr.l0_cat(8))
    print(os.path.exists(nr.l0_cat(ccd_id=8)))
    print(nr.l0_log(8))
    print(os.path.exists(nr.l0_log(ccd_id=8)))
    print(nr.l0_aux(8, "img"))
    print(nr.available_ccd_ids)
    print(nr.l1_sci(25))
