import glob
import re

from astropy.io import fits

from .backbone import VER_SIMS, CCD_FILTER_MAPPING


class CsstMscDataManager:
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

    def __init__(self, ver_sim="C5.1", dir_l0="", dir_l1="", dir_pcref="", path_aux="", force_all_ccds=False):
        assert ver_sim in VER_SIMS

        self.dir_l0 = dir_l0
        self.dir_l1 = dir_l1
        self.dir_pcref = dir_pcref
        self.path_aux = path_aux
        self.ver_sim = ver_sim

        fps_img = self.glob_image(dir_l0, ver_sim=ver_sim)
        fps_cat = self.glob_cat(dir_l0, ver_sim=ver_sim)

        if force_all_ccds:
            assert len(fps_img) == len(CCD_ID_LIST)
        else:
            assert len(fps_img) > 0

        if ver_sim == "C3":
            # get info
            # print(re.split(r"[_.]", fps[0]))
            self._instrument, self._survey, \
                self._exp_start, self._exp_id, \
                _ccd_id, self._l0_suffix, _ext = re.split(r"[_.]", fps_img[0])
            self._cat_id = re.split(r"[_.]", fps_cat[0])[1]

            self._exp_start = int(self._exp_start)
            self._exp_id = int(self._exp_id)

            # available CCD IDs
            self.available_ccd_ids = [int(re.split(r"[_.]", fp)[4]) for fp in fps_img]
            self.available_ccd_ids.sort()

        elif ver_sim == "C5.1":
            # get info
            # print(re.split(r"[_.]", fps[0]))
            self._telescope, self._instrument, self._survey, self._imagetype, \
                self._exp_start, self._exp_stop, self._exp_id, \
                _ccd_id, self._l0_suffix, self._version, _ext = re.split(r"[_.]", fps_img[0])
            self._cat_id = re.split(r"[_.]", fps_cat[0])[1]

            self._exp_start = int(self._exp_start)
            self._exp_stop = int(self._exp_stop)
            self._exp_id = int(self._exp_id)

            # available CCD IDs
            self.available_ccd_ids = [int(re.split(r"[_.]", fp)[7]) for fp in fps_img]
            self.available_ccd_ids.sort()

    @staticmethod
    def glob_image(dir_l0, ver_sim="C5"):
        """ glob files in L0 data directory """
        if ver_sim == "C3":
            pattern = os.path.join(dir_l0, "MSC_MS_*_raw.fits")
        elif ver_sim == "C5.1":
            pattern = os.path.join(dir_l0, "CSST_MSC_MS_SCI_*.fits")
        fps = glob.glob(pattern)
        fps = [os.path.basename(fp) for fp in fps]
        fps.sort()

        print("@DM.glob_dir: {} files found with pattern: {}".format(len(fps), pattern))
        return fps

    @staticmethod
    def glob_cat(dir_l0, ver_sim="C5"):
        """ glob input catalogs in L0 data directory """
        if ver_sim == "C3":
            pattern = os.path.join(dir_l0, "MSC_*.cat")
        elif ver_sim == "C5.1":
            pattern = os.path.join(dir_l0, "MSC_*.cat")
        fps = glob.glob(pattern)
        fps = [os.path.basename(fp) for fp in fps]
        fps.sort()

        print("@DM.glob_dir: {} files found with pattern: {}".format(len(fps), pattern))
        return fps

    def l0_cat(self, ccd_id=6):
        """ the L0 cat file path"""
        if self.ver_sim == "C3":
            fn = "{}_{}_{:07d}_{:02d}.cat".format(
                self._instrument, self._cat_id, self._exp_id - 100000000, ccd_id)
        elif self.ver_sim == "C5.1":
            fn = "{}_{}_chip_{:02d}_filt_{}.cat".format(
                self._instrument, self._exp_id - 90000000, ccd_id, CCD_FILTER_MAPPING[ccd_id])
        return os.path.join(self.dir_l0, fn)

    def l0_log(self, ccd_id=6):
        """ L0 log file path """
        if self.ver_sim == "C5.1":
            fn = "{}_{}_chip_{:02d}_filt_{}.log".format(
                self._instrument, self._exp_id - 90000000, ccd_id, CCD_FILTER_MAPPING[ccd_id])
            return os.path.join(self.dir_l0, fn)

    def l0_sci(self, ccd_id=6):
        """ L0 image file path """
        if self.ver_sim == "C3":
            fn = "{}_{}_{}_{}_{:02d}_raw.fits".format(
                self._instrument, self._survey, self._exp_start, self._exp_id, ccd_id)
        elif self.ver_sim == "C5.1":
            fn = "{}_{}_{}_SCI_{}_{}_{}_{:02d}_L0_1.fits".format(
                self._telescope, self._instrument, self._survey,
                self._exp_start, self._exp_stop, self._exp_id, ccd_id)
        return os.path.join(self.dir_l0, fn)

    def l0_crs(self, ccd_id=6):
        """ L0 cosmic ray file path """
        if self.ver_sim == "C3":
            fn = "{}_CRS_{}_{}_{:02d}_raw.fits".format(
                self._instrument, self._exp_start, self._exp_id, ccd_id)
        elif self.ver_sim == "C5.1":
            fn = "{}_{}_{}_CRS_{}_{}_{}_{:02d}_L0_1.fits".format(
                self._telescope, self._instrument, self._survey,
                self._exp_start, self._exp_stop, self._exp_id, ccd_id)
        return os.path.join(self.dir_l0, fn)

    def l1_sci(self, ccd_id=6, suffix="img_whead", ext="fits"):
        """ generate L1 file path

        Parameters
        ----------
        ccd_id:
            CCD ID
        suffix:
            {"img", "wht", "flg", "img_L1", "wht_L1", "flg_L1", "whead"}
        ext:
            {"fits", "cat", "rcat"}

        Returns
        -------
        L1 file path

        """
        if self.ver_sim == "C3":
            fn = "{}_{}_{}_{}_{:02d}_{}.{}".format(
                self._instrument, self._survey,
                self._exp_start, self._exp_id, ccd_id, suffix, ext)
        elif self.ver_sim == "C5.1":
            fn = "{}_{}_{}_SCI_{}_{}_{}_{:02d}_{}.{}".format(
                self._telescope, self._instrument, self._survey,
                self._exp_start, self._exp_stop, self._exp_id, ccd_id, suffix, ext)
        return os.path.join(self.dir_l1, fn)

    def pc_combined_image(self, suffix="img", ext="fits"):
        """ combined images

        Parameters
        ----------
        suffix:
            {"img", "wht", "flg", "cat"}
        ext:
            {"fits", }

        Returns
        -------
        combined image path

        """
        fn = "combined_" + "{}.{}".format(suffix, ext)
        return os.path.join(self.dir_l1, fn)

    @property
    def pc_ref_cat(self):
        """ reference catalog """
        return os.path.join(self.dir_l1, "ref.cat")

    @property
    def pc_check_fits(self):
        """ check fits """
        return os.path.join(self.dir_l1, "check.fits")

    @property
    def pc_combined_head(self):
        """ combined head """
        return os.path.join(self.dir_l1, "combined_cat.head")

    @property
    def pc_combined_head_fits(self):
        """ combined head (converted to) fits """
        return os.path.join(self.dir_l1, "combined_cat.head.fits")

    @property
    def pc_scamp_coord(self):
        """ SCAMP coord """
        return os.path.join(self.dir_l1, "scamp_coord.txt")

    def get_ccd_ids(self, ccd_ids=None):
        """  """
        if ccd_ids is None:
            # default ccd_ids
            ccd_ids = self.available_ccd_ids
        else:
            try:
                # assert ccd_ids is a subset of available ccd_ids
                assert set(ccd_ids).issubset(set(self.available_ccd_ids))
            except AssertionError as ae:
                print("@DM: available CCD IDs are ", self.available_ccd_ids)
                print("@DM: target CCD IDs are ", ccd_ids)
                raise ae
        return ccd_ids

    def get_bias(self, ccd_id=6):
        fp = glob.glob(self.path_aux.format("CLB", ccd_id))[0]
        return fits.getdata(fp)

    def get_dark(self, ccd_id=6):
        fp = glob.glob(self.path_aux.format("CLD", ccd_id))[0]
        return fits.getdata(fp)

    def get_flat(self, ccd_id=6):
        fp = glob.glob(self.path_aux.format("CLF", ccd_id))[0]
        return fits.getdata(fp)


if __name__ == "__main__":
    # test C3
    import os
    from csst.msc.data_manager import CsstMscDataManager

    dm = CsstMscDataManager(
        ver_sim="C3", dir_l0="/data/L1Pipeline/msc/MSC_0000020", dir_l1="/data/L1Pipeline/msc/work")
    print("----- L0 images -----")
    print(dm.l0_sci(ccd_id=6))
    print(os.path.exists(dm.l0_sci(ccd_id=6)))
    print("----- L0 crs -----")
    print(dm.l0_crs(ccd_id=6))
    print(os.path.exists(dm.l0_sci(ccd_id=8)))
    print("----- L0 input cat -----")
    print(dm.l0_cat(8))
    print(os.path.exists(dm.l0_cat(ccd_id=8)))
    print("----- available ccd_ids -----")
    print(dm.available_ccd_ids)
    print("----- L1 images -----")
    print(dm.l1_sci(25, "img", "fits"))

    # test C5.1
    import os
    from csst.msc.data_manager import CsstMscDataManager
    from csst.msc.backbone import CCD_ID_LIST

    dm = CsstMscDataManager(
        ver_sim="C5.1", dir_l0="/data/sim_data/MSC_0000100", dir_l1="/home/user/L1Pipeline/msc/work")
    print("----- available ccd_ids -----")
    print(dm.available_ccd_ids)
    for ccd_id in dm.available_ccd_ids[:2]:
        print("----- L0 images -----")
        print(dm.l0_sci(ccd_id=ccd_id))
        print(os.path.exists(dm.l0_sci(ccd_id=ccd_id)))
        print("----- L0 crs -----")
        print(dm.l0_crs(ccd_id=ccd_id))
        print(os.path.exists(dm.l0_sci(ccd_id=ccd_id)))
        print("----- L0 input cat -----")
        print(dm.l0_cat(ccd_id=ccd_id))
        print(os.path.exists(dm.l0_cat(ccd_id=ccd_id)))
        print("----- L0 input log -----")
        print(dm.l0_log(ccd_id=ccd_id))
        print(os.path.exists(dm.l0_log(ccd_id=ccd_id)))
        print("----- L1 images -----")
        print(dm.l1_sci(ccd_id, "img", "fits"))
