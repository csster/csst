# not finish yet
import numpy as np
from astropy.io import fits

from ..core.processor import CsstProcessor, CsstProcStatus


class CsstMscRefProc(CsstProcessor):
    _status = CsstProcStatus.empty

    def __init__(self):
        super(CsstMscRefProc).__init__()
        self.bias = None
        self.dark = None
        self.flat = None

    def array_combine(self, ndarray, mode="mean") -> np.ndarray:
        """ Function to combine 3-D data array

        Parameters
        ----------
        ndarray: array, input data cube (3D)
        model: mean, median, sum, mean_clip, median_clip, default is mean
        """
        if mode == "median":
            array = np.median(ndarray, axis=0)
        elif mode == "median_clip":
            ndarray = np.sort(ndarray, axis=0)[1:-1]
            array = np.median(ndarray, axis=0)
        elif mode == "sum":
            array = np.sum(ndarray, axis=0)
        elif mode == "mean":
            array = np.mean(ndarray, axis=0)
        elif mode == "mean_clip":
            ndarray = np.sort(ndarray, axis=0)[1:-1]
            array = np.mean(ndarray, axis=0)
        return array

    def load_bias(self, path: str) -> np.ndarray:
        with fits.open(path) as hdul:
            du = hdul[1].data
        du = du.astype(int)
        return du

    def load_dark(self, path: str) -> np.ndarray:
        with fits.open(path) as hdul:
            du = hdul[1].data
            hu = hdul[0].header
        du = du.astype(int)
        du = du - self.bias
        du = du / hu["EXPTIME"]
        return du

    def load_flat(self, path: str) -> np.ndarray:
        with fits.open(path) as hdul:
            du = hdul[1].data
            hu = hdul[0].header
        du = du.astype(int)
        du = du - self.bias - self.dark * hu["EXPTIME"]
        du = du / hu["EXPTIME"]
        du = du / np.median(du)
        return du

    def combine(self, func, mode: str, path_list, *args) -> np.ndarray:
        du_list = [func(path, *args) for path in path_list]
        du = self.array_combine(du_list, mode)
        return du

    def prepare(self, b_p_lst, d_p_lst, f_p_lst, save_path, mode_list=["median", "median", "median", ]):
        """

        Parameters
        ----------
        b_p_lst:
            List of currently ccd number bias file path
        d_p_lst:
            List of currently ccd number dark file path
        f_p_lst:
            List of currently ccd number flat file path
        save_path:
            as u c
        mode_list:
            [0] bias combine mode
            [1] dark combine mode
            [2] flat combine mode
            mean, median, sum, mean_clip, median_clip
        """
        self.b_p_lst = b_p_lst
        self.d_p_lst = d_p_lst
        self.f_p_lst = f_p_lst
        self.save_path = save_path
        self.mode_list = mode_list

    def run(self):
        self.bias = self.combine(self.load_bias, self.mode_list[0], self.b_p_lst)
        self.dark = self.combine(self.load_dark, self.mode_list[1], self.d_p_lst)
        self.flat = self.combine(self.load_flat, self.mode_list[2], self.f_p_lst)
        return self.bias, self.dark, self.flat

    def cleanup(self):
        self.bias = None
        self.dark = None
        self.flat = None
        pass
