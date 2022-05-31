from pathlib import Path

import numpy as np
from ccdproc import cosmicray_lacosmic
from deepCR import deepCR
from collections import OrderedDict

from ..core.processor import CsstProcessor, CsstProcStatus
from ..msc import CsstMscImgData
from .. import PACKAGE_PATH


DEEPCR_MODEL_PATH = PACKAGE_PATH + "/msc/deepcr_model/CSST_2021-12-30_CCD23_epoch20.pth"


class CsstMscInstrumentProc(CsstProcessor):
    _status = CsstProcStatus.empty
    _switches = {'deepcr': True, 'clean': False}

    def __init__(self):
        super(CsstMscInstrumentProc).__init__()

        self.__img = None
        self.__wht = None
        self.__flg = None

    def _do_fix(self, raw, bias, dark, flat):
        '''仪器效应改正

        将raw扣除本底, 暗场, 平场. 并且避免了除0
        Args:
            raw: 科学图生图
            bias: 本底
            dark: 暗场
            flat: 平场
            exptime: 曝光时长
        '''
        self.__img = np.divide(
            raw.data - bias - dark * raw.exptime, flat,
            out=np.zeros_like(raw.data, float),
            where=(flat != 0),
        )

    def _do_badpix(self, flat):
        '''坏像元标记

        因为探测器本身原因造成不可用于有效科学研究的像元, 像元响应小于0.5中值或大于1.5中值

        Args:
            flat: 平场
        '''
        # TODO: flat-bias?
        med = np.median(flat)
        flg = (flat < 0.5 * med) | (1.5 * med < flat)
        self.__flg = self.__flg | (flg * 1)

    def _do_hot_and_warm_pix(self, dark, exptime, rdnoise):
        '''热像元与暖像元标记

        因为探测器本身原因造成的影响科学研究结果的像元. 像元在曝光时长积分时间内, 
        热像元为: 暗流计数大于探测器平均读出噪声的平方.
        暖像元为: 暗流计数大于0.5倍读出噪声的平方, 但小于1倍读出噪声的平方.

        Args:
            dark: 暗场
            exptime: 曝光时长
            rdnoise: 读出噪声
        '''
        tmp = dark * exptime
        tmp[tmp < 0] = 0
        flg = 1 * rdnoise ** 2 <= tmp  # 不确定是否包含 暂定包含
        self.__flg = self.__flg | (flg * 2)
        flg = (0.5 * rdnoise ** 2 < tmp) & (tmp < 1 * rdnoise ** 2)
        self.__flg = self.__flg | (flg * 4)

    def _do_over_saturation(self, raw):
        '''饱和溢出像元标记

        饱和像元及流量溢出污染的像元.

        Args:
            raw: 科学图生图
        '''
        flg = raw.data == 65535
        self.__flg = self.__flg | (flg * 8)

    def _do_cray(self, gain, rdnoise):
        '''宇宙线像元标记

        宇宙线污染的像元, 用deepCR或lacosmic进行标记
        '''

        if self._switches['deepcr']:
            clean_model = DEEPCR_MODEL_PATH
            inpaint_model = 'ACS-WFC-F606W-2-32'
            model = deepCR(clean_model, inpaint_model, device='CPU', hidden=50)
            masked, cleaned = model.clean(
                self.__img, threshold=0.5, inpaint=True, segment=True, patch=256, parallel=False, n_jobs=self.n_jobs)
        else:
            cleaned, masked = cosmicray_lacosmic(ccd=self.__img,
                                                 sigclip=3.,  # cr_threshold
                                                 sigfrac=0.5,  # neighbor_threshold
                                                 objlim=5.,  # constrast
                                                 gain=gain,
                                                 readnoise=rdnoise,
                                                 satlevel=65535.0,
                                                 pssl=0.0,
                                                 niter=4,
                                                 sepmed=True,
                                                 cleantype='meanmask',
                                                 fsmode='median',
                                                 psfmodel='gauss',
                                                 psffwhm=2.5,
                                                 psfsize=7,
                                                 psfk=None,
                                                 psfbeta=4.765,
                                                 verbose=False,
                                                 gain_apply=True)

        self.__flg = self.__flg | (masked * 16)
        if self._switches['clean']:
            self.__img = cleaned

    def _do_weight(self, bias, gain, rdnoise, exptime):
        '''权重图

        Args:
            bias: 本底
            gain: 增益
            rdnoise: 读出噪声
            exptime: 曝光时长
        '''
        data = self.__img.copy()
        data[self.__img < 0] = 0
        weight_raw = 1. / (gain * data + rdnoise ** 2)
        bias_weight = np.std(bias)
        weight = 1. / (1. / weight_raw + 1. / bias_weight) * exptime ** 2
        weight[self.__flg > 0] = 0
        self.__wht = weight

    def prepare(self, n_jobs=2, **kwargs):
        self.n_jobs = n_jobs
        for name in kwargs:
            self._switches[name] = kwargs[name]

    def run(self, raw: CsstMscImgData, bias: np.ndarray, dark: np.ndarray, flat: np.ndarray):

        assert isinstance(raw, CsstMscImgData)
        self.__img = np.copy(raw.data)
        self.__wht = np.zeros_like(raw.data, dtype=np.float32)
        self.__flg = np.zeros_like(raw.data, dtype=np.uint16)

        exptime = raw.get_keyword('EXPTIME', hdu=0)
        gain = raw.get_keyword('GAIN1', hdu=1)
        rdnoise = raw.get_keyword('RDNOISE1', hdu=1)

        # Flat and bias correction
        self._do_fix(raw, bias, dark, flat)
        self._do_badpix(flat)
        self._do_hot_and_warm_pix(dark, exptime, rdnoise)
        self._do_over_saturation(raw)
        self._do_cray(gain, rdnoise)
        self._do_weight(bias, gain, rdnoise, exptime)

        print('finish the run and save the results back to CsstData')

        img = raw.deepcopy(name="SCI", data=self.__img.astype(np.float32))
        wht = raw.deepcopy(name="WHT", data=self.__wht.astype(np.float32))
        flg = raw.deepcopy(name="FLG", data=self.__flg.astype(np.uint16))
        img.set_keyword("FILENAME", img.get_keyword("FILENAME", hdu=0).replace("_raw", "_img"), hdu=0)
        wht.set_keyword("FILENAME", wht.get_keyword("FILENAME", hdu=0).replace("_raw", "_wht"), hdu=0)
        flg.set_keyword("FILENAME", flg.get_keyword("FILENAME", hdu=0).replace("_raw", "_flg"), hdu=0)

        return img, wht, flg

    def cleanup(self):
        self.__img = None
        self.__wht = None
        self.__flg = None
        pass
