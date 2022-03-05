from collections import OrderedDict
from abc import ABCMeta, abstractmethod
from enum import Enum
import numpy as np
from deepCR import deepCR
from ..core.processor import CsstProcStatus, CsstProcessor


class CsstMscInstrumentProc(CsstProcessor):
    _status = CsstProcStatus.empty
    _switches = {'crosstalk': False, 'nonlinear': False,
                 'deepcr': False, 'cti': False, 'brighterfatter': False}

    def __init__(self):
        pass

    def _do_fix(self, raw, bias, dark, flat, exptime):
        '''仪器效应改正

        将raw扣除本底, 暗场, 平场. 并且避免了除0
        Args:
            raw: 科学图生图
            bias: 本底
            dark: 暗场
            flat: 平场
            exptime: 曝光时长
        '''
        self.__l1img = np.divide(
            raw - bias - dark * exptime, flat,
            out=np.zeros_like(raw, float),
            where=(flat != 0),
        )

    def _do_badpix(self, flat):
        '''坏像元标记

        因为探测器本身原因造成不可用于有效科学研究的像元, 像元响应小于0.5中值或大于1.5中值

        Args:
            flat: 平场
        '''
        med = np.median(flat)
        flg = (flat < 0.5 * med) | (1.5 * med < flat)
        self.__flagimg = self.__flagimg | (flg * 1)

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
        self.__flagimg = self.__flagimg | (flg * 2)
        flg = (0.5 * rdnoise ** 2 < tmp) & (tmp < 1 * rdnoise ** 2)
        self.__flagimg = self.__flagimg | (flg * 4)

    def _do_over_saturation(self, raw):
        '''饱和溢出像元标记

        饱和像元及流量溢出污染的像元.

        Args:
            raw: 科学图生图
        '''
        flg = raw == 65535
        self.__flagimg = self.__flagimg | (flg * 8)

    def _do_cray(self):
        '''宇宙线像元标记

        宇宙线污染的像元, 用deepCR进行标记
        '''


        
        flg, _ = deepCR(self.__l1img)



        self.__flagimg = self.__flagimg | (flg * 16)

    def _do_weight(self, bias, gain, rdnoise, exptime):
        '''权重图

        Args:
            bias: 本底
            gain: 增益
            rdnoise: 读出噪声
            exptime: 曝光时长
        '''
        data = self.__l1img.copy()
        data[self.__l1img < 0] = 0
        weight_raw = 1. / (gain * data + rdnoise ** 2)
        bias_weight = np.std(bias)
        weight = 1. / (1. / weight_raw + 1. / bias_weight) * exptime ** 2
        weight[self.__flagimg > 0] = 0
        self.__weightimg = weight

    def prepare(self, **kwargs):
        for name in kwargs:
            self._switches[name] = kwargs[name]

    def run(self, data):
        if type(data).__name__ == 'CsstMscImgData' or type(data).__name__ == 'CsstMscSlsData':
            raw = data.get_l0data()
            self.__l1img = raw.copy()
            self.__weightimg = np.zeros_like(raw)
            self.__flagimg = np.zeros_like(raw, dtype=np.uint16)

            exptime = data.get_l0keyword('pri', 'EXPTIME')
            gain = data.get_l0keyword('img', 'GAIN1')
            rdnoise = data.get_l0keyword('img' ,'RDNOISE1')
            flat = data.get_flat()
            bias = data.get_bias()
            dark = data.get_dark()
            print('Flat and bias correction')
            self._do_fix(raw, bias, dark, flat, exptime)
            self._do_badpix(flat)
            self._do_hot_and_warm_pix(dark, exptime, rdnoise)
            self._do_over_saturation(raw)
            # self._do_cray()
            self._do_weight(bias, gain, rdnoise, exptime)

            print('fake to finish the run and save the results back to CsstData')

            data.set_l1data('sci', self.__l1img)
            data.set_l1data('weight', self.__weightimg)
            data.set_l1data('flag', self.__flagimg)

            print('Update keywords')
            data.set_l1keyword('SOMEKEY', 'some value',
                               'Test if I can append the header')

            self._status = CsstProcStatus.normal
        else:
            self._status = CsstProcStatus.ioerror
        return self._status

    def cleanup(self):
        pass
