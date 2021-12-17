from collections import OrderedDict
from abc import ABCMeta, abstractmethod
from csst.common.data import CsstData
from enum import Enum
import numpy as np
from ..common.processors import CsstProcStatus, CsstProcessor


class CsstMscInstrumentProc(CsstProcessor):
    _status = CsstProcStatus.empty
    _switches = {'crosstalk': False, 'nonlinear': False, 'deepcr': False, 'cti': False, 'brighterfatter': False}

    def __init__(self):
        pass

    def _do_crosstalk(self):
        if self._switches['crosstalk']:
            print('Crosstalk correction')

    def _do_nonlinear(self):
        if self._switches['nonlinear']:
            print('Nonlinear effect correction')

    def _do_deepcr(self):
        if self._switches['deepcr']:
            print('Deep CR operation')
        else:
            print('Laplace CR correction')

    def _do_cti(self):
        if self._switches['cti']:
            print('CTI effect correction')

    def _do_brighterfatter(self):
        if self._switches['brighterfatter']:
            print('Brighter-Fatter effect correction')

    def prepare(self, **kwargs):
        for name in kwargs:
            self._switches[name] = kwargs[name]

    def run(self, data):
        if type(data).__name__ == 'CsstMscImgData' or type(data).__name__ == 'CsstMscSlsData':
            self.__l1img = data.get_l0data(copy=True)
            self.__weightimg = np.random.uniform(0, 1, (9216, 9232))
            self.__flagimg = np.random.uniform(0, 1, (9216, 9232))

            flat = data.get_flat()
            bias = data.get_bias()
            dark = data.get_dark()
            print('Flat and bias correction')
            self._do_crosstalk()
            self._do_nonlinear()
            self._do_cti()
            self._do_deepcr()
            self._do_brighterfatter()

            print('fake to finish the run and save the results back to CsstData')

            data.set_l1data('sci', self.__l1img)
            data.set_l1data('weight', self.__weightimg)
            data.set_l1data('flag', self.__flagimg)

            print('Update keywords')
            data.set_l1keyword('SOMEKEY', 'some value', 'Test if I can append the header')

            self._status = CsstProcStatus.normal
        else:
            self._status = CsstProcStatus.ioerror
        return self._status

    def cleanup(self):
        pass
