from abc import ABCMeta, abstractmethod
from enum import Enum


class CsstProcStatus(Enum):
    empty = -1
    normal = 0
    ioerror = 1
    runtimeerror = 2


#     self['empty'].info = 'Not run yet.'
#     self['normal'].info = 'This is a normal run.'
#     self['ioerror'].info = 'This run is exceptionally stopped due to IO error.'
#     self['runtimeerror'].info = 'This run is exceptionally stopped due to runtime error.'

class CsstProcessor(metaclass=ABCMeta):

    def __init__(self, **kwargs):
        pass

    @abstractmethod
    def prepare(self, **kwargs):
        pass

    @abstractmethod
    def run(self, kwargs):
        """ """
        return self._status

    @abstractmethod
    def cleanup(self):
        pass
