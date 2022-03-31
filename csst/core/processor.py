from abc import ABC, ABCMeta, abstractmethod
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


class CsstProcessor(ABC):

    def __init__(self, **kwargs):
        # self._status = CsstProcStatus()
        pass

    @abstractmethod
    def prepare(self, **kwargs):
        # do your preparation here
        raise NotImplementedError

    @abstractmethod
    def run(self, **kwargs):
        # run your pipeline
        raise NotImplementedError

    @abstractmethod
    def cleanup(self):
        # clean up environment
        raise NotImplementedError


class CsstDemoProcessor(CsstProcessor):

    def __init__(self, **kwargs):
        super().__init__()

    def some_function(self, **kwargs):
        print("some function")

    def prepare(self):
        print("prepare")

    def run(self):
        print("run")

    def cleanup(self):
        print("clear up")


if __name__ == "__main__":
    cp = CsstDemoProcessor()

