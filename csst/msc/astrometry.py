from ..core.processor import CsstProcessor


class CsstProcMscPositionCalibration(CsstProcessor):
    def prepare(self, **args):
        # prepare the environment
        # for example, if you make use of some third-party software like SEXTRACTOR,
        # do your preparation here.
        pass

    def run(self, data, *args, **kwargs):
        # run your pipeline here
        # make sure that your input data should be a child class instance of CsstData.
        pass

    def cleanup(self, **kwargs):
        # clean up environment
        pass
