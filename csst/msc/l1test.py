import glob
import os
import shutil
from .backbone import VER_SIMS


class CsstMscL1Prep:
    """ This is designed for testing L1 pipeline on dandelion """

    def __init__(self, ver_sim="C5.2", dir_target="/nfsdata/users/cham/L1PipelineTest"):
        """

        Parameters
        ----------
        ver_sim:
            version of simulation, {C3, C5.2}
        dir_target:
            the target directory

        """
        self.dir_target = dir_target
        assert ver_sim in VER_SIMS
        if ver_sim == "C3":
            self.dir_source = "/nfsdata/share/csst_simulation_data/Cycle-3-SimuData/multipleBandsImaging/CSST_shearOFF"
        elif ver_sim == "C5.2":
            self.dir_source = "/nfsdata/share/csst_simulation_data/Cycle-5-SimuData/multipleBandsImaging/NGP_AstrometryON_shearOFF"

        print("globbing directory {}".format(self.dir_source))
        self.dps = glob.glob(os.path.join(self.dir_source, "MSC*"))
        self.dps.sort()
        print(" --- available choices --- ")
        for _ in self.dps:
            print(_)

    def run(self, name="MSC_0000150"):
        src = os.path.join(self.dir_source, name)
        dst = os.path.join(self.dir_target, name)
        if os.path.exists(dst):
            print("removing existing {}".format(dst))
            shutil.rmtree(dst)
        print("copying tree ... it takes ~ 1 min for one exposure ...")
        print("{} -> {}".format(src, dst))
        shutil.copytree(src, dst)
        print("DONE!")
        return os.path.join(self.dir_target, name)


def test():
    """ test copy one exposure"""
    l1prep = CsstMscL1Prep(ver_sim="C5.2", dir_target="/nfsdata/users/cham/L1PipelineTest")
    dst = l1prep.run("MSC_0000150")


if __name__ == "__main__":
    test()
