import os
# os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'  # 如果出现OMP: Error #15 可以使用这行
import glob

from csst.msc.data import CsstMscImgData
from csst.msc.ref_combine import CsstMscRefProc


def test():
    from glob import glob
    root = '/share/simudata/CSSOSDataProductsSims/data/CSSTSimImage_C5/NGP_AstrometryON_shearOFF/'
    bias = glob(root+'*/CSST_MSC_MS_BIAS_*_09_*')
    dark = glob(root+'*/CSST_MSC_MS_DARK_*_09_*')
    flat = glob(root+'*/CSST_MSC_MS_FLAT_*_09_*')
    cmrp = CsstMscRefProc()
    cmrp.prepare(b_p_lst=bias,
                 d_p_lst=dark,
                 f_p_lst=flat,)
    b, d, f = cmrp.run()
    cmrp.cleanup()
    print(b.shape, d.shape, f.shape)

if __name__ == "__main__":
    test()