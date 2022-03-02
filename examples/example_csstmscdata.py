fp = "/Users/cham/projects/csst/examples/MSC_MS_210527171000_100000279_16_raw.fits"
fp = "/Users/cham/projects/csst/examples/bias/MSC_CLB_210525120000_100000000_06_raw.fits"
fp = "/Users/cham/projects/csst/examples/lam/MSC_CLB_210525120000_100000000_06_raw.fits"

import os
os.chdir("/Users/cham/projects/csst/test/")
fp = "MSC_MS_210527171000_100000279_16_raw.fits"

from csst.msc import CsstMscImgData
data = CsstMscImgData.read(fp)

print("data: ", data)
print("instrument: ", data.get_l0keyword("pri", "INSTRUME"))
print("object: ", data.get_l0keyword("pri", "OBJECT"))
