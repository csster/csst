#%% code


fp = "/Users/cham/projects/csst/test/MSC_MS_210527171000_100000279_16_raw.fits"

from csst.common.factory import CsstDataFactory

data = CsstDataFactory.createData(fp)

print(data)

print(data.get_l0keyword("pri", "INSTRUME"))
