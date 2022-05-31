import os
import time

import numpy as np

from csst_photometry import do_phot

ccds = np.arange(6, 26, 1)
indir = '/line17/zouhu/csst/simulation_new'
outdir = '/line17/zouhu/csst/simulation_new/cat_tan/'

for i in range(len(ccds)):
    iccdstr = '%2.2d' % ccds[i]
    start = time.time()
    ifits = os.path.join(indir, 'MSC_MS_210525121500_100000001_' + iccdstr + '_img.fits')
    print(ifits)
    if not os.path.isfile(ifits): continue
    do_phot(ifits, outdir=outdir, stage="phot")
