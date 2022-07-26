## how to test msc pipeline

last modified: 2022-07-26

aims: to test position calibration module

Steps:
1. log in:
   - `ssh csstpipeline@10.3.10.30`  passwd:pipeline@123
2. `cd ~/L1Pipeline/csst/`
3. edit your code
4. install your code: `sh install.sh`
5. test your code:
   - open a new python concole and run the code below
```python
from csst.msc.pipeline import do_one_exposure, CONFIG_ANYUSER_DDL_C52
print(CONFIG_ANYUSER_DDL_C52)
do_one_exposure(runproc=(0, 1, 0, 0), **CONFIG_ANYUSER_DDL_C52)
```
6. check the performance, outputs are here: `~/L1Pipeline/msc/work/`
7. Loop Steps 3-4-5-6 if necessary.

## tips on updates of `DataManager`
The `csst.msc.data_manager.CsstMscDataManager` has several updates.
Currently, 

```python
from csst.msc.data_manager import CsstMscDataManager
dm = CsstMscDataManager(...)

# L0 file paths:

dm.l0_cat(ccd_id=6)  # the input catalog of CCD 6
dm.l0_crs(ccd_id=6)  # the cosmic ray image of CCD 6
dm.l0_log(ccd_id=6)  # the log file of CCD 6
dm.l0_ccd(ccd_id=6)  # the raw image of CCD 6

# L1 synthetic file paths:

# the L1 image for CCD 6
# e.g., CSST_MSC_MS_SCI_20270810081950_20270810082220_100000100_06_img_L1.fits
dm.l1_ccd(ccd_id=6, post="img_L1.fits")
# in "this function", you want a combined image file
dm.l1_file(name="combined_image.fits", comment="this function")

```