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
do_one_exposure(runproc=(1, 1, 0, 0), **CONFIG_ANYUSER_DDL_C52)
```
6. check the performance, outputs are here: `~/L1Pipeline/msc/work/`
7. Loop Steps 3-4-5-6 if necessary.