## how to test msc pipeline

1. **log in**:
   1. `ssh csstpipeline@10.3.10.30`  passwd:pipeline@123
2. `cd ~/L1Pipeline/csst/`
3. edit your code
4. install your code: `sh install.sh`
5. test your code:
   1. open a new python concole
   2. run the code below
```python
from csst.msc.pipeline import do_one_exposure, CONFIG_ANYUSER_DDL_C52
print(CONFIG_ANYUSER_DDL_C52)
do_one_exposure(runproc=(1, 1, 0, 0), **CONFIG_ANYUSER_DDL_C52)
```
Go to Step 3-4-5 again if necessary.