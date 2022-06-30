## Multi-Color Imaging (MCI)

### 1. instrumental effect correction
   1. bias, dark, flat correction
   2. mark bad pixels
   3. hot and warm pixels
   4. over-saturated pixels
   5. cosmic rays
   6. evaluate weight map
### 2. position calibration
   1. *join_data*: combine img/wht/flg data
      - `combined_img/wht/flg.fits`
      - 18*2 HDUs ("", "SCI")
   2. *run_sextractor*: extract sources for each image
      - `csst_realtime.param` PARAMETERS
      - `csst_realtime.conv` FILTER
      - `csst_realtime.nnw` STARNNW
      - `{ccd_id}_img.fits`-->`{ccd_id}_img.acat`
   3. *combine_catalog*: combine sex catalogs into one
      - `{ccd_id}_img.acat` --> `combined_acat.fits`
      - 18*3 HDUs ("", "LDAC_IMHEAD", "LDAC_OBJECTS")
   4. *get_refcat*: get Gaia eDR3 reference catalog
      - `gaia.fits`
      - 2 HDUs ("Primary", "")
      - `gaia_ldac.fits` (LDAC format)
      - 3 HDUs ("Primary", "LDAC_IMHEAD", "LDAC_OBJECTS")
   5. *run_scamp*: run scamp to calibrate position
      - `default2.scamp` CONFIG
      - `combined_acat.fits` -> `combined_acat.head` (auto-generated) WCS head file from scamp
      - `merged.cat` MERGEOUTCAT
      - `full.cat` FULLOUTCAT
      - `scamp_coord.txt`: ra-dec pairs with scamp WCS
   7. *check_astrometry*: collect sex & scamp results
      - *rewrite_wcs_head*
        - `combined_acat.head` -> `combined_acat.head.fits`
      - merge `{ccd_id}_img.acat` into `combined_acat.change.fits` with WCS header `combined_acat.head.fits`
      - `combined_acat.change.fits` has 18*3 HDUs but with updated RA&Dec relative to `combined_acat.fits`
   8. `write_headers` : write (WCS) history to headers
      - merge scamp WCS `combined_acat.head.fits` into `{ccd_id}_img.head`
      - mode=update / overwrite???
      - actually it has no effect...

### 3. flux calibration
For each image, do (`calib` function):
   1. rewrite_sex_catalog

      split `combined_acat.fits` and write down image head
    
   2. filenames
      - `psname=flux_calib.ps`
      - `ckf=checkwcs.ls`
      - `ckim=checkim.ls`
   3. rewrite_sex_cat:
      
      split `combined_acat`
   4. rewrite_wcs_head
      
      convert `combined_acat.head` to `combined_acat.head.fits`
   5. split_wcs_head


### 4. Photometry



## Slit-Less Spectra (SLS)
