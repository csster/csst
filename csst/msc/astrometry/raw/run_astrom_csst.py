#!/usr/bin/env python3.8
#code for astrometry
#written by Jundan NIE
#history
#2021.6 not use missfits

import os
import numpy as np
from astropy.wcs import WCS
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
import healpy as hp
from astropy import table
from astropy.io import fits
from astropy.time import Time
import time
import matplotlib.pyplot as plt
from subprocess import Popen, PIPE
from multiprocessing import Pool
from functools import partial
from astropy import units as u
import re

def prepare_data(fitslist,image_prefix,path_config,output_dir):
    """
    Prepare data for running scamp; Combine all image data, weight files, flag files to their one frame.

    Parameters
    ----------
    fitslist: the image name, e.g.,/line17/Jundan/csst/combine/MSC_210304093000_0000000_06_img.fits 
    image_prefix: the prefix of the combined image fits, e.g.,MSC_210304093000_0000000
    path_config:path of the configure files
    output_dir: directory of output

    Outputs
    ----------
    The joined multi-extension file(not stacked), weight,flag files. e.g., MSC_210304093000_0000000_img.fits,MSC_210304093000_0000000_wht.fits, MSC_210304093000_0000000_flg.fits.
    """
     
    print('################## use missfits to combine images ###################')

    fn=output_dir+image_prefix
    miss_fits1='missfits -c '+path_config+'default.missfits -OUTFILE_TYPE MULTI '+fn+' -SPLIT_SUFFIX _%02d_img.fits -SAVE_TYPE NEW -NEW_SUFFIX _img'
    miss_fits2='missfits -c '+path_config+'default.missfits -OUTFILE_TYPE MULTI '+fn+' -SPLIT_SUFFIX _%02d_flg.fits -SAVE_TYPE NEW -NEW_SUFFIX _flg'
    miss_fits3='missfits -c '+path_config+'default.missfits -OUTFILE_TYPE MULTI '+fn+' -SPLIT_SUFFIX _%02d_wht.fits -SAVE_TYPE NEW -NEW_SUFFIX _wht'
    p=Popen(miss_fits1,shell=True)
    p.wait()
    p=Popen(miss_fits2,shell=True)
    p.wait()
    p=Popen(miss_fits3,shell=True)
    p.wait()

    print('################## missfits done ########################')

def run_sextractor(fitslist,path_config,output_dir):
    """
    Run sextractor
    Parameters 
    ----------
    input file: file name, e.g.,/line17/Jundan/csst/combine/MSC_210304093000_0000000_06_img.fits 
    Outputs
    ----------
    The catalog, with position and flux.
    """

    #run sextractor under specific sex setting.
    path_sex=path_config+'new_csst_realtime.no.weight.sex'
    image_name=fitslist
    fn=fitslist[fitslist.rfind('MSC'):-5] #MSC_210304093000_0000000_06_img
    fn_wt=fitslist[:-8]+'wht.fits'
    fn_flg=fitslist[:-8]+'flg.fits'
    hdul=fits.open(fitslist)
    hdr=hdul[1].header
    gain=hdr['GAIN1']
    sex_comd1='sex -c '+path_sex+' '
    sex_comd2=fitslist+' -CATALOG_NAME '+output_dir+fn+'.acat'
    sex_comd3=' -GAIN '+str(gain) +' -GAIN_KEY '+ 'abcdefg '
    sex_comd4='-PARAMETERS_NAME '+path_config+'csst_realtime.param '+'-FILTER_NAME '+path_config+'csst_realtime.conv '+'-STARNNW_NAME '+path_config+'csst_realtime.nnw'
    sex_comd=sex_comd1+sex_comd2+sex_comd3+sex_comd4
    print(sex_comd)
    p=Popen(sex_comd,shell=True)
    p.wait()
    
def check_data(fitslist,path_config,image_prefix,output_dir):
    """
    Check the catalog of the sextractor. Extract a good one for running scamp.
    Parameters 
    ----------
    input file: file name, e.g., /line17/Jundan/csst/combine/MSC_210304093000_0000000_06_img.fits
    Outputs
    ----------
    The joined catalog, e.g.,MSC_210304093000_0000000.acat.fits
    """

    #Join single sextractor catalogs together. not use missfits 
    fn=output_dir+image_prefix
    output_catnm=str(fn+'.acat.fits')
    hdul=fits.HDUList()
    if (len(fitslist)==18):
    	for i in range(0,len(fitslist)):
            cat_nm=output_dir+fitslist[i][fitslist[i].rfind('MSC'):-5]+'.acat' #MSC*_img.acat
            cat_i=fits.open(cat_nm)
            hdul.append(cat_i[0])
            hdul.append(cat_i[1])
            hdul.append(cat_i[2])
    	hdul.writeto(output_catnm,overwrite=True)
    else:
        print('the length of fitslist is not equal to 18,needs to check')	
    
    
    #check the number of stars of sextractor catalog.
    c = output_catnm.index('.acat.fits')
    filenm=output_catnm[c-c:c+10] #the file MSC_210304093000_0000000.acat.fits. 10 is the length of '.acat.fits'
    hdulist = fits.open(filenm)
    cols = 0
    for n in range(0,len(fitslist)):  #top header +9ccds of csst(image header +data). need to modify for every survey.
         tdata=hdulist[3*n+2].data
         cols+=len(tdata)
    print ('#detections for '+filenm+':',cols)
    if (cols >200):
         print('number of stars >200,good')
    else:
         print('number of stars <200,bad')
         stop
    print('###############check_data done#######################')
    

def convert_hdu_to_ldac(hdu):
    """
    Convert an hdu table to a fits_ldac table (format used by astromatic suite)

    Parameters
    ----------
    hdu: `astropy.io.fits.BinTableHDU` or `astropy.io.fits.TableHDU`
        HDUList to convert to fits_ldac HDUList

    Returns
    -------
    tbl1: `astropy.io.fits.BinTableHDU`
        Header info for fits table (LDAC_IMHEAD)
    tbl2: `astropy.io.fits.BinTableHDU`
        Data table (LDAC_OBJECTS)
    """
    tblhdr = np.array([hdu[1].header.tostring()])
    col1 = fits.Column(name='Field Header Card', array=tblhdr, format='13200A')
    cols = fits.ColDefs([col1])
    tbl1 = fits.BinTableHDU.from_columns(cols)
    tbl1.header['TDIM1'] = '(80, {0})'.format(len(hdu[1].header))
    tbl1.header['EXTNAME'] = 'LDAC_IMHEAD'

    dcol = fits.ColDefs(hdu[1].data)
    tbl2 = fits.BinTableHDU.from_columns(dcol)
    tbl2.header['EXTNAME'] = 'LDAC_OBJECTS'
    return (tbl1, tbl2)

def get_refcat(image,output_dir,search_radius,path_gaia,silent=True):
    """
    Get reference catalog for scamp. The reference cat is GAIA DR2.
    Parameters
    ----------
    image: a image to get its reference catalog, e.g.,MSC_210304093000_0000000_img.fits. Usually the center of the image is the wcs parameters CRVAL1,CRVAL1.
    search_radius: circle radius for searching, units: degree. e.g., 2 degree for a 1x1 deg^2 image. For large ccd size, use larger radius. csst, r=3 deg.
    path: directory of the reference catalog.
    Returns
    -------
    outcat: filename of the cross matched catalog. This catalog is used as a reference catalog for running scamp. e.g.,MSC_000001_r_210103150723gaia.fits
    """    
    
    print('############## getting reference catalog....################')    
    line = output_dir+image
    c = line.index('_img.fits')
    fname=line[0:c]+'_img.fits'
    gaianame=line[0:c]+'.gaia.fits'
    outcat=gaianame
    gaialacnm=line[0:c]+'.gaialac.fits'
    hdu = fits.open(fname)
    header1=hdu[0].header
    header2=hdu[1].header
    ################proper motion ##########
    #oday=header1['DATE-OBS']
    #otime=header1['TIME-OBS']
    #exptime=header1['EXPTIME']
    #odaytime=header1['DATE-OBS']+'T'+header1['TIME-OBS']
    #print('calculating proper motion...',oday, odaytime)
    #t = Time(oday, format='isot', scale='utc')
    #print('t=', t)
    #t.format = 'decimalyear'
    #print(t.format)
    #t1=Time(2015.5, format='decimalyear', scale='utc')
    #print('t1, t.value =', t1, t.value)
    #deltatime= t.value -2015.5
    #print('(Time - 2015.5) = deltatime =', deltatime)
    #print(deltatime, 'deltatime')
    deltatime=0.0
    ra =  float(header2['CRVAL1'])
    dec = float(header2['CRVAL2'])
    c = SkyCoord(ra,dec, unit=(u.deg, u.deg))
    print('ra, dec ra.deg dec.deg= ', ra, dec, c.ra.deg, c.dec.deg)
    ra = c.ra.deg
    dec = c.dec.deg
    c = SkyCoord(ra, dec, unit=(u.deg, u.deg))
    phi = c.ra.deg / (180. / np.pi)
    theta = (90. - c.dec.deg) / (180 / np.pi)
    vec = hp.ang2vec(theta, phi) # Cartesian coordinates of center of circle: xyz
    #find the indices of all the pixels within 2 degrees of that point
    # Array of indices of pixels inside circle. Use nside=32,because zhouzhimin gaia catalog is 12*32*32=12288
    list1 = hp.query_disc(nside=32, vec=vec, radius=np.radians(search_radius))
    if -1 in list1:
        list1.remove(-1)
    ipring = np.array(list1)
    # uniq and sort
    pix = np.unique(ipring) #pix =00000,03061....
    npix = pix.size
    print(ipring, 'ipring', pix, 'pix', npix, 'npix')  #check if there are reduplicates
    dt = np.dtype([('ra', '>f8'), ('dec', '>f8'), ('ra_error', '>f8'), ('dec_error', '>f8'),('phot_g_mean_mag', '>f8'),
                   ('pmra', '>f8'), ('pmra_error', '>f8'), ('pmdec', '>f8'),  ('pmdec_error', '>f8')])
    refcat = table.Table(dtype=dt)
    for i in pix:
        print('i= %5.5d' % i, path_gaia)
        fname = path_gaia + 'chunk-' + '%5.5d' % i + '.fits'
        print('fname=', fname)
        if not silent: print('Reading ', fname)
        d = table.Table.read(fname)
        refcat = [refcat, d]
        refcat = table.vstack(refcat, join_type='inner')
    refcat.rename_column('ra', 'X_WORLD')
    refcat.rename_column('dec', 'Y_WORLD')
    # print(refcat) above are just getting gaia information for ra,dec 2 deg ring.
    print('delta_time between obs_cat and ref_cat:',deltatime)

    mask = (refcat['pmdec']!=refcat['pmdec'])
    refcat['pmdec'][mask]=0
    mask = (refcat['pmra']!=refcat['pmra'])
    refcat['pmra'][mask]=0
    refcat['X_WORLD']=refcat['X_WORLD']+deltatime*refcat['pmra']/3600.0/1000.0/np.cos(refcat['Y_WORLD']*np.pi/180.)
    refcat['Y_WORLD']=refcat['Y_WORLD']+deltatime*refcat['pmdec']/3600.0/1000.0
    refcat['ra_error']=refcat['ra_error']/1000.0/3600.0
    refcat['dec_error']=refcat['dec_error']/1000.0/3600.0
    refcat.rename_column('ra_error', 'ERRA_WORLD')
    refcat.rename_column('dec_error', 'ERRB_WORLD')
    refcat.rename_column('phot_g_mean_mag', 'MAG')
    if outcat: refcat.write(outcat, format='fits',  overwrite=True)

    if os.path.isfile(gaianame):
        print('exist')
        hdu = fits.open(gaianame)
        hdu1= convert_hdu_to_ldac(hdu)
        hdup = fits.PrimaryHDU()
        hdu=hdu1[0]
        tbhdu=hdu1[1]
        thdulist = fits.HDUList([hdup,hdu, tbhdu])
        if os.path.isfile(gaialacnm): os.remove(gaialacnm)
        thdulist.writeto(gaialacnm)

    print('##################### end #####################')
    return (outcat)

def rewrite_wcs_head(head):
    #rewrite the WCS head from Scamp to the standard fits header.
    wcshead=head+'.fits'
    f = open(head, 'r')
    f1 = open(wcshead, 'w') 
    a=''
    i=0
    for u in f.readlines(): 
        v = re.sub("é", "e", u)
        sp=''
        asp=''
        i+=1
        if len(v)<=81:
          sp=' '*(81-len(v))
        if 'END' in v:
            asp=' '*80*(36-i%36)
            i=i+(36-i%36)
            #print(i)
        a=a+v+sp+asp
    f1.write(a.replace('\n','')) 
    f1.close()
    f.close()
    return wcshead


def check_astrometry(output_dir,path_config,image_prefix,fitslist,run_twice=False):
    #update the header of the image fits
    print('############## check the astrometry quality ################')
    r1=[]
    d1=[]
    fn1_1=output_dir+image_prefix
    
    ###################################
    wcshead=rewrite_wcs_head(fn1_1+'.acat.head')#MSC_MS_210525121500_100000001.acat.head->MSC*.acat.head.fits$
    acat=fits.open(fn1_1+'.acat.fits')
    acat_change=str(fn1_1+'.acat.change.fits')
    cat_suffix='.acat'

    if run_twice:
        wcshead=rewrite_wcs_head(fn1_1+'.bcat.head')#MSC_MS_210525121500_100000001.bcat.head->MSC*.bcat.head.fits$
        acat=fits.open(fn1_1+'.bcat.fits')
        acat_change=str(fn1_1+'.bcat.change.fits')
        #cat_suffix='.bcat'
    ####################################

    hdul=fits.HDUList()
    if (len(fitslist)==18):
        for i in range(0,len(fitslist)):
            print(fitslist[i])
            wcshdr=fits.getheader(wcshead,i)# read headers and change to RA---TPV,DEC--TPV for wcs_transfer package
            wcshdr['CTYPE1']='RA---TPV'
            wcshdr['CTYPE2']='DEC--TPV'
            w=WCS(wcshdr)
            #print(wcshdr)
            cat_nm=output_dir+fitslist[i][fitslist[i].rfind('MSC'):-5]+cat_suffix
            cat_i=fits.open(cat_nm)
            sexcat=cat_i[2].data
            ra_sex=sexcat['ALPHA_J2000']
            dec_sex=sexcat['DELTA_J2000']
            x=sexcat['XWIN_IMAGE']
            y=sexcat['YWIN_IMAGE']
            r,d=w.all_pix2world(x,y,0) #convert xwin,ywin to ra,de
            sexcat['ALPHA_J2000']=r
            sexcat['DELTA_J2000']=d
            cat_i[2].data=sexcat
            hdul.append(cat_i[0])
            hdul.append(cat_i[1])
            hdul.append(cat_i[2])
            r1=np.hstack((r1,r))
            d1=np.hstack((d1,d)) 
        obsc=SkyCoord(ra=r1*u.degree,dec=d1*u.degree)
        tmp_cat=np.zeros((len(obsc),2))
        tmp_cat[:,0]=obsc.ra
        tmp_cat[:,1]=obsc.dec
        np.savetxt(output_dir+'scamp_coord.txt',tmp_cat,fmt="%.10f %.10f",delimiter="\n")
        hdul.writeto(acat_change,overwrite=True)#update the cat with new ra,dec (from 1st scamp wcs.)
    else:
        print('the length of fitslist is not equal to 18,needs to check')

    gaia_cat=Table.read(fn1_1+'.gaia.fits') 
    gaia_ra=gaia_cat['X_WORLD']
    gaia_dec=gaia_cat['Y_WORLD']
    refc=SkyCoord(ra=gaia_ra*u.degree,dec=gaia_dec*u.degree)
    idx, d2d, d3d = obsc.match_to_catalog_sky(refc)
    ref_uid=np.unique(idx)
    obs_uid=np.full_like(ref_uid,-1)
    tmpj=-1
    ccdraoff_med=ccddecoff_med=ccdra_rms=ccddec_rms=-1
    for i in ref_uid:
        tmpj=tmpj+1
        iid=(idx == i)
        iiid = (d2d.deg[iid] == d2d.deg[iid].min())
        obs_uid[tmpj]=iid.nonzero()[0][iiid.nonzero()[0]][0]

    uidlim=d2d[obs_uid].arcsecond <1. # set match radius=1 arcsec
    if uidlim.sum()>0:
       obs_uidlim=obs_uid[uidlim]
       ref_uidlim=ref_uid[uidlim]
       ccdraoff=(obsc[obs_uidlim].ra- refc[ref_uidlim].ra).arcsec*np.cos(obsc[obs_uidlim].dec.deg*np.pi/180.)
       ccdraoff_med=np.median(ccdraoff)
       ccdra_rms=np.std(ccdraoff)
       ccddecoff=(obsc[obs_uidlim].dec- refc[ref_uidlim].dec).arcsec
       ccddec_rms=np.std(ccddecoff)
       ccddecoff_med=np.median(ccddecoff)
       match_num=len(ccdraoff)
       print('################# astrometry result: ##############')
       if (match_num<100):
           print('### bad astrometry ###')
    print('median ra_off, dec_off (mas):',ccdraoff_med*1000.,ccddecoff_med*1000.)  
    print('rms ra_off, dec_off (mas):',ccdra_rms*1000.,ccddec_rms*1000.)  
    print('############################################')    
    return ccdraoff,ccddecoff,ccdraoff_med,ccddecoff_med,ccdra_rms,ccddec_rms


def write_headers(fitslist,output_dir,run_twice=False):
    """
    Wrtie history to header
    """
    if run_twice:
        head_suffix=output_dir+fitslist[0][fitslist[0].rfind('MSC'):-12]+'.bcat.head.fits'
    else:
        head_suffix=output_dir+fitslist[0][fitslist[0].rfind('MSC'):-12]+'.acat.head.fits'
    print(head_suffix)
    hdul2=fits.open(head_suffix)
    if (len(fitslist)==18):
        for i in range(0,len(fitslist)):
            fits_nm=output_dir+fitslist[i][fitslist[i].rfind('MSC'):-5]+'.head'
            hdul1=fits.open(fits_nm,mode='update')
            hdr=hdul1[0].header
            hdr2=hdul2[i].header
            hdr.extend(hdr2,unique=True, update=True)
            WCS_S=0
            WCS_V='2.0.4'
            WCS_P='default.scamp'
            WCS_TOL=time.strftime('%Y-%m-%d %H:%M:%S %p')
            hdr.set('WCS_S','0','0=done')
            hdr.set('WCS_V',WCS_V,'Version of WCS calibration')
            hdr.set('WCS_P',WCS_P,'Configure file name of WCS')
            hdr.set('WCS_TOL',WCS_TOL,'Time of last wcs calibration')
            hdul1.flush()
            hdul1.close()
    else:
        print('total numner of the fits files are not 18.')



def do_astrometry(fitslist,search_radius,path_gaia,path_config,output_dir,plot=False,remove=False,run_twice=False):
    fn1_1=fitslist[0][fitslist[0].rfind('MSC'):-12] #MSC_210304093000_0000000
    fn1_2=fitslist[0][fitslist[0].rfind('MSC'):-5] #MSC_210304093000_0000000_06_img
    image_prefix=fn1_1
    #######################preparing data########################
    #print('############## preparing data ##############')
    prepare_data(fitslist,image_prefix,path_config,output_dir)  

    print('############## start running sextractor, in multi-processing ##############')
    p=Pool()
    prod_x=partial(run_sextractor, path_config=path_config,output_dir=output_dir)#multi-processing with multiple arguments.
    result=p.map(prod_x,fitslist)
    p.close()
    p.join()
    print('############## sextractore done ##############')

   # ###############combine and check the sextractor_catalog,extract a good one for running scamp############
    print('############## check sextractor catalog ##############')
    check_data(fitslist,path_config,image_prefix,output_dir)
   
    ###############match with reference catalog#####################
    image=image_prefix+'_img.fits'
    get_refcat(image,output_dir,search_radius,path_gaia,silent=True)
    Popen('cp '+output_dir+image_prefix+'.gaialac.fits '+output_dir+'ref.cat',shell=True)
    
    ###############run scamp################
    print('############## start running scamp for the first time ##############')
    scamp_comd='scamp '+output_dir+image_prefix+'.acat.fits -ASTREFCAT_NAME= '+output_dir+'ref.cat\
    -MERGEDOUTCAT_NAME '+output_dir+'merged.cat -FULLOUTCAT_NAME '+output_dir+'full.cat\
    -c '+path_config+'default.scamp'
    print(scamp_comd)
    p=Popen(scamp_comd,shell=True)
    p.wait()
    print('############## scamp done ##############')

    #######################check the astrometry quality#########################
    ccdoff=check_astrometry(output_dir,path_config,image_prefix,fitslist,run_twice=run_twice)

    #########################write headers#########################
    print('############### updating headers......###############')
    write_headers(fitslist,output_dir,run_twice=run_twice)
    ############################################# 

    ########################## make some analyze plot.#######################
    if plot:
           print('##### Analyzing the scampe result, making some pltos.... ####')
           plt.figure(figsize=(11,5))
           ax1=plt.subplot(121)
           bin=0.05
           plt.grid(color='grey',ls='--')
           plt.plot(ccdoff[0],ccdoff[1],'ko',markersize=3,alpha=0.3)
           plt.xlabel(r'$\Delta$ RA (arcsec)',fontsize=12)
           plt.ylabel(r'$\Delta$ Dec (arcsec)',fontsize=12)
           ax2=plt.subplot(122)
           plt.grid(color='grey',ls='--')
           plt.hist(ccdoff[0],bins=np.arange(-1,1, bin),histtype="step",color="r",label=r'$\Delta$RA (arcsec)')
           plt.hist(ccdoff[1],bins=np.arange(-1,1, bin),histtype="step",color="b",label=r'$\Delta$Dec (arcsec)')
           plt.legend()
           a=str(float(ccdoff[2]))
           b=str(float(ccdoff[4]))
           c=str(float(ccdoff[3]))
           d=str(float(ccdoff[5]))
           plt.text(-0.95,45,r'$\mu$='+a[0:6]+r',  $\sigma$='+b[0:5]+' (arcsec)',color='red')
           plt.text(-0.95,35,r'$\mu$='+c[0:6]+r',  $\sigma$='+d[0:5]+' (arcsec)',color='blue')
           plt.xlabel('coord_diff (arcsec)',fontsize=12)

           plt.savefig(output_dir+image_prefix+'_radec_off.png',dpi=300) 


###################### remove redundant files########################
    if remove:
        print('#################remove tmp files......#################')
        if os.path.isfile(output_dir+image_prefix+'_??_img.acat'): os.remove(output_dir+image_prefix+'_??_img.acat')
        if os.path.isfile(output_dir+image_prefix+'.gaia.fits'): os.remove(output_dir+image_prefix+'.gaia.fits')
        if os.path.isfile(output_dir+image_prefix+'gaialac.fits'): os.remove(output_dir+image_prefix+'gaialac.fits')
        if os.path.isfile(output_dir+'scamp.xml'): os.remove(output_dir+'scamp.xml')
        if os.path.isfile(output_dir+image_prefix+'_img.fits.back'): os.remove(output_dir+image_prefix+'_img.fits.back')
        if os.path.isfile(output_dir+image_prefix+'_wht.fits'): os.remove(output_dir+image_prefix+'_wht.fits')
        if os.path.isfile(output_dir+image_prefix+'_flg.fits'): os.remove(output_dir+image_prefix+'_flg.fits')

    print('###################### all steps are done #################')


################ main program ####################
def main():
    import argparse
    import glob
    import time
    from fnmatch import fnmatch, fnmatchcase
    parser=argparse.ArgumentParser(description='Do astrometry.')
    parser.add_argument('fitsfile', metavar='FITSFILE', type=np.str,help='input file for the CCD frame')	
    parser.add_argument('-r','--search_radius',type=float,help='search radius for getting reference catalog.',default=2.0)
    parser.add_argument('-p','--path_gaia',type=str,help='direcotory of the reference catalog',default='/line12/gaia/chunks-gaia-dr2-astrom/')
    parser.add_argument('-s','--path_config',type=str,help='path of the config file',default='/line17/Jundan/csst/c3_data/pro/')
    parser.add_argument('-o','--path_output',type=str,help='path of the output files',default='/line17/Jundan/csst/c3_data/out_put/')
    parser.add_argument("-f", "--figure",action="store_true",help='plot the astrometry result.')
    parser.add_argument("-d", "--delete",action="store_true",help='delete redundant files.')
    parser.add_argument("-t", "--twice",action="store_true",help='run sextractor and scampe for the second time to get more accurate wcs parameters.')
    args = parser.parse_args()
    fns=args.fitsfile
    fitslist_tmp0=sorted(glob.glob(fns))
    search_radius=args.search_radius
    path_gaia=args.path_gaia
    path_config=args.path_config
    outdir=args.path_output
    path_origin=fitslist_tmp0[0][:fitslist_tmp0[0].rfind('MSC')] #e.g.,/line17/Jundan/csst/combine/
    fitslist_tmp1=fitslist_tmp0[0][fitslist_tmp0[0].rfind('MSC'):-8]#the first file,e.g.,MSC_210304093000_0000000_06_img.fits
    fitslist_tmp2=outdir+fitslist_tmp1+'img.fits'
    fitslist_tmp3=outdir+fitslist_tmp1+'flg.fits'
    fitslist_tmp4=outdir+fitslist_tmp1+'wht.fits'
    fitslist_tmp5=outdir+fitslist_tmp1+'img.head'
    #check if the linked files are exist. If not, make a soft link.
    if os.path.isfile(fitslist_tmp2):
        print('original MSC*_img.fits already exist')
    else:
        print('make a soft link to the original data...')
        p=Popen('ln -s '+path_origin+'MSC*_img.fits '+outdir,shell=True)#make a soft link to the original data.
        p.wait()
    if os.path.isfile(fitslist_tmp3):
        print('original MSC*_flg.fits already exist')
    else:
        print('make a soft link to the original data...')
        p=Popen('ln -s '+path_origin+'MSC*_flg.fits '+outdir,shell=True)#make a soft link to the original da                                p.wait()
    if os.path.isfile(fitslist_tmp4):
        print('original MSC*_wht.fits already exist')                
    else:                   
        print('make a soft link to the original data...')                        
        p=Popen('ln -s '+path_origin+'MSC*_wht.fits '+outdir,shell=True)#make a soft link to the original da                                p.wait()
    if os.path.isfile(fitslist_tmp5):
        print('original MSC*_img.head already exist')                
    else:                    
        print('original MSC*_img.head does not exist, make a soft link to the original data...')                        
        p=Popen('ln -s '+path_origin+'MSC*_img.head '+outdir,shell=True)#make a soft link to the original da                                p.wait()

    fitslist1=[fitslist1 for fitslist1 in fitslist_tmp0 if fnmatchcase(fitslist1,'/*_??_img.fits')]
    fitslist2=fitslist1
    for i in range(0,len(fitslist1)):
        fitslist2[i]=outdir+fitslist1[i][fitslist1[i].rfind('MSC'):]
    fitslist=fitslist2
    print('fitslist for reducing:',fitslist)
    do_astrometry(fitslist,search_radius=search_radius,path_gaia=path_gaia,path_config=path_config,output_dir=outdir,plot=args.figure,remove=args.delete,run_twice=args.twice)

if (__name__=="__main__"):
    main()
