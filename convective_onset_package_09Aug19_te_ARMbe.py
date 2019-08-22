##### Convective Onset Diagnostics Package, for DOE ARM Diagnostics Package
##### [in progress - last update 21 Mar 2019 by K. Schiro]
##### Original scripts by Kathleen Schiro (MATLAB)
##### Translations (version 3.5.2) and additions made by Fiaz Ahmed and Baird Langenbrunner

import scipy
import scipy.io
import numpy as np
from glob import glob
import os
import re
from netCDF4 import Dataset
from getSites import getSites
from trmm_process import trmm_cwv_match
##### User-Specified Input [paths, default values, package options, ARM site] #####

## SPECIFY PATHS 
working_directory = os.getcwd()+'/'


## DEFAULT VALUES
number_of_bins = 26 #(default = 28)
cwv_max = 67 #(default = 70)
cwv_min = 28 #(default = 28)
bin_width = 1.5 #(default = 1.5)
precip_threshold = 0.5 #(default = 0.5)
window_length = 6 #(for MWR and precip interpolation, in hrs; default = 6)
avg_interval = 1 # (in hours) match to model, if applicable (default = 1/12 = 5 mins)

## PACKAGE OPTIONS
## ["BASIC" is Fig. 1d-f of Schiro et al. (2016) with MWR CWV and rain gauge precip at 15 min resolution]
## ["BASIC_model" is same as "BASIC" but with model output overlain]
## ["BASIC_time" is same as "BASIC+model" with MWR CWV and rain gauge precip averaged at specified interval]
## ["BASIC_space" is same as "BASIC+model" with TRMM 3B42 precip spatially averaged over ARM site]   
## ["ADV_sonde" is Fig. 1a-c of Schiro et al. (2016) with radiosonde CWV and rain gauge precip 1-hr avg surrounding radiosonde]
## ["ADV_temp" adds temperature dependence to Fig. 1d of Schiro et al. (2016); MWR CWV, rain gauge precip, interpolated radiosonde T (?)]
## ["ADV_time" compares temporal averages of in situ gauge precip binned by MWR]
  
BASIC_model = 1 
BASIC_space = 0
ADV_sonde   = 0
ADV_temp    = 0
ADV_time    = 0


use_trmm = False
## CHOOSE ARM SITE

arm_site    = 5
#Amazon     = 1
#Nauru      = 2
#SGP        = 3
#Manus_Isl  = 4
#Darwin     = 5
#Niamey     = 6
#Barrow     = 7


model       = 2
#GISS       = 1
#GFDL       = 2 
#NCAR       = 3  is currently broken
#E3SM       = 4
#IPSL       = 5
#CanAm4     = 6
#HadGem2    = 7


#Specify Paths from model/ARM site selection
# BE is for use of ARM best estimate data
be =    True
season = 'annual'







if arm_site==1:
    cwv_data_input_directory    = working_directory + 'goamazon_data/raw_data/MWRRET/'
    precip_data_input_directory = working_directory + 'goamazon_data/raw_data/AOSMET/'
    hours_minus_GMT             = 4
    location                    = 'GOAmazon'        #location is specified here for title of plots
    cwv_varname_mwr             = 'vap'             #in some instances
    cwv_varname_mwrret          = 'be_pwv'
    precip_varname              = 'rain_intensity'
    file_time_resolution        = 'daily'
    number_of_bins              = 26 #(default = 28)
    cwv_max                     = 67 #(default = 70)
    cwv_min                     = 28 #(default = 28)

elif arm_site==2:
    if be:
        cwv_data_input_directory    = working_directory + 'ARM_Nauru/nauru_be_prw/'
        precip_data_input_directory = working_directory + 'ARM_Nauru/nauru_met/'
    
    else:
        cwv_data_input_directory    = working_directory + 'ARM_Nauru/nauru_mwrret/'
        precip_data_input_directory = working_directory + 'ARM_Nauru/nauru_met/'
        
    
    hours_minus_GMT             = -12
    cwv_varname_mwr             = 'vap'
    cwv_varname_mwrret          = 'be_pwv'
    precip_varname              = 'org_precip_rate_mean'
    location                    = 'Nauru'
    file_time_resolution        = 'daily'
    number_of_bins              = 30 #(default = 28)
    cwv_max                     = 90 #(default = 70)
    cwv_min                     = 28 #(default = 28)
    
elif arm_site==3:
    if be:
        cwv_data_input_directory    = working_directory + 'ARM_SGP/sgp_be_prw/'
        precip_data_input_directory = working_directory + 'ARM_SGP/sgp_be_precip/'
    
    else: 
        cwv_data_input_directory    = working_directory + 'ARM_SGP/sgp_mwrret/'
        precip_data_input_directory = working_directory + 'ARM_SGP/sgp_be_precip/'
    
    hours_minus_GMT             = -12
    cwv_varname_mwr             = 'vap'
    cwv_varname_mwrret          = 'be_pwv'
    precip_varname              = 'prec_sfc'
    location                    = 'SGP'
    file_time_resolution        = 'annual'
    number_of_bins              = 26 #(default = 28)
    cwv_max                     = 67 #(default = 70)
    cwv_min                     = 20 #(default = 28)
        

elif arm_site==4:
    
    if be:    
        cwv_data_input_directory    = working_directory + 'ARM_Manus_Island/manus_be_prw/'
        precip_data_input_directory = working_directory + 'ARM_Manus_Island/manus_be_precip_1998_2010/'
    
    else:
        cwv_data_input_directory    = working_directory + 'ARM_Manus_Island/manus_mwrret_c1_cwv_1998_2010/'
        precip_data_input_directory = working_directory + 'ARM_Manus_Island/manus_be_precip_1998_2010/'
    
    hours_minus_GMT             = -12
    cwv_varname_mwr             = 'vap'
    cwv_varname_mwrret          = 'be_pwv'
    precip_varname              = 'prec_sfc'
    location                    = 'Manus Island'
    file_time_resolution        = 'annual'
    number_of_bins              = 28 #(default = 28)
    cwv_max                     = 85 #(default = 70)
    cwv_min                     = 28 #(default = 28)
        
elif arm_site==5:
    if be: 
        cwv_data_input_directory    = working_directory + 'ARM_darwin/darwin_be_prw/'
        precip_data_input_directory = working_directory + 'ARM_darwin/darwin_pr/'     
    
    else:
        cwv_data_input_directory    = working_directory + 'ARM_darwin/darwin_cwv/'
        precip_data_input_directory = working_directory + 'ARM_darwin/darwin_pr/'
    
    hours_minus_GMT             = -12
    cwv_varname_mwr             = 'vap'
    cwv_varname_mwrret          = 'be_pwv'
    precip_varname              = 'prec_sfc'
    location                    = 'Darwin'
    file_time_resolution        = 'annual'
    number_of_bins              = 30 #(default = 28)
    cwv_max                     = 85 #(default = 70)
    cwv_min                     = 28 #(default = 28)
        
elif arm_site==6:
    cwv_data_input_directory    = working_directory + 'ARM_niamey/niamey_mwrret/'
    precip_data_input_directory = working_directory + 'ARM_niamey/niamey_precip/'
    hours_minus_GMT             = -12
    cwv_varname_mwr             = 'vap'
    cwv_varname_mwrret          = 'be_pwv'
    precip_varname              = 'precip_rate_mean'
    location                    = 'Niamey'
    file_time_resolution        = 'daily'
    number_of_bins              = 13 #(default = 28)
    cwv_max                     = 67 #(default = 70)
    cwv_min                     = 28 #(default = 28)
    
    

if model ==1:
    model_input_directory       = working_directory + 'GISS/'
    model_name                  = 'GISS'
    avg_interval                = 6
    
elif model == 2:
    model_input_directory       = working_directory + 'example_model_output_GFDL/'
    model_name                  = 'GFDL'
    avg_interval                = 1
        
elif model==3:
    model_input_directory       = working_directory + 'example_model_output_NCAR/'
    model_name                  = 'NCAR'
    avg_interval                = 6
    
elif model ==4:
    model_input_directory       = working_directory + 'e3sm_timeseries/'
    model_name                  = 'E3SM'
    avg_interval                = 3
    
elif model==5: 
    model_input_directory       = working_directory +'IPSL_model_output/'
    model_name                  = 'IPSL'
    avg_interval                = 0.5
    
elif model==6:
    model_input_directory       = working_directory+'CanAm4_model/'
    model_name                  = 'CanAm4'
    avg_interval                = 0.5
    
elif model==7:
    model_input_directory       = working_directory+'HadGem2_model_output/'
    model_name                  = 'Had-Gem2'
    avg_interval                = 1.0
    
    
    
    
data_output_directory = working_directory + '/processed_data/'
plot_output_directory = working_directory + '/plots/'




avg_interval = 1
##### Run BASIC Package #####

## begin BASIC

#if Amazon == 1:
	
#if Nauru == 1:
#	hours_minus_GMT = -12
	#cwv_varname_mwr = 'vap'
	#cwv_varname_mwrret = 'be_pwv'
	#precip_varname = 'rain_intensity'
#if Manus_Isl == 1:
#	hours_minus_GMT = -10
	#cwv_varname_mwr = 'vap'
	#cwv_varname_mwrret = 'be_pwv'
	#precip_varname = 'rain_intensity'

## Load data and compute mean MWR timeseries 
## [window_length, hours_minus_GMT, paths: specified above] 	
# Read filename prefixes to determine which script to use for pre-processing

os.chdir(cwv_data_input_directory)
fnames = glob(cwv_data_input_directory + '*.cdf')
filename_prefix = fnames[0].split('/')[-1][3:9]
os.chdir(working_directory)
# Load data and compute mean MWR timeseries 
# [window_length, hours_minus_GMT, paths: specified above] 
if np.logical_and(np.logical_and(np.logical_or(np.logical_or(arm_site==3,arm_site==4),arm_site==5),be),avg_interval==1):
       cwv_timeseries_obs = []
       fprw = sorted(glob(cwv_data_input_directory+'/*.cdf'))
       for a,i in enumerate(fprw):
           cwv_timeseries_obs = np.append(cwv_timeseries_obs,Dataset(i,'r').variables['pwv'][:]*10)
           cwv_timeseries_obs[cwv_timeseries_obs<0] = np.nan
else:
    if filename_prefix=='mwrlos':
    	print('Calculating for '+filename_prefix)
    	from pre_processing_scripts import calculate_MWR_interp
    	cwv_timeseries_obs = calculate_MWR_interp(avg_interval,cwv_varname_mwr,window_length,hours_minus_GMT,season,cwv_data_input_directory)
    if filename_prefix=='mwrret':
    	print('Calculating for '+filename_prefix)
    	from pre_processing_scripts import calculate_MWRRET_interp
    	[cwv_timeseries_obs,cwv_time_bnds] = calculate_MWRRET_interp(avg_interval,\
    	window_length,season,cwv_data_input_directory)
        
    
#if arm_site == 2:
    #    cwv_timeseries_obs = 

## Load and compute mean precip timeseries
## [avg_interval, paths: specified above]
if use_trmm:
    cfmip_precip    = np.load('cfmip_trmm.npy')
    trmm_time       = np.load('trmm_time.npy')
    local           = True
    if local:
        suffix = ' local '
        if arm_site == 1:
            precip_cfmip = cfmip_precip[118][5,5,:]
    
        elif arm_site ==2:
            precip_cfmip = cfmip_precip[30][5,5,:]
    
        elif arm_site ==3:
            precip_cfmip = cfmip_precip[36][5,5,:]
    
        elif arm_site ==4:
            precip_cfmip = cfmip_precip[32][5,5,:]
    
        elif arm_site ==5:
            precip_cfmip = cfmip_precip[35][5,5,:]
    
        elif arm_site ==6:
            precip_cfmip = cfmip_precip[109][5,5,:]
    
    else:
        suffix = ' 1 deg '
        if arm_site == 1:
            precip_cfmip = np.nanmean(cfmip_precip[30],axis=(0,1))
    
        elif arm_site ==2:
            precip_cfmip = np.nanmean(cfmip_precip[30],axis=(0,1))
    
        elif arm_site ==3:
            precip_cfmip = np.nanmean(cfmip_precip[36],axis=(0,1))
    
        elif arm_site ==4:
            precip_cfmip = np.nanmean(cfmip_precip[32],axis=(0,1))
    
        elif arm_site ==5:
            precip_cfmip = np.nanmean(cfmip_precip[35],axis=(0,1))
    
        elif arm_site ==6:
            precip_cfmip = np.nanmean(cfmip_precip[109],axis=(0,1))

    precip_timeseries_obs = trmm_cwv_match(cwv_time_bnds,precip_cfmip,trmm_time)
    filename_prefix = 'TRMM' + suffix
        
            
            
            
            
            
            
    
else:
    os.chdir(working_directory)
    if np.logical_and(np.logical_and(np.logical_or(np.logical_or(arm_site==3,arm_site==4),arm_site==5),be),avg_interval==1):
       precip_timeseries_obs = []
       fpr = sorted(glob(precip_data_input_directory+'/*.cdf'))
       for a,i in enumerate(fpr):
           precip_timeseries_obs = np.append(precip_timeseries_obs,Dataset(i,'r').variables['prec_sfc'][:]) 
           precip_timeseries_obs[precip_timeseries_obs<0]= np.nan
        
    else:
        from pre_processing_scripts import calculate_precip_timeseries
        [precip_timeseries_obs,precip_time_bnds] = calculate_precip_timeseries(avg_interval,precip_varname,season,file_time_resolution,precip_data_input_directory)
    ## Compute Convective Onset Statistics 

os.chdir(working_directory)
from convective_onset_data import convective_onset_data
[cwv_bin_centers,precip_mean,errorbars_precip,probability_precip,\
freq_cwv,freq_precip_points,pdf_cwv,pdf_precip_points,\
errorbar_precip_points,errorbar_precip_binom] = convective_onset_data(number_of_bins,cwv_max,cwv_min,\
precip_threshold,cwv_timeseries_obs,precip_timeseries_obs)

## Plot Convective Onset Statistics
os.chdir(working_directory)
from convective_onset_plot import convective_onset_plot
convective_onset_plot(number_of_bins,cwv_bin_centers,precip_mean,errorbars_precip,probability_precip,\
freq_cwv,freq_precip_points,errorbar_precip_points,errorbar_precip_binom,plot_output_directory,location,filename_prefix,season,avg_interval)

    

## end BASIC
print('End Basic Package')

##### Run BASIC_model #####
if BASIC_model == 1:
    os.chdir(working_directory)
    from convective_onset_plot_model import convective_onset_plot_model
    os.chdir(model_input_directory)
    if arm_site == 1:
        print('Start Basic_model')
    	## Pre-process model output here
        ## [For now, read in a few sample files]
        
    	
        if model==1:
            cwv_data = Dataset('cwv_GISS_E3dev.nc')
            cwv_timeseries_model = cwv_data.variables['cwv_GISS_E3dev'][:]
            rain_data = Dataset('pr6_GISS_E3dev.nc')
            precip_model = rain_data.variables['pr6_GISS_E3dev'][:]
            #convert precip from kg m-2 s-1 to mm hr-1
            precip_timeseries_model = precip_model; 
	
        elif model==2:
            cwv_timeseries_model = scipy.io.loadmat('cwv_Manaus.mat')['cwv_Manaus']
            precip_timeseries_model = scipy.io.loadmat('rain_Manaus.mat')['rain_Manaus']
            
        elif model==3:
            cwv_timeseries_model = scipy.io.loadmat('cwv_Manaus.mat')['cwv_Manaus']
            precip_timeseries_model = scipy.io.loadmat('rain_Manaus.mat')['rain_Manaus']
            
      
                            
                

    elif arm_site == 2:
        
        if model==2:
            cwv_timeseries_model = getSites('/home/toddemmen/Documents/shared_ARM_diag/cmip6-cfsites-locations.txt','atmos.2009010100-2010123123.PRW.nc','PRW',0)[30,:]
            precip_timeseries_model = getSites('/home/toddemmen/Documents/shared_ARM_diag/cmip6-cfsites-locations.txt','atmos.2009010100-2010123123.pr.nc','pr',0) [30,:]*3600
            # if you are looking at how this index (30) is different from the rest for this site, this is because the model above uses a code to extract the timeseries based on the nearest 
            # point in the mode to the site. The output from the 'getSites' is an array based off the indices of the sites for the upcoming CMIP6 run. The index below (48) is
            # is the index of the site from the CMIP5 run, which is not consistent with the CMIP6 indexing of sites
        elif model==3:
            cwv_timeseries_model = getSites('/home/toddemmen/Documents/shared_ARM_diag/cmip6-cfsites-locations.txt','NCAR-CAM5.atmos.1990010100-1991123123.prw.nc','prw',0)[30,:]
            precip_timeseries_model = getSites('/home/toddemmen/Documents/shared_ARM_diag/cmip6-cfsites-locations.txt','NCAR-CAM5.atmos.1990010100-1991123123.pr.nc','pr',0)[30,:]*3.6e9
        
        elif model==4:
            moisture_data = Dataset('twpc2_e3sm_v1_0001-01-01_0007-01-01.nc','r')
            precip_timeseries_model = moisture_data.variables['PRECT'][:,0,1]*3.6e6
            cwv_timeseries_model = moisture_data.variables['TMQ'][:,0,1]
            
        elif model==5:
            from glob import glob
            precip_timeseries_model = []
            fpr = glob(model_input_directory+'/ipsl_pr/*.nc')
            for a,i in enumerate(fpr):
                precip_timeseries_model = np.append(precip_timeseries_model,Dataset(i,'r').variables['pr'][:,48]*3600)
                
            cwv_timeseries_model = []
            fcwv = glob(model_input_directory+'/ipsl_prw/*.nc')
            for a,i in enumerate(fcwv):
                cwv_timeseries_model = np.append(cwv_timeseries_model,Dataset(i,'r').variables['prw'][:,48])
                
        elif model==6:
            cwv_timeseries_model = Dataset('prw_cfSites_CanAM4_amip_r1i1p1_19790101003000-20100101000000.nc','r').variables['prw'][:,48]
            precip_timeseries_model = Dataset('pr_cfSites_CanAM4_amip_r1i1p1_19790101003000-20100101000000.nc','r').variables['pr'][:,48]*3600    
        
        elif model==7:
            from glob import glob
            precip_timeseries_model = []
            fpr = glob(model_input_directory+'/hadgem2_pr/*.nc')
            for a,i in enumerate(fpr):
                precip_timeseries_model = np.append(precip_timeseries_model,Dataset(i,'r').variables['pr'][:,30]*3600)
                
            cwv_timeseries_model = []
            fcwv = glob(model_input_directory+'/hadgem2_cwv/*.nc')
            for a,i in enumerate(fcwv):
                cwv_timeseries_model = np.append(cwv_timeseries_model,Dataset(i,'r').variables['prw'][:,30])
                
                
            
    elif arm_site==3:
        
        if model==2:
            cwv_timeseries_model = getSites('/home/toddemmen/Documents/shared_ARM_diag/cmip6-cfsites-locations.txt','atmos.2009010100-2010123123.PRW.nc','PRW',0)[36,:]
            precip_timeseries_model = getSites('/home/toddemmen/Documents/shared_ARM_diag/cmip6-cfsites-locations.txt','atmos.2009010100-2010123123.pr.nc','pr',0)[36,:]*3600
        
        elif model==3:
            cwv_timeseries_model = getSites('/home/toddemmen/Documents/shared_ARM_diag/cmip6-cfsites-locations.txt','NCAR-CAM5.atmos.1990010100-1991123123.prw.nc','prw',0)[36,:]
            precip_timeseries_model = getSites('/home/toddemmen/Documents/shared_ARM_diag/cmip6-cfsites-locations.txt','NCAR-CAM5.atmos.1990010100-1991123123.pr.nc','pr',0) [36,:]*3.6e9
        
        elif model==4:
            moisture_data = Dataset('sgp_e3sm_v1_0001-01-01_0007-01-01.nc','r')
            precip_timeseries_model = moisture_data.variables['PRECT'][:,0,1]*3.6e6
            cwv_timeseries_model = moisture_data.variables['TMQ'][:,0,1]
            
        elif model==5:
            from glob import glob
            precip_timeseries_model = []
            fpr = glob(model_input_directory+'/ipsl_pr/*.nc')
            for a,i in enumerate(fpr):
                precip_timeseries_model = np.append(precip_timeseries_model,Dataset(i,'r').variables['pr'][:,64]*3600)
                
            cwv_timeseries_model = []
            fcwv = glob(model_input_directory+'/ipsl_prw/*.nc')
            for a,i in enumerate(fcwv):
                cwv_timeseries_model = np.append(cwv_timeseries_model,Dataset(i,'r').variables['prw'][:,64])
                
                
        elif model==6:
            cwv_timeseries_model = Dataset('prw_cfSites_CanAM4_amip_r1i1p1_19790101003000-20100101000000.nc','r').variables['prw'][:,64]
            precip_timeseries_model = Dataset('pr_cfSites_CanAM4_amip_r1i1p1_19790101003000-20100101000000.nc','r').variables['pr'][:,64]*3600  
            
            
        elif model==7:
            from glob import glob
            precip_timeseries_model = []
            fpr = glob(model_input_directory+'/hadgem2_pr/*.nc')
            for a,i in enumerate(fpr):
                precip_timeseries_model = np.append(precip_timeseries_model,Dataset(i,'r').variables['pr'][:,64]*3600)
                
            cwv_timeseries_model = []
            fcwv = glob(model_input_directory+'/hadgem2_cwv/*.nc')
            for a,i in enumerate(fcwv):
                cwv_timeseries_model = np.append(cwv_timeseries_model,Dataset(i,'r').variables['prw'][:,64])
                
                
            
                

    #Manus Island        
    elif arm_site==4:
        
        if model==2:
            cwv_timeseries_model = getSites('/home/toddemmen/Documents/shared_ARM_diag/cmip6-cfsites-locations.txt','atmos.2009010100-2010123123.PRW.nc','PRW',0)[32,:]
            precip_timeseries_model = getSites('/home/toddemmen/Documents/shared_ARM_diag/cmip6-cfsites-locations.txt','atmos.2009010100-2010123123.pr.nc','pr',0) [32,:]*3600
        
        elif model==3:
            cwv_timeseries_model = getSites('/home/toddemmen/Documents/shared_ARM_diag/cmip6-cfsites-locations.txt','NCAR-CAM5.atmos.1990010100-1991123123.prw.nc','prw',0)[32,:]
            precip_timeseries_model = getSites('/home/toddemmen/Documents/shared_ARM_diag/cmip6-cfsites-locations.txt','NCAR-CAM5.atmos.1990010100-1991123123.pr.nc','pr',0) [32,:]*3.6e9
        
        
        elif model==4:
            moisture_data = Dataset('twpc1_e3sm_v1_0001-01-01_0007-01-01.nc','r')
            precip_timeseries_model = moisture_data.variables['PRECT'][:,0,1]*3.6e6
            cwv_timeseries_model = moisture_data.variables['TMQ'][:,0,1]
            
        elif model==5:
            from glob import glob
            precip_timeseries_model = []
            fpr = glob(model_input_directory+'/ipsl_pr/*.nc')
            for a,i in enumerate(fpr):
                precip_timeseries_model = np.append(precip_timeseries_model,Dataset(i,'r').variables['pr'][:,57]*3600)
                
            cwv_timeseries_model = []
            fcwv = glob(model_input_directory+'/ipsl_prw/*.nc')
            for a,i in enumerate(fcwv):
                cwv_timeseries_model = np.append(cwv_timeseries_model,Dataset(i,'r').variables['prw'][:,57])

        elif model==6:
            cwv_timeseries_model = Dataset('prw_cfSites_CanAM4_amip_r1i1p1_19790101003000-20100101000000.nc','r').variables['prw'][:,57]
            precip_timeseries_model = Dataset('pr_cfSites_CanAM4_amip_r1i1p1_19790101003000-20100101000000.nc','r').variables['pr'][:,57]*3600
            
        elif model==7:
            from glob import glob
            precip_timeseries_model = []
            fpr = glob(model_input_directory+'/hadgem2_pr/*.nc')
            for a,i in enumerate(fpr):
                precip_timeseries_model = np.append(precip_timeseries_model,Dataset(i,'r').variables['pr'][:,32]*3600)
                
            cwv_timeseries_model = []
            fcwv = glob(model_input_directory+'/hadgem2_cwv/*.nc')
            for a,i in enumerate(fcwv):
                cwv_timeseries_model = np.append(cwv_timeseries_model,Dataset(i,'r').variables['prw'][:,32])
                
                
            
    elif arm_site==5:
        
        if model==2:
            cwv_timeseries_model = getSites('/Users/temmen/Documents/shared_ARM_diag/cmip6-cfsites-locations.txt','atmos.2009010100-2010123123.PRW.nc','PRW',0)[35,:]
            precip_timeseries_model = getSites('/Users/temmen/Documents/shared_ARM_diag/cmip6-cfsites-locations.txt','atmos.2009010100-2010123123.pr.nc','pr',0) [35,:]*3600
        
        elif model==3:
            cwv_timeseries_model = getSites('/Users/temmen/Documents/shared_ARM_diag/cmip6-cfsites-locations.txt','NCAR-CAM5.atmos.1990010100-1991123123.prw.nc','prw',100)[35,:]
            precip_timeseries_model = getSites('/Users/temmen/Documents/shared_ARM_diag/cmip6-cfsites-locations.txt','NCAR-CAM5.atmos.1990010100-1991123123.pr.nc','pr',100) [35,:]*3.6e9
        
        
        
        elif model==4:
            moisture_data = Dataset('twpc3_e3sm_v1_0001-01-01_0007-01-01.nc','r')
            precip_timeseries_model = moisture_data.variables['PRECT'][:,0,1]*3.6e6
            cwv_timeseries_model = moisture_data.variables['TMQ'][:,0,1]

        elif model==5:
            from glob import glob
            precip_timeseries_model = []
            fpr = glob(model_input_directory+'/ipsl_pr/*.nc')
            for a,i in enumerate(fpr):
                precip_timeseries_model = np.append(precip_timeseries_model,Dataset(i,'r').variables['pr'][:,63]*3600)
                
            cwv_timeseries_model = []
            fcwv = glob(model_input_directory+'/ipsl_prw/*.nc')
            for a,i in enumerate(fcwv):
                cwv_timeseries_model = np.append(cwv_timeseries_model,Dataset(i,'r').variables['prw'][:,63])

        elif model==6:
            cwv_timeseries_model = Dataset('prw_cfSites_CanAM4_amip_r1i1p1_19790101003000-20100101000000.nc','r').variables['prw'][:,63]
            precip_timeseries_model = Dataset('pr_cfSites_CanAM4_amip_r1i1p1_19790101003000-20100101000000.nc','r').variables['pr'][:,63]*3600
               
        elif model==7:
            from glob import glob
            precip_timeseries_model = []
            fpr = glob(model_input_directory+'/hadgem2_pr/*.nc')
            for a,i in enumerate(fpr):
                precip_timeseries_model = np.append(precip_timeseries_model,Dataset(i,'r').variables['pr'][:,35]*3600)
                
            cwv_timeseries_model = []
            fcwv = glob(model_input_directory+'/hadgem2_cwv/*.nc')
            for a,i in enumerate(fcwv):
                cwv_timeseries_model = np.append(cwv_timeseries_model,Dataset(i,'r').variables['prw'][:,35])
                
                

               
    elif arm_site==6:
        
        if model==2:
            cwv_timeseries_model = getSites('/Users/temmen/Documents/shared_ARM_diag/cmip6-cfsites-locations.txt','atmos.2009010100-2010123123.PRW.nc','PRW',0)[109,:]
            precip_timeseries_model = getSites('/Users/temmen/Documents/shared_ARM_diag/cmip6-cfsites-locations.txt','atmos.2009010100-2010123123.pr.nc','pr',0) [109,:]*3600
        
        elif model==3:
            cwv_timeseries_model = getSites('/Users/temmen/Documents/shared_ARM_diag/cmip6-cfsites-locations.txt','NCAR-CAM5.atmos.1990010100-1991123123.prw.nc','prw',0)[109,:]
            precip_timeseries_model = getSites('/Users/temmen/Documents/shared_ARM_diag/cmip6-cfsites-locations.txt','NCAR-CAM5.atmos.1990010100-1991123123.pr.nc','pr',0) [109,:]*3.6e9
        
        
        
        elif model==4:
            moisture_data = Dataset('twpc3_e3sm_v1_0001-01-01_0007-01-01.nc','r')
            precip_timeseries_model = moisture_data.variables['PRECT'][:,0,1]*3.6e6
            cwv_timeseries_model = moisture_data.variables['TMQ'][:,0,1]

        elif model==5:
            from glob import glob
            precip_timeseries_model = []
            fpr = glob(model_input_directory+'/ipsl_pr/*.nc')
            for a,i in enumerate(fpr):
                precip_timeseries_model = np.append(precip_timeseries_model,Dataset(i,'r').variables['pr'][:,56]*3600)
                
            cwv_timeseries_model = []
            fcwv = glob(model_input_directory+'/ipsl_prw/*.nc')
            for a,i in enumerate(fcwv):
                cwv_timeseries_model = np.append(cwv_timeseries_model,Dataset(i,'r').variables['prw'][:,56])

        elif model==6:
            cwv_timeseries_model = Dataset('prw_cfSites_CanAM4_amip_r1i1p1_19790101003000-20100101000000.nc','r').variables['prw'][:,56]
            precip_timeseries_model = Dataset('pr_cfSites_CanAM4_amip_r1i1p1_19790101003000-20100101000000.nc','r').variables['pr'][:,56]*3600
               
        elif model==7:
            from glob import glob
            precip_timeseries_model = []
            fpr = glob(model_input_directory+'/hadgem2_pr/*.nc')
            for a,i in enumerate(fpr):
                precip_timeseries_model = np.append(precip_timeseries_model,Dataset(i,'r').variables['pr'][:,106]*3600)
                
            cwv_timeseries_model = []
            fcwv = glob(model_input_directory+'/hadgem2_cwv/*.nc')
            for a,i in enumerate(fcwv):
                cwv_timeseries_model = np.append(cwv_timeseries_model,Dataset(i,'r').variables['prw'][:,106])
                
                

               
    
    
    [cwv_bin_centers,precip_mean_model,errorbars_precip_model,probability_precip_model,\
	freq_cwv_model,freq_precip_points_model,pdf_cwv_model,pdf_precip_points_model,\
	errorbar_precip_points_model,errorbar_precip_binom_model] = convective_onset_data(number_of_bins,cwv_max,cwv_min,\
	precip_threshold,cwv_timeseries_model,precip_timeseries_model)

    convective_onset_plot_model(number_of_bins,cwv_bin_centers,precip_mean,precip_mean_model,errorbars_precip,\
	errorbars_precip_model,probability_precip,probability_precip_model,pdf_cwv,\
	pdf_cwv_model,pdf_precip_points,pdf_precip_points_model,errorbar_precip_points,\
	errorbar_precip_points_model,errorbar_precip_binom_model,plot_output_directory,model_name,location)

## end BASIC_model

##### Run BASIC_space #####

#if BASIC_space == 1:

##### Run ADV_time #####
## end computation of convective onset statistics
