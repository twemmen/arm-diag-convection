'''
PURPOSE : To pre-process the CWV (from MWR and MWRRET) and precip for use in convective onset statistics
          There are three functions in this file: (i) calculate_MWR_interp, (ii) calculate_MWRRET_interp 
          and (iii)calculate_precip_timeseries. With the appropriate input paramters, i & ii produce
          5-min averaged, interpolated, equinox-corrected CWV data; iii produces 5-min averaged precipitation data.

AUTHORS : Kathleen Schiro & Fiaz Ahmed

KAS [Matlab versions - 2016]
FIAZ [Pythonization 2/8/2017]
KAS [Updates 3/28/2019]
TE []

'''

## Import require modules
import numpy as np
import datetime as dt
from dateutil.relativedelta import relativedelta
from netCDF4 import Dataset
from glob import glob
from operator import itemgetter
from itertools import groupby

###### CALCULATE MWR INTERP ##########

def calculate_MWR_interp(avg_interval,cwv_varname_mwr,window_length,hours_minus_GMT,season,cwv_data_input_directory):
    
    '''
    This function takes raw MWR column integrated water vapor (CWV) data.
    It calculates an average according to avg_interval and linearly
    interpolates over missing data periods less than those specified by
    window_length (hours)
    '''
    
    # timesteps_in_day
    timesteps_in_day = (1./avg_interval)*24;

    # window_length is the length (in hours) of the window over which one
    # wishes to linearly interpolate the CWV data

    window = timesteps_in_day/(24./window_length);

    #Load each individual file
    f=glob(cwv_data_input_directory+'*.cdf') # Get file list
    length_fnames = len(f)

    # Create date file from filenames
    dts=[]  # Create an empty list to hold dates
    dts_yymmdd=[]
#     f=iter(f) # Make list iterable
    
    
    for i in f:
        if i[-10:-4] == 'custom':
            dts.append(dt.datetime.strptime(i[-26:-11],'%Y%m%d.%H%M%S')) # Extract dates from filename
            dts_yymmdd.append(dt.datetime.strptime(i[-26:-18],'%Y%m%d')) # Extract dates from filename in yymmdd format
        else:
            dts.append(dt.datetime.strptime(i[-19:-4],'%Y%m%d.%H%M%S')) # Extract dates from filename
            dts_yymmdd.append(dt.datetime.strptime(i[-19:-11],'%Y%m%d')) # Extract dates from filename in yymmdd format
    
    
    
    ## Ensure that there no multiple files on the same day
    for i,j in enumerate(dts_yymmdd[:-1]):
        td=dts_yymmdd[i+1]-dts_yymmdd[i]
        if td==0:
            print(dts[i])
            print('There are multiple files for the above date. Need to delete or combine files.')
            return
    
    # Get start and end dates
    strt_date=dts_yymmdd[0]
    temp_date=strt_date
    end_date=dts_yymmdd[-1]
    
    # Find missing days    
    missing_dates=[]  # Storing missing dates in this list
    continuous_dates=[]
    ctr=0
    while temp_date<=end_date:
        continuous_dates.append(temp_date)
        if (temp_date not in dts_yymmdd):
            missing_dates.append(temp_date)
        temp_date+=dt.timedelta(days=1)
        ctr+=1
    
    #Define days to ignore afternoon data due to problems around equinox
    #Ignore data in the 30 days surrounding equinox
       
    yy_strt=dts_yymmdd[0].timetuple().tm_year
    yy_end=dts_yymmdd[-1].timetuple().tm_year
    year_diff=yy_end-yy_strt
    equinox=[]
    
    while yy_strt<=yy_end:

        d01,d02=dt.datetime(yy_strt,3,6),dt.datetime(yy_strt,4,4)  # Spring equinox
        d11,d12=dt.datetime(yy_strt,9,6),dt.datetime(yy_strt,10,5) # Autumn equinox
        
        while d01<=d02:
            equinox.append(d01)
            d01+=dt.timedelta(days=1)
        
        while d11<=d12:
            equinox.append(d11)
            d11+=dt.timedelta(days=1)
        yy_strt+=1
        
    H=np.zeros((int(len(f)),int(timesteps_in_day)))
    H[:]=np.nan # Initialize to nan
    
    #Read in time and cwv 
    for a,i in enumerate(f):
        
        f1=Dataset(i,'r')
        cwv=f1.variables[cwv_varname_mwr][:] # column vapor in cm
        time=f1.variables['time_offset'][:] # seconds since yy-mm-dd 00:00:00 0:00
        TB23=f1.variables['tbsky23'][:] # sky brightness temperature in K
        step=86400./timesteps_in_day
        array_loop=np.arange(step,86400+step,step)
        
        for j,k in enumerate(array_loop):
    
            if k==step:
                ind=np.where(time<k+1)
            
            elif (k>step):
                ind=np.where(np.logical_and(time>k+1-step,time<k+1))
                
            Q=cwv[ind]
            TB23_Q=TB23[ind]

            #Ridding of erroneous data 
            #Rids of -9999 (missing data)
            #Q[Q==-9999]=np.nan
            #Q[Q==0]=np.nan
            #Q[TB23_Q>100]=np.nan
            
            
            #Rids of erroneous data that occur surrounding solstices and
            #equinoxes from 11 am to 2 pm local time

        
            # If the day is within +/- 30 days of equinox        
            if dts_yymmdd[a] in equinox:
            # Get local time 
                dates_day=[dts_yymmdd[a]+dt.timedelta(seconds=z)-dt.timedelta(hours=hours_minus_GMT) for z in time[ind]]
                ind_eq=np.int_([y for y,z in enumerate(dates_day) if (z.timetuple().tm_hour>=11) & (z.timetuple().tm_hour<=14)])
                Q[ind_eq]=np.nan
            
            #H[a,j]=np.nanmean(Q)

            if H[a,j]<=0:
                H[a,j]=np.nan
    
    Hnew=np.zeros((int(len(f)+len(missing_dates)),int(timesteps_in_day))) # Create new array to fill in missing dates
    Hnew[:]=np.nan # Initialize to nan
        
    #Fill in missing days with NaN by defining new array Hnew (Hnew = H + missing days)

    if len(missing_dates)>0:
        for i,j in enumerate(continuous_dates):
            if j in dts_yymmdd:
                ind=np.int_([a for a,b in enumerate(dts_yymmdd) if b == j])[0]
                Hnew[i,:]=H[ind,:]
            else:
                Hnew[i,:]=np.nan

    elif len(missing_dates)==0:
        Hnew=np.copy(H)
    
    #Reshape Hnew into timeseries with length = timesteps_in_day*(length_fnames+length(missing_days)
    cwv_timeseries_raw=Hnew.flatten()
    
    #Do linear interpolation over pre-defined window of NaN (defined in the
    #first part of the code)
    nan_idx=np.isnan(cwv_timeseries_raw)
    nan_ind=np.where(np.isnan(cwv_timeseries_raw))[0] # 
    nan_window=[]
    
    # Get the lengths of NaN blocks. Ignore if length > specified window length
    for l, m in groupby(enumerate(nan_ind), lambda x: x[1]-x[0]):
        temp=(list(map(itemgetter(1), m)))
        
        if len(temp)<=window_length/avg_interval:
            nan_window+=temp
    
    nan_window=(np.asarray(nan_window,dtype=int))
    cwv_timeseries_interp=np.copy(cwv_timeseries_raw)
    indices=np.arange(0,cwv_timeseries_raw.size)
    not_nan=np.logical_not(np.isnan(cwv_timeseries_raw))
    cwv_timeseries_interp[nan_window]=np.interp(indices[nan_window],indices[not_nan],cwv_timeseries_raw[not_nan])    

    
    return cwv_timeseries_interp


###### CALCULATE MWRRET INTERP ##########


def calculate_MWRRET_interp(avg_interval,window_length,season,cwv_data_input_directory):
    
    '''
    This function takes raw MWRRET column integrated water vapor (CWV) data.
    It calculates an average according to avg_interval and linearly
    interpolates over missing data periods less than those specified by
    window_length (hours)
    '''

    # timesteps_in_day
    timesteps_in_day = (1./avg_interval)*24;

    # window_length is the length (in hours) of the window over which one
    # wishes to linearly interpolate the CWV data
    window = timesteps_in_day/(24./window_length);

    #Load each individual file
    f=sorted(glob(cwv_data_input_directory+'*.cdf'))# Get file list
    length_fnames = len(f)

    # Create date file from filenames
    dts=[]          # Create an empty list to hold dates
    dts_yymmdd=[]
    
    # Extract datetimes from the filenames
    for i in f:
        if i[-10:-4] == 'custom':
            dts.append(dt.datetime.strptime(i[-26:-11],'%Y%m%d.%H%M%S')) # Extract dates from filename
            dts_yymmdd.append(dt.datetime.strptime(i[-26:-18],'%Y%m%d')) # Extract dates from filename in yymmdd format
        else:
            dts.append(dt.datetime.strptime(i[-19:-4],'%Y%m%d.%H%M%S')) # Extract dates from filename
            dts_yymmdd.append(dt.datetime.strptime(i[-19:-11],'%Y%m%d')) # Extract dates from filename in yymmdd format
    
    
    # Ensure that there no multiple files on the same day
    for i,j in enumerate(dts_yymmdd[:-1]):
        td=dts_yymmdd[i+1]-dts_yymmdd[i]
        if td==0:
            print(dts[i])
            print('There are multiple files for the above date. Need to delete or combine files.')
            return
    
    # Get start and end dates
    strt_date=dts_yymmdd[0]
    temp_date=strt_date
    end_date=dts_yymmdd[-1]
    
    # Find missing days. This includes days from the beginning of the year (if the series does not start at Jan 1st)
    # and days at the end of the year (if the series does not end on Dec. 31st)
    missing_dates       =[]  # Storing missing dates in this list
    continuous_dates    =[]  # Storing continuous dates in this list (dates available + missing dates)
    ctr                 =0   # counter for debugging reasons
    
    # Fill in days beginning at Jan 1st to the start date
    if (temp_date.month,temp_date.day) != (1,1):                
        pretemp_date = dt.datetime(temp_date.year,1,1,0,0)
        while pretemp_date<temp_date:
            continuous_dates.append(pretemp_date)
            if (pretemp_date not in dts_yymmdd):
                missing_dates.append(pretemp_date)
            pretemp_date+=dt.timedelta(days=1)
            ctr+=1
    
    # Make continuous dates and find missing dates from what is included in the files
    while temp_date<=end_date:                              
        continuous_dates.append(temp_date)
        if (temp_date not in dts_yymmdd):
            missing_dates.append(temp_date)
        temp_date+=dt.timedelta(days=1)
        ctr+=1
    
    # Fill in days from end date to Dec 31st    
    if (end_date.month,end_date.day) != (12,31):
        post_end_date = dt.datetime(end_date.year,12,31,0,0)
        end_temp_date = end_date 
        while end_temp_date<post_end_date:
            continuous_dates.append(end_temp_date)   
            missing_dates.append(end_temp_date)
            end_temp_date+=dt.timedelta(days=1)
            ctr+=1
    
    
    # define some matrix H, which we will fill in with the cwv data 
    H=np.zeros((int(len(f)),int(timesteps_in_day)))
    H[:]=np.nan                                     # Initialize to nan
    
    #Read in time and cwv 
    for a,i in enumerate(f):
        
        f1=Dataset(i,'r')                   # read in file data
        cwv=f1.variables['be_pwv'][:]       # column vapor in cm
        time=f1.variables['time_offset'][:] # seconds since yy-mm-dd 00:00:00 0:00
        
        step=86400./timesteps_in_day                # timestep
        array_loop=np.arange(step,86400+step,step)  # the timestep intervals in the day
        
        # now loop through and average data between points in the days. Ex: if you choose 3 hr avg interval,
        # we look for data points between time bounds evenly spaced by 3 hr, starting at 00 hr, so points
        # between 0 and 3 hrs will be averaged as will 3 to 6 hrs, etc.
        for j,k in enumerate(array_loop):
            
            if k==step:
                ind=np.where(time<k+1)
                
            
            elif (k>step):
                ind=np.where(np.logical_and(time>k+1-step,time<k+1))
           
            ############ QUALITY CONTROL SHOULD BE PLACED HERE ################
            cwv = np.array(cwv)     # cwv as array
            Q1 = cwv[ind]           # indices defined above, points between the timesteps
            Q1 = np.array(Q1)       
            Q1[Q1<=0]=np.nan        # no negative cwv, nan it	            
            H[a,j]=np.nanmean(Q1)   # average the points together and store them in H

    Hnew=np.zeros((int(len(f)+len(missing_dates)),int(timesteps_in_day))) # Create new array to fill in missing dates
    Hnew[:]=np.nan # Initialize to nan
        

    # if we filter by season    
    if season != 'annual':
        season_ind=np.zeros(H.shape[0],dtype='bool') # logical array for season indices
        # define seasons
        if season == 'DJF':
            season_month = [12,1,2]
        elif season == 'MAM':
            season_month = [3,4,5]
        elif season == 'JJA':
            season_month = [6,7,8]
        elif season == 'SON':
            season_month = [9,10,11]

        
        # loop through the dts_yymmdd and if the date is within the selected season, the logical array value for that index becomes true
        for i in np.arange(len(dts_yymmdd)):
            if dts_yymmdd[i].month in season_month:
                season_ind[i] = True
                
        H[season_ind==False,:] =np.nan # change the matrix to nan where the values are not within the season
        # loop through as before    
    
    if len(missing_dates)>0:
        for i,j in enumerate(continuous_dates): # loop through continuous dates

            if j in dts_yymmdd:                 # if continuous date is in a data file:
                ind         =np.int_([a for a,b in enumerate(dts_yymmdd) if b == j])[0] # find index of where j is in dts_yymmdd
                Hnew[i,:]   =H[ind,:]           # fill in the new matrix with the value        
            else:
                Hnew[i,:]   =np.nan             # else the new matrix has a nan value where the continuous date does not exist in the files

    elif len(missing_dates)==0:
        Hnew=np.copy(H)
    
    
    cwv_timeseries_raw=Hnew.flatten()
    
    #Do linear interpolation over pre-defined window of NaN (defined in the
    #first part of the code)
    nan_idx=np.isnan(cwv_timeseries_raw)
    nan_ind=np.where(np.isnan(cwv_timeseries_raw))[0] # 
    nan_window=[]
    
    # Get the lengths of NaN blocks. Ignore if length > specified window length
    for l, m in groupby(enumerate(nan_ind), lambda x: x[1]-x[0]):
        temp=(list(map(itemgetter(1), m)))
        
        if len(temp)<=window_length/avg_interval:
            nan_window+=temp
            
    nan_window=(np.asarray(nan_window,dtype=int))
    cwv_timeseries_interp=np.copy(cwv_timeseries_raw)
    indices=np.arange(0,cwv_timeseries_raw.size)
    not_nan=np.logical_not(np.isnan(cwv_timeseries_raw))
    cwv_timeseries_interp[nan_window]=np.interp(indices[nan_window],indices[not_nan],cwv_timeseries_raw[not_nan])    
    cwv_timeseries_interp = cwv_timeseries_interp*10
    
    
    # create time bnds
    date_list           = [continuous_dates[0] + dt.timedelta(hours=x) for x in np.arange(0, timesteps_in_day*len(continuous_dates)*avg_interval,avg_interval)]
    time_uppr_bnd       = np.roll(date_list,-1)
    time_uppr_bnd[-1]   = time_uppr_bnd[-2] + dt.timedelta(seconds=step) 
    time_bnds           = np.stack((np.asarray(time_uppr_bnd),np.asarray(date_list)))
    
    
    
    return cwv_timeseries_interp, time_bnds


###### CALCULATE PRECIP INTERP ##########

def calculate_precip_timeseries(avg_interval,precip_varname,season,file_time_resolution,precip_data_input_directory):

    '''
    This function takes raw AOSMET precip data from GoAmazon site.
    It currently calculates an average according to avg_interval.
    The extract dates from filename section will need to be ammended based on the raw data files being extracted.
    Need to find a way to automate this more seamlessly among different ARM sites/datasets
    '''
    # timesteps_in_day
    timesteps_in_day = (1./avg_interval)*24;

    #Load each individual file
    f               = sorted(glob(precip_data_input_directory+'*.cdf'))# Get file list
    length_fnames   = len(f)

    # Create date file from filenames
    dts         =[]          # Create an empty list to hold dates
    dts_yymmdd  =[]
    
    # Extract datetimes from the filenames
    for i in f:
        if i[-10:-4] == 'custom':
            dts.append(dt.datetime.strptime(i[-26:-11],'%Y%m%d.%H%M%S')) # Extract dates from filename
            dts_yymmdd.append(dt.datetime.strptime(i[-26:-18],'%Y%m%d')) # Extract dates from filename in yymmdd format
        else:
            dts.append(dt.datetime.strptime(i[-19:-4],'%Y%m%d.%H%M%S')) # Extract dates from filename
            dts_yymmdd.append(dt.datetime.strptime(i[-19:-11],'%Y%m%d')) # Extract dates from filename in yymmdd format
    
    
    # Ensure that there no multiple files on the same day
    for i,j in enumerate(dts_yymmdd[:-1]):
        td=dts_yymmdd[i+1]-dts_yymmdd[i]
        if td==0:
            print(dts[i])
            print('There are multiple files for the above date. Need to delete or combine files.')
            return
    
    # Get start and end dates
    strt_date=dts_yymmdd[0]
    temp_date=strt_date
    end_date=dts_yymmdd[-1]
    
    # Find missing days. This includes days from the beginning of the year (if the series does not start at Jan 1st)
    # and days at the end of the year (if the series does not end on Dec. 31st)
    missing_dates       =[]  # Storing missing dates in this list
    continuous_dates    =[]  # Storing continuous dates in this list (dates available + missing dates)
    ctr                 =0   # counter for debugging reasons
    
    # Fill in days beginning at Jan 1st to the start date
    if (temp_date.month,temp_date.day) != (1,1):                
        pretemp_date = dt.datetime(temp_date.year,1,1,0,0)
        while pretemp_date<temp_date:
            continuous_dates.append(pretemp_date)
            if (pretemp_date not in dts_yymmdd):
                missing_dates.append(pretemp_date)
            pretemp_date+=dt.timedelta(days=1)
            ctr+=1
    
    # Make continuous dates and find missing dates from what is included in the files
    while temp_date<=end_date:                              
        continuous_dates.append(temp_date)
        if (temp_date not in dts_yymmdd):
            missing_dates.append(temp_date)
        temp_date+=dt.timedelta(days=1)
        ctr+=1
    
    # Fill in days from end date to Dec 31st    
    if (end_date.month,end_date.day) != (12,31):
        post_end_date = dt.datetime(end_date.year,12,31,0,0)
        end_temp_date = end_date 
        while end_temp_date<post_end_date:
            continuous_dates.append(end_temp_date)   
            missing_dates.append(end_temp_date)
            end_temp_date+=dt.timedelta(days=1)
            ctr+=1
     
    # file_time_resolution is inluced to indicate whether the .cdf file is for daily data or for an enitre year
    # we define a time factor here to help construct the timeseries
    if file_time_resolution == 'daily':
        time_factor = 1
        daily = True
    else:
        time_factor = 366   # 366 in case of leap yr
        daily = False
    
    # create array, H, that is either of the form (day, timestep) or (year, timestep), in the latter case
    # the array will be much longer       
    
    H=np.zeros((int(len(f)),int(timesteps_in_day)*time_factor))-1
    
    # now loop through and average data between points in the days. Ex: if you choose 3 hr avg interval,
    # we look for data points between time bounds evenly spaced by 3 hr, starting at 00 hr, so points
    # between 0 and 3 hrs will be averaged as will 3 to 6 hrs, etc.
    for a,i in enumerate(f):
        
        f1=Dataset(i,'r')           # read in data
        
        if daily==False:
            tot_seconds = f1.variables['time'][-1] 
            time_factor = np.ceil(tot_seconds/60/60/24)          # this will return eith either 365 or 366
        
        precip=f1.variables[precip_varname][:] # precipitation in mm hr-1
        time=f1.variables['time_offset'][:] # seconds since yy-mm-dd 00:00:00 0:00
        step=86400./timesteps_in_day    
        array_loop=np.arange(step,86400*time_factor+step,step)
        
        # now loop through and average data between points in the days. Ex: if you choose 3 hr avg interval,
        # we look for data points between time bounds evenly spaced by 3 hr, starting at 00 hr, so points
        # between 0 and 3 hrs will be averaged as will 3 to 6 hrs, etc.
        for j,k in enumerate(array_loop): 
            if k==step:
                ind=np.where(time<k+1)
            
            elif (k>step):
                ind=np.where(np.logical_and(time>k+1-step,time<k+1))
            
            
            ############ QUALITY CONTROL SHOULD BE PLACED HERE ################
            precip = np.array(precip)
            Q2 = precip[ind]                # find data points between the timestpes to average
            Q2[np.logical_or(Q2<0,Q2>100)]=np.nan
            H[a,j]=np.nanmean(Q2)           # store average value in the array
    
    
    

     
    #if season == 'annual':
    if (daily==True):    
        Hnew=np.zeros((int(len(continuous_dates)),int(timesteps_in_day)*time_factor)) # Create new array to fill in missing dates
        Hnew[:]=np.nan # Initialize to nan
    
        #Fill in missing days with NaN by defining new array Hnew (Hnew = H + missing days)
        if len(missing_dates)>0:
            for i,j in enumerate(continuous_dates): # loop through continuous dates
                if j in dts_yymmdd:
                    ind=np.int_([a for a,b in enumerate(dts_yymmdd) if b == j])[0] # find index of where j is in dts_yymmdd
                    Hnew[i,:]=H[ind,:]  
                else:
                    Hnew[i,:]=np.nan
    
        elif len(missing_dates)==0:
            Hnew=np.copy(H)
            

        precip_timeseries=Hnew.flatten()
    else:
        precip_timeseries = []              # store precip in here
        for i in np.arange(H.shape[0]):     
            if np.sum(H[i,-int(timesteps_in_day):])<0: # we initialized H 
                # to have values of -1, if those values are still -1, then 
                #the yr is not a leap year, so we store the first 365*steps_in_day points
                precip_timeseries = np.append(precip_timeseries,H[i,:-int(timesteps_in_day)])
                
            else: # else this is a leap year and we store all points 
                precip_timeseries = np.append(precip_timeseries,H[i,:])
    
    

#    else: 
#        if (daily==True): 
#        
#            season_ind=np.zeros(H.shape[0],dtype='bool') # logical array for season indices
#            # define seasons
#            if season == 'DJF':
#                season_month = [12,1,2]
#            elif season == 'MAM':
#                season_month = [3,4,5]
#            elif season == 'JJA':
#                season_month = [6,7,8]
#            elif season == 'SON':
#                season_month = [9,10,11]
#
#        
#            # loop through the dts_yymmdd and if the date is within the selected season, the logical array value for that index becomes true
#            for i in np.arange(len(dts_yymmdd)):
#                if dts_yymmdd[i].month in season_month:
#                    season_ind[i] = True
#                    
#            H[season_ind==False,:] =np.nan # change the matrix to nan where the values are not within the season
#            # loop through as before    
#                
#                
#                
#            Hnew=np.zeros((int(len(f)+len(missing_dates)),int(timesteps_in_day)*time_factor)) # Create new array to fill in missing dates
#            Hnew[:]=np.nan # Initialize to nan
#                
#                  
#            #Fill in missing days with NaN by defining new array Hnew (Hnew = H + missing days)
#            if len(missing_dates)>0:
#                for i,j in enumerate(continuous_dates):
#        
#                    if j in dts_yymmdd:
#                        ind=np.int_([a for a,b in enumerate(dts_yymmdd) if b == j])[0]
#                        Hnew[i,:]=H[ind,:]
#                    else:
#                        Hnew[i,:]=np.nan
#        
#            elif len(missing_dates)==0:
#                Hnew=np.copy(H)
#        
#            #Reshape Hnew into timeseries with length = timesteps_in_day*(length_fnames+length(missing_days)
#            
#
#            precip_timeseries=Hnew.flatten()
#    
#        else:     
#            #seasonal filter in days
#            if season == 'DJF':
#                strt_ind = 0
#                end_ind = 31+31+28-1
#            elif season == 'MAM':
#                strt_ind = 31+31+28-1
#                end_ind = strt_ind + 31+30+31
#            elif season == 'JJA':
#                strt_ind = 31+31+28+31+30+31-1
#                end_ind = strt_ind + 30+31+31
#            elif season == 'SON':
#                strt_ind = 31+31+28+31+30+31+30+31+31-1
#                end_ind = strt_ind + 30+31+30
#                
#            strt_ind    =int(timesteps_in_day*strt_ind)
#            end_ind     =int(timesteps_in_day*end_ind)
#      
#            
#           #seasonal filter, since rows of H are years, we use the indices above to filter each column, then store it 
#            
#            precip_timeseries = []
#            for i in np.arange(H.shape[0]):
#                if np.sum(H[i,-int(timesteps_in_day):])<0:
#                    if season == 'DJF':
#                        x = np.zeros(366*int(timesteps_in_day),dtype='bool')
#                        x[strt_ind:end_ind+1] = True
#                        H[i,x==False] = np.nan
#                        precip_timeseries = np.append(precip_timeseries,H[i,:])
#                    else:
#                        x = np.zeros(366*int(timesteps_in_day),dtype='bool')
#                        x[strt_ind+1:end_ind] = True
#                        H[i,x==False] = np.nan
#                        precip_timeseries = np.append(precip_timeseries,H[i,:])
#                    
#                else:
#                    x = np.zeros(366*int(timesteps_in_day),dtype='bool')
#                    x[strt_ind:end_ind] = True
#                    H[i,x==False] = np.nan
#                    precip_timeseries = np.append(precip_timeseries,H[i,:-int(timesteps_in_day)])
            
            
    # create time bounds       
    date_list           = [continuous_dates[0] + dt.timedelta(hours=x) for x in np.arange(0, timesteps_in_day*len(continuous_dates)*avg_interval,avg_interval)]
    time_uppr_bnd       = np.roll(date_list,-1)
    time_uppr_bnd[-1]   = time_uppr_bnd[-2] + dt.timedelta(seconds=step) 
    time_bnds           = np.stack((np.asarray(time_uppr_bnd),np.asarray(date_list)))
    
        
    if season!='annual':
        season_ind=np.zeros(precip_timeseries.shape[0],dtype='bool') # logical array for season indices
        # define seasons
        if season == 'DJF':
            season_month = [12,1,2]
        elif season == 'MAM':
            season_month = [3,4,5]
        elif season == 'JJA':
            season_month = [6,7,8]
        elif season == 'SON':
            season_month = [9,10,11]

    
        # loop through the dts_yymmdd and if the date is within the selected season, the logical array value for that index becomes true
        for i in np.arange(len(date_list)):
            if date_list[i].month in season_month:
                season_ind[i] = True
        precip_timeseries[season_ind==False] = np.nan
        
    
    
    
    return precip_timeseries,time_bnds
 
    


 



