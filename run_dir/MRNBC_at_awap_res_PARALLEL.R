# Still need to test whether I need to remove variables
rm(list=ls())

#Read some arguments from the command line.
#args=(commandArgs(trailingOnly=TRUE))
args=commandArgs(trailingOnly=TRUE)


#Assign the input values to the variables
#For operational running
gcm    <- args[1]
scn    <- args[2]
lonval <- args[3]
lonj   <- args[4]
latval <- args[5]
latj   <- args[6]

## For debugging
#gcm    <- "ACCESS1-0"
#scn    <- "rcp45"
###scn    <- "rcp85"
#lonval <- "112.0"
#lonj   <- "0"
###latval <- "-20.0"
###latj   <- "201"
##lonval <- "140.0"
##lonj   <- "560"
#latval <- "-10.0"
#latj   <- "0"

# For debugging - Cairns
## Cairns - check with Fiona and Raj
#gcm <- "ACCESS1-0"
#scn <- "rcp85"
#lonval <- "145.8"
#lonj <- "676"
#latval <- "-16.90"
#latj <- "138"



# Set the lonsliceidx variable based on the lonj index
#if (0 <= lonj   & lonj <= 209){lonsliceidx="1"}
#if (210 <= lonj & lonj <= 419){lonsliceidx="2"}
#if (420 <= lonj & lonj <= 629){lonsliceidx="3"}
#if (630 <= lonj & lonj <= 840){lonsliceidx="4"}
# Set the increment to that the input AWAP and GCM data is split into
# Use below for other models
#lon_inc    <- 0.5 # AWAP is 0.05 so this is 10 longitudes
# Use below for ACCESS
lon_inc    <- 0.05 # AWAP is 0.05 so this is 1 longitudes
first_lon  <- 112.0
last_lon   <- 154.0
delta_awap <- 0.05 # The native AWAP grid increment
# Evaluate the first longitude that is in each netcdf file longitude strip
b_lons     <- seq(from=first_lon,to=(last_lon-lon_inc),by=lon_inc)
b_lons[2:length(b_lons)]<- b_lons[2:length(b_lons)]+delta_awap

# Determine the lonslice number from the element in the b_lons array
lonsliceidx <- max(which(lonval >= b_lons))


##### LIBRARIES #####
suppressMessages(library(ncdf.tools))
suppressMessages(library(ncdf4))
suppressMessages(library(zoo))
suppressMessages(library(filesstrings))

# Information on variables to be bias corrected
#vars <- c("pr","rsds","tasmax","tasmin","sfcWind")
vars <- c("pr","tasmax","tasmin","sfcWind","rsds")
#vars <- c("pr","tasmax","tasmin","sfcWind")
#vars <- c("pr","rsds","tasmax")
#vars <- c("pr","tasmax","tasmin","sfcWind")
nvar <- length(vars)

#variablelowerlimit <- c(0,-20,-20,0,0) # Physical lower limit for variable i.e. values less than this are not sensible
variablelowerlimit <- c(0,-50,-50,0,0) # Physical lower limit for variable i.e. values less than this are not sensible

#variableupperlimit <- c(800,50,50,50) # Physical upper limit for variable i.e. values greater that this are not sensible
#variableupperlimit <- c(800,60,60,50,40) # Physical upper limit for variable i.e. values greater that this are not sensible
variableupperlimit <- c(1000,65,65,50,50) # Physical upper limit for variable i.e. values greater that this are not sensible
#variableupperlimit <- c(1000,70,70,43,50) # Physical upper limit for variable i.e. values greater that this are not sensible

#Need to check whether sfcWind and rsds should be aggregated
variableaggregation <- c(1,0,0,0,0) # Should variable be aggregated or averaged when calculations at longer time scales are performed - 0 for average, 1 for sum (applies to both  monthly and annual level calculations)

variablethresholdTRUE <- c(1,0,0,0,0) # Should a threshold be applied for occurence/non occurence of climate variable

# Note the individual days will be sensitive to these values, however, the
# means of the time series will not.
variablethrehold <- c(1,1,1,1,1) # set value for threshold if used
#variablethrehold <- c(1.0,1,1,1,1) # set value for threshold if used
#variablethrehold <- c(0.5,1,1,1,1) # set value for threshold if used
# Used below to test against Raj's values
#variablethrehold <- c(0.3,1,1,1,1) # set value for threshold if used

#variabledetails <- cbind(variablelowerlimit,variableupperlimit,variableaggregation,variableaggregation,variablethresholdTRUE,variablethrehold)
# Needed to remove duplication of variable aggregation
variabledetails <- cbind(variablelowerlimit,variableupperlimit,variableaggregation,variablethresholdTRUE,variablethrehold)

#This is what should be done for final product
nyh <- 30   # number of years of observed data
nsh <- 1976 # Start year for observed data

# For first iteration of MRNBC - incorrect!
#nyc <- 46   # number of years of gcm historical data
#nsc <- 1960 # start year for GCM historical data
# For second iteration of MRNBC - correct!
#nyc <- 30   # number of years of gcm historical data
#nsc <- 1976 # start year for GCM historical data
# For calculating historical data for 1960-1975
# The MRNBC has to run in time slices equal to the reference period (AWAP
# data)
#nyc <- 30   # number of years of gcm historical data
#nsc <- 1960 # start year for GCM historical data
nyc <- 30   # number of years of gcm historical data
nsc <- 1976 # start year for GCM historical data - same as AWAP data


#nyf <- 95   # number of gcm future data
#nyf <- 94   # number of gcm future data
nyf <- 30   # number of gcm future data
#nsf <- 2006 # start year for GCM future data
#nsf <- 2036 # start year for GCM future data
#nsf <- 2066 # start year for GCM future data
#nsf <- 2070 # start year for GCM future data

# Use below for building complete time series
#nsf <- c(2006,2036,2066,2070)
# Use below to calculate the 1960-1976 results
nsf <- c(2006)
#nsf <- c(2036)
#nsf <- c(2066)
#nsf <- c(2070)
#nsf <- 2010 # start year for GCM future data
leapyear = TRUE

# To reference the appropriate file names
awap_syear <- 1976 #Same as nsh
awap_eyear <- 2005

gcmh_syear <- 1960 #Same as nsc
#gcmh_syear <- 1976 #Same as nsc
gcmh_eyear <- 2005

gcmf_syear <- nsf #Same as nsf
#gcmf_eyear <- 2099
gcmf_eyear <- nsf+nyf-1

##### GENERAL SET UP #####
# File names to write out netcdf data into text file as input to MRNBC
#obsfile <-"awap.dat"
#gcmcurfile <- "gcmc.dat"
#gcmfutfile <- "gcmf.dat"

#obsfile <- paste(output_dir,"/","awap.dat",sep="")
#gcmcurfile <- paste(output_dir,"/","gcmc.dat",sep="")
#gcmfutfile <- paste(output_dir,"/","gcmf.dat",sep="")


#maindir <- "C:/Work/OneDrive - UNSW/Projects/BOM_Downscaling/"
#maindir <- "/g/data1a/er4/jp0715/HydroProj/code/unsw/run_mrnbc/"
maindir <- "/g/data/er4/jp0715/HydroProj/code/unsw/run_mrnbc/"

# Set the input data directories
# data_dir <- "/g/data/eg3/jp0715/HydroProj/data/unsw/"
#data_dir <- "/scratch/eg3/jp0715/HydroProj/data/unsw/"
awap_data_dir <- "/g/data/eg3/jp0715/HydroProj/data/unsw/awap/awap_merged_data/"
gcm_data_dir  <- "/g/data/eg3/jp0715/HydroProj/data/unsw/gcm/gridded_to_awap/"

# Set the output directory
#out_dir <- "/g/data/eg3/jp0715/HydroProj/data/unsw/mrnbc_output/"
#out_dir <- "/g/data/eg3/jp0715/HydroProj/data/unsw/mrnbc_output/test_mpi_run/"

#out_dir <- "/scratch/eg3/jp0715/HydroProj/data/unsw/mrnbc_output/"
#out_dir <- "/scratch/eg3/jp0715/HydroProj/data/unsw/mrnbc_output/awap_res/"
# Usually use one of the below
#out_dir <- "/g/data/eg3/jp0715/HydroProj/data/unsw/mrnbc_output/awap_res/"
out_dir <- "/scratch/eg3/jp0715/HydroProj/data/unsw/mrnbc_output/awap_res/"
#out_dir <- "/scratch/er4/jp0715/HydroProj/data/unsw/mrnbc_output/awap_res/"

#print("Here")

# Loop through time slices
for (itime in 1:length(nsf)){
    #print(itime)
    #Following worked before I added facility for scenarios
    #output_dir <- paste(out_dir,gcm,"/",lonj,"_",latj,sep="")
    #output_dir <- paste(out_dir,gcm,"/",scn,"/",lonj,"_",latj,sep="")
    #output_dir <- paste(out_dir,gcm,"/",scn,"/",gcmf_syear,"-",gcmf_eyear,"/",lonj,"_",latj,sep="")
    output_dir <- paste(out_dir,gcm,"/",scn,"/",gcmf_syear[itime],"-",gcmf_eyear[itime],"/",lonj,"_",latj,sep="")
    #print(output_dir)
    
    # Create a directory to move the files to after completion of MRNBC. This is
    # one directory up from the output directory
    #mv_dir <- paste(out_dir,gcm,"/",scn,"/",sep="")
    mv_dir <- paste(out_dir,gcm,"/",scn,"/",gcmf_syear[itime],"-",gcmf_eyear[itime],"/",sep="")
    
    # Set the names of the output files
    obsfile <- paste(output_dir,"/","awap.dat",sep="")
    gcmcurfile <- paste(output_dir,"/","gcmc.dat",sep="")
    gcmfutfile <- paste(output_dir,"/","gcmf.dat",sep="")
    
    locname <- paste(lonj,"_",latj,sep="")
    bc_curfile <- paste(output_dir,"/","bc_cur","_",locname,".dat",sep="")
    bc_futfile <- paste(output_dir,"/","bc_fut","_",locname,".dat",sep="")

    # Clean up any directories that may have been left half-finished in the
    # previous run first
    if (dir.exists(output_dir)){
       unlink(output_dir,recursive=TRUE)
    }

    
    # Stop execution of the script if the bias correction has already been
    # perfomed for this grid cell
    #if (file.exists(bc_curfile) & file.exists(bc_futfile)){
    if (scn == "rcp45"){
        if (file.exists(paste(mv_dir,basename(bc_futfile),sep=""))){
           #stop("Bias correction already completed: exiting!")
           next
           #stop()
           #quit(save="no")
        }
    }
    
    if (scn == "rcp85"){
        #if (file.exists(paste(mv_dir,basename(bc_curfile),sep="")) & 
        #    file.exists(paste(mv_dir,basename(bc_futfile),sep=""))){
        # USe below if running complete time slices
        #if (file.exists(paste(mv_dir,basename(bc_futfile),sep=""))){
        # Use below if only running 1960-1975 slice
        if (file.exists(paste(mv_dir,basename(bc_curfile),sep=""))){
           #stop("Bias correction already completed: exiting!")
           next
           #quit(save="no")
        }
    }
    
    # Clean up any directories that may have been left half-finished in the
    # previous run first
    #if (dir.exists(output_dir)){
    #   unlink(output_dir,recursive=TRUE)
    #   stop("Bias correction already completed: exiting!")
    #   #quit(save="no")
    #}
    # Create the output directory if it doesn't exist
    if (!dir.exists(output_dir)){
        #system(paste("mkdir -p ",output_dir,sep=""))
        dir.create(output_dir,recursive=TRUE)
    }
    
    #setwd(paste(maindir,"UNSW",sep=""))
    #setwd(paste(maindir,sep=""))
    setwd(output_dir)
    
    #source(paste(maindir,"RFUN_basicdat.R",sep=""))
    source(paste(maindir,"RFUN_basicdat_PARALLEL.R",sep=""))
    #file.copy("./MRNBC_PARALLEL.R",paste0(output_dir,"/","MRNBC_PARALLEL.R"))
    
    
    
    
    #### LOOP THROUGH LOCATIONS ####
    # File names to save outputs from MRNBC
    ##loc <- c("Katherine","Mildura")
    #loc <- c("katherine","mildura")
    ##loc <- c("katherine")
    ##loc <- c("mildura")
    ##for (iloc in 1:2)
    
    #for (iloc in seq_along(loc))
    #{
      #locname <- loc[iloc]
      #bc_curfile <- paste("bc_cur",locname,".dat",sep="")
      #bc_futfile <- paste("bc_fut",locname,".dat",sep="")
      #locname <- paste(lonj,"_",latj,sep="")
      #bc_curfile <- paste(output_dir,"/","bc_cur","_",locname,".dat",sep="")
      #bc_futfile <- paste(output_dir,"/","bc_fut","_",locname,".dat",sep="")
    
    #browser()
     # Stop execution of the script if the bias correction has already been
     # perfomed for this grid cell
     #if (file.exists(bc_curfile) & file.exists(bc_futfile)){
     #if (scn == "rcp45"){
     #    if (file.exists(paste(mv_dir,basename(bc_futfile),sep=""))){
     #       stop("Bias correction already completed: exiting!")
     #       quit(save="no")
     #    }
     #}
    
     #if (scn == "rcp85"){
     #    if (file.exists(paste(mv_dir,basename(bc_curfile),sep="")) & 
     #        file.exists(paste(mv_dir,basename(bc_futfile),sep=""))){
     #       stop("Bias correction already completed: exiting!")
     #       quit(save="no")
     #    }
     #}
    
     
    
      #### WRITE OUT BASIC.DAT ####
      #writebasicdat(obsnyr = nyh,obsstart = nsh,obsdatafile = obsfile,gcmcurnyr = nyc,gcmcurstart = nsc,gcmcurdatafile = gcmcurfile,bcgcmcurfile = bc_curfile,gcmfutnyr = nyf,gcmfutstart = nsf,gcmfutdatafile = gcmfutfile,bcgcmfutfile = bc_futfile,leap = leapyear,nv = nvar,vardet = variabledetails,tscale = "Daily",missval = -9000.0,windwidth = 11)
      writebasicdat(obsnyr = nyh,obsstart = nsh,obsdatafile = obsfile,gcmcurnyr = nyc,gcmcurstart = nsc,gcmcurdatafile = gcmcurfile,bcgcmcurfile = bc_curfile,gcmfutnyr = nyf,gcmfutstart = nsf[itime],gcmfutdatafile = gcmfutfile,bcgcmfutfile = bc_futfile,leap = leapyear,nv = nvar,vardet = variabledetails,tscale = "Daily",missval = -9000.0,windwidth = 11)
      
      #### READ IN DATA AND PROCESS TO MRNBC FORMAT ####
      for (ivar in 1:nvar)
      #for (ivar in 1:(nvar-1)) #Without rsds for now - NEED TO REMOVE FOR PROPER PROCESSING!
      {
        ##### Read in Observed Data ####
        varname <- vars[ivar]
    
        #ncfilename <- paste(awap_data_dir,varname,"/","merged_awap_",awap_syear,"-",awap_eyear,"_", varname,"_","chunked.nc",sep="")
        #fnm <- paste(awap_data_dir,varname,"/","merged_awap_",awap_syear,"-",awap_eyear,"_", varname,"_","chunked.nc",sep="")
        #Used below for ACCESS1-0
        fnm <- paste(awap_data_dir,varname,"/split/","merged_awap_",awap_syear,"-",awap_eyear,"_", varname,"_","chunked_lonslice_",lonsliceidx,".nc",sep="")
        # Use below for other models
        #fnm <- paste(awap_data_dir,varname,"/split/84_slices/","merged_awap_",awap_syear,"-",awap_eyear,"_", varname,"_","chunked_lonslice_",lonsliceidx,".nc",sep="")
    
        # Now perform extraction on the fly
        sys_cmd <- paste("ncks -O -d longitude,",lonval," -d latitude,",latval," ",fnm," ",output_dir,"/","temp_awap_ext.nc",sep="") 
        # Invoke the system command
        system(sys_cmd)
        ncfilename <- paste(output_dir,"/","temp_awap_ext.nc",sep="")
    
        #print(ncfilename)
        clim <- readNcdf(ncfilename)
        clim_coords <- readNcdfCoordinates(ncfilename)
        nlat <- length(clim_coords$lat)
        nlon <- length(clim_coords$lon)
        awapdates <- convertDateNcdf2R(ncfilename)
        #browser()
        awapp <- matrix(NA, nrow = length(clim), ncol = 4) # Save dates and data from netcdf file
        awapp[,1] <- substr(as.character(awapdates),1,4)
        awapp[,2] <- substr(as.character(awapdates),6,7)
        awapp[,3] <- substr(as.character(awapdates),9,10)
        awapp[,4] <- clim
        #browser()
        start <- min(which(awapp[,1]==as.character(nsh)))  # specified start year for observed data
        end <- max(which(awapp[,1]==as.character(nsh+nyh-1)))  # calculate end year based on number of years from start
        temp <- awapp[start:end,]  # Create subset of data to match starting date and number of years specified (ensures that all variables have the same data availability)
        #browser()
        if (ivar == 1) # For the first variable, set up matrix to save all data into
        {
          obs <- matrix(NA, nrow = dim(temp)[1], ncol = 3+nvar)
          obs[,1] <-temp[,1]  # Year
          obs[,2] <-temp[,2]  # Month
          obs[,3] <-temp[,3]  # Day
        }
        nn <- 3+ivar # Climate variable values
        obs[,nn] <- temp[,4]
        #obs[,nn] <- round(as.numeric(temp[,4]),digits=2)
        #browser()
      }
      # MRNBC needs missing values as negative number - replace NA by -999
      nas <- which(obs == 'NA')
      nans <- which(obs == 'NaN')
      obs[nas] <- -999
      obs[nans] <- -999
      obs[is.na(obs)] <- -999
      obs[is.nan(obs)] <- -999
      #browser()
    
      #write output to file in format required for MRNBC
      write.table(obs,file=obsfile,row.names=F,col.names=F,quote=F)
      flush.console()
    
      # Remove the temporary netCDF file created
      #system(paste("rm -f temp_awap_ext.nc",sep=""))
      #gc()
    
      ##### Read in Historical GCM Data ####
      for (ivar in 1:nvar)
      {
        varname = vars[ivar]
    
        #fnm <- paste(gcm_data_dir,gcm,"/",varname,"_day_",gcm,"_historical_r1i1p1_19600101-20051231_chunked.nc",sep="")
        fnm <- paste(gcm_data_dir,gcm,"/split/historical/",varname,"/",varname,"_day_",gcm,"_historical_r1i1p1_19600101-20051231_chunked_lonslice_",lonsliceidx,".nc",sep="")
    
        # Now perform extraction on the fly
        sys_cmd <- paste("ncks -O -d longitude,",lonval," -d latitude,",latval," ",fnm," ",output_dir,"/","temp_gcm_hist_ext.nc",sep="")
        # Invoke the system command
        system(sys_cmd)
        #Sys.sleep(1)
        ncfilename <- paste(output_dir,"/","temp_gcm_hist_ext.nc",sep="")
    
        #print(ncfilename)
        clim <- readNcdf(ncfilename,var.name = varname)
        clim_coords <- readNcdfCoordinates(ncfilename)
        nlat <- length(clim_coords$lat)
        nlon <- length(clim_coords$lon)
        #gcmdates <- convertDateNcdf2R(ncfilename) 
        # convertDateNcdf2R doesn't read dates correctly - not sure why, use ncdf4
        f1<-nc_open(ncfilename)
        time<-ncvar_get(f1,"time")
        tunits<-ncatt_get(f1,"time",attname="units")
        tustr<-strsplit(tunits$value, " ")
        gcmdates<-as.Date(time,origin=unlist(tustr)[3])
    
        if (varname == "pr")
        {
          clim <- clim*86400 # Convert pr to mm/day from kg m-2 s-1
        }
        if (varname =="rsds")
        {
          clim <- clim*86400/1e6 # convert rsds from W/m-2 to MJ/m-2 as AWAP data is solar exposure
          #clim <- clim # convert rsds from W/m-2 to MJ/m-2 as AWAP data is solar exposure
        }
        if (varname == "tasmax" | varname == "tasmin")
        {
          clim <- clim-273.15 # Convert degree K to degree C to match AWAP data format for temperature
        }
        gcmd <- matrix(NA, nrow = length(clim), ncol = 4) # Save dates and data from netcdf file
        gcmd[,1] <- substr(as.character(gcmdates),1,4)
        gcmd[,2] <- substr(as.character(gcmdates),6,7)
        gcmd[,3] <- substr(as.character(gcmdates),9,10)
        gcmd[,4] <- clim
        start <- min(which(gcmd[,1]==as.character(nsc)))  # specified start year for observed data
        end <- max(which(gcmd[,1]==as.character(nsc+nyc-1)))  # calculate end year based on number of years from start
        temp <- gcmd[start:end,]  # Create subset of data to match starting date and number of years specified (ensures that all variables have the same data availability)
        if (ivar == 1) # For the first variable, set up matrix to save all data into
        {
          hist <- matrix(NA, nrow = dim(temp)[1], ncol = 3+nvar)
          hist[,1] <-temp[,1]  # Year
          hist[,2] <-temp[,2]  # Month
          hist[,3] <-temp[,3]  # Day
        }
        nn <- 3+ivar # Climate variable values
        hist[,nn] <- temp[,4]
        #hist[,nn] <- as.numeric(temp[,4])
        #hist[,nn] <- round(as.numeric(temp[,4]),digits=2)
        #browser()
      }
    
      # MRNBC needs missing values as negative number - replace NA by -999
      #hist[is.na(hist)] <- -999
      #hist[is.nan(hist)] <- -999
      nas <- which(hist == 'NA')
      nans <- which(hist == 'NaN')
      hist[nas] <- -999
      hist[nans] <- -999
    
      
      #write output to file in format required for MRNBC
      write.table(hist,file=gcmcurfile,row.names=F,col.names=F,quote=F)
      flush.console()
    
      # Remove the temporary historical GCM data
      #system(paste("rm -f temp_gcm_hist_ext.nc",sep=""))
      #gc()
    
    
      ##### Read in Future GCM Data ####
      for (ivar in 1:nvar)
      {
        varname = vars[ivar]
    
        #fnm <- paste(gcm_data_dir,gcm,"/",varname,"_day_",gcm,"_",scn,"_r1i1p1_20060101-20991231_chunked.nc",sep="")
        fnm <- paste(gcm_data_dir,gcm,"/split/",scn,"/",varname,"/",varname,"_day_",gcm,"_",scn,"_r1i1p1_20060101-20991231_chunked_lonslice_",lonsliceidx,".nc",sep="")
    
        # Now perform extraction on the fly
        sys_cmd <- paste("ncks -O -d longitude,",lonval," -d latitude,",latval," ",fnm," ",output_dir,"/","temp_gcm_fut_ext.nc",sep="")
        # Invoke the system command
        system(sys_cmd)
        ncfilename <- paste(output_dir,"/","temp_gcm_fut_ext.nc",sep="")
    
    
        #print(ncfilename)
        clim <- readNcdf(ncfilename,var.name = varname)
        clim_coords <- readNcdfCoordinates(ncfilename)
        nlat <- length(clim_coords$lat)
        nlon <- length(clim_coords$lon)
        #gcmdates <- convertDateNcdf2R(ncfilename)
        # convertDateNcdf2R doesn't read dates correctly - not sure why, use ncdf4
        f1<-nc_open(ncfilename)
        time<-ncvar_get(f1,"time")
        tunits<-ncatt_get(f1,"time",attname="units")
        tustr<-strsplit(tunits$value, " ")
        gcmdates<-as.Date(time,origin=unlist(tustr)[3])
    
        if (varname == "pr")
        {
          clim <- clim*86400 # Convert pr to mm/day from kg m-2 s-1
        }
        if (varname =="rsds")
        {
          clim <- clim*86400/1e6 # convert rsds from W/m-2 to MJ/m-2 as AWAP data is solar exposure
          #clim <- clim # convert rsds from W/m-2 to MJ/m-2 as AWAP data is solar exposure
        }
        if (varname == "tasmax" | varname == "tasmin")
        {
          clim <- clim-273.15 # Convert degree K to degree C to match AWAP data format for temperature
        }
        gcmd <- matrix(NA, nrow = length(clim), ncol = 4) # Save dates and data from netcdf file
        gcmd[,1] <- substr(as.character(gcmdates),1,4)
        gcmd[,2] <- substr(as.character(gcmdates),6,7)
        gcmd[,3] <- substr(as.character(gcmdates),9,10)
        gcmd[,4] <- clim
        #start <- min(which(gcmd[,1]==as.character(nsf)))  # specified start year for observed data
        #end <- max(which(gcmd[,1]==as.character(nsf+nyf-1)))  # calculate end year based on number of years from start
        start <- min(which(gcmd[,1]==as.character(nsf[itime])))  # specified start year for observed data
        end <- max(which(gcmd[,1]==as.character(nsf[itime]+nyf-1)))  # calculate end year based on number of years from start
        temp <- gcmd[start:end,]  # Create subset of data to match starting date and number of years specified (ensures that all variables have the same data availability)
        if (ivar == 1) # For the first variable, set up matrix to save all data into
        {
          fut <- matrix(NA, nrow = dim(temp)[1], ncol = 3+nvar)
          fut[,1] <-temp[,1]  # Year
          fut[,2] <-temp[,2]  # Month
          fut[,3] <-temp[,3]  # Day
        }
        nn <- 3+ivar # Climate variable values
        fut[,nn] <- temp[,4]
        #fut[,nn] <- as.numeric(temp[,4])
      }
      # MRNBC needs missing values as negative number - replace NA by -999
      #fut[is.na(fut)] <- -999
      #fut[is.nan(fut)] <- -999
      nas <- which(fut == 'NA')
      nans <- which(fut == 'NaN')
      fut[nas] <- -999
      fut[nans] <- -999
    
      
      #write output to file in format required for MRNBC
      write.table(fut,file=gcmfutfile,row.names=F,col.names=F,quote=F)
      flush.console()
    
      # Remove the temporary future GCM data
      #system(paste("rm -f temp_gcm_fut_ext.nc",sep=""))
      #gc()
    
      
      ##### Run MRNBC ####
      # Removed below line to reduce output
      #message("Running bias correction model")
    # system(shQuote("MRNBC.exe",type="cmd"))
    # system("mrnbc")
    
      #To run in parallel we need to copy the executable to the appropriate
      #directory and run from there
     
      #input_exe <- "/g/data1a/er4/jp0715/HydroProj/code/unsw/run_mrnbc/mrnbc"
      input_exe <- "/g/data1a/er4/jp0715/HydroProj/code/unsw/run_mrnbc/mrnbc8"
      #file.copy("mrnbc",paste0(output_dir,"/","mrnbc"))
      #file.copy(input_exe,paste0(output_dir,"/","mrnbc"))
      file.copy(input_exe,paste0(output_dir,"/","mrnbc8"))
      #mrnbc_exe <- paste0(output_dir,"/","mrnbc")
      mrnbc_exe <- paste0(output_dir,"/","mrnbc8")
      #system(paste0(output_dir,"/","mrnbc"))
      system(mrnbc_exe)
    
    #  # Do some cleaning up
    #  # The historical bc'd data is produced for both scenarios so remove it for
    #  # one the scenarios and move the bc file up one level
    #browser()
       if (scn == "rcp45"){
           #system(paste("rm -f ",bc_curfile,sep=""))
           #file.remove(bc_curfile)
           file.move(c(bc_futfile),c(mv_dir))
       }
    #
    #  # If it is the rcp85 scenario, then move the bc'd historical and future
    #  # files
       if (scn == "rcp85"){
          if (nsf[itime] == "2006"){
             #file.move(c(bc_curfile,bc_futfile),mv_dir)
             file.move(c(bc_curfile),mv_dir)
             #Remove line below for 1960 runs
             #file.move(c(bc_futfile),mv_dir)

             #file.copy(c(bc_curfile),mv_dir)
             #file.copy(c(bc_futfile),mv_dir)
          } else {
             file.move(c(bc_futfile),mv_dir)
          }
          #file.move(c(bc_curfile,bc_futfile),mv_dir)
       }
    #      
    #  # Move the output files up on directory
    #  #system(paste("mv -f ","bc_* ", "../",sep=""))
    #  #file.move(c(bc_curfile,bc_futfile),mv_dir)
    #
    #  # Remove the basic.dat file after completion of MRNBC
    #  #system(paste("rm -f basic.dat",sep=""))
    #  #system(paste("rm -f awap.dat",sep=""))
    #  #system(paste("rm -f gcmc.dat",sep=""))
    #  #system(paste("rm -f gcmf.dat",sep=""))
    #  #system(paste("rm -f mrnbc",sep=""))
    #
    #  #system(paste("rm -rf ",output_dir,sep=""))
       # Set the wd to mv_dir worked when not looping through time slices
       #setwd(mv_dir)
       #unlink(output_dir,recursive=TRUE)
       unlink(output_dir,recursive=TRUE)
       setwd("/g/data/er4/jp0715/HydroProj/code/unsw/run_mrnbc/")
    #
    
    #}
} #End loop over time slices
