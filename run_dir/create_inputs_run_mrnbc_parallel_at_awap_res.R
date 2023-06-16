# Script to create the input files to the parallel MRNBC jobs at AWAP
# resolution. Based on create_inputs_run_mrnbc_parallel.R

# Justin Peter, BoM Research, 20 May 2020

rm(list=ls())

#models <- c('ACCESS1-0','CNRM-CM5','GFDL-ESM2M','MIROC5')
#models <- c('ACCESS1-0')
#models <- c('CNRM-CM5')
models <- c('GFDL-ESM2M')
#models <- c('MIROC5')
#models <- c('ACCESS1-0','CNRM-CM5')
#models <- c('GFDL-ESM2M','MIROC5')
#models <- c('ACCESS1-0','CNRM-CM5','GFDL-ESM2M')

scenarios <- c("rcp45","rcp85")
#scenarios <- ("rcp45")
#scenarios <- ("rcp85")

#t_slice_start <- c("2006","2036","2066","2070")
#t_slice_start <- c("2006")

#The length of the lon and lat grids
xsize <- c()
ysize <- c()

# AWAP grid size specification
# Unlike the MRNBC performed at GCM res, when it is done at AWAP res the
# grid specified is that of the AWAP grid and not the host GCM
xfirst  = 112.0
xinc    = 0.05
yfirst  = -10.0
yinc    = -0.05
xsize   = 841
ysize   = 681
# Make the grid
# Subtract 1 since at grid cell centres
xlast  = xfirst+(xsize-1.)*xinc 
ylast  = yfirst+(ysize-1.)*yinc
lonvals = seq(xfirst,xlast,by=xinc)
latvals = seq(yfirst,ylast,by=yinc)

#for (tstart in t_slice_start){
for (gcm in models){ 
    for (scn in scenarios){
        fname <- paste("inputs_at_awap_res_",gcm,"_",scn,".txt",sep="")
        # For running in reverse
        #fname <- paste("inputs_at_awap_res_reversed_",gcm,"_",scn,".txt",sep="")
        # fname <- paste("inputs_at_awap_res_time_slices_",gcm,"_",scn,".txt",sep="")
        # For splitting by time chunks as well
        sink(fname)
        #for (tstart in t_slice_start){
            #for (loni in seq(1,xsize)){
            # For reverse
            for (loni in seq(xsize,1)){
                #for (lati in seq(1,ysize)){
                # For reverse
                for (lati in seq(ysize,1)){
                    cat("Rscript --vanilla /g/data1a/er4/jp0715/HydroProj/code/unsw/run_mrnbc/MRNBC_at_awap_res_PARALLEL.R",gcm,scn,sprintf(lonvals[loni],fmt='%#.2f'),loni-1,sprintf(latvals[lati],fmt='%#.2f'),lati-1,'\n',sep=" ")
                    #cat("Rscript --vanilla /g/data1a/er4/jp0715/HydroProj/code/unsw/run_mrnbc/MRNBC_at_awap_res_time_slices_PARALLEL.R",gcm,scn,sprintf(lonvals[loni],fmt='%#.2f'),loni-1,sprintf(latvals[lati],fmt='%#.2f'),lati-1,'\n',sep=" ")
                    #cat("Rscript --vanilla /g/data1a/er4/jp0715/HydroProj/code/unsw/run_mrnbc/MRNBC_at_awap_res_time_slices_PARALLEL.R",gcm,scn,sprintf(lonvals[loni],fmt='%#.2f'),loni-1,sprintf(latvals[lati],fmt='%#.2f'),lati-1,tstart,'\n',sep=" ")
                    #cat("Rscript --vanilla /g/data1a/er4/jp0715/HydroProj/code/unsw/run_mrnbc/MRNBC_at_awap_res_time_slices_PARALLEL.R",gcm,scn,sprintf(lonvals[loni],fmt='%#.2f'),loni-1,sprintf(latvals[lati],fmt='%#.2f'),lati-1,tstart,"$PBS_JOBFS",'\n',sep=" ")
                }
            }
        #Close the file connection
        #sink()
        #}
        #Close the file connection
        sink()
    }
}
