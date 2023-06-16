#!bin/bash

# Script to create and run the PBS jobs for the bias correction step of the
# MRNBC. The input command are created using the script:
#      "create_inputs.R".
# The bias correction is submitted as separate jobs for each 

# Justin Peter
# BMRP - 06 Mar 2020

# Set the working directory to the current directory
JOBDIR=$PWD

# Set the job name
JOBNAME="run_mrnbc_awapres"

# Set the template PBC submission file
TEMPLATE_FILE="run_mrnbc_parallel_at_awap_res.pbs_template"

# Set the output for the PBS jobs and create the directory if it doesn't exist
PBSDIR=${JOBDIR}/pbs_jobs
mkdir -p ${PBSDIR}

#declare -a models=("ACCESS1-0" "CNRM-CM5" "GFDL-ESM2M" "MIROC5")
#declare -a models=("CNRM-CM5" "GFDL-ESM2M" "MIROC5")
declare -a models=("ACCESS1-0")
#declare -a models=("CNRM-CM5")
#declare -a models=("GFDL-ESM2M")
#declare -a models=('MIROC5')
#declare -a models=('GFDL-ESM2M' 'MIROC5')

# Don't submit a historical scenario job as it is produced as a by-product
# when on of rcp45 or rcp85 are run. In the current setup the histrocial
# output is deleted during rcp45 run and only output for rcp85
#declare -a scenarios=("historical" "rcp45" "rcp85")
#declare -a scenarios=("historical")
#declare -a scenarios=("rcp45" "rcp85")
#declare -a scenarios=("rcp45")
declare -a scenarios=("rcp85")

# Change to the job directory as a precaution
cd ${JOBDIR}

# Loop over all the GCMs, scenarios and models
for gcm in ${models[@]}; do
    for scn in ${scenarios[@]}; do
        echo 'Submitting job: ' ${gcm} ${scn}
        # Copy the template file to the job submission script
        cp ${JOBDIR}/${TEMPLATE_FILE} ${JOBDIR}/run_mrnbc_at_awap_res_${gcm}_${scn}.pbs
        # Change the values of the command line input file
        sed -i.bak -e "s/JOBNAME/${JOBNAME}/g" -e "s/GCM/${gcm}/g" -e "s/SCN/${scn}/g" run_mrnbc_at_awap_res_${gcm}_${scn}.pbs
        wait

        # Submit the job
        qsub run_mrnbc_at_awap_res_${gcm}_${scn}.pbs
        #qsub -kod run_mrnbc_at_awap_res_${gcm}_${scn}.pbs

        wait #Make sure it has been submitted before moving the PBS script

                # Move the *.pbs and *.pbs_bak file to the PBSDIR
        mv -f run_mrnbc_at_awap_res_${gcm}_${scn}.pbs ${PBSDIR}
        mv -f run_mrnbc_at_awap_res_${gcm}_${scn}.pbs.bak ${PBSDIR}
    done
done

