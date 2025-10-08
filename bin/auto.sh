#!/bin/bash -l
#SBATCH --export=NONE
#SBATCH --partition=mwa
#SBATCH --account=mwasci
#SBATCH --ntasks=16
#SBATCH --mem=64GB
#SBATCH --time=20:00:00
#SBATCH --mail-type FAIL,TIME_LIMIT
#SBATCH --mail-user sirmcmissile47@gmail.com

## load modules
module load singularity/4.1.0-slurm
shopt -s expand_aliases
source /scratch/mwasci/sprabu/moonshine/aliases

set -x
{

obsnum=OBSNUM
base=BASE
model=CAL
scratchID=SCRATCHID

cd ${base}/processing

## create folder
if [ -d "${obsnum}" ]; then
  echo "${obsnum} folder does exist."
else
  mkdir ${obsnum}
fi

cd ${obsnum}

## copy data using scratch id
if [ -d "${obsnum}.ms" ]; then
  echo "measurement set already exists."
else
  cp -r /scratch/mwasci/asvo/${scratchID}/* ${base}/processing/${obsnum}
fi

## step 1) run aoflagger
aoflagger ${obsnum}.ms

if [ -f "round1.bin" ]; then
  echo "calibration file already exits"
else
  ## step 2) calibrate using source model
  calibrate -d ${obsnum}.ms ${obsnum}.metafits -s ../../models/model-${model}-*_withalpha.txt \
    -o round1.bin --uvw-min 100m --uvw-max 2000m \
    --beam-file /scratch/mwasci/sprabu/moonshine/containers/mwa_full_embedded_element_pattern.h5
fi

if [ -f "round1*.png" ]; then
  echo "calibration plots already exits"
else
  ## step 3) plot calibration solution
  plotsolution round1.bin
fi

if [ -d "calibrated.ms" ]; then
  echo "calibrated ms exits"
else
  ## step 4) apply solutions
  applysolution --data ${obsnum}.metafits ${obsnum}.ms \
  -s round1.bin --outputs calibrated.ms
fi

## step 5) imaging
wsclean -name ${obsnum}-img-narrowband -size 2000 2000 -scale 40asec -weight briggs 1\
 -niter 10000 -mgain 0.8 -auto-threshold 1.3 -pol I -apply-primary-beam \
 -mwa-path /scratch/mwasci/sprabu/moonshine/containers -channels-out 24 -maxuvw-m 2500 -circular-beam \
  calibrated.ms/ 

# wsclean -name ${obsnum}-img-wideband -size 1400 1400 -scale 40asec -weight briggs 1\
#  -niter 50000 -mgain 0.8 -auto-threshold 1.3 -pol I -apply-primary-beam
#  -mwa-path /scratch/mwasci/sprabu/moonshine/containers -circular-beam
#   calibrated.ms/ 

}





