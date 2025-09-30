#!/bin/bash -l
#SBATCH --export=NONE
#SBATCH --partition=mwa
#SBATCH --account=mwasci
#SBATCH --ntasks=28
#SBATCH --mem=124GB

## load modules
module load hyperdrive/0.6.1
module load singularity/4.1.0-slurm

## step 1
cp -r * /scratch/mwasci/sprabu/moonshine/processing/1384513256

## step 2) ## aoflagger
singularity exec /scratch/mwasci/sprabu/moonshine/containers/wsclean_2.10.0-build-1.sif aoflagger 1384513256.ms

## calibrate ## 
singularity exec /scratch/mwasci/sprabu/moonshine/containers/hyperdrive_main.sif \
hyperdrive di-calibrate -d 1384513256.ms 1384513256.metafits \
-s ../../models/model-3C444-10comp_withalpha.txt -o round1.bin \
--uvw-min 100m --uvw-max 2000m --beam-file /scratch/mwasci/sprabu/moonshine/containers/mwa_full_embedded_element_pattern.h5

##  make calibration solutions
singularity exec /scratch/mwasci/sprabu/moonshine/containers/hyperdrive_main.sif\
 hyperdrive solutions-plot round1.bin

## apply solutions
singularity exec /scratch/mwasci/sprabu/moonshine/containers/hyperdrive_main.sif 
    hyperdrive solutions-apply --data round1.bin 1384513256.ms --outputs calibrated.ms

## imaging
singularity exec /scratch/mwasci/sprabu/moonshine/containers/wsclean_2.10.0-build-1.sif \
wsclean -name test -size 1400 1400 -scale 20asec -weight natural\
 -niter 1000 -mgain 0.8 calibrated.ms/



