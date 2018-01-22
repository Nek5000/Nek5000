#!/bin/bash
#SBATCH -A SNIC2017-1-126
#SBATCH -n 1
#SBATCH --gres=gpu:k80:1
##SBATCH --ntasks-per-node=28
#SBATCH -t 1:30:00

module PGI/16.9-GCC-5.4.0-2.26


./nek5000 > eddy_dir443.log.gpu 2>&1 
