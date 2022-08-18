#!/bin/bash
#
#SBATCH --job-name=gk_detr
#SBATCH --time=01-10:00:00
#SBATCH --output=job_output_%j.out
#SBATCH --error=job_output_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=grahamking@cusacklab.org


python K_detrend.py

