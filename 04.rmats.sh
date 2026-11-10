#!/bin/bash
#SBATCH -J rmats
#SBATCH -p AMD_9654
#SBATCH -N 1
#SBATCH -n 10
#SBATCH -o rmats.out
#SBATCH -e rmats.err

cd workpath
python ~/mambaforge/envs/rmats/bin/rmats.py --b1 b1.txt --b2 b2.txt --gtf /ref/GRCm38/gencode.vM25.annotation.gtf --od AS_out --tmp tmp_AS_out -t paired --readLength 150 --cstat 0.0001 --nthread 10
