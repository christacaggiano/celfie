#!/bin/bash
#$ -cwd
#$ -j yes
#$ -l h_data=5G
#$ -l h_rt=72:00:00
#$ -l highp
#$ -t 1-10
#### UGE parameters used, parallel job with 10 replicates 

conda activate py37

# run em on pregnancy data 
python EM/em.py data/sample_data.txt EM/sample_output 15 1000 1 1 0.001 1 