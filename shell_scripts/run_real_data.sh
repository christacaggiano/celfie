#!/bin/bash
#$ -cwd
#$ -j yes
#$ -l h_data=5G
#$ -l h_rt=72:00:00
#$ -l highp
#$ -t 1-10
#### UGE parameters used, parallel job with 10 restarts 

conda activate py37

# run em on pregnancy data 
# requires mixing function to be specificed in EM/mix_function.py 
python em.py ../data/pregnancy_sample_data.txt pregnancy_sample 1000 1000 preg $SGE_TASK_ID 0.0001
