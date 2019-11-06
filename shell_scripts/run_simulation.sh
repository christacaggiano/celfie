#!/bin/bash
#$ -cwd
#$ -j yes
#$ -l h_data=5G
#$ -l h_rt=72:00:00
#$ -l highp
#$ -t 1-50
#### UGE parameters used #######

# run 50 replicates in a parallel job  

conda activate py36  

# simulates data at fixed percents for 1 person and 1 unknown 
# python vary_percent.py OUTPUT_DIR REPLICATE_NUM TISSUES CPGS PEOPLE CFDNA_DEPTH REF_DEPTH 
python vary_percent.py simple_vary_perc_10x_1per_1unk $SGE_TASK_ID 10 1000 1 10 10 


# simulates mixtures of roadmap/encode WGBS data 
# python roadmap INPUT OUTPUT_DIR ITERATIONS CPGS REPLICATE_NUM CONVERGANCE CRITERIA MIXING_PROP
python roadmap_complex_mix.py ../data/encode_input_and_ref_TIMs.txt chip-tims-complex_same10_nomissing 1000 1000 $SGE_TASK_ID 0.0001 complex_same_10.pkl
