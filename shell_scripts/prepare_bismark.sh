#!/bin/bash

conda activate celfie

bismark_file="../data/sample.bismark.cov"
sites="../data/reference_file_regions.txt"
reference="../data/reference_file_tims.txt"


# make bismark files bed file in format as CpG chrom, start, end, methylation, depth
awk '{print $1 "\t" $2 "\t" $3+1 "\t" $5 "\t" $5+$6}' $bismark_file > $bismark_file".bed"

# find TIM sites and sum +/-250 bp around TIM 
# python script to make a file _tims that has the same number of sites as TIM region file. 
# if the site is missing in the input file, indicate this missing data as 0.0, 0.0
bedtools intersect -a $bismark_file".bed" -b $sites > $bismark_file"_500.txt"
python sum_by_list.py $sites $bismark_file"_500.txt" $bismark_file"_tims.txt" 1
rm $bismark_file"_500.txt"

# get rid of mac line endings UGH 
tr -d '\r' < $bismark_file"_tims.txt"  > $bismark_file"_tims_line_ending.txt" 
mv $bismark_file"_tims_line_ending.txt" $bismark_file"_tims.txt"

# need to adapt for multiple input samples
paste $bismark_file"_tims.txt" $reference > "../data/sample_reference_file_tims.txt"