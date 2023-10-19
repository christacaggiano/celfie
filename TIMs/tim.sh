#!/bin/bash
################## for SGE jobs #######################
####### comment out if you aren't running on SGE#######
#$ -cwd
#$ -j y
#$ -l h_data=10G
#$ -l h_rt=8:00:00
#######################################################

#### optional conda environment
# conda activate py37

############### define parameters here ################

input_file=sample_input.txt
output_file=sample_tims.txt
summed_file=sample_tims_summed.txt
sum_window_size=500
number_tims=100
number_tissues=19
depth_filter=15
na_filter=2

#######################################################

#### run tim caller
### python tim.py <input file> <output file> <# of tim/tissue> <# of tissues> <depth filter> <max # of missing>
python tim.py $input_file $output_file $number_tims $number_tissues $depth_filter $na_filter


#### sum reference data into window

window_size=$(($sum_window_size/2))  # define window size

sed '1d' $output_file > $output_file"_no_header"  # remove header
cut -f 1-3 $output_file"_no_header" > $output_file"_sites"  # keep only position information
awk -v window="$window_size" '{print $1 "\t" $2-window "\t" $3+window}' $output_file"_sites" > $output_file"_"$window_size

### speed up the computation over large files by using bedtools to select Cpgs already in region
bedtools sort -i $output_file"_"$window_size > $output_file"_sorted"
bedtools intersect -a $input_file -b $output_file"_sorted" > $input_file"_tims"

### sum all CpGs in a region together 
python sum_by_list.py $output_file"_sorted" $input_file"_tims" $summed_file $number_tissues

##### remove intermediate files
rm $input_file"_tims"
rm $output_file"_"$window_size
rm $output_file"_sites"
rm $output_file"_no_header"
rm $output_file"_sorted"
