import csv
import numpy as np
import pandas as pd
import sys
import collections


def get_region_dict(file, tissues): 

	regions_dict = collections.defaultdict(list)

	with open(file, "r") as input_file: 
		list_file = csv.reader(input_file, delimiter="\t")

		for line in list_file: 
			chrom, start, end = line[0], line[1], line[2]
			regions_dict[chrom].append({"start":int(start), "end":int(end), "meth":np.zeros(tissues), "depth":np.zeros(tissues)})

	return regions_dict


def get_methylation_counts(file, regions_dict): 
	with open(file, "r") as input_file: 
		cpg_file = csv.reader(input_file, delimiter="\t")

		for line in cpg_file: 
			chrom, start, end = line[0], line[1], line[2]
			meth = np.array(line[3::2], dtype=np.float64)
			depth = np.array(line[4::2], dtype=np.float64) 
			
			if chrom in regions_dict: 
				for region in regions_dict[chrom]: 
					if int(start) >= region["start"] and int(end) <= region["end"]: 
						region["meth"] += meth
						region["depth"] += depth
						
	return regions_dict		


def write_bed_file(output_file, regions_dict): 
	with open(output_file, "w") as output: 
		bed_file = csv.writer(output, delimiter="\t")

		for chrom in regions_dict: 
			for region in regions_dict[chrom]:
				values = [] 
				for m, d in zip(region["meth"], region["depth"]): 
					values.append(m)
					values.append(d)
				bed_file.writerow([chrom] + [region["start"]] + [region["end"]] + values)
		
if __name__ == "__main__": 

	list_file = sys.argv[1] 
	tissue_cpg_file = sys.argv[2]
	output_file_name = sys.argv[3]
	tissues = int(sys.argv[4])
	
	regions = get_region_dict(list_file, tissues)
	get_methylation_counts(tissue_cpg_file, regions)
	write_bed_file(output_file_name, regions)

