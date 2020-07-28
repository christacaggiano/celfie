# CELFIE (CELl Free dna decomposItion Expectation maximization)

## Overview
Expectation maximization algorithm to decompose complex mixtures of cfDNA into the tissues generating the cfDNA fragments. Required input is methylated cfDNA, either WGBS or Illumina 450K with downsampled read counts, and a reference panel of tissues to estimate their contribution to the cfDNA. CelFiE can estimate an arbitrary number of unknown or missing tissues from your reference that are truly in the cfDNA mixtures.


## Preparing Data

CelFiE expects the methylation data for a cfDNA individual or reference cell type is in the form of # of methylated reads, # of total reads. For example for one sample, the file could like like:
```
# CHR START END METH  DEPTH   
chr1	186917271	186917772	446.0	630.0
chr15	92708070	92708571	71.0	133.0
chr14	55296905	55297406	89.0	115.0
```

If using the provided TIMs (in `data/reference_file_regions.txt`) the chrom, start, end would be a 500bp region around the TIM. To obtain the summed reads per sample, the file `sum_by_list.py` is provided. First, all CpGs within the specified TIM region can be found using bedtools. Then, all reads can be summed using `sum_by_list.py`. If a region has no coverage, meaning that region is 'missing', `sum_by_list.py` will return 0 methylated, 0 depth, for that region.
```
bedtools intersect -a <sample file> -b <data/reference_file_regions.txt" > <output> "
python sum_by_list.py <data/reference_file_regions.txt" >  <output>  <output summed> 1  # run for one sample
```

All cfDNA and reference data should then be compiled into one input tab separated file for CelFiE, with one set of bed columns before the sample data, and one set before the reference data. For example, for 2 input individuals and 4 reference cell types would look like

```
chr1	186917271	186917772	446.0	630.0	230.0	304.0	chr1	186917271	186917772	1156.062	1196.0	3.968	224.0	1172.14	1234.0	782.018	852.0
chr1	23291168	23291669	71.0	133.0	68.0	94.0 chr1	23291168	23291669	341.02	352.0	2.994	87.0	332.996	359.0	295.99	314.0
chr1	23291168	23291669	89.0	115.0	74.0	83.0  chr1	23291168	23291669	168.987	173.0	1.0	39.0	182.98	196.0	94.00	119.0
```

To prepare CelFiE files from Bismark output, see `prepare_bismark.sh`

## Code

### Installation 

To install CelFiE, clone or fork this repository using `git clone https://github.com/christacaggiano/celfie.git`. All required packages can be installed using [Anaconda](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) using the environment file specified in `celfie_conda_env.yml`. Run `conda env create -f celfie_conda_env.yml -n celfie_env`

Full implementation of the EM model at `EM/em.py`. Code to generate simulations can be found in `EM/simulations`.

### EM Script 

After preparing data as above, the EM script as follows:

``` python EM/em.py <input_file> <output_directory> <num of cfDNA samples> <max EM iterations> <num of unknown categories> <parallel job ID> <convergence> <num of random restarts per replicate> ```

Since there is some stochasticity in the EM intialization > 10 random restarts is desirable for a final data run. 

### Sample Code 

``` python EM/em.py data/sample_data.txt EM/sample_output 15 1000 1 1 0.001 10 ```

Currently, the estimated methylation proportions for the reference and the estimated cell type proportions are output in pickled python numpy arrays. These can be read back into python for further analysis by the following

``` python 
import pickle as pkl
cell_proportions = pkl.load(open("EM/sample_output/1_alpha.pkl", "rb"))
methylation_proportions = pkl.load(open("EM/sample_output/1_gamma.pkl", "rb"))
```

### Parallelization 

To run many parallel replicates on a SGE or UGE cluster configuration, see run_real_data.sh
`qsub run_real_data.sh`

### Figures 

Jupyter notebooks to reproduce figures and statistical analyses for the final version of this manuscript can be found in `paper_figures` directory. 

## Contact 
For any questions with this code, please contact christa@g.ucla.edu 

## Citation
For more details on CelFiE, see:

Christa Caggiano, Barbara Celona, Fleur Garton, Joel Mefford, Brian Black, Catherine Lomen-Hoerth, Andrew Dahl, Noah Zaitlen, *"Estimating the rate of cell type degeneration from epigenetic sequencing of cell-free DNA"*, BioRxiv, Jan 2020, [DOI]( https://doi.org/10.1101/2020.01.15.907022)
