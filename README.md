
# CelFiE

[![DOI](https://zenodo.org/badge/217374131.svg)](https://zenodo.org/badge/latestdoi/217374131)


## Overview
CelFiE (CELl Free dna decomposItion Expectation maximization) is an expectation maximization algorithm that takes as input, a reference panel and cell-free DNA methylation of several individuals. From this data, CelFiE will estimate the contribution of the reference tissues to the cfDNA of each individual, along with an arbitrary number of missing tissues not contained in the reference data.

For more details, please see our paper. 


## Installation

To install CelFiE, clone or fork this repository:

`git clone https://github.com/christacaggiano/celfie.git`.

All required packages can be installed using [Anaconda](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) using the environment file specified in `celfie_conda_env.yml`. Run:

`conda env create -f celfie_conda_env.yml -n celfie_env`

CelFiE was developed in Python 3.7.

## TL;DR Running CelFiE

To run CelFiE:

```bash
python EM/em.py <input_path> <output_directory> <num_samples> <--max_iterations> <--unknowns> <--parallel_job_id <--convergence> <--random_restarts>
```
For a detailed description of the parameters, see below. To run a test run of CelFiE with the default parameters run: 
```bash
python EM/em.py celfie_demo/sample_data.txt celfie_demo/sample_output 15 
```


#### CelFiE demo

Sample data is provided in `celfie_demo/`, along with a sample Jupyter Notebook for analyzing the output `demo.ipynb`. 



## Details of Running CelFiE

### Preparing Data

CelFiE expects the methylation data to be in the form # of methylated reads, # of total reads. For example it could look like:

```
CHR   START END METH DEPTH
chr1	10	11	44.0	63.0
chr1	50	51	71.0	133.0
chr1	60	61	89.0	115.0
```

CelFiE should work, in theory, on Illumina Chip data, if you estimate the read depth of each of the sites. However, we do not officially recommend this.

The input of CelFiE is a single txt file including both the reference data and the cfDNA, with a header indicating sample names (see `celfie_demo/sample_data.txt`). Essentially the file is the reference and cfDNA sample bed files combined. This data should look something like this:


```
CHROM START END SAMPLE1_METH SAMPLE1_DEPTH CHROM START END TISSUE1_METH TISSUE1_DEPTH
chr1	10	11	44.0	63.0  chr1	10	11	25.0	29.0
chr1	50	51	71.0	133.0 chr1	50	51	85.0	99.0
chr1	60	61	89.0	115.0 chr1	60	61	92.0	117.0
```


## Code

### EM Script

After preparing data as above, you can run EM script as follows:

```bash
python EM/em.py <input_path> <output_directory> <num_samples> <--max_iterations> <--unknowns> <--parallel_job_id <--convergence> <--random_restarts>
```

CelFiE takes several parameters. `Input_path`, `output_directory,` and `num_samples` are the only mandatory parameters. 

```bash
usage: em.py [-h] [-m MAX_ITERATIONS] [-u UNKNOWNS] [-p PARALLEL_JOB_ID]
             [-c CONVERGENCE] [-r RANDOM_RESTARTS]
             input_path output_directory num_samples

CelFiE - Cell-free DNA decomposition. CelFie estimated the cell type of origin
proportions of a cell-free DNA sample.

positional arguments:
  input_path            The path to the input file
  output_directory      The path to the output directory
  num_samples           Number of cfdna samples

optional arguments:
  -h, --help            show this help message and exit
  -m MAX_ITERATIONS, --max_iterations MAX_ITERATIONS
                        How long the EM should iterate before stopping, unless
                        convergence criteria is met. Default 1000.
  -u UNKNOWNS, --unknowns UNKNOWNS
                        Number of unknown categories to be estimated along
                        with the reference data. Default 1. Can be increased to 2+ for large samples. 
  -p PARALLEL_JOB_ID, --parallel_job_id PARALLEL_JOB_ID
                        Replicate number in a simulation experiment. Default
                        1.
  -c CONVERGENCE, --convergence CONVERGENCE
                        Convergence criteria for EM. Default 0.001.
  -r RANDOM_RESTARTS, --random_restarts RANDOM_RESTARTS
                        CelFiE will perform several random restarts and select
                        the one with the highest log-likelihood. Default 10.
```

### Output

CelFiE will output the tissue estimates for each sample in your input - i.e. the proportion of each tissue in the reference making up the cfDNA sample. See `celfie_demo/sample_output/1_tissue_proportions.txt` for an example of this output.

```
        tissue1 tissue2 .... unknown
sample1 0.05 0.08 .... 0.1
sample2 0.7 0.12 .... 0.2

```

CelFiE also outputs the methylation proportions for each of the tissues plus however many unknowns were estimated. This output will look like this:

```   
      tissue1  tissue2 ... unknown
CpG1  0.99 1.0 ... 0.3
CpG2  0.45 0.88 ... 0.1
```

Sample code for processing both of these outputs can be seen in `demo.ipynb`.

### L1 projection method

We also developed a method to project estimates onto the L1 ball, based on Duchi et al 2008. The code for this method is available at `EM/projection.py`. It can be ran as

```python
python projection.py <output_dir> <replicate> <number of tissues> <number of sites> <number of individuals> <input depth> <reference depth> <tissue_proportions.pkl>
```

Sample tissue proportions are included at `EM/simulations/unknown_sim_0201_10people.pkl`.

## Tissue Informative Markers

In our paper, we identified a set of tissue informative markers (TIMs). We claim that these are a good set of CpGs to use for decomposition.

#### Pre-selected TIMs

TIMs are available at `TIMs/sample_tims.txt` for individual CpG TIMs, and `TIMs/sample_tims_summed.txt` for reads summed +/-250bp around a TIM. We recommend using the `TIMs/sample_tims_summed.txt` for improved decomposition performance.

The TIMs represent markers for the following tissues:

- dendritic cells
- endothelial cells
- eosinophils
- erythroblasts
- macrophages
- monocytes
- neutrophils
- placenta
- T-cells
- adipose
- brain
- fibroblasts
- heart left ventricle
- hepatocytes
- lung
- mammary gland
- megakaryocytes
- skeletal muscle myoblasts
- small intestine

Data was retrieved from the [ENCODE](https://www.encodeproject.org/search/?type=Experiment&assay_slims=DNA+methylation&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&assay_title=WGBS) and [Blueprint](http://dcc.blueprint-epigenome.eu/#/files) data portals. When available, two biological replicates per tissue were combined into one sample. The TIMs were then calculated on the combined sample.

Please note all data was converted to hg38 and all CpGs are reported as (Chrom, start, end), where the end position indicates the C in the CpG dinucleotide.  

#### Selecting TIMs

Code to find TIMs is located at `TIMs/tim.py`. This code takes a reference bedfile of all the tissues you would like to calculate TIMs for as input. See `TIMs/sample_input.txt.`

The TIM code can be run as:

```bash
python tim.py <input file> <output file> <num of tim/tissue> <num of tissues> <depth filter> <nan filter>
```

The number of TIMs per tissue can be adjusted, but note that as the number of TIMs approaches the number of CpGs, the less informative that TIM will be for that tissue.

The **depth filter** only will consider CpGs that have a median depth across all tissues greater than a user specified value. This is to ensure that low-coverage CpGs do not get selected as TIMs. The **NaN filter** will only consider CpGs that have less than a user specified number of missing values. This is to ensure a TIM isn't selected for a tissue because it is one of the few tissues with data at that location. The **number of tims/tissue** can vary. We find that 100 is a good number, and note that as the number of TIMs increase, the lower quality the TIMs will be, since we are selecting the top most informative CpGs/tissue (in other words, the top 100 most informative CpGs for pancreas will by definition, be "better" than the top 500).

For the sample data provided, we suggest:

```bash
python tim.py sample_input.txt tim.txt 100 19 15 2
```

#### Combining Reads

In our paper, we found that summing all reads +/-250bp offered improved performance when decomposing. To do this for TIMs generated as output of `tim.py`, we provide a shell script `TIMs/tim.sh` to call TIMs and sum data.

This script can be updated to change the following parameters:

```bash
input_file=sample_input.txt
output_file=sample_tims.txt
summed_file=sample_tims_summed.txt
sum_window_size=500
number_tims=100
number_tissues=19
depth_filter=15
na_filter=2
```

The pipeline can then be ran as
```bash
./tim.sh
```

## Figures

Jupyter notebooks to reproduce figures and statistical analyses for the final version of this manuscript can be found in `paper_figures` directory.

## Acknowledgements

Thanks to Arya Boudaie for help with writing and reviewing this code and to Antoine Passemiers for their tremendous help in speeding up the EM calculation.

## Contact
For any questions with this code, please contact christa@g.ucla.edu. I am more than happy to help and very open to any suggestions!

## Citation

Christa Caggiano, Barbara Celona, Fleur Garton, Joel Mefford, Brian Black, Catherine Lomen-Hoerth, Andrew Dahl, Noah Zaitlen, *"Estimating the rate of cell type degeneration from epigenetic sequencing of cell-free DNA"*, BioRxiv, Jan 2020, [DOI]( https://doi.org/10.1101/2020.01.15.907022)
