# CELFIE (CELl Free dna decomposItion Expectation maximization)

## Overview
Expectation maximization algorithm to decompose complex mixtures of cfDNA into the tissues generating the cfDNA fragments. Required input is methylated cfDNA, either WGBS or Illumina 450K with downsampled read counts, and a reference panel of tissues to estimate their contribution to the cfDNA. CelFiE can estimate an arbitrary number of unknown or missing tissues from your reference that are truly in the cfDNA mixtures.

## Code
Anaconda environment file specified in `celfie_conda_env.yml`

Full implementation of the EM model at `EM/em.py`. Code to generate simulations can be found in `EM/simulations`.

To run, modify run_real_data.sh
`qsub run_real_data.sh`

## Citation
Pending
