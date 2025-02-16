# radseq_variance_analysis.py: RAD-seq Variant Analysis

This script performs **RAD-seq variant analysis** to locate sequencing sites, detect variations within those sites, and calculate relevant statistics. It supports two modes: **SingleRad** and **ddRad**, which are different methods for sequencing based on restriction enzyme digestion.

## Requirements

Before running the script, ensure you have the following dependencies installed:

- Python 3.x
- Required Python packages (install with `conda/mamba`):
  - `pyvcf`
  - `pyfaidx`
  - `Bio.Restriction.Restriction_Dictionary`
 

Example of running the script using example data:
- python radseq_variance_analysis.py Mzebra_GT3_Chr1.fasta YH_MC_samples_Chr1.vcf SingleRad EcoRI 
- python radseq_variance_analysis.py Mzebra_GT3_Chr1.fasta YH_MC_samples_Chr1.vcf ddRad EcoRI --RE2 HindIII


## Running the script
```
python radseq_analysis.py --GenomeFile /path/to/genome.fasta --Mode SingleRad --RE1 EcoRI --VCFFile /path/to/variants.vcf
```


