# lost-in-space-public
public repository to recreate analyses in allelic correlation + noise transmission paper

## General
This code will generate all simulation data to replicate the analysis in the paper. However, despite using seeds for random number generation, small differences will likely arise from the data used in the paper. The generated data used in the paper is available at the dropbox link provided in the key resources table or upon request from the authors.

The limited data from Larsson et al. 2019 that was re-analyzed for the paper is provided in this repository and so does not need to be downloaded separately.

## Usage
All scripts use relative path names and are made to be run from their directory.

To generate simulation data:
Run /Data/GenerateDataScripts/allele_constant_signal/Scripts/make* scripts
then /Data/GenerateDataScripts/allele_constant_signal/Scripts/save* scripts
then /Data/GenerateDataScripts/allele_constant_signal/Scripts/simulate* scripts

Warning, this will generate ~300GB of data and took ~4 weeks of computational time on a 2019 iMac with a 3GHz 6 core CPU and 40GB of RAM.

To analyze the simulated data, either run the above scripts or download the provided data from Dropbox and place into /Data/GenerateDataScripts/allele_constant_signal/Data/

Then run /Data/GenerateDataScripts/allele_constant_signal/Scripts/analyze* scripts

Finally, to recreate final analysis and plots, please run R scripts provided in 
/Paper/plotScripts

## Citation
Please consider citing https://doi.org/10.1101/2021.11.26.470134
