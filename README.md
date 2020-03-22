# CoroNet

[![CoroNet Video Explaining the method](http://img.youtube.com/vi/D_jKHQ0AHQA/0.jpg)](http://www.youtube.com/watch?v=D_jKHQ0AHQA)

## Installation
You need python and R installed. You will need the following python packages:
* matplotlib
* basemap
* geopy
* Biopython

and the following R packages
* ape
* ggsci

Additionally you will need the following commandline tools

MAFFT for constructing multiple sequence alignments (https://mafft.cbrc.jp/alignment/software/)
RAxML-NG for building phylogenetic trees (https://academic.oup.com/bioinformatics/article/35/21/4453/5487384)


## Usage
The [notebook](investigate_spread.ipynb) analyzes the Covid-19 spread through sequence evolution. The number of sequences in the maps does not reflect the number of cases because there iss only a very(!) small subset of patients for which viral sequences are obtained and those sequences are also publically available.
For the connections between places, we use the genomic distances (number of basepair substitutions). This way, we can track the virus spread back in time.
