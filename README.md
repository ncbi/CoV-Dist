# CoV-Dist 

## Overview
A python script to measure the dissimilarity distance between SARS-CoV-2 metagenome samples.

It will 
1) Calculate the Yue & Clayton dissimilarity index of two samples based on the base-pair proportions at each genomic position;
2) Sum all the dissimilarity indexes across the genome;
3) Output the sum as a distance metric.

## Prerequisites
#### TBA
+ [Biopython](https://biopython.org/)
+ [pysam](https://github.com/pysam-developers/pysam)(>=0.15.3)
+ [scikit-bio](https://github.com/biocore/scikit-bio) (>=0.5.6)
+ [seaborn](https://github.com/mwaskom/seaborn)

```
# set conda environment
conda create --name cov-dist python=3.7 biopython seaborn scikit-bio pysam
conda activate cov-dist
```


## Tutorial

```
Usage: covdist.py [OPTIONS] COMMAND [ARGS]...

Options:
  -h, --help  Show this message and exit.

Commands:
  dist  Compute pairwise metagenome distance for SARS-CoV-2 samples.
  plot  Perform ordination analysis and visualize results.
```

### dist command
```
Usage: covdist.py dist [OPTIONS]

  Compute pairwise metagenome distance for SARS-CoV-2 samples.

Options:
  -p, --prefix_list PATH  A file that lists the prefix for vcf and
                          corresponidng depth files.  [required]
  -r, --ref PATH          Input reference file. [default:
                          data/NC_045512.2.fasta]
  -c, --cov FLOAT         Only samples with reads mapped to equal to or
                          greater than this fraction of the genome will be
                          used for PcoA analysis [default: 0.5].
  -d, --depth INTEGER     Depth cutoff to include positions when calculating
                          coverage [default: 10].
  -t, --threads INTEGER   Number of threads used for computing [default: all
                          available cpus].
  -o, --out TEXT          Output folder name  [required]
  -h, --help              Show this message and exit.
```

### dist input
* prefix_list: A prefix list file (a perfix per line). For each prefix, a `<prefix>.vcf` and a `<prefix>.depth` are supposed to be found. If a depth file is not provided, default to set all position depths are valid.
* ref: A reference fasta file which is used for the vcf files indicated in prefix_list.
* cov: A threshold for coverage to filter out disqualified vcf (and depth) files.
* depth: A depth threshold to include positions when calculating coverage.
* threads: Number of cpus used for computing in parallel.
* out: Output folder name.


### plot command
```
Usage: covdist.py plot [OPTIONS]

  Perform ordination analysis and visualize results.

Options:
  -d, --distance_file PATH  Distance file obtained from cov-dist "dist"
                            command.  [required]
  -m, --meta PATH           Meta data for samples.  [required]
  -c, --column TEXT         The column name in the meta data to color the
                            samples.  [required]
  -o, --out TEXT            Output folder name  [required]
  -h, --help                Show this message and exit.
```
### plot input
* distance_file: The distance matrix resulting from "dist" command or any other distance matrix as long as its rownames and colnames are the same to the basename of the prefices and in the same order as the those in prefix file (e.g. if prefix is "a/b/c/xyz", the corresponding rowname and colname are "xyz").
* meta: A metadata file (tab-delimited) but it must has a "sample" column.
* column: A column name in the meta data used to color samples in the pcoa plot.
* out: Output folder name.


### Output summary
After running dist and plot, in the output directory, 3 files will be generated as listed below:
```
data/test/testout
├── distance.txt
├── ordination_result.txt
└── pcoa_plot.html
```

* Calculated distance matrix: `distance.txt`. It could be used as input for other dimensionality reduction methods.
* Output from the principal coordinates analysis (PcoA): `ordination_result.txt`. It could be used as input for your preferred data visualization tool to plot the PcoA.
* Plotly html for the PcoA, user can interacte with it in the broswer and download a perferable snapshot as PNG: `pcoa_plot.html`

## Notes
* If there is zero coverage for a position in the alignment, we assume the dissimilarity in that position is zero, which could lead to an underestimation of the distance between samples with low coverage.
* This method in theory could be measure dissimilarity between metagenomes and isolate genomes.
