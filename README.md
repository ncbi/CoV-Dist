![CoV-Dist ](./img.svg)

## Overview
A python script to measure the dissimilarity distance between SARS-CoV-2 metagenome samples.

It will 
1) Calculate the Yue & Clayton dissimilarity index of two samples based on the base-pair proportions at each genomic position;
2) Sum all the dissimilarity indexes across the genome;
3) Output the sum as a distance metric;
4) Visualize the pairwise distance matrix among samples using PCoA, MDS or t-SNE method.

## Prerequisites
```
# set conda environment
conda config --add channels plotly
conda create --name cov-dist python=3.8 biopython scikit-bio click pyvcf plotly=5.10.0 scipy=1.8.1 # scipy=1.8.1 is to make scipy compatible with scikit-bio
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
  -p, --prefix_list PATH  A file that lists the prefix path for vcf and
                          corresponding depth files.  [required]
  -r, --ref PATH          Input reference file. [default:
                          data/NC_045512.2.fasta]
  -c, --cov FLOAT         Only samples with reads mapped to equal to or
                          greater than this fraction of the genome will be
                          used for ordination analysis [default: 0.5].
  -d, --depth INTEGER     Depth cutoff to include positions when calculating
                          coverage [default: 10].
  -t, --threads INTEGER   Number of threads used for computing [default: all
                          available cpus].
  -v, --voc TEXT          VOC (Variant Of Concern) folder name containing VOC
                          vcf files, no depth files are needed (default is
                          None).
  -o, --outpath TEXT      Output file path.  [required]
  -h, --help              Show this message and exit.
```

### dist input
* prefix_list: A prefix list file (a prefix per line). For each prefix, a `<prefix>.vcf` and a `<prefix>.depth` are supposed to be found. If a depth file is not provided, default to set all position depths are valid.
* ref: A reference fasta file which is used for the vcf files indicated in prefix_list.
* cov: A threshold for coverage to filter out disqualified vcf (and depth) files.
* depth: A depth threshold to include positions when calculating coverage.
* threads: Number of cpus used for computing in parallel.
* voc: A VOC (Variant Of Concern) folder path containing VOC vcf files, no depth files are needed. It can be used to investigate how close custom samples are to VOC lineages.
* outpath: Output file path.

> prefix list file example:
```
data/test/3194
data/test/3270
data/test/B199
data/test/B319
...
```

> VOC folder structure example:
```
voc_vcf/
├── AY.10.vcf
├── AY.11.vcf
├── AY.133.vcf
├── ...
├── AY.9.2.vcf
├── AY.9.vcf
├── B.1.1.529.vcf
└── B.1.617.2.vcf
```

### plot command
```
Usage: covdist.py plot [OPTIONS]

  Perform ordination analysis and visualize results.

Options:
  -d, --distance_file PATH        Distance file obtained from cov-dist "dist"
                                  command.  [required]
  -m, --meta PATH                 Tab-delimited metadata for samples. Must
                                  have at least three columns - "sample",
                                  "collection_date" and "collection_site".
                                  [required]
  -c, --column TEXT               The column name in the meta data to color
                                  the samples. Default collection_date and
                                  collection_site are used to plot figures.
  -a, --analysis [PCoA|MDS|TSNE]  Method to show samples. Can be choosen from
                                  PCoA (default), MDS, TSNE. Case is not
                                  sensitive.
  -v, --voc_meta PATH             Tab-delimited metadata file for VOC. At
                                  least two columns are needed: "sample" (VOC
                                  name, e.g: BA.1, AY.4 and etc.) and "group"
                                  (lineage group, e.g. Omicron, Delta and
                                  etc.).
  -o, --out TEXT                  Output folder name  [required]
  -h, --help                      Show this message and exit.
```

### plot input
* distance_file: The distance matrix resulting from "dist" command or any other distance matrix as long as its rownames and colnames are the same to the basename of the prefices and in the same order as the those in prefix file (e.g. if prefix is "a/b/c/xyz", the corresponding rowname and colname are "xyz").
* meta: A metadata file (tab-delimited) and it must has "sample", "collection_date" and "collection_site" columns.
* column: A column name in the meta data used to color samples in the pcoa plot. By default, "collection_date" and "collection_site" will be used to generate a special scatter plots. Otherwise, the program will use a provided column.
* analysis: Method to show samples. Can be choosen from PCoA (default), MDS, TSNE. Since case is not sensitive, uppercase and lowercase letters are both allowed.
* voc_meta: A metadata file for VOC (tab-delimited). It must has "sample" and "group" columns.
* out: Output folder name.

> metadata example:
```
sample  collection_site collection_date sample_group
2367    Arizona 2021-07-22      Batch1
2370    Arizona 2021-07-22      Batch1
4010    Arizona 2022-06-06      Batch10
3939    Arizona 2022-05-23      Batch10
4011    Arizona 2022-06-06      Batch10
4012    Arizona 2022-06-06      Batch10
4013    Arizona 2022-06-06      Batch10
...
```

> metadata for voc example:
```
sample  group
BA.1.14.2       Omicron
BA.5.7  Omicron
BF.21   Omicron
BA.5.6.1        Omicron
BA.1.14 Omicron
BA.1.6  Omicron
AY.107  Delta
AY.30   Delta
...
```

### Output summary
After running dist and plot, a total of 4 files will be generated as listed below:
```
1. distance.txt
2. <PCoA|MDS|TSNE>_ordination_result.txt
3. <pcoa|mds|tsne>_plot_2d.html
4. <pcoa|mds|tsne>_plot_3d.html
```

* Calculated distance matrix: `distance.txt`. It could be used as input for other dimensionality reduction methods.
* Output from the PCoA, MDS or t-SNE: `<PCoA|MDS|TSNE>_ordination_result.txt`. It could be used as input data for your preferred data visualization tool to plot the results.
* Plotly html (including both 3d and 2d scatter plots) with which user can interacte in the broswer and download a perferable snapshot as PNG: `<pcoa|mds|tsne>_plot_<3d|2d>.html`

## Notes
* You can use (UShER)[https://github.com/yatisht/usher] to retrive and update SARS-CoV-2 lineage information and convert that into `vcf` format
  ```
  wget http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/public-latest.all.masked.pb.gz
  matUtils extract -i public-latest.all.masked.pb.gz -C lineagePaths.txt
  mkdir -p voc_vcf
  python ./usher_lineage2vcf.py lineagePaths.txt voc_vcf
  ```
* If there is zero coverage for a position in the alignment, we assume the dissimilarity in that position is zero, which could lead to an underestimation of the distance between samples with low coverage.
