# RNAseq-kallisto-sleuth
Using for Kallisto and Sleuth scripts used for RNA-seq trimmed (right now) data.

## Installing
- Install Conda or Miniconda
- Activate Bioconda channel
- Install Kallisto
- Install Bioconductor Manager
- Via BiocManager install biomaRt and rhdf5
- Install devtools R package
- Via devtools install Sleuth with `devtools::install_github("pachterlab/sleuth")`

## Using [UNDER CONSTRUCTION]
- Run `python3 main.py -h` or `python3 main.py --help` for a help message

### Arguments
- -c or --complement is the complement for paired-ended file names, if read 1 is always sample_R1.fq.gz and read 2 is sample_R2.fq.gz use `-c _R1 _R2` or `--complement _R1 _R2` so the code will iterate over samples with this complementary name.
- -s or --samples is the list of samples used to integrate with complement and iterate in the directory, e.g. `-s sample1 sample2 sample3` or `--sample sample1 sample2 sample3` the program will iterate as `sample1_R1.fq.gz` and `sample1_R2.fq.gz` as paired-ended.
- -o or --output is the directory used to create a folder and save the output of every Kallisto pseudo-alignment, creating a folder for every sample passed.
- -d or --dir is the directory where all `.fq.gz` files are.