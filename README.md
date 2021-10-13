# MinPipe <a href="https://github.com/ThomazGR/MinPipe"> <img align="right" src="./img/minpipe.png"> </a>
MinPipe is a minimal pipeline or workflow to be used for a series of RNA-seq data analysis from raw data to differentially expressed genes by using a mixture of pure Python, R and Shell script mainly for workspace management.

### Installing
- Install Conda or Miniconda by using the `conda.sh` inside `src` folder
	- Run `chmod +x conda.sh` and then `bash conda.sh`, and follow the next steps
- Activate Bioconda channel by running `packages.sh` inside `src` folder
	- Run `chmod +x packages.sh` and then `bash packages.sh`, and follow the next steps
- Install Bioconductor Manager
	- Open an R session and run:
	```{r}
	if (!requireNamespace("BiocManager", quietly = TRUE))
		install.packages("BiocManager")
	```
	- And then 
	```{r}
	BiocManager::install(version = "3.13")
	```
- Via BiocManager install all packages used
	- Run 
	```{r}
	BiocManager::install(c("biomaRt", "rhdf5", "ggplot2","clusterProfiler",
	"dplyr","org.Hs.eg.db","org.Mm.eg.db"))
	```
- Install devtools R package and docopt for CLI R interface for args
	- Run
	```{r}
	install.packages(c("devtools", "docopt"))
	```
- Via devtools install Sleuth
	- Install via this command 
	```{r}
	devtools::install_github("pachterlab/sleuth")
	```

### Using the main pipeline
- Run `python3 main.py -h` or `python3 main.py --help` for a help message
- Example of paired-ended analysis without pre-built index: 
	- `python3 main.py -c _read1 _read2 -s AT1_G1 AT2_G1 AT1_G2 AT2_G2 -t MMU_39.fastq.gz --threads 4 -b 100`
- Example of single-ended analysis with pre-built index: 
	- `python3 main.py -s AT1_G1 AT2_G1 AT1_G2 AT2_G2 -i MMU_39.idx --threads 4 -b 100`
#### Arguments
- -c or --complement is the complement for paired-ended file names, if read 1 is always sample_R1.fq.gz and read 2 is sample_R2.fq.gz use `-c _R1 _R2` or `--complement _R1 _R2` so the code will iterate over samples with this complementary name.
- -s or --samples is the list of samples used to integrate with complement and iterate in the directory, e.g. `-s sample1 sample2 sample3` or `--sample sample1 sample2 sample3` the program will iterate as `sample1_R1.fq.gz` and `sample1_R2.fq.gz` as paired-ended.
- -i or --index is the name of the index file to be used in pseudoalignment. If no file has been passed it has to have a UNIQUE .idx file in `index` folder or it will raise an error, or it has to be passed the -t/--transcript argument to build a new index. It is an optional argument.
- -t or --transcript is the name of the transcript file to be indexed if no index has been passed. Also an optional argument.
- --threads refers to the number of threads to be used in quantification for Kallisto. Default: 1.
- -b or --bootstrap is the number of bootstrap samples. Default: 100
- --single is the flag to indicate single-ended quantification without complements. An optional argument.

### How to work with Kallisto results using Sleuth R package
- Run `Rscript sleuth_results.R [arguments]`

#### Arguments [UNDER CONSTRUCTION]
