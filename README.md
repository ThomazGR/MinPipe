# MinPipe <a href="https://github.com/ThomazGR/MinPipe"> <img align="right" src="./img/minpipe.png"> </a>
MinPipe is a minimal but fully logged pipeline or workflow to be used for a series of RNA-seq data analysis from raw data to differentially expressed genes by using a mixture of Python, R and Shell script mainly for workspace management.

### Installing
Rewritten version of installers include:
- env.yml: The yaml file that sets specification for environment building used in miniconda/mamba. It lists all needed packages to run every part of MinPipe.
- installers/conda.sh: A shell script to download and install miniconda, mamba and create the needed environment. Go to the project root directory and run `bash installers/conda.sh`. 

### Using the main pipeline
- Run `python3 minpipe.py -h` or `python3 minpipe.py --help` for a help message
- Example of paired-ended analysis without pre-built index: 
	- `python3 minpipe.py -c _read1 _read2 -s AT1_G1 AT2_G1 AT1_G2 AT2_G2 -t MMU_39.fastq.gz --threads 4 -b 100`
- Example of single-ended analysis with pre-built index: 
	- `python3 minpipe.py -s AT1_G1 AT2_G1 AT1_G2 AT2_G2 -i MMU_39.idx --threads 4 -b 100`
#### Arguments
- -c or --complement is the complement for paired-ended file names, if read 1 is always sample_R1.fq.gz and read 2 is sample_R2.fq.gz use `-c _R1 _R2` or `--complement _R1 _R2` so the code will iterate over samples with this complementary name.
- -s or --samples is the list of samples used to integrate with complement and iterate in the directory, e.g. `-s sample1 sample2 sample3` or `--sample sample1 sample2 sample3` the program will iterate as `sample1_R1.fq.gz` and `sample1_R2.fq.gz` as paired-ended.
- -i or --index is the Name of the index file to be used in pseudoalignment. Either `index` or `transcript` has to be passed.
- -t or --transcript is the Name of the transcript file to be indexed. `mmu` or `hsa` can be passed so the transcript will be downloaded automatically and index will be built.
- --threads refers to the number of threads to be used in quantification for Kallisto. Default: 1.
- -b or --bootstrap is the number of bootstrap samples. Default: 100
- --single is the flag to indicate single-ended quantification without complements. An optional argument.
- --ext-qc is a flag to indicate that will have extensive QC. **MAY NEED MORE FILES**
- --json pass the Json file name that has to be located inside the input folder. The user can create separated folders inside the input, e.g. input/params/parameters.json.
- --yaml pass the YAML/YML file name that has to be located inside the input folder. The user can do the same as the Json file creating folders, e.g. input/params/parameters.yml.

### How to work with Kallisto results using Sleuth R package
- Run `Rscript minpipe.R [arguments]`
- Pass `Rscript minpipe.R -h` to see the help text
#### Arguments
- -f or --file is the argument needed for the metadata.txt, passing as a path/to/metadata.txt
- -o or --organism is the name of the organism to be used for gene annotation, either `mmu` or `hsa`, others will be supported lately
- -p or --path is the path/to/kallisto/results where it would have SAMPLE_FOLDER/abundance.h5 files.
- -r or --results is the name of the path/to/save/results for tables and visualization.
- -s or --separator string used as a separator for metadata file. Default is ;
- --no-volcano is a flag that will force no volcano image creation

