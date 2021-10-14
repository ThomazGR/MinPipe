# Append new channel to download conda packages
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda

# Download all packages needed to run RNA-seq workflow (star not included)
conda install -c anaconda git wget --yes
conda install -c bioconda fastqc trim-galore cutadapt kallisto, picard --yes

echo "All packages have been installed from Bioconda"
echo "Packages included: fastqc, trim galore, cutadapt, kallisto and picard."