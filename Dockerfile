FROM continuumio/miniconda3:master
RUN conda config --add channels conda-forge defaults r bioconda
RUN conda install r-base r-essentials
RUN conda install -c anaconda git wget --yes
RUN conda install -c bioconda fastqc trim-galore cutadapt kallisto picard --yes
RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')"
RUN R -e "BiocManager::install(version = '3.13')"
RUN R -e "BiocManager::install(c('biomaRt', 'rhdf5', 'ggplot2','dplyr','clusterProfiler', 'org.Hs.eg.db','org.Mm.eg.db'))"
RUN R -e "install.packages(c('devtools', 'docopt'))"
RUN R -e "devtools::install_github('pachterlab/sleuth')"
RUN git clone https://github.com/ThomazGR/MinPipe.git
# Copy files
# Transcript if needed, index if needed
# Samples is obligatory
# YAML or Json file if needed
COPY source dest
# How to pass parameters from docker run to main.py
CMD [ "python", "main.py" , "parameters_needed"]