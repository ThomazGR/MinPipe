FROM ubuntu

SHELL ["/bin/bash", "-c"]
ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=America/Sao_Paulo

# Installs miniconda and removes R base to be installed outisde of miniconda
RUN apt update -qq
RUN apt upgrade --yes
RUN apt install wget git --yes
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
RUN bash ~/miniconda.sh -b -p miniconda
RUN export PATH="miniconda/bin:$PATH"
RUN source /root/.bashrc
RUN /miniconda/bin/conda config --add channels conda-forge
RUN /miniconda/bin/conda config --add channels defaults
RUN /miniconda/bin/conda config --add channels bioconda
RUN /miniconda/bin/conda install -c bioconda fastqc trim-galore cutadapt kallisto picard --yes
RUN /miniconda/bin/conda remove r-base --yes

# Install R and its libraries
RUN apt update -qq
RUN apt install --no-install-recommends software-properties-common dirmngr --yes
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
RUN apt install --no-install-recommends r-base r-base-dev --yes
RUN apt install libxml2-dev libcurl4-openssl-dev libssl-dev --yes
RUN R -e "if (requireNamespace('BiocManager', quietly = TRUE) == FALSE) install.packages('BiocManager')"
RUN R -e "BiocManager::install(version = '3.13', ask = FALSE)"
RUN R -e "BiocManager::install(c('biomaRt', 'rhdf5', 'ggplot2','dplyr','clusterProfiler', 'org.Hs.eg.db','org.Mm.eg.db'))"
RUN R -e "install.packages(c('devtools', 'docopt'))"
RUN R -e "devtools::install_github('pachterlab/sleuth')"

# Clone the project to be used
RUN git clone https://github.com/ThomazGR/MinPipe.git