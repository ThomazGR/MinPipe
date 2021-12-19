#!/usr/bin/bash

RINST=`which R`

if [ "$RINST" == "" ]; then
    echo "You do not have R installed. Please install it before installing the packages. 
        Check the MinPipe documention to know how."
    exit
else 
    R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')"
    R -e "BiocManager::install(version = '3.13')"
    R -e "BiocManager::install(c('biomaRt', 'rhdf5', 'ggplot2','dplyr','clusterProfiler',
        'org.Hs.eg.db','org.Mm.eg.db', 'cowplot'))"
    R -e "install.packages(c('devtools', 'docopt'))"
    R -e "devtools::install_github('pachterlab/sleuth')"

    echo "All packages have been installed for R execution"
    echo "Packages included: Bioconductor Manager, biomaRt, rhdf5, ggplot2, dplyr, DOSE, hsa & mmu organisms,
        devtools, docopt and sleuth."
    exit
fi;