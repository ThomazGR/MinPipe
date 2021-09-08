#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.13")
#BiocManager::install("biomaRt", "rhdf5")
#install.packages("devtools")
#devtools::install_github("pachterlab/sleuth")

suppressMessages({
  library("sleuth")
  library("biomaRt")
  library("cowplot")
})

getwd()

# Start metadata and build paths to abundance.h5 files
metadata <- read.table("resultados_kallisto_sleuth/metadata.txt", sep=";", header=T, stringsAsFactors = F)
metadata <- dplyr::mutate(metadata, 
                          path=file.path("resultados_kallisto_sleuth", "results", "hipotalamo",
                                         Run_s, "abundance.h5"))
metadata <- dplyr::rename(metadata, sample=Run_s)

# Select base (control) group that comaprison will be made with
metadata$treatment <- as.factor(metadata$treatment)
metadata$treatment <- relevel(metadata$treatment, "PO")
factor(metadata$treatment)

head(metadata)

# Start biomaRt annotation genes from ensembl
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "mmusculus_gene_ensembl",
                         host = "ensembl.org")
ttg <- biomaRt::getBM(
  attributes = c("ensembl_transcript_id", "transcript_version",
                 "ensembl_gene_id", "external_gene_name", "description",
                 "transcript_biotype"),
  mart = mart)
ttg <- dplyr::rename(ttg, target_id = ensembl_transcript_id,
                    ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
ttg <- dplyr::select(ttg, c('target_id', 'ens_gene', 'ext_gene'))
head(ttg)

# Build Sleuth object and use transformation function so beta will be log2FC
r <- sleuth_prep(metadata, ~ treatment, 
                 transformation_function = function(x) log2(x + 0.5))

# Fit the object with full or reduced formula
r <- sleuth_fit(r)
r <- sleuth_fit(r, ~1, 'reduced')

# Check possible models and tests
models(r)
tests(r)

# Run likelihood-ratio test comparing the two models for more than 2 groups
r <- sleuth_lrt(r, 'reduced', 'full')
res_lrt <- sleuth_results(r, 'reduced:full', test_type = 'lrt')
res_lrt_de <- dplyr::filter(res_lrt, pval <= 0.05)
res_lrt_de_ann <- dplyr::left_join(res_lrt_de, ttg, by=c("target_id"))
head(dplyr::filter(res_lrt_de_ann, !is.na(ext_gene)))
nrow(dplyr::filter(res_lrt_de_ann, !is.na(ext_gene)))

# Run Wald test comparing two groups
r <- sleuth_wt(r, "treatmentINTPO")
res_wt <- sleuth_results(r, 'treatmentINTPO')
res_wt_de <- dplyr::filter(res_wt, pval <= 0.05) # pval <= 0.05) qval <= 0.05)
res_wt_de_ann <- dplyr::left_join(res_wt_de, ttg, by=c("target_id"))
head(dplyr::filter(res_wt_de_ann, !is.na(ext_gene)))
nrow(dplyr::filter(res_wt_de_ann, !is.na(ext_gene)))
table <- dplyr::filter(res_wt_de_ann, !is.na(ext_gene))
table <- dplyr::filter(table, b <= -1 | b >= 1)
write.table(x = as.data.frame(dplyr::select(table, c("target_id", "ens_gene", "ext_gene", "pval", "qval", "b", "se_b"))), 
            file = "~/resultados_kallisto_sleuth/pairwise_POvsINTPO_annotated.txt", 
            sep = '\t', 
            quote = F,
            row.names = F)
