#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.13")
#BiocManager::install("biomaRt")

suppressMessages({
  library("sleuth")
  library("biomaRt")
  library("cowplot")
})

getwd()
metadata <- read.table("kallisto/metadata.txt", sep=";", header=T, stringsAsFactors = F)
metadata <- dplyr::mutate(metadata, 
                          path=file.path("kallisto", "results", "hipotalamo",
                                         Run_s, "abundance.h5"))
metadata <- dplyr::rename(metadata, sample=Run_s)
head(metadata)

mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "mmusculus_gene_ensembl",
                         host = "ensembl.org")
# host = "dec2015.archive.ensembl.org")
# host = "ensembl.org")
ttg <- biomaRt::getBM(
  attributes = c("ensembl_transcript_id", "transcript_version",
                 "ensembl_gene_id", "external_gene_name", "description",
                 "transcript_biotype"),
  mart = mart)
ttg <- dplyr::rename(ttg, target_id = ensembl_transcript_id,
                    ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
ttg <- dplyr::select(ttg, c('target_id', 'ens_gene', 'ext_gene'))
head(ttg)


r <- sleuth_prep(metadata, ~ treatment, 
                 transformation_function = function(x) log2(x + 0.5))
r <- sleuth_fit(r)
r <- sleuth_fit(r, ~1, 'reduced')

r <- sleuth_lrt(r, 'reduced', 'full')
r <- sleuth_wt(r, paste0("treatmentPO"))	100

res_wt <- sleuth_results(r, 'treatmentPO')
res_wt_de <- dplyr::filter(res_wt, pval <= 0.05) # pval <= 0.05) qval <= 0.05)
res_wt_de_ann <- dplyr::left_join(res_wt_de, ttg, by=c("target_id"))
#nrow(res_wt_de_ann)
#nrow(dplyr::filter(res_wt_de_ann, !is.na(ext_gene)))

res_lrt <- sleuth_results(r, 'reduced:full', test_type = 'lrt')
res_lrt_de <- dplyr::filter(res_lrt, pval <= 0.05)
res_lrt_de_ann <- dplyr::left_join(res_lrt_de, ttg, by=c("target_id"))
#head(res_lrt_de_ann)

#head(kallisto_table(r))
#head(r)
#tests(r)

#res <- merge(res_lrt, res_wt[, c('target_id', 'b', 'se_b', 'mean_obs')], on = 'target_id', sort = FALSE)

#sleuth_live(r)

# so <- sleuth_prep(metadata, target_mapping = ttg, 
#                   aggregation_column = 'ens_gene', extra_bootstrap_summary=T)
# 
# so <- sleuth_fit(so, ~treatment, 'full')
# 
# so <- sleuth_wt(so, which_beta = "treatmentPO", which_model = 'full')
# gene_table <- sleuth_results(so, test="treatmentPO", test_type = "wt", which_model = "full", show_all = F)
# tests(so)
# sleuth_live(so)
# 
# head(so)
