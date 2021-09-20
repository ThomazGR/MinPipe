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
  library("dplyr")
  library("ggplot2")
  library("docopt")
})

create.args <- function(){
  "Sleuth implementation from Kallisto pseudoalignment
  
  Usage:
    sleuth_results.R (-f --file) [-g --groups] [--no-volcano] [-o --organism] [-p --path]
  
  Options:
    -h --help     Show this screen.
    -f <file> --file <file>     path/for/metadata.txt
    -g <groups>... --groups <groups>...  Gets groups.
    -o <org> --organism <org> Organism name to be used in biomaRt, either mmu or hsa
    -p <path> --path <path> path/to/save/volcano_and_tables
    --no-volcano  Does not build volcano.
  " -> doc
  
  opt <- docopt::docopt(doc)
  
  return(opt)
}

start.metadata <- function(args){
  # Start metadata and build paths to abundance.h5 files
  metadata <- read.table(args[["--file"]], sep=";", header=T, stringsAsFactors = F)
  metadata <- dplyr::mutate(metadata, 
                            path=file.path("results", "figado",
                                            Run_s, "abundance.h5"))
  metadata <- dplyr::rename(metadata, sample=Run_s)
  
  if(is.null(args[["--groups"]][[1]])){
    groups <- unique(metadata$treatment)
  } else{
    groups <- strsplit(args[["--groups"]], ",")[[1]]
  }
  volcano <- !args[["--no-volcano"]]
  
  return(list(metadata=metadata, groups=groups, volcano=volcano))
}

change.base <- function(metadata, group){
  # Select base (control) group that comaprison will be made with
  metadata$treatment <- as.factor(metadata$treatment)
  metadata$treatment <- relevel(metadata$treatment, group)
  
  return(metadata)
}

get.annotation <- function(args){
  # Start biomaRt annotation genes from ensembl
  if(args[["--organism"]] == "mmu"){
    ds <- "mmusculus_gene_ensembl"
  } else if(args[["--organism"]] == "hsa"){
    ds <- "hsapiens_gene_ensembl"
  } else{
    stop("-o/--organism needs to be passed, either mmu or hsa.")
  }
  mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                           dataset = ds,
                           host = "ensembl.org")
  ttg <- biomaRt::getBM(
    attributes = c("ensembl_transcript_id", "transcript_version",
                   "ensembl_gene_id", "external_gene_name", "description",
                   "transcript_biotype"),
    mart = mart)
  ttg <- dplyr::rename(ttg, target_id = ensembl_transcript_id,
                       ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
  ttg <- dplyr::select(ttg, c('target_id', 'ens_gene', 'ext_gene'))
  
  return(ttg)
}

prepare.sleuth <- function(metadata){
  # Build Sleuth object and use transformation function so beta will be log2FC
  so <- sleuth_prep(metadata, ~ treatment, 
                          transformation_function = function(x) log2(x + 0.5))
  
  # Fit the object with full or reduced formula
  so <- sleuth_fit(so)
  so <- sleuth_fit(so, ~1, 'reduced')
  
  return(so)
}

run.lrt <- function(so, ttg){
  # Run likelihood-ratio test comparing the two models for more than 2 groups
  so <- sleuth_lrt(so, 'reduced', 'full')
  res_lrt <- sleuth_results(so, 'reduced:full', test_type = 'lrt')
  res_lrt_de <- dplyr::filter(res_lrt, pval <= 0.05)
  res_lrt_de_ann <- dplyr::left_join(res_lrt_de, ttg, by=c("target_id"))
  res_lrt_de_ann_f <- dplyr::filter(res_lrt_de_ann, !is.na(ext_gene))
  res_lrt_de_ann_f <- dplyr::filter(res_lrt_de_ann_f, b <= -1 | b >= 1)
  head(dplyr::filter(res_lrt_de_ann_f, !is.na(ext_gene)))
  nrow(dplyr::filter(res_lrt_de_ann_f, !is.na(ext_gene)))
  
  return(res_lrt_de_ann_f)
}

run.wald <- function(so, ttg, comp){
  # Run Wald test comparing two groups
  so <- sleuth_wt(so, comp)
  res_wt <- sleuth_results(so, comp)
  res_wt_de_ann <- dplyr::left_join(res_wt, ttg, by=c("target_id"))
  res_wt_de_ann_f <- dplyr::filter(res_wt_de_ann, !is.na(ext_gene))

  return(res_wt_de_ann_f)
}

run.volcano <- function(table){
  data <- data.frame(pval = -log10(table$pval), 
                            lfc = table$b, 
                            row.names = table$target_id)
  
  # Remove any rows that have NA as an entry
  data <- na.omit(data)
  
  # Color the points which are up or down
  ## If fold-change > 0 and pvalue > 1.3 (Increased significant)
  ## If fold-change < 0 and pvalue > 1.3 (Decreased significant)
  data <- mutate(data, color = case_when(data$lfc > 1 & data$pval > 1.3 ~ "Increased",
                                         data$lfc < -1 & data$pval > 1.3 ~ "Decreased",
                                         data$pval < 1.3 ~ "nonsignificant"))
  # Make a basic ggplot2 object with x-y values
  volc <- ggplot(data, aes(x = lfc, y = pval, color = color))
  
  # Add ggplot2 layers
  volcano_plot <- volc +   
    ggtitle(label = "Volcano Plot", subtitle = "Colored by fold-change direction") +
    geom_point(size = 2.5, alpha = 0.8, na.rm = T) +
    scale_color_manual(name = "Directionality",
                       values = c(Increased = "#008B00", Decreased = "#CD4F39", nonsignificant = "darkgray")) +
    theme_bw(base_size = 14) + # change overall theme
    theme(legend.position = "right") + # change the legend
    xlab("logFC") + # Change X-Axis label
    ylab(expression(-log[10]("adjusted p-value"))) + # Change Y-Axis label
    geom_hline(yintercept = 1.3, colour = "darkgrey") + # Add p-adj value cutoff line
    scale_y_continuous(trans = "log1p") # Scale yaxis due to large p-values
  
  return(volcano_plot)
}

write <- function(table, path){
  write.table(x = as.data.frame(dplyr::select(table, c("target_id", "ens_gene", "ext_gene", "pval", "qval", "b", "se_b"))), 
              file = path, 
              sep = '\t', 
              quote = F,
              row.names = F)
}

# Main paths
arg <- create.args()
# Example args
# arg <- list(
#   "--groups"="PO,HFD,INTPO,INTHFD",
#   "--no-volcano"=FALSE,
#   "--file"="metadata_figado.txt",
#   "--organism"="mmu"
# )

`%notin%` <- Negate(`%in%`)

if(substr(arg[["--file"]], nchar(arg[["--file"]])-4+1, nchar(arg[["--file"]])) != ".txt"){
  stop("Metadata has to be a metadata.txt file. Take a look at the example in the Github repo.")
}
if(arg$organism %notin% c("mmu", "hsa")){stop("-o/--organism not supported, please use only for hsa or mmu.")}
if(length(strsplit(arg$groups, ",")[[1]]) <= 1){stop("-g/--groups argument has only 1 value (group) which is not right.")}
if(substr(arg$path, nchar(arg$path), nchar(arg$path)) == "/"){
  arg$path <- substr(arg$path, 1, nchar(arg$path)-1)
  arg[["--path"]] <- substr(arg$path, 1, nchar(arg$path)-1)
} else{}

attach(start.metadata(args = arg))

cat(paste("\nThe base (or control) group is", groups[1]))
cat("\nPlease make sure this is right\n")

comb <- combn(groups, 2)

ann <- get.annotation(arg)

for(i in 1:ncol(comb)){
  gs <- comb[,i]
  if(levels(as.factor(metadata$treatment))[1] == gs[1]){
  } else{
    metadata <- change.base(metadata = metadata, group = gs[1])
  }
  print(paste("Base (control) group:", levels(as.factor(metadata$treatment))[1]))
  print(paste("Compared group:", gs[2]))
  
  if(levels(as.factor(metadata$treatment))[1] == gs[1]){
    names <- paste0(gs[1], "vs", gs[2])
  } else{
    stop("Not the same base group")
  }
  
  so <- prepare.sleuth(metadata)
  
  table <- run.wald(so, ann, paste0("treatment",gs[2]))
  
  if(!args[["--no-volcano"]]){
    volcano <- run.volcano(table)
    ggsave(paste0(arg$path, "/volcano_", names, ".png"), plot=volcano, device="png")
    print(paste("Saved volcano for", names))
  } else{}
  
  res <- dplyr::filter(table, pval <= 0.05)
  res <- dplyr::filter(res, b <= -1 | b >= 1)
  write(res, path=paste0(arg$path, "/tabela_annotated_", names, ".txt"))
  print(paste("Saved table for", names))
}

#models(so)
#tests(so)
#sleuth_live(so)
