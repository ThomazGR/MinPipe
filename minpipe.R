# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Motivation: Automated way of finding differentially expressed genes for every match
# Author: Thomaz Guadagnini
# Github: @tzgr
# Email: thzgr@tuta.io
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

suppressMessages({
    requireNamespace("tximport") # import Kallisto results .h5 to DESeq2 matrix-way
    requireNamespace("DESeq2") # statistical analysis by gene counts from Kallisto
    requireNamespace("stringr") # needed to deal with string manipulation
    requireNamespace("biomaRt") # retrieve annotation data for gene list
    requireNamespace("ggplot2") # used for plotting kegg/go results in a bar plot
    requireNamespace("cowplot") # used for plotting kegg/go results in a bar plot
    requireNamespace("pheatmap")
    requireNamespace("SummarizedExperiment")
    requireNamespace("clusterProfiler") # enrichment for KEGG and GO from gene list
    library("dplyr") # import if needed to manage deeper a dataframe
})

create.args <- function() {
    "DEGs discovery using Kallisto abundance results, tximport and DESeq2 for statistical analysis 

    Usage:
      minpipe.R (-f --file) [--no-volcano] [-o --organism] [-p --path] (-r --results)

    Options:
        -h --help     Show this screen.

        -f <file> --file <file>     path/for/metadata.txt in a format of te following:
            |-------|-----------|
            | Run_s | treatment |
            |-------|-----------|
            | Samp1 |  Control  |
            | Samp2 |  Control  |
            | Samp3 |  Treated  |
            | Samp4 |  Treated  |
            |-------|-----------|

        -o <org> --organism <org> Organism name to be used in biomaRt, either mmu or hsa.

        -p <path> --path <path> path/to/kallisto/results where it would have SAMPLE_FOLDER/abundance.h5 files.

        -r <results> --results <results> path/to/save/results for tables and visualization.
    
        -s <separator> --separator <separator> String used as a separator for metadata file. Default is ;

        --no-volcano  Does not build volcano.
    " -> doc

    opt <- docopt::docopt(doc)

    return(opt)
}

start.metadata <- function(file = NULL, separator = NULL) {
    # Start metadata and build paths to abundance.h5 files
    if (!is.character(file)) {quit("File is not a string to be used. Check --help.")}
    if (grepl("https://docs.google.com/spreadsheets", file, fixed = T) &
        grepl("output=csv", file, fixed = T)) {
       tryCatch(
           expr = {
               metadata <- read.table(file, sep = ",", header = T, stringsAsFactors = F)
           },
           error = function(e) {
               print(paste("ERROR: \n", e))
               quit()
           },
           warning = function(w) {
               print(paste("WARNING: \n", w))
               quit()
           }
       )
    } else {
       metadata <- read.table(file, sep = separator, header = T, stringsAsFactors = F)
    }
    return(metadata)
}

import.kallisto.tx <- function(base_path = NULL, metadata = NULL) {
    if (!is.data.frame(metadata)) {
        quit("Metadata passed is not a data.frame.")
    }
    if (!is.character(base_path)) {
        quit("Base path is not a string. Check --help.")
    }
    
    needed_fields <- c("Run_s", "treatment")
    if (!all(needed_fields %in% colnames(metadata))) {
        quit("Metadata file needs to have both Run_s and treatment field. Check --help.")
    }

    files <- file.path(base_path, metadata$Run_s, "abundance.h5")
    
    txi.kallisto <- tximport::tximport(files, 
                                       type="kallisto",
                                       txOut=TRUE)
    
    sampleTable <- data.frame(condition=metadata$treatment)
    rownames(sampleTable) <- colnames(txi.kallisto$counts)
    
    imported_matrix <- DESeq2::DESeqDataSetFromTximport(txi.kallisto, 
                                                        sampleTable, 
                                                        ~condition)

    return(imported_matrix)
}

get.annotation <- function(organism = NULL) {
    if (organism %in% c("mmu", "mmusculus_gene_ensembl")) {
        org_ensembl <- "mmusculus_gene_ensembl"
    } else if (organism %in% c("hsa", "hsapiens_gene_ensembl")) {
        org_ensembl <- "hsapiens_gene_ensembl"
    } else {
        stop("Organism needs to be passed as either mmu or hsa.")
    }

    biomart_dataset <- biomaRt::useMart(biomart = "ensembl",
                         dataset = org_ensembl,
                         host = "https://www.ensembl.org")
    attributes_BM <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "transcript_version",
                                                "ensembl_gene_id", "external_gene_name", "description",
                                                "transcript_biotype", "entrezgene_id"),
                                    mart = biomart_dataset)
    attributes_BM <- dplyr::rename(attributes_BM, 
                                    target_id = ensembl_transcript_id,ens_gene = ensembl_gene_id,
                                    ext_gene = external_gene_name, entrez_id = entrezgene_id)
    attributes_BM <- dplyr::select(attributes_BM, 
                                    c('target_id', 'ens_gene', 'ext_gene', 'entrez_id'))
    
    return(attributes_BM)
}

run.volcano <- function(pval = NULL, lfc = NULL, names = NULL, lfc_cutoff = 1, pval_cutoff = 0.05) {
    if (all(c(is.null(pval), is.null(lfc), is.null(names)))) {quit("pval, lfc or names is null.")}

    pval_cutoff <- -log10(pval_cutoff)

    data <- data.frame(pval = -log10(pval),
                       lfc = lfc,
                       row.names = names)

    # Remove any rows that have NA as an entry
    data <- na.omit(data)

    # Color the points which are up or down
    ## If fold-change > 0 and pvalue > 1.3 (Increased significant)
    ## If fold-change < 0 and pvalue > 1.3 (Decreased significant)
    data <- dplyr::mutate(data, color = case_when(data$lfc > lfc_cutoff & data$pval > pval_cutoff ~ "Increased",
                                           data$lfc < -(lfc_cutoff) & data$pval > pval_cutoff ~ "Decreased",
                                           data$pval < pval_cutoff ~ "Nonsignificant"))
    # Make a basic ggplot2 object with x-y values
    volc <- ggplot2::ggplot(data, ggplot2::aes(x = lfc, y = pval, color = color))

    # Add ggplot2 layers
    volcano_plot <- volc +
        ggplot2::ggtitle(label = "Volcano Plot", subtitle = "Colored by fold-change direction") +
        ggplot2::geom_point(size = 2.5, alpha = 0.8, na.rm = T) +
        ggplot2::scale_color_manual(name = "Directionality",
                               values = c(Increased = "#008B00",
                                          Decreased = "#CD4F39",
                                          nonsignificant = "darkgray")) +
        ggplot2::theme_bw(base_size = 14) + # change overall theme
        ggplot2::theme(legend.position = "right") + # change the legend
        ggplot2::xlab("logFC") + # Change X-Axis label
        ggplot2::ylab(expression(-log[10]("adjusted p-value"))) + # Change Y-Axis label
        ggplot2::geom_hline(yintercept = 1.3, colour = "darkgrey") + # Add p-adj value cutoff line
        ggplot2::scale_y_continuous(trans = "log1p") # Scale yaxis due to large p-values

    return(volcano_plot)
}

pathway.enrichment <- function(gene_list = NULL, organism = NULL, path = NULL, kegg_file = NULL, go_file = NULL) {
    # - - - - - - - - - - - - -
    # Enrich with KEGG database
    # - - - - - - - - - - - - -
    if (is.null(gene_list)) { stop("No gene list has been passed.") }
    if (is.null(path)) { stop("No path has been passed.") }
    if (is.null(organism)) { stop("No organism passed. Needs to be either mmu or hsa.") }

    kegg_enrich <- clusterProfiler::enrichKEGG(gene = gene_list,
                                organism = organism,
                                pvalueCutoff = 0.05)

    print("Generate KEGG pathway results")
    # Get table of results
    kegg_table <- as.data.frame(kegg_enrich) %>% 
                    dplyr::arrange(desc(-log10(pvalue)))
    kegg_table <- head(kegg_table, n = 20)
            

    # KEGG plot
    kegg_bar <- ggplot2::ggplot(kegg_table, 
                                ggplot2::aes(x = reorder(Description, 
                                                         -log10(pvalue)), 
                                             y = -log10(pvalue))) +
      ggplot2::geom_bar(stat = 'identity', fill = "#6B2525") +
      ggplot2::geom_col(width = 0.7) +
      ggplot2::labs(title = "KEGG Enrichment Pathways", 
                    x = "KEGG Terms") +
      ggplot2::coord_flip()

    kegg_full <- as.data.frame(kegg_enrich) %>% 
      dplyr::arrange(desc(-log10(pvalue)))
    
    write.table(kegg_full, 
                file.path(path, paste0(kegg_file, "__kegg_full_table.txt")), 
                sep = "\t", 
                quote = F, 
                row.names = F)
    
    ggplot2::ggsave(file.path(path, paste0(kegg_file, "__kegg_bar.png")), 
           plot = kegg_bar, 
           device = "png")
    
    print("Saved KEGG full table and bar plot")

    # - - - - - - - - - - - - -
    # Enrich with GO
    # - - - - - - - - - - - - -
    print("Generate GO barplot")
    if (is.null(organism)) {
        stop("No organism passed. Needs to be either mmu or hsa.")
    } else if (organism == "mmu") {
        go_enrich <- clusterProfiler::enrichGO(gene = gene_list,
                              OrgDb = "org.Mm.eg.db",
                              ont = "BP",
                              pvalueCutoff = 0.05)
    } else if (organism == "hsa") {
        go_enrich <- clusterProfiler::enrichGO(gene = gene_list,
                              OrgDb = "org.Hs.eg.db",
                              ont = "BP",
                              pvalueCutoff = 0.05)
    }

    # Get table of results
    go_table <- as.data.frame(go_enrich) %>%
                dplyr::arrange(desc(-log10(pvalue)))
    go_table <- head(go_table, n = 20)

    # Plot results
    go_bar <- ggplot2::ggplot(go_table, 
                              ggplot2::aes(x = reorder(Description, 
                                                       -log10(pvalue)), 
                                           y = -log10(pvalue))) +
      ggplot2::geom_bar(stat = 'identity', fill = "#157296") +
      ggplot2::geom_col(width = 0.7) +
      ggplot2::labs(title = "GO Biological Process", x = "Go Terms") +
      ggplot2::coord_flip()

    go_full <- as.data.frame(go_enrich) %>% 
      dplyr::arrange(desc(-log10(pvalue)))
    
    write.table(go_full, 
                file.path(path, paste0(go_file, "__go_full_table.txt")), 
                sep = "\t", 
                quote = F, 
                row.names = F)
    
    ggplot2::ggsave(file.path(path, paste0(go_file, "__go_bar.png")), 
           plot = go_bar, 
           device = "png")
    
    print("Saved GO full table and bar plot")
}

JAS.results <- function(result_df = NULL, annotation_df = NULL) {
  final_data_frame <- result_df %>%
    dplyr::left_join(annotation_df, by="target_id") %>%
    dplyr::arrange(padj) %>%
    dplyr::select(c("target_id", "ext_gene", "entrez_id",
                    "pvalue", "padj", "log2FoldChange", 
                    "lfcSE", "stat", "baseMean"))
  
  return(final_data_frame)
}

calculateEnricment <- function(x, y){
    xs <- base::eval(
        base::parse(text = x)
    )
    
    ys <- base::eval(
        base::parse(text = y)
    )
    
    return(xs / ys)
}

generateBubble <- function(data = NULL, save_path = NULL){
    data$foldEnrichmentCalc <- base::mapply(calculateEnricment, 
                                            data$GeneRatio,
                                            data$BgRatio)
    
    dh <- arrange(data, desc(pvalue))
    dh <- head(dh, n=10)

    stopifnot(
        all(c("ID", "FoldEnrichmentCalc", "Count", "pvalue") %in% colnames(data))
        )
    
    bubble_enrich <- ggplot2::ggplot(dh, 
                            mapping = ggplot2::aes(x=foldEnrichmentCalc, 
                                          y=ID, 
                                          size=Count, 
                                          color=-log10(pvalue))) +
        ggplot2::geom_point(alpha=1) + 
        ggplot2::scale_x_continuous(limits = c(0, (base::max(dh$foldEnrichmentCalc) + 5))) + 
        ggplot2::scale_size(range = c(base::min(dh$Count), 
                                      base::max(dh$Count)),
                   name = "Gene counts") + 
        ggplot2::scale_colour_gradient(low = "white",
                              high = "red") + 
        ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
              panel.grid.minor = ggplot2::element_blank(),
              panel.background = ggplot2::element_rect(fill = "#D3D3D3"),
              axis.line = ggplot2::element_line(colour = "black")
        ) + 
        ggplot2::ylab("Biological process (Gene Ontology)") + 
        ggplot2::xlab("Fold enrichment (GeneRatio / BgRatio)")
    
    ggplot2::ggsave(save_path,
                    plot = bubble_enrich,
                    device = "png")
}

# Main paths
arg <- create.args()

`%notin%` <- Negate(`%in%`)

# # # # # # # # # # # # # # # #
# Check if results dir exists #
# # # # # # # # # # # # # # # #
if (!dir.exists(arg$results)) {
    print(paste(arg$results, "folder for results does not exists."))
    quit()
}

# # # # # # # # # # # # # # #
# Check supported organisms #
# # # # # # # # # # # # # # #
if (arg$organism %notin% c("mmu", "hsa")) { stop("-o/--organism not supported, please use only for hsa or mmu.") }

# # # # # # # # # # # # # # # # # # # # # # # #
# Check if metadata is a google spreadsheets  #
# # # # # # # # # # # # # # # # # # # # # # # #
if (substr(arg$file, nchar(arg$file) - 4 + 1, nchar(arg$file)) != ".txt") {
    stop("Metadata has to be a metadata.txt file. Take a look at the example in the Github repo.")
}

# # # # # # # # # # # # # # # # #
# Remove / from path if needed  #
# # # # # # # # # # # # # # # # #
if (substr(arg$path, nchar(arg$path), nchar(arg$path)) == "/") {
    arg$path <- substr(arg$path, 1, nchar(arg$path) - 1)
} else { }

metadata_df <- start.metadata(file = arg$file, separator=arg$separator)

import_matrix <- import.kallisto.tx(base_path = arg$path, metadata = metadata_df)

deseq_matrix <- DESeq2::DESeq(import_matrix)

SummEx_df <- as.data.frame(
    SummarizedExperiment::assay(
        DESeq2::vst(deseq_matrix)
        )
    )
SummEx_df$Gene <- (
    stringr::str_split_fixed(
        row.names(SummEx_df),
        stringr::fixed("."), 2))[,1]

column_select_names <- list(import_matrix@colData[,c("condition")])[[1]]

combinations <- combn(unique(metadata_df$treatment), 2)

annotation <- get.annotation(organism = arg$organism)

for(col in 1:ncol(combinations)) {
    groups <- combinations[, col]
    
    print(paste0("Pair-wise comparison using Wald test: ", 
                 groups[1], " vs ", groups[2], 
                 " running with BH (or false discovery rate, FDR) for p-value adjusted method."))

    general_results <- as.data.frame(DESeq2::results(deseq_matrix, 
                                     contrast=c("condition", 
                                                groups[1], 
                                                groups[2])
                                     ))
    general_results$target_id <- (
      stringr::str_split_fixed(
        row.names(general_results),
        stringr::fixed("."), 2))[,1]
    
    general_annotated <- JAS.results(result_df = general_results,
                                     annotation_df = annotation)
    
    results_filtered <- as.data.frame(
                      base::subset(general_results, 
                                    padj <= 0.05 & 
                                      (log2FoldChange <= -1 | 
                                         log2FoldChange >= 1)
                                    ))
    
    filtered_annotated <- JAS.results(result_df = results_filtered,
                                      annotation_df = annotation)
    
    results_filtered_01 <- as.data.frame(
        base::subset(general_results,
                     padj <= 0.1 &
                         (log2FoldChange <= -1 |
                              log2FoldChange >= 1)
                     )
    )
    
    filtered_annotated_01 <- JAS.results(result_df = results_filtered_01,
                                         annotation_df = annotation)

    FOLDER_GROUP <- paste0(groups[1], "_", groups[2])
    dir.create(file.path(arg$results, FOLDER_GROUP))
    
    write.table(general_annotated,
                file=file.path(arg$results,
                               FOLDER_GROUP,
                               "general_no_filter.tsv"),
                sep="\t",
                quote = F,
                row.names = F)
    
    write.table(filtered_annotated, 
                file=file.path(arg$results, 
                               FOLDER_GROUP, 
                               "ann_de_filt_padj005_logfc1.tsv"), 
                sep="\t",
                quote = F,
                row.names = F)
    
    write.table(filtered_annotated_01,
                file=file.path(arg$results,
                               FOLDER_GROUP,
                               "ann_de_filt_padj01_logfc1.tsv"),
                sep="\t",
                quote = F,
                row.names = F)

    if (!arg$no_volcano) {
        volc_df <- dplyr::distinct(na.omit(general_annotated), 
                                   ext_gene, 
                                   .keep_all = TRUE)
        
        volcano <- run.volcano(pval = volc_df$padj,
                    lfc = volc_df$log2FoldChange,
                    names = volc_df$ext_gene)
        
        ggplot2::ggsave(file.path(arg$results, 
                                  FOLDER_GROUP, 
                                  "volcano_plot.png"),
                        plot = volcano, 
                        device = "png")
    }
    
    matrix_hmp <- na.omit(
        dplyr::inner_join(SummEx_df, 
                          filtered_annotated, 
                          by = c("Gene" = "target_id"), 
                          copy = TRUE)
        )
    matrix_hmp <- matrix_hmp %>% dplyr::arrange(padj)
    rownames(matrix_hmp) <- make.names(matrix_hmp$ext_gene, unique = T)
    matrix_hmp_sel <- matrix_hmp[,1:length(column_select_names)]
    colnames(matrix_hmp_sel) <- column_select_names
    older_columns <- sort(colnames(matrix_hmp_sel))
    matrix_hmp_sel <- matrix_hmp_sel[, order(names(matrix_hmp_sel))]
    colnames(matrix_hmp_sel) <- older_columns
    
    pheatmap::pheatmap(mat=head(matrix_hmp_sel, 50),
                       scale = "row",
                       cluster_rows = T,
                       show_rownames = T,
                       cluster_cols = F,
                       filename = file.path(arg$results,
                                            FOLDER_GROUP,
                                            "heatmap_de_genes_top_50.png"))
    
    pathway.enrichment(gene_list = na.omit(filtered_annotated$entrez_id),
                       organism = arg$organism,
                       path = file.path(arg$results, FOLDER_GROUP),
                       kegg_file = "filtered_005",
                       go_file = "filtered_005")
    
    pathway.enrichment(gene_list = na.omit(filtered_annotated_01$entrez_id),
                       organism = arg$organism,
                       path = file.path(arg$results, FOLDER_GROUP),
                       kegg_file = "filtered_01",
                       go_file = "filtered_01")
}
