# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Motivation: Automated way of finding differentially expressed genes for every match
# Author: Thomaz Guadagnini
# Github: @ThomazGR
# Email: thomaz@tutanota.de
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

requireNamespace("tximport") # import Kallisto results .h5 to DESeq2 matrix-way
requireNamespace("DESeq2") # statistical analysis by gene counts from Kallisto
requireNamespace("dplyr") # import if needed to manage deeper a dataframe
requireNamespace("stringr") # needed to deal with string manipulation
requireNamespace("biomaRt") # retrieve annotation data for gene list
requireNamespace("ggplot2") # used for plotting kegg/go results in a bar plot
requireNamespace("cowplot") # used for plotting kegg/go results in a bar plot
#requireNamespace("clusterProfiler") # enrichment for KEGG and GO from gene list


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

        --no-volcano  Does not build volcano.
    " -> doc

    opt <- docopt::docopt(doc)

    return(opt)
}

start.metadata <- function(file = NULL) {
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
       metadata <- read.table(file, sep = ";", header = T, stringsAsFactors = F)
    }
    #setwd(path.expand("~/"))

    # metadata <- dplyr::mutate(metadata,
    #                           path = file.path(args$results, "3_kallisto_results",
    #                                            Run_s, "abundance.h5"))
    # metadata <- dplyr::rename(metadata, sample = Run_s)

    # if (is.null(args$groups[[1]])) {
    #     groups <- unique(metadata$treatment)
    # } else {
    #     groups <- strsplit(args$groups, ",")[[1]]
    # }
    # volcano <- !args[["--no-volcano"]]

    #return(list(metadata = metadata, groups = groups, volcano = volcano))
    return(metadata)
}

import.kallisto.tx <- function(base_path = NULL, metadata = NULL) {
    if (!is.data.frame(metadata)) {quit("Metadata passed is not a data.frame.")}
    if (!is.character(base_path)) {quit("Base path is not a string. Check --help.")}
    needed_fields <- c("Run_s", "treatment")
    if (!all(needed_fields %in% colnames(metadata))) {quit("Metadata file needs to have both Run_s and treatment field. Check --help.")}

    files <- file.path(base_path, metadata$Run_s, "abundance.h5")
    txi.kallisto <- tximport::tximport(files, type="kallisto", txOut=TRUE)
    sampleTable <- data.frame(condition=metadata$treatment)
    rownames(sampleTable) <- colnames(txi.kallisto$counts)
    imported_matrix <- DESeq2::DESeqDataSetFromTximport(txi.kallisto, sampleTable, ~condition)

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
    
    return(ttg)
}

run.volcano <- function(pval = NULL, lfc = NULL, names = NULL) {
    if (all(list(is.null(pval), is.null(lfc), is.null(names)))) {quit("pval, lfc or names is null.")}
    data <- data.frame(pval = -log10(pval),
                       lfc = lfc,
                       row.names = names)

    # Remove any rows that have NA as an entry
    data <- na.omit(data)

    # Color the points which are up or down
    ## If fold-change > 0 and pvalue > 1.3 (Increased significant)
    ## If fold-change < 0 and pvalue > 1.3 (Decreased significant)
    data <- dplyr::mutate(data, color = case_when(data$lfc > 1 & data$pval > 1.3 ~ "Increased",
                                           data$lfc < -1 & data$pval > 1.3 ~ "Decreased",
                                           data$pval < 1.3 ~ "Nonsignificant"))
    # Make a basic ggplot2 object with x-y values
    volc <- ggplot2::ggplot(data, aes(x = lfc, y = pval, color = color))

    # Add ggplot2 layers
    volcano_plot <- volc +
            ggtitle(label = "Volcano Plot", subtitle = "Colored by fold-change direction") +
            geom_point(size = 2.5, alpha = 0.8, na.rm = T) +
            scale_color_manual(name = "Directionality",
                               values = c(Increased = "#008B00",
                                          Decreased = "#CD4F39",
                                          nonsignificant = "darkgray")) +
            theme_bw(base_size = 14) + # change overall theme
            theme(legend.position = "right") + # change the legend
            xlab("logFC") + # Change X-Axis label
            ylab(expression(-log[10]("adjusted p-value"))) + # Change Y-Axis label
            geom_hline(yintercept = 1.3, colour = "darkgrey") + # Add p-adj value cutoff line
            scale_y_continuous(trans = "log1p") # Scale yaxis due to large p-values

    return(volcano_plot)
}

pathway.enrichment <- function(gene_list = NULL, organism = NULL, path = NULL) {
    # - - - - - - - - - - - - -
    # Enrich with KEGG database
    # - - - - - - - - - - - - -
    if (is.null(gene_list)) { stop("No gene list has been passed.") }
    if (is.null(path)) { stop("No path has been passed.") }
    if (is.null(organism)) {
        stop("No organism passed. Needs to be either mmu or hsa.")
    } else if (organism == "mmu") {
        kegg_enrich <- enrichKEGG(gene = gene_list,
                                  organism = "mmu",
                                  pvalueCutoff = 0.05)
    } else if (organism == "hsa") {
        kegg_enrich <- enrichKEGG(gene = gene_list,
                                  organism = "hsa",
                                  pvalueCutoff = 0.05)
    }

    print("Generate KEGG pathway results")
    # Get table of results
    kegg_table <- head(as.data.frame(kegg_enrich), n = 10) %>%
            arrange(desc(-log10(pvalue)))

    # KEGG plot
    kegg_bar <- ggplot(kegg_table, aes(x = reorder(Description, -log10(pvalue)), y = -log10(pvalue))) +
            geom_bar(stat = 'identity', fill = "#6B2525") +
            geom_col(width = 0.7) +
            labs(title = "KEGG Enrichment Pathways", x = "Termos de KEGG") +
            coord_flip()

    kegg_full <- as.data.frame(kegg_enrich) %>% arrange(desc(-log10(pvalue)))
    write.table(kegg_full, paste0(path, "/kegg_full_table.txt"), sep = "\t", quote = F, row.names = F)
    ggsave(paste0(path, "/kegg_bar.png"), plot = kegg_bar, device = "png")
    print("Saved KEGG full table and bar plot")

    # - - - - - - - - - - - - -
    # Enrich with GO
    # - - - - - - - - - - - - -
    print("Generate GO barplot")
    if (is.null(organism)) {
        stop("No organism passed. Needs to be either mmu or hsa.")
    } else if (organism == "mmu") {
        go_enrich <- enrichGO(gene = df$entrez,
                              OrgDb = "org.Mm.eg.db",
                              ont = "BP",
                              pvalueCutoff = 0.05)
    } else if (organism == "hsa") {
        go_enrich <- enrichGO(gene = df$entrez,
                              OrgDb = "org.Hs.eg.db",
                              ont = "BP",
                              pvalueCutoff = 0.05)
    }

    # Get table of results
    go_table <- head(as.data.frame(go_enrich), n = 10) %>%
            arrange(desc(-log10(pvalue)))

    # Plot results
    go_bar <- ggplot(go_table, aes(x = reorder(Description, -log10(pvalue)), y = -log10(pvalue))) +
            geom_bar(stat = 'identity', fill = "#157296") +
            geom_col(width = 0.7) +
            labs(title = "GO Biological Pathways", x = "Termos de GO") +
            coord_flip()

    go_full <- as.data.frame(go_enrich) %>% arrange(desc(-log10(pvalue)))
    write.table(go_full, paste0(path, "/go_full_table.txt"), sep = "\t", quote = F, row.names = F)
    ggsave(paste0(path, "/go_bar.png"), plot = go_bar, device = "png")
    print("Saved GO full table and bar plot")
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

# # # # # # # # # # # # # # # #
# Check if results dir exists #
# # # # # # # # # # # # # # # #
if (!dir.exists(arg$results)) {
    print(paste(arg$results, "folder for results does not exists."))
    quit()
}

# # # # # # # # # # # # # # # # # # # # # # # #
# Check if metadata is a google spreadsheets  #
# # # # # # # # # # # # # # # # # # # # # # # #
if (!grepl("https://docs.google.com/spreadsheets", args$file, fixed = T) &
    !grepl("output=csv", args$file, fixed = T) &
    substr(arg[["--file"]], nchar(arg[["--file"]]) - 4 + 1, nchar(arg[["--file"]])) != ".txt") {
    stop("Metadata has to be a metadata.txt file. Take a look at the example in the Github repo.")
}

# # # # # # # # # # # # # # #
# Check supported organisms #
# # # # # # # # # # # # # # #
if (arg$organism %notin% c("mmu", "hsa")) { stop("-o/--organism not supported, please use only for hsa or mmu.") }

# # # # # # # # # # # # # # # # #
# Remove / from path if needed  #
# # # # # # # # # # # # # # # # #
if (substr(arg$path, nchar(arg$path), nchar(arg$path)) == "/") {
    arg$path <- substr(arg$path, 1, nchar(arg$path) - 1)
    arg[["--path"]] <- substr(arg$path, 1, nchar(arg$path) - 1)
} else { }

metadata_df <- start.metadata(file = arg$file)

import_matrix <- import.kallisto.tx(base_path = arg$path, metadata = metadata_df)

deseq_matrix <- DESeq2::DESeq(import_matrix)

combinations <- combn(unique(metadata_df$treatment), 2)

annotation <- get.annotation(organism = arg$organism)

for(col in 1:ncol(combinations)) {
    groups <- combinations[, col]
    print(paste0("Pair-wise comparison using Wald test: ", 
                 groups[1], " vs ", groups[2], 
                 " running with BH (or false discovery rate, FDR) for p-value adjusted method."))

    final_results <- DESeq2::results(deseq_matrix, contrast=c("condition", groups[1], groups[2]))
    ss_padj_logfc <- base::subset(final_results, padj <= 0.05 & (log2FoldChange <= -1 | log2FoldChange >= 1))
    filt_mres_df <- as.data.frame(ss_padj_logfc)
    filt_mres_df$target_id <- (stringr::str_split_fixed(row.names(filt_mres_df), fixed("."), 2))[,1]
    
    gene_name <- filt_mres_df %>%
        dplyr::inner_join(annotation, by="target_id") %>%
        dplyr::arrange(padj) %>%
        dplyr::select(c("target_id", "ext_gene", "pvalue", 
                        "padj", "log2FoldChange", 
                        "lfcSE", "stat", "baseMean"))
    print(head(gene_name))

    FOLDER_GROUP <- paste0(groups[1], "_", groups[2])
    dir.create(file.path(arg$results, FOLDER_GROUP))
    write.table(gene_name, 
                file=file.path(arg$results, FOLDER_GROUP, "annotated_DE_genes.tsv"), 
                sep="\t",
                quote = F,
                row.names = F)

    if (!arg$no_volcano) {
        volcano <- run.volcano(pval = gene_name$padj,
                    lfc = gene_name$log2FoldChange,
                    names = gene_name$ext_gene)
        ggplot2::ggsave(file.path(arg$results, FOLDER_GROUP, "volcano_plot.png"), plot = volcano, device = "png")
    }    

}