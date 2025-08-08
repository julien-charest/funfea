###########################################################################
# FunFEA: A Tool for Enrichment Analysis of Fungal Functional Annotations #
#                                                                         #
# Written by Julien Charest & Paul Loebenstein                            #
#                                                                         #
# Version 1.2.0                                                           #
# 23.04.2025                                                              #
###########################################################################

#############
# Utilities #
#############

###########################
# Load GTF/GFF Annotation #
###########################
#' Load GTF/GFF Annotation
#'
#' @description
#' This function loads a GTF/GFF/GFF3 annotation into a dataframe.
#'
#' @param path Path to GTF/GFF/GFF3 annotation
#' @return A dataframe with GTF/GFF/GFF3 annotation information
#' @examples
#' gtf_df <- load_gtf_annotation(path/to/gtf/annotation)
#' @export

load_gtf_annotation <- function(path){

  all_lines <- readLines(path)
  gtf_lines <- all_lines[!grepl("^#", all_lines)]
  gtf_data <- read.delim(textConnection(gtf_lines), header = FALSE, sep = "\t", quote = "")
  colnames(gtf_data) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
  rownames(gtf_data) <- NULL

  return(gtf_data)
}

###########################################################
# Create Transcript ID to Protein ID Conversion Dataframe #
###########################################################
#' Create Transcript ID to Protein ID Conversion Dataframe
#'
#' @description
#' This function parses a GTF/GFF/GFF3 annotation dataframe to generate a transcript ID to protein ID conversion dataframe.
#'
#' @param gtf_df A dataframe with GTF/GFF/GFF3 annotation information (see load_gtf_annotation)
#'
#' @return A dataframe with gene ID, gene name, transcript ID and protein ID information
#' @examples
#' transcript2protein_id_df <- create_transcript2protein_id_df(gtf_df)
#' @export

create_transcript2protein_id_df <- function(gtf_df) {

  print("Parsing Gene Annotation for Gene IDs, Gene Names, Transcript IDs, Protein IDs & Products...")

  gene_ids <- c()
  gene_names <- c()
  transcript_ids <- c()
  protein_ids <- c()
  products <- c()

  if ("gene" %in% gtf_df$feature){

    cds_df <- gtf_df[gtf_df$feature == "CDS",]
    total <- nrow(cds_df)
    pb <- txtProgressBar(min = 0, max = total, style = 3)
    progress <- 0

    for (i in 1:nrow(cds_df)){

      setTxtProgressBar(pb, progress)
      attribute_i <- trimws(strsplit(cds_df[i, ]$attribute, ";")[[1]])

      gene_id <- NA
      gene_name <- NA
      transcript_id <- NA
      protein_id <- NA
      product <- NA

      for (j in attribute_i){

        # Gene ID
        if (grepl("^locus_tag=", j)){
          gene_id <- substr(j,11,nchar(j))
        }
        # Gene Name
        if (grepl("^gene=", j)){
          gene_name <- substr(j,6,nchar(j))
        }
        # Transcript ID
        if (grepl("^Parent=", j)){
          transcript_id <- substr(j,8,nchar(j))
          transcriptidfound <- TRUE
        }
        # Protein ID
        if (grepl("^protein_id=", j)){
          protein_id <- substr(j,12,nchar(j))
        }
        # Product
        if (grepl("^product=", j)){
          product <- substr(j,9,nchar(j))
        }
      }

      gene_ids <- c(gene_ids, gene_id)
      gene_names <- c(gene_names, gene_name)
      transcript_ids <- c(transcript_ids, transcript_id)
      protein_ids <- c(protein_ids, protein_id)
      products <- c(products, product)

      progress <- progress + 1

    }

    close(pb)

    if (all(is.na(protein_ids)) | all(is.na(transcript_ids))){

      print("Parsing for Transcript IDs...")

      gene_ids <- c()
      gene_names <- c()
      transcript_ids <- c()
      protein_ids <- c()
      products <- c()

      mrna_df <- gtf_df[gtf_df$feature == "mRNA",]
      total <- nrow(mrna_df)
      pb <- txtProgressBar(min = 0, max = total, style = 3)
      progress <- 0

      for (i in 1:nrow(mrna_df)){

        setTxtProgressBar(pb, progress)
        attribute_i <- trimws(strsplit(mrna_df[i, ]$attribute, ";")[[1]])

        gene_id <- NA
        gene_name <- NA
        transcript_id <- NA
        protein_id <- NA
        product <- NA

        for (j in attribute_i){

          # Gene ID
          if (grepl("^Parent=", j)){
            gene_id <- substr(j,8,nchar(j))
          }
          # Gene Name
          if (grepl("^gene=", j)){
            gene_name <- substr(j,6,nchar(j))
          }
          # Transcript ID
          if (grepl("^ID=", j)){
            transcript_id <- substr(j,4,nchar(j))
            transcriptidfound <- TRUE
          }
          # Protein ID
          if (grepl("^protein_id=", j)){
            protein_id <- substr(j,12,nchar(j))
          }
          else if (grepl("^ID=", j)){
            protein_id <- substr(j,4,nchar(j))
          }
          # Product
          if (grepl("^product=", j)){
            product <- substr(j,9,nchar(j))
          }
        }

        gene_ids <- c(gene_ids, gene_id)
        gene_names <- c(gene_names, gene_name)
        transcript_ids <- c(transcript_ids, transcript_id)
        protein_ids <- c(protein_ids, protein_id)
        products <- c(products, product)

        progress <- progress + 1

      }

      close(pb)

    }

    transcript2protein_id_df <- data.frame(gene_id = gene_ids,
                                           gene_name = gene_names,
                                           transcript_id = transcript_ids,
                                           protein_id = protein_ids,
                                           product = products,
                                           stringsAsFactors = FALSE)

    transcript2protein_id_df <- unique(transcript2protein_id_df)
    rownames(transcript2protein_id_df) <- NULL
    close(pb)

    if (all(is.na(transcript2protein_id_df$gene_name))){

      print("Parsing for Gene Names...")

      gene_names <- c()
      gene_df <- gtf_df[gtf_df$feature == "gene",]
      total <- nrow(gene_df)
      pb <- txtProgressBar(min = 0, max = total, style = 3)
      progress <- 0

      for (k in 1:nrow(transcript2protein_id_df)){

        setTxtProgressBar(pb, progress)
        gene_name <- NA
        attribute_k <- trimws(strsplit(gene_df[grepl(transcript2protein_id_df[k,]$gene_id, gene_df$attribute), ]$attribute, ";")[[1]])

        for (j in attribute_k){
          # Gene Name
          if (grepl("^Name=", j)){
            gene_name <- substr(j,6,nchar(j))
          }
        }

        gene_names <- c(gene_names, gene_name)
        progress <- progress + 1

      }

      transcript2protein_id_df$gene_name <- gene_names
      close(pb)

    }

  }

  else {

    cds_df <- gtf_df[gtf_df$feature == "CDS",]
    exon_df <- gtf_df[gtf_df$feature == "exon",]
    total <- nrow(cds_df) + nrow(exon_df)
    pb <- txtProgressBar(min = 0, max = total, style = 3)
    progress <- 0

    gene_ids <- c()
    gene_names <- c()
    transcript_ids <- c()
    protein_ids <- c()
    products <- c()

    for (i in 1:nrow(cds_df)){

      setTxtProgressBar(pb, progress)
      attribute_i <- trimws(strsplit(cds_df[i, ]$attribute, ";")[[1]])
      protein_id <- NA
      gene_id <- NA

      for (j in attribute_i){

        # Protein ID
        if (grepl("^protein_id ", j)){
          protein_id <- substr(j,12,nchar(j))
        }
        else if (grepl("^proteinId ", j)){
          protein_id <- substr(j,11,nchar(j))
        }
        else if (grepl("^proteinId=", j)){
          protein_id <- substr(j,11,nchar(j))
        }
        # Gene Name
        if (grepl("^name ", j)){
          gene_id <- gsub('"', '', substr(j,6,nchar(j)))
        }
      }

      protein_ids <- c(protein_ids, protein_id)
      gene_ids <- c(gene_ids, gene_id)
      progress <- progress + 1

    }

    transcript2protein_id_df <- data.frame(gene_id = gene_ids,
                                           gene_name = rep(NA, length(protein_ids)),
                                           transcript_id = rep(NA, length(protein_ids)),
                                           protein_id = protein_ids,
                                           product = rep(NA, length(protein_ids)),
                                           stringsAsFactors = FALSE)

    transcript2protein_id_df <- unique(transcript2protein_id_df)
    rownames(transcript2protein_id_df) <- NULL

    for (k in 1:nrow(transcript2protein_id_df)){

      setTxtProgressBar(pb, progress)
      gene_id <- transcript2protein_id_df[k,]$gene_id
      gene_name <- NA
      transcript_id <- NA
      product <- NA

      attribute_k <- trimws(strsplit(exon_df[grepl(gene_id, exon_df$attribute), ]$attribute, ";")[[1]])

      for (j in attribute_k){

        # Transcript ID
        if (grepl("^transcriptId ", j)){
          transcript_id <- substr(j,14,nchar(j))
          break
        }
      }

      gene_names <- c(gene_names, gene_name)
      transcript_ids <- c(transcript_ids, transcript_id)
      products <- c(products, product)
      progress <- progress + 1

    }

    progress <- total
    setTxtProgressBar(pb, progress)
    close(pb)

    transcript2protein_id_df$gene_name <- gene_names
    transcript2protein_id_df$transcript_id <- transcript_ids
    transcript2protein_id_df$product <- products

  }
  return(transcript2protein_id_df)
}

#####################################################
# Convert a Vector of Transcript IDs to Protein IDs #
#####################################################
#' Convert a Vector of Transcript IDs to Protein IDs
#'
#' @description
#' This function converts a vector of transcript IDs to protein IDs, using a transcript2protein_id dataframe constructed from GTF/GFF/GFF3 annotation.
#'
#' @param transcript_ids A vector of transcript IDs
#' @param transcript2protein_id_df A dataframe with gene ID, gene name, transcript ID and protein ID information (see create_transcript2protein_id_df)
#'
#' @return A vector of protein IDs
#'
#' @examples
#' protein_ids <- transcript2protein_id(transcript_ids, transcript2protein_id_df)
#' @export

transcript2protein_id <- function(transcript_ids, transcript2protein_id_df){

  mapping_vector <- setNames(transcript2protein_id_df$protein_ids, transcript2protein_id_df$transcript_ids)

  return(mapping_vector[transcript_ids])
}


###############################################
# Convert a Vector of Gene IDs to Protein IDs #
###############################################
#' Convert a Vector of Gene IDs to Protein IDs
#'
#' @description
#' This function converts a vector of gene IDs and/or gene names to protein IDs, using a transcript2protein_id dataframe constructed from GTF/GFF/GFF3 annotation.
#'
#' @param gene_ids A vector of gene IDs and/or gene names
#' @param transcript2protein_id_df A dataframe with gene ID, gene name, transcript ID and protein ID information (see create_transcript2protein_id_df)
#'
#' @return A vector of protein IDs
#'
#' @examples
#' protein_ids <- gene2protein_id(gene_ids, transcript2protein_id_df)
#' @export

gene2protein_id <- function(gene_ids, transcript2protein_id_df) {
  genes_upper <- toupper(gene_ids)

  # Make a copy of the dataframe with uppercase keys for case-insensitive matching
  df <- transcript2protein_id_df
  df$gene_ids <- toupper(df$gene_ids)
  df$gene_names <- toupper(df$gene_names)

  # Create lookup tables
  gene_id_map <- setNames(df$protein_ids, df$gene_ids)
  gene_name_map <- setNames(df$protein_ids, df$gene_names)

  # Look up each gene in the maps, preferring gene_ids over gene_names
  protein_ids <- sapply(genes_upper, function(gene) {
    if (!is.na(gene_id_map[gene])) return(gene_id_map[gene])
    if (!is.na(gene_name_map[gene])) return(gene_name_map[gene])
    return(NA)
  })

  return(protein_ids)
}


###############################
# COG/KOG Enrichment Analysis #
###############################

###########################
# Load COG/KOG Annotation #
###########################
#' Load COG/KOG Annotation
#'
#' @description
#' This function loads a COG/KOG (Clusters of Orthologous Genes) annotation file into a dataframe.
#'
#' @param path Path to COG/KOG annotation
#'
#' @return A dataframe with COG/KOG annotation information
#'
#' @examples
#' kog_annotation_df <- load_kog_annotation(path/to/kog/annotation)
#' @export

load_kog_annotation <- function(path){

  annotation_data <- read.delim(path, header = TRUE, sep = "\t")

  if (dim(annotation_data)[2] != 6) {
    stop("KOG annotation file must have six columns: transcriptId,	proteinId,	kogid,	kogdefline,	kogClass & kogGroup")
  }

  colnames(annotation_data) <- c("transcriptId",	"proteinId",	"kogid",	"kogdefline",	"kogClass",	"kogGroup")

  annotation_data$transcriptId <- trimws(annotation_data$transcriptId)
  annotation_data$proteinId <- trimws(annotation_data$proteinId)
  annotation_data$kogid <- trimws(annotation_data$kogid)
  annotation_data$kogdefline <- trimws(annotation_data$kogdefline)
  annotation_data$kogClass <- trimws(annotation_data$kogClass)
  annotation_data$kogGroup <- trimws(annotation_data$kogGroup)

  return(annotation_data)
}

########################
# Create COG/KOG Model #
########################
#' Create COG/KOG Model
#'
#' @description
#' This function generates a background frequency model for COG/KOG categories enrichment analysis from a COG/KOG annotation.
#'
#' @param kog_annotation A dataframe with COG/KOG annotation information (see load_kog_annotation)
#'
#' @return A dataframe with COG/KOG categories background frequencies
#'
#' @examples
#' kog_model_df <- create_kog_model(kog_annotation_df)
#' @export

create_kog_model <- function(kog_annotation){

  data(cog_classes, envir = environment())

  groups <- c("CELLULAR PROCESSES AND SIGNALING",
              "Cell wall/membrane/envelope biogenesis",
              "Cell motility",
              "Posttranslational modification, protein turnover, chaperones",
              "Signal transduction mechanisms",
              "Intracellular trafficking, secretion, and vesicular transport",
              "Defense mechanisms",
              "Extracellular structures",
              "Nuclear structure",
              "Cytoskeleton",
              "INFORMATION STORAGE AND PROCESSING",
              "RNA processing and modification",
              "Chromatin structure and dynamics",
              "Translation, ribosomal structure and biogenesis",
              "Transcription",
              "Replication, recombination and repair",
              "METABOLISM",
              "Energy production and conversion",
              "Cell cycle control, cell division, chromosome partitioning",
              "Amino acid transport and metabolism",
              "Nucleotide transport and metabolism",
              "Carbohydrate transport and metabolism",
              "Coenzyme transport and metabolism",
              "Lipid transport and metabolism",
              "Inorganic ion transport and metabolism",
              "Secondary metabolites biosynthesis, transport and catabolism",
              "POORLY CHARACTERIZED",
              "General function prediction only",
              "Function unknown",
              "ANNOTATED",
              "UNANNOTATED")

  kog_model <- data.frame(kogClass = groups)
  kog_model$kogModel <- 0
  kog_model$proteinIds <- NA

  for (i in cog_classes_df$class){
    i_ids <- kog_annotation[kog_annotation$kogClass == i,]$proteinId

    if (length(i_ids) >0){
      kog_model[kog_model$kogClass == i,]$proteinIds <- paste(sort(unique(kog_annotation[kog_annotation$kogClass == i,]$proteinId)), collapse = ", ")
    }
  }

  kog_model[kog_model$kogClass == "CELLULAR PROCESSES AND SIGNALING",]$proteinIds <- paste(sort(unique(unlist(strsplit(paste(kog_model[2:10,]$proteinIds[!is.na(kog_model[2:10,]$proteinIds)], collapse = ", "), split = ", ")))), collapse = ", ")
  kog_model[kog_model$kogClass =="INFORMATION STORAGE AND PROCESSING",]$proteinIds <- paste(sort(unique(unlist(strsplit(paste(kog_model[12:16,]$proteinIds[!is.na(kog_model[12:16,]$proteinIds)], collapse = ", "), split = ", ")))), collapse = ", ")
  kog_model[kog_model$kogClass =="METABOLISM",]$proteinIds <- paste(sort(unique(unlist(strsplit(paste(kog_model[18:26,]$proteinIds[!is.na(kog_model[18:26,]$proteinIds)], collapse = ", "), split = ", ")))), collapse = ", ")
  kog_model[kog_model$kogClass =="POORLY CHARACTERIZED",]$proteinIds <- paste(sort(unique(unlist(strsplit(paste(kog_model[28:29,]$proteinIds[!is.na(kog_model[28:29,]$proteinIds)], collapse = ", "), split = ", ")))), collapse = ", ")
  kog_model[kog_model$kogClass =="ANNOTATED",]$proteinIds <- paste(sort(unique(unlist(strsplit(paste(kog_model[1:29,]$proteinIds[!is.na(kog_model[1:29,]$proteinIds)], collapse = ", "), split = ", ")))), collapse = ", ")

  for (i in kog_model$kogClass){
    if (length(na.omit(kog_model[kog_model$kogClass == i,]$proteinIds)) > 0){
      kog_model[kog_model$kogClass == i,]$kogModel <- length(unique(unlist(strsplit(kog_model[kog_model$kogClass == i,]$proteinIds, split = ", "))))
    }
  }

  kog_model <- kog_model[, c("kogClass", "kogModel", "proteinIds")]
  colnames(kog_model) <- c("kogClass", "kogModel", "proteinIdsModel")

  return(kog_model)
}

################################
# Calculate COG/KOG Enrichment #
################################
#' Calculate COG/KOG Enrichment
#'
#' @description
#' This function generates a dataframe containing background/sample frequencies and enrichment statistics per COG/KOG category from a COG/KOG model and a vector of protein IDs. Protein IDs must match the Protein IDs from the annotation used to build the model.
#'
#' @param kog_model A dataframe with COG/KOG categories background frequencies (see create_kog_model)
#' @param protein_ids A vector of protein IDs. Must match the Protein IDs from the annotation used to build the model.
#' @param test A statistical method for testing the null of independence of rows and columns in a contingency table (default = "fisher"; options = c("fisher", "chisq))
#' @param p.adjust.method A statistical method to adjust p-values for multiple comparisons (default = "BH"; options = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")(see p.adjust.methods))
#'
#' @return A dataframe with COG/KOG categories enrichment statistics
#'
#' @examples
#' kog_enrichment_df <- kog_enrichment(kog_model_df, protein_ids)
#' kog_enrichment_df <- kog_enrichment(kog_model_df, protein_ids, test = "fisher", p.adjust.method = "BH")
#' kog_enrichment_df <- kog_enrichment(kog_model_df, protein_ids, test = "chisq", p.adjust.method = "bonferroni")
#' @export

kog_enrichment <- function(kog_model, protein_ids, test= "fisher", p.adjust.method= "BH"){

  kog_count_df <- kog_model
  kog_count_df$proteinCount <- 0
  kog_count_df$proteinIds <- NA
  protein_ids <- sort(unique(protein_ids))

  for (i in kog_count_df$kogClass){

    i_proteins <- sort(unlist(strsplit(kog_count_df[kog_count_df$kogClass == i,]$proteinIdsModel, split = ", ")))
    intersect_protein_ids <- intersect(protein_ids, i_proteins)

    if (length(intersect_protein_ids) > 0){
      kog_count_df[kog_count_df$kogClass == i,]$proteinIds <- paste(sort(unique(intersect_protein_ids)), collapse = ", ")
      kog_count_df[kog_count_df$kogClass == i,]$proteinCount <- length(unique(intersect_protein_ids))
    }
  }

  unannotated_proteins <- sort(setdiff(protein_ids, (unique(unlist(strsplit(paste(kog_count_df[kog_count_df$kogClass == "ANNOTATED",]$proteinIds, collapse = ", "), split = ", "))))))
  kog_count_df[kog_count_df$kogClass =="UNANNOTATED",]$proteinCount <- length(unannotated_proteins)
  kog_count_df[kog_count_df$kogClass =="UNANNOTATED",]$proteinIds <- paste(unannotated_proteins, collapse = ", ")

  enrichment_df <- kog_count_df
  enrichment_df$enrichment <- (enrichment_df$proteinCount/enrichment_df[enrichment_df$kogClass == "ANNOTATED",]$proteinCount)/(enrichment_df$kogModel/enrichment_df[enrichment_df$kogClass == "ANNOTATED",]$kogModel)
  enrichment_df$enrichment[is.infinite(enrichment_df$enrichment)] <- NaN
  enrichment_df$enrichment[enrichment_df$kogClass == "UNANNOTATED"] <- NaN
  enrichment_df$enrichment[enrichment_df$kogClass == "ANNOTATED"] <- NaN

  stat_results <- c()

  for (i in enrichment_df$kogClass){
    if (is.nan(enrichment_df[enrichment_df$kogClass == i,]$enrichment)){
      stat_results <- append(stat_results, NaN)
      next
    }

    contingency_table <- matrix(c(enrichment_df[enrichment_df$kogClass == i,]$proteinCount,
                                  enrichment_df[enrichment_df$kogClass == "ANNOTATED",]$proteinCount - enrichment_df[enrichment_df$kogClass == i,]$proteinCount,
                                  enrichment_df[enrichment_df$kogClass == i,]$kogModel,
                                  enrichment_df[enrichment_df$kogClass == "ANNOTATED",]$kogModel - enrichment_df[enrichment_df$kogClass == i,]$kogModel),2,2)
    if (test == "fisher"){
      stat_results <- append(stat_results, fisher.test(contingency_table, alternative='greater')$p.value)
    }

    else if (test == "chisq"){
      stat_results <- append(stat_results, chisq.test(contingency_table)$p.value)
    }
  }

  enrichment_df$p.value <- stat_results
  enrichment_df$adj.p.value <- p.adjust(enrichment_df$p.value, method = p.adjust.method)

  enrichment_df <- enrichment_df[, c("kogClass", "kogModel", "proteinCount", "enrichment", "p.value", "adj.p.value", "proteinIds")]

  return(enrichment_df)
}

####################################
# Generate COG/KOG Enrichment Plot #
####################################
#' Generate COG/KOG Enrichment Plot
#'
#' @description
#' This function generates a bar plot for visualization of COG/KOG categories enrichment statistics.
#'
#' @param enrichment_df A dataframe with COG/KOG enrichment statistics (see kog_enrichment)
#' @param significant A color assignment for statistically significant COG/KOG categories (default = "#0B775E")
#' @param plot_title A plot title (default = NA)
#' @param plot_type A plot type option (default = "bar"; options = c("bar", "lollipop"))
#'
#' @return A ggplot2 bar plot object
#'
#' @examples
#' generate_kog_plot(kog_enrichment_df)
#' generate_kog_plot(kog_enrichment_df, plot_type = "lollipop")
#' generate_kog_plot(kog_enrichment_df, plot_title = "KOG Enrichment", significant = "blue")
#' @export

generate_kog_plot <- function(enrichment_df, significant = "#0B775E", plot_title = NA, plot_type = "bar"){

  fig_tab <- enrichment_df[1:29,]

  fig_tab$enrichment[is.nan(fig_tab$enrichment)] <- 0
  fig_tab$p.value[is.nan(fig_tab$p.value)] <- 1
  fig_tab$adj.p.value[is.nan(fig_tab$adj.p.value)] <- 1

  fig_tab$color <- "grey"

  for (i in 1:nrow(fig_tab)){
    if (fig_tab[i,]$adj.p.value <= 0.05 & fig_tab[i,]$enrichment > 1){
      fig_tab[i,]$color <- significant
    }
  }

  fig_tab$kogClass <- factor(fig_tab$kogClass, levels = rev(fig_tab$kogClass))

  if (plot_type == "bar"){

    p <- ggplot2::ggplot(data = fig_tab, ggplot2::aes(x = kogClass, y = enrichment)) + ggplot2::geom_bar(stat = "identity", width = 0.75, fill = fig_tab$color) +
      ggplot2::coord_flip() + ggplot2::xlab("") + ggplot2::ylab("Enrichment") + ggplot2::ylim(0, (max(fig_tab$enrichment, na.rm = TRUE) + 0.03)) +
      ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "#B40F20") + ggplot2::theme_minimal()

    if (!is.na(plot_title)){
      p <- p + ggplot2::ggtitle(plot_title)
    }

    for (i in fig_tab$kogClass){
      if (fig_tab$adj.p.value[fig_tab$kogClass == i] <= 0.05 & fig_tab$enrichment[fig_tab$kogClass == i] > 1){
        if (fig_tab$adj.p.value[fig_tab$kogClass == i] <= 0.001){
          p <- p + ggplot2::geom_text(x=dim(fig_tab)[1]-which(fig_tab$kogClass == i) + 1,
                            y=fig_tab$enrichment[fig_tab$kogClass == i] + 0.05,
                            label="***", angle = 90)
        }

        else if (fig_tab$adj.p.value[fig_tab$kogClass == i] <= 0.01){
          p <- p + ggplot2::geom_text(x=dim(fig_tab)[1]-which(fig_tab$kogClass == i) + 1,
                            y=fig_tab$enrichment[fig_tab$kogClass == i] + 0.05,
                            label="**", angle = 90)
        }

        else {
          p <- p + ggplot2::geom_text(x=dim(fig_tab)[1]-which(fig_tab$kogClass == i) + 1,
                           y=fig_tab$enrichment[fig_tab$kogClass == i] + 0.05,
                           label="*" , , angle = 90)
        }

      }
    }
  }

  else if (plot_type == "lollipop"){

    p <- ggplot2::ggplot(data = fig_tab, ggplot2::aes(x = kogClass, y = enrichment)) +
      ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "#B40F20") +
      ggplot2::geom_segment(ggplot2::aes(xend = kogClass, yend = 0), color = "gray", size = 0.5) +
      ggplot2::geom_point(ggplot2::aes(size = proteinCount, color = -log10(adj.p.value)), alpha = 1) +
      ggplot2::coord_flip() +
      ggplot2::scale_color_gradient2(low = "gray99", mid = "gray85", high = significant,
                                     midpoint = 0.001, limits = c(0, 5), oob = scales::squish, name = "-log10(adj.p.value)") +
      ggplot2::scale_size(range = c(3, 8), name = "Protein Count") +
      ggplot2::labs(title = "", x = "", y = "Enrichment") +
      ggplot2::theme_minimal()

    if (!is.na(plot_title)){
      p <- p + ggplot2::ggtitle(plot_title)
    }
  }

  return(p)
}


##########################
# GO Enrichment Analysis #
##########################

######################
# Load GO Annotation #
######################
#' Load GO Annotation
#'
#' @description
#' This function loads a GO (Gene Ontology) annotation into a dataframe.
#'
#' @param path Path to GO annotation
#'
#' @return A dataframe with GO annotation information
#'
#' @examples
#' go_annotation_df <- load_go_annotation(path/to/go/annotation)
#' @export

load_go_annotation <- function(path){

  annotation_data <- read.delim(path, header = TRUE, sep = "\t")

  if (dim(annotation_data)[2] != 5) {
    stop("GO annotation file must have five columns: proteinId, gotermId, goName, gotermType & goAcc")
  }

  colnames(annotation_data) <- c("proteinId", "gotermId", "goName", "gotermType", "goAcc")

  return(annotation_data)
}

###################
# Create GO Model #
###################
#' Create GO Model
#'
#' @description
#' This function generates a background frequency model for GO term enrichment analysis from a GO annotation.
#'
#' @param go_annotation A dataframe with GO annotation information (see load_go_annotation)
#'
#' @return A list of dataframes with GO terms background frequencies by gotermType (biological_processes, cellular_component, molecular_function)
#'
#' @examples
#' go_model_df <- create_go_model(go_annotation_df)
#' go_model_biological_processes_df <- go_model_df$biological_processes
#' go_model_cellular_component_df <- go_model_df$cellular_component
#' go_model_molecular_function_df <- go_model_df$molecular_function
#' @export

create_go_model <- function(go_annotation){

  data(go_terms, envir = environment())
  go_model <- go_terms_df
  go_model <- go_model[, c("goAcc", "alt_id", "is_a", "goName", "gotermType")]
  go_model$goModel <- 0
  go_model$proteinIds <- NA

  go_annotated <- unique(go_annotation$proteinId)

  go_graph_df <- go_terms_df[,c("goAcc", "is_a")]
  go_graph_df$is_a <- lapply(go_graph_df$is_a, function(x) ifelse(x == "", NA, x))
  go_graph_df <- tidyr::separate_rows(go_graph_df, is_a, sep = ",")
  go_graph_df <- as.data.frame(go_graph_df)
  colnames(go_graph_df) <- c("accession", "parent")
  go_graph_df$protein_id <- NA

  alt_id_mapping <- strsplit(go_terms_df$alt_id, split = ",")
  names(alt_id_mapping) <- go_terms_df$goAcc
  protein_id_mapping <- split(go_annotation$proteinId, go_annotation$goAcc)

  print("Parsing Gene Annotation for Protein IDs...")
  total <- length(unique(go_graph_df$accession))
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  progress <- 0

  for (i in unique(go_graph_df$accession)) {
    setTxtProgressBar(pb, progress)

    go_acc_alt_ids <- unique(unlist(alt_id_mapping[[i]]))
    go_model[go_model$goAcc == i, "alt_id"] <- paste(sort(go_acc_alt_ids), collapse = ", ")
    go_acc_all_ids <- unique(c(i, go_acc_alt_ids))
    all_protein_ids <- unique(unlist(protein_id_mapping[go_acc_all_ids]))

    if (length(all_protein_ids) > 0) {
      go_graph_df[go_graph_df$accession == i, "protein_id"] <- paste(sort(all_protein_ids), collapse = ", ")
    }

    progress <- progress + 1
  }

  close(pb)

  print("Propagating Protein IDs across the Gene Ontology Hierarchy...")

  iteration_complete <- FALSE
  total <- 10
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  progress <- 0

  edges <- go_graph_df[!is.na(go_graph_df$parent), c("accession", "parent")]
  g <- igraph::graph_from_data_frame(edges, directed = TRUE)

  setTxtProgressBar(pb, progress)

  while (!iteration_complete){
    protein_mapping <- split(go_graph_df$protein_id, go_graph_df$accession)

    order <- igraph::topo_sort(g, mode = "out")

    for (node in rev(order$name)) {
      parents <- igraph::neighbors(g, node, mode = "out")

      for (parent in parents$name) {
        protein_mapping[[parent]] <- unique(c(protein_mapping[[parent]], protein_mapping[[node]]))
      }
    }

    result <- data.frame(
      accession = names(protein_mapping),
      protein_ids = sapply(protein_mapping, function(ids) paste(sort(unique(unlist(strsplit(ids, split = ", ")))), collapse = ", ")))

    result <- result[match(go_graph_df$accession, result$accession), ]

    if(identical(go_graph_df$protein_id, result$protein_ids)){
      iteration_complete <- TRUE
    }

    else{
      go_graph_df$protein_id <- result$protein_ids
    }

    progress <- progress + 1
    setTxtProgressBar(pb, progress)
  }

  progress <- 10
  setTxtProgressBar(pb, progress)
  close(pb)

  print("Propagation complete.")

  go_graph_df$protein_id <- lapply(go_graph_df$protein_id, function(x) ifelse(x == "", NA, x))
  go_graph_df$protein_id <- unlist(go_graph_df$protein_id)

  print("Finalizing Gene Ontology Model...")
  total <- length(unique(go_graph_df$accession))
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  progress <- 0

  accession_to_protein_ids <- split(go_graph_df$protein_id, go_graph_df$accession)

  for (i in unique(go_graph_df$accession)) {
    setTxtProgressBar(pb, progress)

    protein_ids <- accession_to_protein_ids[[i]]
    protein_ids <- unique(unlist(strsplit(na.omit(protein_ids), split = ", ")))
    protein_ids <- sort(protein_ids)

    if (length(protein_ids) > 0) {
      go_model[go_model$goAcc == i, "goModel"] <- length(protein_ids)
      go_model[go_model$goAcc == i, "proteinIds"] <- paste(protein_ids, collapse = ", ")
    }

    progress <- progress + 1
  }

  close(pb)

  go_model <- go_model[go_model$goModel > 0,]
  go_model <- go_model[, c("goAcc", "alt_id", "is_a", "gotermType", "goName", "goModel", "proteinIds")]

  annotated <- c(NA,	NA, NA, NA, "ANNOTATED", 0, NA)
  unannotated <- c(NA,	NA, NA, NA, "UNANNOTATED", 0, NA)

  biological_process <- go_model[go_model$gotermType == "biological_process",]
  biological_process <- rbind(biological_process, annotated)
  biological_process_annotated <- sort(unique(unlist(strsplit(paste(biological_process$proteinIds[!is.na(biological_process$proteinIds)], collapse = ", "), split = ", "))))
  biological_process[biological_process$goName == "ANNOTATED",]$proteinIds <- paste(biological_process_annotated, collapse = ", ")
  biological_process[biological_process$goName == "ANNOTATED",]$goModel <- length(biological_process_annotated)
  biological_process <- rbind(biological_process, unannotated)
  biological_process_unannotated <- setdiff(go_annotated, biological_process_annotated)
  biological_process[biological_process$goName == "UNANNOTATED",]$proteinIds <- paste(biological_process_unannotated, collapse = ", ")
  biological_process[biological_process$goName == "UNANNOTATED",]$goModel <- length(biological_process_unannotated)
  biological_process <- biological_process[, c("goAcc", "alt_id", "is_a", "gotermType", "goName", "goModel", "proteinIds")]
  colnames(biological_process) <- c("goAcc", "alt_id", "is_a", "gotermType", "goName", "goModel", "proteinIdsModel")
  rownames(biological_process) <- NULL

  cellular_component <- go_model[go_model$gotermType == "cellular_component",]
  cellular_component <- rbind(cellular_component, annotated)
  cellular_component_annotated <- sort(unique(unlist(strsplit(paste(cellular_component$proteinIds[!is.na(cellular_component$proteinIds)], collapse = ", "), split = ", "))))
  cellular_component[cellular_component$goName == "ANNOTATED",]$proteinIds <- paste(cellular_component_annotated, collapse = ", ")
  cellular_component[cellular_component$goName == "ANNOTATED",]$goModel <- length(cellular_component_annotated)
  cellular_component <- rbind(cellular_component, unannotated)
  cellular_component_unannotated <- setdiff(go_annotated, cellular_component_annotated)
  cellular_component[cellular_component$goName == "UNANNOTATED",]$proteinIds <- paste(cellular_component_unannotated, collapse = ", ")
  cellular_component[cellular_component$goName == "UNANNOTATED",]$goModel <- length(cellular_component_unannotated)
  cellular_component <- cellular_component[, c("goAcc", "alt_id", "is_a", "gotermType", "goName", "goModel", "proteinIds")]
  colnames(cellular_component) <- c("goAcc", "alt_id", "is_a", "gotermType", "goName", "goModel", "proteinIdsModel")
  rownames(cellular_component) <- NULL

  molecular_function <- go_model[go_model$gotermType == "molecular_function",]
  molecular_function <- rbind(molecular_function, annotated)
  molecular_function_annotated <- sort(unique(unlist(strsplit(paste(molecular_function$proteinIds[!is.na(molecular_function$proteinIds)], collapse = ", "), split = ", "))))
  molecular_function[molecular_function$goName == "ANNOTATED",]$proteinIds <- paste(molecular_function_annotated, collapse = ", ")
  molecular_function[molecular_function$goName == "ANNOTATED",]$goModel <- length(molecular_function_annotated)
  molecular_function <- rbind(molecular_function, unannotated)
  molecular_function_unannotated <- setdiff(go_annotated, molecular_function_annotated)
  molecular_function[molecular_function$goName == "UNANNOTATED",]$proteinIds <- paste(molecular_function_unannotated, collapse = ", ")
  molecular_function[molecular_function$goName == "UNANNOTATED",]$goModel <- length(molecular_function_unannotated)
  molecular_function <- molecular_function[, c("goAcc", "alt_id", "is_a", "gotermType", "goName", "goModel", "proteinIds")]
  colnames(molecular_function) <- c("goAcc", "alt_id", "is_a", "gotermType", "goName", "goModel", "proteinIdsModel")
  rownames(molecular_function) <- NULL

  go_model <- list(biological_process = biological_process, cellular_component = cellular_component, molecular_function = molecular_function)

  print("Done!")

  return(go_model)
}

###########################
# Calculate GO Enrichment #
###########################
#' Calculate GO Enrichment
#'
#' @description
#' This function generates a dataframe containing background/sample frequencies and enrichment statistics per GO term from a GO model and a vector of protein IDs. Protein IDs must match the Protein IDs from the annotation used to build the model.
#'
#' @param go_model 	A list of dataframes with GO terms background frequencies per GO term type (see create_go_model)
#' @param protein_ids A vector of protein IDs. Must match the Protein IDs from the annotation used to build the model.
#' @param test A statistical method for testing the null of independence of rows and columns in a contingency table (default = "fisher"; options = c("fisher", "chisq))
#' @param p.adjust.method A statistical method to adjust p-values for multiple comparisons (default = "BH"; options = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")(see p.adjust.methods))
#' @param reduce_redundancy_by_ancestor A method to remove enriched GO terms that are ancestors of other enriched terms to reduce redundancy and retain the most specific, informative terms. (default = FALSE; options = c(TRUE, FALSE))
#' @param max_ancestor_depth Maximum distance (in terms of graph depth) allowed between a GO term and its descendant when applying redundancy reduction; higher values allow removal of more distant ancestor terms. (default = 5)
#'
#' @return A list of dataframes with GO terms enrichment statistics per GO term type
#'
#' @examples
#' go_enrichment_df <- go_enrichment(go_model_df, protein_ids)
#' go_enrichment_df <- go_enrichment(go_model_df, protein_ids, test = "fisher", p.adjust.method = "BH")
#' go_enrichment_df <- go_enrichment(go_model_df, protein_ids, test = "chisq", p.adjust.method = "bonferroni")
#' go_enrichment_df <- go_enrichment(go_model_df, protein_ids, test = "fisher", p.adjust.method = "BH", reduce_redundancy_by_ancestor = TRUE, max_ancestor_depth = 5)
#' @export

go_enrichment <- function(go_model, protein_ids, test= "fisher", p.adjust.method= "BH", reduce_redundancy_by_ancestor = FALSE, max_ancestor_depth = 5){

  go_count_df <- go_model
  go_groups <- c("biological_process", "cellular_component", "molecular_function")
  protein_ids <- sort(unique(protein_ids))

  for (group in go_groups){
    go_count_df[[group]]$proteinCount <- 0
    go_count_df[[group]]$proteinIds <- NA

    annotated_df <- go_count_df[[group]][go_count_df[[group]]$goName == "ANNOTATED",]
    unannotated_df <- go_count_df[[group]][go_count_df[[group]]$goName == "UNANNOTATED",]

    go_count_df[[group]] <- go_count_df[[group]][!is.na(go_count_df[[group]]$goAcc),]

    for (i in go_count_df[[group]]$goAcc){

      i_proteins <- sort(unlist(strsplit(go_count_df[[group]][go_count_df[[group]]$goAcc == i,]$proteinIdsModel, split = ", ")))
      intersect_protein_ids <- intersect(protein_ids, i_proteins)

      if (length(intersect_protein_ids) > 0){
        go_count_df[[group]][go_count_df[[group]]$goAcc == i,]$proteinIds <- paste(sort(unique(intersect_protein_ids)), collapse = ", ")
        go_count_df[[group]][go_count_df[[group]]$goAcc == i,]$proteinCount <- length(unique(intersect_protein_ids))
      }
    }

    group_annotated <- sort(unique(unlist(strsplit(paste(go_count_df[[group]]$proteinIds[!is.na(go_count_df[[group]]$proteinIds)], collapse = ", "), split = ", "))))
    group_unannotated <- sort(setdiff(protein_ids, group_annotated))

    go_count_df[[group]] <- rbind(go_count_df[[group]], annotated_df)
    go_count_df[[group]][go_count_df[[group]]$goName == "ANNOTATED",]$proteinIds <- paste(group_annotated, collapse = ",")
    go_count_df[[group]][go_count_df[[group]]$goName == "ANNOTATED",]$proteinCount <- length(group_annotated)

    go_count_df[[group]] <- rbind(go_count_df[[group]], unannotated_df)
    go_count_df[[group]][go_count_df[[group]]$goName == "UNANNOTATED",]$proteinIds <- paste(group_unannotated, collapse = ",")
    go_count_df[[group]][go_count_df[[group]]$goName == "UNANNOTATED",]$proteinCount <- length(group_unannotated)

    go_count_df[[group]] <- go_count_df[[group]][go_count_df[[group]]$proteinCount > 0,]

  }

  enrichment_df <- go_count_df

  for (group in go_groups){
    annotated_model <- as.numeric(enrichment_df[[group]][enrichment_df[[group]]$goName == "ANNOTATED",]$goModel)
    annotated_proteins <- as.numeric(enrichment_df[[group]][enrichment_df[[group]]$goName == "ANNOTATED",]$proteinCount)

    enrichment_df[[group]]$enrichment <- (as.numeric(enrichment_df[[group]]$proteinCount)/annotated_proteins)/(as.numeric(enrichment_df[[group]]$goModel)/annotated_model)

    enrichment_df[[group]]$enrichment[is.infinite(enrichment_df[[group]]$enrichment)] <- NaN
    enrichment_df[[group]]$enrichment[enrichment_df[[group]]$goName == "UNANNOTATED"] <- NaN
    enrichment_df[[group]]$enrichment[enrichment_df[[group]]$goName == "ANNOTATED"] <- NaN

    stat_results <- c()

    for (i in enrichment_df[[group]]$goName){
      if (is.nan(enrichment_df[[group]][enrichment_df[[group]]$goName == i,]$enrichment)){
        stat_results <- append(stat_results, NaN)
        next
    }

        contingency_table <- matrix(c(as.numeric(enrichment_df[[group]][enrichment_df[[group]]$goName == i,]$proteinCount),
                                      annotated_proteins - as.numeric(enrichment_df[[group]][enrichment_df[[group]]$goName == i,]$proteinCount),
                                  as.numeric(enrichment_df[[group]][enrichment_df[[group]]$goName == i,]$goModel),
                                  annotated_model - as.numeric(enrichment_df[[group]][enrichment_df[[group]]$goName == i,]$goModel)),2,2)
        if (test == "fisher"){
          stat_results <- append(stat_results, fisher.test(contingency_table, alternative='greater')$p.value)
    }

        else if (test == "chisq"){
          stat_results <- append(stat_results, chisq.test(contingency_table)$p.value)
    }
    }

    enrichment_df[[group]]$p.value <- stat_results
    enrichment_df[[group]] <- enrichment_df[[group]][order(enrichment_df[[group]]$p.value),]

    ##############################################
    # Redundancy Reduction

    if (reduce_redundancy_by_ancestor){

      go_graph_df <- go_terms_df[,c("goAcc", "is_a")]
      go_graph_df$is_a <- lapply(go_graph_df$is_a, function(x) ifelse(x == "", NA, x))
      go_graph_df <- tidyr::separate_rows(go_graph_df, is_a, sep = ",")
      go_graph_df <- as.data.frame(go_graph_df)
      colnames(go_graph_df) <- c("accession", "parent")
      edges <- go_graph_df[!is.na(go_graph_df$parent), c("accession", "parent")]
      g <- igraph::graph_from_data_frame(edges, directed = TRUE)

      all_graph_terms <- igraph::V(g)$name
      enriched_terms <- enrichment_df[[group]]$goAcc
      valid_terms <- enriched_terms[enriched_terms %in% all_graph_terms]
      dropped_terms <- setdiff(enriched_terms, valid_terms)
      enrichment_df[[group]] <- enrichment_df[[group]][enrichment_df[[group]]$goAcc %in% valid_terms, ]

      get_ancestors <- function(go_id, max_depth) {
        node <- igraph::V(g)[name == go_id]
        if (length(node) == 0){
          return(character(0))
        }
        ancestors <- igraph::ego(g, order = max_ancestor_depth, nodes = node, mode = "in", mindist = 1)[[1]]
        return(igraph::V(g)[ancestors]$name)
      }

      ancestor_list <- setNames(lapply(valid_terms, get_ancestors, max_depth), valid_terms)
      ancestor_to_terms <- list()

      for (term in valid_terms) {
        for (ancestor in ancestor_list[[term]]) {
          ancestor_to_terms[[ancestor]] <- c(ancestor_to_terms[[ancestor]], term)
        }
      }

      redundant_groups <- Filter(function(x) length(x) > 1, ancestor_to_terms)

      get_representative <- function(term_group) {
        group_df <- enrichment_df[[group]][enrichment_df[[group]]$goAcc %in% term_group, ]
        return(group_df$goAcc[which.min(group_df$p.value)])
      }

      reps <- unique(unlist(lapply(redundant_groups, get_representative)))
      redundant_members <- unique(unlist(redundant_groups))
      nonredundant_terms <- setdiff(valid_terms, redundant_members)
      final_terms <- unique(c(nonredundant_terms, reps))
      enrichment_df[[group]] <- enrichment_df[[group]][enrichment_df[[group]]$goAcc %in% final_terms, ]

    }

    ##############################################

  row.names(enrichment_df[[group]]) <- NULL

  enrichment_df[[group]]$adj.p.value <- p.adjust(enrichment_df[[group]]$p.value, method = p.adjust.method)

  enrichment_df[[group]] <- enrichment_df[[group]][, c("goAcc", "alt_id", "gotermType", "goName", "goModel", "proteinCount", "enrichment", "p.value", "adj.p.value", "proteinIds")]

  enrichment_df[[group]]$goModel <- as.numeric(enrichment_df[[group]]$goModel)
  enrichment_df[[group]]$proteinCount <- as.numeric(enrichment_df[[group]]$proteinCount)
  enrichment_df[[group]]$p.value <- as.numeric(enrichment_df[[group]]$p.value)
  enrichment_df[[group]]$adj.p.value <- as.numeric(enrichment_df[[group]]$adj.p.value)

  }

  return(enrichment_df)
}

#####################################
# Generate GO Terms Enrichment Plot #
#####################################
#' Generate GO Terms Enrichment Plot
#'
#' @description
#' This function generates a bar plot for visualization of GO terms enrichment statistics.
#'
#'
#' @param enrichment_df A list of dataframes with GO terms enrichment statistics per GO term type (see go_enrichment)
#' @param gotermType A GO term type for plotting (default = "biological_process"; options = c("biological_process", "molecular_function", "cellular_component"))
#' @param n A number of top-ranking GO terms to include in the plot (default = 10)
#' @param significant A color assignment for statistically significant GO terms (default = "#E58601")
#' @param plot_title A plot title (default = NA)
#' @param plot_type A plot type option (default = "bar"; options = c("bar", "lollipop"))
#'
#' @return A ggplot2 bar plot object
#'
#' @examples
#' generate_go_plot(go_enrichment_df)
#' generate_go_plot(go_enrichment_df, gotermType = "biological_process", n = 15, plot_title = "GO Enrichment: Biological Process", plot_type = "lollipop")
#' generate_go_plot(go_enrichment_df, gotermType = "molecular_function")
#' generate_go_plot(go_enrichment_df, gotermType = "cellular_component")
#' @export

generate_go_plot <- function(enrichment_df, gotermType = "biological_process", n = 10, significant = "#E58601", plot_title = NA, plot_type = "bar"){

  fig_tab <- enrichment_df[[gotermType]][1:n,]
  fig_tab$color <- "grey"

  for (i in 1:n){
    if (fig_tab[i,]$adj.p.value <= 0.05 & fig_tab[i,]$enrichment > 1){
      fig_tab[i,]$color <- significant
    }
  }

  fig_tab$goName <- factor(fig_tab$goName, levels = rev(fig_tab$goName))

  if (plot_type == "bar"){

    p <- ggplot2::ggplot(data = fig_tab, ggplot2::aes(x = goName, y = enrichment)) + ggplot2::geom_bar(stat = "identity", width = 0.75, fill = fig_tab$color) +
      ggplot2::coord_flip() + ggplot2::xlab("") + ggplot2::ylab("Enrichment") + ggplot2::ylim(0, (max(fig_tab$enrichment, na.rm = TRUE) + 0.03)) +
      ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "#B40F20") + ggplot2::theme_minimal()

    if (!is.na(plot_title)){
      p <- p + ggplot2::ggtitle(plot_title)
    }

    for (i in fig_tab$goName){
      if (fig_tab$adj.p.value[fig_tab$goName == i] <= 0.05 & fig_tab$enrichment[fig_tab$goName == i] > 1){
        if (fig_tab$adj.p.value[fig_tab$goName == i] <= 0.001){
          p <- p + ggplot2::geom_text(x=dim(fig_tab)[1]-which(fig_tab$goName == i) + 1,
                           y=fig_tab$enrichment[fig_tab$goName == i] + 0.05,
                           label="***", angle = 90)
        }

        else if (fig_tab$adj.p.value[fig_tab$goName == i] <= 0.01){
          p <- p + ggplot2::geom_text(x=dim(fig_tab)[1]-which(fig_tab$goName == i) + 1,
                           y=fig_tab$enrichment[fig_tab$goName == i] + 0.05,
                           label="**", angle = 90)
        }

        else {
          p <- p + ggplot2::geom_text(x=dim(fig_tab)[1]-which(fig_tab$goName == i) + 1,
                           y=fig_tab$enrichment[fig_tab$goName == i] + 0.05,
                           label="*" , , angle = 90)
        }

      }
    }
  }

  else if (plot_type == "lollipop"){

    p <- ggplot2::ggplot(data = fig_tab, ggplot2::aes(x = goName, y = enrichment)) +
      ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "#B40F20") +
      ggplot2::geom_segment(ggplot2::aes(xend = goName, yend = 0), color = "gray", size = 0.5) +
      ggplot2::geom_point(ggplot2::aes(size = proteinCount, color = -log10(adj.p.value)), alpha = 1) +
      ggplot2::coord_flip() +
      ggplot2::scale_color_gradient2(low = "gray99", mid = "gray85", high = significant,
                                     midpoint = 0.001, limits = c(0, 5), oob = scales::squish, name = "-log10(adj.p.value)") +
      ggplot2::scale_size(range = c(3, 8), name = "Protein Count") +
      ggplot2::labs(x = "", y = "Enrichment") +
      ggplot2::theme_minimal()

    if (!is.na(plot_title)){
      p <- p + ggplot2::ggtitle(plot_title)
    }
  }

  return(p)
}


#################################
# Load EggNOG-Mapper Annotation #
#################################
#' Load EggNOG-Mapper Annotation
#'
#' @description
#' This functions loads a eggNOG-mapper annotation into a dataframe.
#'
#' @param path Path to eggNOG-mapper annotation file
#'
#' @return A dataframe with eggNOG-mapper annotation information
#' @export
#'
#' @examples
#' eggnog_annotation_df <- load_eggnog_annotation(path/to/eggnog-mapper/annotation)

load_eggnog_annotation <- function(path){

  annotation_data <- read.delim(path, header = FALSE, sep = "\t")

  if (dim(annotation_data)[2] != 21) {
    stop("Eggnog-mapper annotation file must have 21 columns: query,	seed_ortholog,	evalue,	score,	eggNOG_OGs,	max_annot_lvl,	COG_category,
                                 Description,	Preferred_name,	GOs,	EC,	KEGG_ko,	KEGG_Pathway,	KEGG_Module,	KEGG_Reaction,
                                 KEGG_rclass,	BRITE,	KEGG_TC,	CAZy,	BiGG_Reaction &	PFAMs")
  }

  colnames(annotation_data) <- c("query",	"seed_ortholog",	"evalue",	"score",	"eggNOG_OGs",	"max_annot_lvl",	"COG_category",
                                 "Description",	"Preferred_name",	"GOs",	"EC",	"KEGG_ko",	"KEGG_Pathway",	"KEGG_Module",	"KEGG_Reaction",
                                 "KEGG_rclass",	"BRITE",	"KEGG_TC",	"CAZy",	"BiGG_Reaction",	"PFAMs")
  annotation_data <- annotation_data[!grepl("^#", annotation_data$query),]
  rownames(annotation_data) <- NULL

  return(annotation_data)
}


######################################################
# Create COG/KOG Model from eggNOG-mapper Annotation #
######################################################
#' Create COG/KOG Model from eggNOG-mapper Annotation
#'
#' @description
#' This function generates a background frequency model for COG/KOG categories enrichment analysis based on eggNOG-mapper annotation.
#'
#' @param eggnog_annotation A dataframe with eggNOG-mapper annotation information (see load_eggnog_annotation)
#'
#' @return A dataframe with COG/KOG categories background frequencies
#' @export
#'
#' @examples
#' kog_model_df <- create_kog_model_eggnog(eggnog_annotation_df)

create_kog_model_eggnog <- function(eggnog_annotation){

  data(cog_classes, envir = environment())

  kog_annotation <- eggnog_annotation
  kog_annotation <- kog_annotation[kog_annotation$COG_category != "-",c("query", "COG_category")]
  kog_annotation <- tidyr::separate_rows(kog_annotation, COG_category, sep = "")
  kog_annotation <- kog_annotation[kog_annotation$COG_category != "",]
  colnames(kog_annotation) <- c("proteinId", "kogClass")
  mapping_vector <- setNames(cog_classes_df$class, cog_classes_df$code)
  kog_annotation$kogClass <- mapping_vector[kog_annotation$kogClass]

  return(create_kog_model(kog_annotation))
}

#################################################
# Create GO Model from eggNOG-mapper Annotation #
#################################################
#' Create GO Model from eggNOG-mapper Annotation
#'
#' @description
#' This function generates a background frequency model for GO term enrichment analysis based on eggNOG-mapper annotation.
#'
#' @param eggnog_annotation A dataframe with eggNOG-mapper annotation information (see load_eggnog_annotation)
#'
#' @return A list of dataframes with GO terms background frequencies by gotermType (biological_processes, cellular_component, molecular_function)
#' @export
#'
#' @examples
#' go_model_df <- create_go_model_eggnog(eggnog_annotation_df)
#' go_model_biological_processes_df <- go_model_df$biological_processes
#' go_model_cellular_component_df <- go_model_df$cellular_component
#' go_model_molecular_function_df <- go_model_df$molecular_function

create_go_model_eggnog <- function(eggnog_annotation){

  go_annotation <- eggnog_annotation[, c("query", "GOs")]
  go_annotation <- go_annotation[go_annotation$GOs != "-",]
  colnames(go_annotation) <- c("proteinId", "goAcc")
  go_annotation <- tidyr::separate_rows(go_annotation, goAcc, sep = ",")

  return(create_go_model(go_annotation))
}

####################################
# KEGG Pathway Enrichment Analysis #
####################################

########################
# Load KEGG Annotation #
########################
#' Load KEGG Annotation
#'
#' @description
#' This function loads a KEGG annotation into a dataframe.
#'
#'
#' @param path 	Path to KEGG annotation
#'
#' @return A dataframe with KEGG annotation information
#'
#' @examples
#' kegg_annotation_df <- load_kegg_annotation(path/to/kegg/annotation)
#' @export

load_kegg_annotation <- function(path){

  annotation_data <- read.delim(path, header = FALSE, sep = "\t")

  firstup <- function(x) {
    x <- paste(substr(x, 1, 1), tolower(substr(x, 2, nchar(x))), sep = "")
    return(x)
  }

  if (dim(annotation_data)[2] != 2) {
    stop("KEGG annotation file must have two columns: transcript_id,	KEGG_ko")
  }

  colnames(annotation_data) <- c("proteinId",	"KEGG_ko")
  annotation_data$proteinId <- sub(".*:", "", annotation_data$proteinId)
  annotation_data$KEGG_ko <- sub("ko:", "", annotation_data$KEGG_ko)
  annotation_data$KEGG_ko <- sub("ko", "", annotation_data$KEGG_ko)

  return(annotation_data)
}


#############################
# Create KEGG Pathway Model #
#############################
#' Create KEGG Pathway Model
#'
#' @description
#' This function generates a background frequency model for KEGG pathway enrichment analysis.
#'
#' @param kegg_annotation A dataframe with KEGG annotation information (see load_kegg_annotation)
#'
#' @return A list of dataframes with KEGG pathway background frequencies per KEGG category (pathway_type, pathway_class, pathway_name, definition, complete)
#'
#' @examples
#' kegg_model_df <- create_kegg_model(kegg_annotation_df)
#' kegg_model_type_df <- kegg_model_df$pathway_type
#' kegg_model_class_df <- kegg_model_df$pathway_class
#' kegg_model_name_df <- kegg_model_df$pathway_name
#' kegg_model_definition_df <- kegg_model_df$definition
#' kegg_model_complete_df <- kegg_model_df$complete
#' @export

create_kegg_model <- function(kegg_annotation){

  protein_ids <- sort(unique(kegg_annotation$proteinId))
  data(kegg_pathways, envir = environment())
  data(kegg_pathways_filter, envir = environment())
  pathways_df <- kegg_pathways

  pathways_df$keggModel <- 0
  pathways_df$proteinIdsModel <- NA

  KEGG_ko <- sort(unique(pathways_df$KEGG_ko))
  filter = "All"

  for (i in KEGG_ko){

    proteins_i <- sort(unique(na.omit(kegg_annotation[kegg_annotation$KEGG_ko == i,]$proteinId)))

    if (length(proteins_i) > 0) {
      pathways_df[pathways_df$KEGG_ko == i,]$keggModel <- length(proteins_i)
      pathways_df[pathways_df$KEGG_ko == i,]$proteinIdsModel <- paste(proteins_i, collapse = ", ")
    }
  }

  pathways_df <- pathways_df[!is.na(pathways_df$proteinIdsModel),]

  if (filter == "All"){
    pathways_df <- pathways_df
  }

  else if (filter != "All"){

    if (filter %in% kegg_pathways_filter$kingdom){
      filtered_kos <- sort(unique(kegg_pathways_filter[kegg_pathways_filter$kingdom == filter,]$pathway))
    }

    else if (filter %in% kegg_pathways_filter$phylum){
      filtered_kos <- sort(unique(kegg_pathways_filter[kegg_pathways_filter$phylum == filter,]$pathway))
    }

    else if (filter %in% kegg_pathways_filter$class){
      filtered_kos <- sort(unique(kegg_pathways_filter[kegg_pathways_filter$class == filter,]$pathway))
    }

    else if (filter %in% kegg_pathways_filter$specie){
      filtered_kos <- sort(unique(kegg_pathways_filter[kegg_pathways_filter$specie == filter,]$pathway))
    }

    else {
      filters <- unique(c(kegg_pathways_filter$kingdom, kegg_pathways_filter$phylum, kegg_pathways_filter$class, kegg_pathways_filter$specie))
      stop(sprintf("Invalid filter option. Supported filter options are: '%s'", paste(filters, collapse = ", ")))
    }

    pathways_df <- pathways_df[pathways_df$pathway %in% filtered_kos,]

  }

  complete_df <- as.data.frame(pathways_df)
  complete_annotated_proteins <- sort(unique(unlist(strsplit(paste(na.omit(complete_df$proteinIdsModel), collapse = ", "), split = ", "))))
  complete_unannotated_proteins <- sort(setdiff(protein_ids, complete_annotated_proteins))
  complete_annotated <- c(NA, NA, NA, NA, NA, "ANNOTATED", NA , 0, NA)
  complete_unannotated <- c(NA, NA, NA, NA, NA, "UNANNOTATED", NA , 0, NA)
  complete_df <- rbind(complete_df, complete_annotated, complete_unannotated)
  complete_df$keggModel <- as.numeric(complete_df$keggModel)
  complete_df[complete_df$definition == "ANNOTATED",]$keggModel <- length(complete_annotated_proteins)
  complete_df[complete_df$definition == "UNANNOTATED",]$keggModel <- length(complete_unannotated_proteins)
  complete_df[complete_df$definition == "ANNOTATED",]$proteinIdsModel <- paste(complete_annotated_proteins, collapse = ", ")
  complete_df[complete_df$definition == "UNANNOTATED",]$proteinIdsModel <- paste(complete_unannotated_proteins, collapse = ", ")

  definition_df <- dplyr::summarize(dplyr::group_by(pathways_df, definition), proteinIds = paste(na.omit(proteinIdsModel), collapse = ", "))
  definition_df$proteinIdsModel <- NA
  definition_df$keggModel <- 0

  for (i in definition_df$definition){

    proteins_i <- sort(unique(unlist(strsplit(na.omit(definition_df[definition_df$definition == i,]$proteinIds), split = ", "))))

    if (length(proteins_i) > 0){
      definition_df[definition_df$definition == i,]$proteinIdsModel <- paste(proteins_i, collapse = ", ")
      definition_df[definition_df$definition == i,]$keggModel <- length(proteins_i)
    }
  }

  definition_df <- as.data.frame(definition_df[,c("definition", "keggModel", "proteinIdsModel")])
  definition_annotated <- c("ANNOTATED", 0, NA)
  definition_unannotated <- c("UNANNOTATED", 0, NA)
  definition_df <- rbind(definition_df, definition_annotated, definition_unannotated)
  definition_df$keggModel <- as.numeric(definition_df$keggModel)
  definition_df[definition_df$definition == "ANNOTATED",]$keggModel <- length(complete_annotated_proteins)
  definition_df[definition_df$definition == "UNANNOTATED",]$keggModel <- length(complete_unannotated_proteins)
  definition_df[definition_df$definition == "ANNOTATED",]$proteinIdsModel <- paste(complete_annotated_proteins, collapse = ", ")
  definition_df[definition_df$definition == "UNANNOTATED",]$proteinIdsModel <- paste(complete_unannotated_proteins, collapse = ", ")

  pathway_name_df <- dplyr::summarize(dplyr::group_by(pathways_df, pathway_type, pathway_class, pathway_name), proteinIds = paste(na.omit(proteinIdsModel), collapse = ", "))
  pathway_name_df$proteinIdsModel <- NA
  pathway_name_df$keggModel <- 0

  for (i in pathway_name_df$pathway_name){

    proteins_i <- sort(unique(unlist(strsplit(na.omit(pathway_name_df[pathway_name_df$pathway_name == i,]$proteinIds), split = ", "))))

    if (length(proteins_i) > 0){
      pathway_name_df[pathway_name_df$pathway_name == i,]$proteinIdsModel <- paste(proteins_i, collapse = ", ")
      pathway_name_df[pathway_name_df$pathway_name == i,]$keggModel <- length(proteins_i)
    }
  }

  pathway_name_df <- as.data.frame(pathway_name_df[,c("pathway_type", "pathway_class", "pathway_name", "keggModel", "proteinIdsModel")])
  pathway_name_df$keggModel <- as.numeric(pathway_name_df$keggModel)
  pathway_name_annotated <- data.frame(pathway_type = c(NA, NA), pathway_class = c(NA, NA), pathway_name = c("ANNOTATED", "UNANNOTATED"), keggModel = c(0, 0), proteinIdsModel = c(NA, NA))
  pathway_name_df <- rbind(pathway_name_df, pathway_name_annotated)
  pathway_name_df[pathway_name_df$pathway_name == "ANNOTATED",]$keggModel <- length(complete_annotated_proteins)
  pathway_name_df[pathway_name_df$pathway_name == "UNANNOTATED",]$keggModel <- length(complete_unannotated_proteins)
  pathway_name_df[pathway_name_df$pathway_name == "ANNOTATED",]$proteinIdsModel <- paste(complete_annotated_proteins, collapse = ", ")
  pathway_name_df[pathway_name_df$pathway_name == "UNANNOTATED",]$proteinIdsModel <- paste(complete_unannotated_proteins, collapse = ", ")

  pathway_class_df <- dplyr::summarize(dplyr::group_by(pathways_df, pathway_type, pathway_class), proteinIds = paste(na.omit(proteinIdsModel), collapse = ", "))
  pathway_class_df$proteinIdsModel <- NA
  pathway_class_df$keggModel <- 0

  for (i in pathway_class_df$pathway_class){

    proteins_i <- sort(unique(unlist(strsplit(na.omit(pathway_class_df[pathway_class_df$pathway_class == i,]$proteinIds), split = ", "))))

    if (length(proteins_i) > 0){
      pathway_class_df[pathway_class_df$pathway_class == i,]$proteinIdsModel <- paste(proteins_i, collapse = ", ")
      pathway_class_df[pathway_class_df$pathway_class == i,]$keggModel <- length(proteins_i)
    }
  }

  pathway_class_df <- as.data.frame(pathway_class_df[,c("pathway_type", "pathway_class", "keggModel", "proteinIdsModel")])
  pathway_class_df$keggModel <- as.numeric(pathway_class_df$keggModel)
  pathway_class_annotated <- data.frame(pathway_type = c(NA, NA), pathway_class = c("ANNOTATED", "UNANNOTATED"), keggModel = c(0, 0), proteinIdsModel = c(NA, NA))
  pathway_class_df <- rbind(pathway_class_df, pathway_class_annotated)
  pathway_class_df[pathway_class_df$pathway_class == "ANNOTATED",]$keggModel <- length(complete_annotated_proteins)
  pathway_class_df[pathway_class_df$pathway_class == "UNANNOTATED",]$keggModel <- length(complete_unannotated_proteins)
  pathway_class_df[pathway_class_df$pathway_class == "ANNOTATED",]$proteinIdsModel <- paste(complete_annotated_proteins, collapse = ", ")
  pathway_class_df[pathway_class_df$pathway_class == "UNANNOTATED",]$proteinIdsModel <- paste(complete_unannotated_proteins, collapse = ", ")

  pathway_type_df <- dplyr::summarize(dplyr::group_by(pathways_df, pathway_type), proteinIds = paste(na.omit(proteinIdsModel), collapse = ", "))
  pathway_type_df$proteinIdsModel <- NA
  pathway_type_df$keggModel <- 0

  for (i in pathway_type_df$pathway_type){

    proteins_i <- sort(unique(unlist(strsplit(na.omit(pathway_type_df[pathway_type_df$pathway_type == i,]$proteinIds), split = ", "))))

    if (length(proteins_i) > 0){
      pathway_type_df[pathway_type_df$pathway_type == i,]$proteinIdsModel <- paste(proteins_i, collapse = ", ")
      pathway_type_df[pathway_type_df$pathway_type == i,]$keggModel <- length(proteins_i)
    }
  }

  pathway_type_df <- as.data.frame(pathway_type_df[,c("pathway_type", "keggModel", "proteinIdsModel")])
  pathway_type_df$keggModel <- as.numeric(pathway_type_df$keggModel)
  pathway_type_annotated <- data.frame(pathway_type = c("ANNOTATED", "UNANNOTATED"), keggModel = c(0, 0), proteinIdsModel = c(NA, NA))
  pathway_type_df <- rbind(pathway_type_df, pathway_type_annotated)
  pathway_type_df[pathway_type_df$pathway_type == "ANNOTATED",]$keggModel <- length(complete_annotated_proteins)
  pathway_type_df[pathway_type_df$pathway_type == "UNANNOTATED",]$keggModel <- length(complete_unannotated_proteins)
  pathway_type_df[pathway_type_df$pathway_type == "ANNOTATED",]$proteinIdsModel <- paste(complete_annotated_proteins, collapse = ", ")
  pathway_type_df[pathway_type_df$pathway_type == "UNANNOTATED",]$proteinIdsModel <- paste(complete_unannotated_proteins, collapse = ", ")

  kegg_model_df <- list(complete = complete_df, definition = definition_df, pathway_name = pathway_name_df, pathway_class = pathway_class_df,
                        pathway_type = pathway_type_df)

  return(kegg_model_df)
}


#####################################
# Calculate KEGG Pathway Enrichment #
#####################################
#' Calculate KEGG Pathway Enrichment
#'
#' @description
#' This function generates a list of dataframes containing background/sample frequencies and enrichment statistics per KEGG category.
#'
#' @param kegg_model 	A list of dataframes with KEGG pathway background frequencies per KEGG category (pathway_type, pathway_class, pathway_name, definition, complete) (see create_kegg_model)
#' @param protein_ids A vector of protein IDs. Must match the Protein IDs from the annotation used to build the model.
#' @param test A statistical method for testing the null of independence of rows and columns in a contingency table (default = "fisher"; options = c("fisher", "chisq))
#' @param p.adjust.method A statistical method to adjust p-values for multiple comparisons (default = "BH"; options = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")(see p.adjust.methods))
#'
#' @return A list of dataframes with KEGG pathways enrichment statistics per KEGG category (pathway_type, pathway_class, pathway_name, definition, complete)
#'
#' @examples
#' kegg_enrichment_df <- kegg_enrichment(kegg_model_df, protein_ids)
#' kegg_enrichment_df <- kegg_enrichment(kegg_model_df, protein_ids, test = "fisher", p.adjust.method = "BH")
#' kegg_enrichment_df <- kegg_enrichment(kegg_model_df, protein_ids, test = "chisq", p.adjust.method = "bonferroni")
#' kegg_enrichment_type_df <- kegg_enrichment_df$pathway_type
#' kegg_enrichment_class_df <- kegg_enrichment_df$pathway_class
#' kegg_enrichment_name_df <- kegg_enrichment_df$pathway_name
#' kegg_enrichment_definition_df <- kegg_enrichment_df$definition
#' kegg_enrichment_complete_df <- kegg_enrichment_df$complete
#' @export

kegg_enrichment <- function(kegg_model, protein_ids, test= "fisher", p.adjust.method= "BH"){

  kegg_count_df <- kegg_model
  kegg_groups <- c("complete", "definition", "pathway_name", "pathway_class", "pathway_type")

  protein_ids <- sort(unique(protein_ids))

  for (group in kegg_groups){
    kegg_count_df[[group]]$proteinCount <- 0
    kegg_count_df[[group]]$proteinIds <- NA

    annotated_df <- kegg_count_df[[group]][nrow(kegg_count_df[[group]]) - 1,]
    unannotated_df <- kegg_count_df[[group]][nrow(kegg_count_df[[group]]),]

    kegg_count_df[[group]] <- kegg_count_df[[group]][1:(nrow(kegg_count_df[[group]]) - 2),]

    for (i in 1:nrow(kegg_count_df[[group]])){

      i_proteins <- sort(unlist(strsplit(kegg_count_df[[group]][i,]$proteinIdsModel, split = ", ")))
      intersect_protein_ids <- intersect(protein_ids, i_proteins)

      if (length(intersect_protein_ids) > 0){
        kegg_count_df[[group]][i,]$proteinIds <- paste(sort(unique(intersect_protein_ids)), collapse = ", ")
        kegg_count_df[[group]][i,]$proteinCount <- length(unique(intersect_protein_ids))
      }
    }

    group_annotated <- sort(unique(unlist(strsplit(paste(kegg_count_df[[group]]$proteinIds[!is.na(kegg_count_df[[group]]$proteinIds)], collapse = ", "), split = ", "))))
    group_unannotated <- sort(setdiff(protein_ids, group_annotated))

    kegg_count_df[[group]] <- rbind(kegg_count_df[[group]], annotated_df)
    kegg_count_df[[group]][nrow(kegg_count_df[[group]]),]$proteinIds <- paste(group_annotated, collapse = ",")
    kegg_count_df[[group]][nrow(kegg_count_df[[group]]),]$proteinCount <- length(group_annotated)

    kegg_count_df[[group]] <- rbind(kegg_count_df[[group]], unannotated_df)
    kegg_count_df[[group]][nrow(kegg_count_df[[group]]),]$proteinIds <- paste(group_unannotated, collapse = ",")
    kegg_count_df[[group]][nrow(kegg_count_df[[group]]),]$proteinCount <- length(group_unannotated)

    kegg_count_df[[group]] <- kegg_count_df[[group]][kegg_count_df[[group]]$proteinCount > 0,]
    kegg_count_df[[group]] <- subset(kegg_count_df[[group]], select = -proteinIdsModel)

  }

  enrichment_df <- kegg_count_df

  for (group in kegg_groups){

    protein_ids_group <- enrichment_df[[group]]$proteinIds
    enrichment_df[[group]] <- subset(enrichment_df[[group]], select = -proteinIds)

    annotated_model <- as.numeric(kegg_count_df[[group]][nrow(kegg_count_df[[group]]) - 1,]$keggModel)
    annotated_proteins <- as.numeric(kegg_count_df[[group]][nrow(kegg_count_df[[group]]) - 1,]$proteinCount)

    enrichment_df[[group]]$enrichment <- (as.numeric(enrichment_df[[group]]$proteinCount)/annotated_proteins)/(as.numeric(enrichment_df[[group]]$keggModel)/annotated_model)

    enrichment_df[[group]]$enrichment[is.infinite(enrichment_df[[group]]$enrichment)] <- NaN
    enrichment_df[[group]][nrow(kegg_count_df[[group]]),]$enrichment <- NaN
    enrichment_df[[group]][(nrow(kegg_count_df[[group]]) - 1),]$enrichment <- NaN

    stat_results <- c()

    for (i in 1:nrow(enrichment_df[[group]])){
      if (is.nan(enrichment_df[[group]][i,]$enrichment)){
        stat_results <- append(stat_results, NaN)
        next
      }

      contingency_table <- matrix(c(as.numeric(enrichment_df[[group]][i,]$proteinCount),
                                    annotated_proteins - as.numeric(enrichment_df[[group]][i,]$proteinCount),
                                    as.numeric(enrichment_df[[group]][i,]$keggModel),
                                    annotated_model - as.numeric(enrichment_df[[group]][i,]$keggModel)),2,2)
      if (test == "fisher"){
        stat_results <- append(stat_results, fisher.test(contingency_table, alternative='greater')$p.value)
      }

      else if (test == "chisq"){
        stat_results <- append(stat_results, chisq.test(contingency_table)$p.value)
      }
    }

    enrichment_df[[group]]$p.value <- stat_results

    enrichment_df[[group]] <- enrichment_df[[group]][order(enrichment_df[[group]]$p.value),]
    row.names(enrichment_df[[group]]) <- NULL

    enrichment_df[[group]]$adj.p.value <- p.adjust(enrichment_df[[group]]$p.value, method = p.adjust.method)

    enrichment_df[[group]]$proteinIds <-protein_ids_group

    enrichment_df[[group]] <- as.data.frame(enrichment_df[[group]])

    enrichment_df[[group]]$keggModel <- as.numeric(enrichment_df[[group]]$keggModel)
    enrichment_df[[group]]$proteinCount <- as.numeric(enrichment_df[[group]]$proteinCount)
    enrichment_df[[group]]$p.value <- as.numeric(enrichment_df[[group]]$p.value)
    enrichment_df[[group]]$adj.p.value <- as.numeric(enrichment_df[[group]]$adj.p.value)

  }

  return(enrichment_df)
}

###########################################################
# Create KEGG Pathway Model from eggNOG-mapper Annotation #
###########################################################
#' Create KEGG Pathway Model from eggNOG-mapper Annotation
#'
#' @description
#' This function generates a background frequency model for KEGG pathway enrichment analysis based on eggNOG-mapper annotation.
#'
#' @param eggnog_annotation A dataframe with eggNOG-mapper annotation information (see load_eggnog_annotation)
#'
#' @return A list of dataframes with KEGG pathway background frequencies per KEGG category (pathway_type, pathway_class, pathway_name, definition, complete).
#' @export
#'
#' @examples
#' kegg_model_df <- create_kegg_model_eggnog(eggnog_annotation_df)

create_kegg_model_eggnog <- function(eggnog_annotation){

  filter = "All"
  ec_annotation <- eggnog_annotation[, c("query", "KEGG_ko", "EC")]
  ec_annotation <- tidyr::separate_rows(ec_annotation, KEGG_ko, sep = ",")
  ec_annotation$KEGG_ko <- gsub("ko:", "", ec_annotation$KEGG_ko)
  ec_annotation <- tidyr::separate_rows(ec_annotation, EC, sep = ",")
  ec_annotation[ec_annotation$KEGG_ko == "-",]$KEGG_ko <- NA
  ec_annotation[ec_annotation$EC == "-",]$EC <- NA
  colnames(ec_annotation) <- c("proteinId", "KEGG_ko", "ecNum")

  return(create_kegg_model(ec_annotation))
}


#########################################
# Generate KEGG Pathway Enrichment Plot #
#########################################
#' Generate KEGG Pathway Enrichment Plot
#'
#' @description
#' This function generates a bar plot for visualization of KEGG pathway enrichment statistics.
#'
#' @param enrichment_df A list of dataframes with KEGG pathway enrichment statistics per KEGG category (pathway_type, pathway_class, pathway_name, definition, complete) (see kegg_enrichment)
#' @param pathwayType A KEGG category for plotting (default = "pathway_class"; options = c("definition", "pathway_name", "pathway_class", "pathway_type"))
#' @param n A number of top-ranking KEGG pathways to include in the plot (default = 10)
#' @param significant A color assignment for statistically significant GO terms (default = "#9986A5")
#' @param plot_title A plot title (default = NA)
#' @param plot_type A plot type option (default = "bar"; options = c("bar", "lollipop"))
#'
#' @return A ggplot2 plot object
#'
#' @examples
#' generate_kegg_plot(kegg_enrichment_df)
#' generate_kegg_plot(kegg_enrichment_df, model = "pathway_class", plot_title = "KEGG Pathway Enrichment: Class")
#' generate_kegg_plot(kegg_enrichment_df, model = "pathway_name", significant = "blue")
#' generate_kegg_plot(kegg_enrichment_df, model = "pathway_name", n = 5, plot_type = "lollipop")
#' @export

generate_kegg_plot <- function(enrichment_df, model = "pathway_class", n = 10, significant = "#9986A5", plot_title = NA, plot_type = "bar"){

  allowed_models <- c("definition", "pathway_name", "pathway_class", "pathway_type")

  if (!model %in% allowed_models) {
    stop(sprintf(
      "Invalid model parameter: '%s'. Allowed values are: %s",
      model,
      paste(allowed_models, collapse = ", ")
    ))
  }

  if (model == "pathway_type"){
    n <- nrow(enrichment_df[[model]]) - 2
  }
  else if (n > nrow(enrichment_df[[model]])){
    n <- nrow(enrichment_df[[model]]) - 2
  }

  fig_tab <- enrichment_df[[model]][1:n,]

  if (model == "pathway_name"){
    fig_tab <- fig_tab[, c("pathway_name", "keggModel", "proteinCount", "enrichment", "p.value", "adj.p.value")]
    colnames(fig_tab) <- c("keggPathway", "keggModel", "proteinCount", "enrichment", "p.value", "adj.p.value")
    fig_tab$keggPathway <- factor(fig_tab$keggPathway, levels = rev(fig_tab$keggPathway))
  }

  else if (model == "definition"){
    fig_tab <- fig_tab[, c("definition", "keggModel", "proteinCount", "enrichment", "p.value", "adj.p.value")]
    colnames(fig_tab) <- c("keggPathway", "keggModel", "proteinCount", "enrichment", "p.value", "adj.p.value")
    fig_tab$keggPathway <- factor(fig_tab$keggPathway, levels = rev(fig_tab$keggPathway))
  }

  else if (model == "pathway_class"){
    fig_tab <- fig_tab[, c("pathway_class", "keggModel", "proteinCount", "enrichment", "p.value", "adj.p.value")]
    colnames(fig_tab) <- c("keggPathway", "keggModel", "proteinCount", "enrichment", "p.value", "adj.p.value")
    fig_tab$keggPathway <- factor(fig_tab$keggPathway, levels = rev(fig_tab$keggPathway))
  }

  else if (model == "pathway_type"){
    fig_tab <- fig_tab[, c("pathway_type", "keggModel", "proteinCount", "enrichment", "p.value", "adj.p.value")]
    colnames(fig_tab) <- c("keggPathway", "keggModel", "proteinCount", "enrichment", "p.value", "adj.p.value")
    fig_tab$keggPathway <- factor(fig_tab$keggPathway, levels = rev(fig_tab$keggPathway)) # Correct this for order in KEGG
  }

  fig_tab$color <- "grey"

  for (i in 1:n){
    if (fig_tab[i,]$adj.p.value <= 0.05 & fig_tab[i,]$enrichment > 1){
      fig_tab[i,]$color <- significant
    }
  }

  if (plot_type == "bar"){

    p <- ggplot2::ggplot(data = fig_tab, ggplot2::aes(x = keggPathway, y = enrichment)) + ggplot2::geom_bar(stat = "identity", width = 0.75, fill = fig_tab$color) +
      ggplot2::coord_flip() + ggplot2::xlab("") + ggplot2::ylab("Enrichment") + ggplot2::ylim(0, (max(fig_tab$enrichment, na.rm = TRUE) + 0.03)) +
      ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "#B40F20") +ggplot2::theme_minimal()

    if (!is.na(plot_title)){
      p <- p + ggplot2::ggtitle(plot_title)
    }

    for (i in fig_tab$keggPathway){
      if (fig_tab$adj.p.value[fig_tab$keggPathway == i] <= 0.05 & fig_tab$enrichment[fig_tab$keggPathway == i] > 1){
        if (fig_tab$adj.p.value[fig_tab$keggPathway == i] <= 0.001){
          p <- p + ggplot2::geom_text(x=dim(fig_tab)[1]-which(fig_tab$keggPathway == i) + 1,
                                    y=fig_tab$enrichment[fig_tab$keggPathway == i] + 0.05,
                                    label="***", angle = 90)
        }

        else if (fig_tab$adj.p.value[fig_tab$keggPathway == i] <= 0.01){
          p <- p + ggplot2::geom_text(x=dim(fig_tab)[1]-which(fig_tab$keggPathway == i) + 1,
                                    y=fig_tab$enrichment[fig_tab$keggPathway == i] + 0.05,
                                    label="**", angle = 90)
        }

        else {
          p <- p + ggplot2::geom_text(x=dim(fig_tab)[1]-which(fig_tab$keggPathway == i) + 1,
                                    y=fig_tab$enrichment[fig_tab$keggPathway == i] + 0.05,
                                    label="*" , , angle = 90)
        }

      }
    }
  }

  else if (plot_type == "lollipop"){

    p <- ggplot2::ggplot(data = fig_tab, ggplot2::aes(x = keggPathway, y = enrichment)) +
      ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "#B40F20") +
      ggplot2::geom_segment(ggplot2::aes(xend = keggPathway, yend = 0), color = "gray", size = 0.5) +
      ggplot2::geom_point(ggplot2::aes(size = proteinCount, color = -log10(adj.p.value)), alpha = 1) +
      ggplot2::coord_flip() +
      ggplot2::scale_color_gradient2(low = "gray99", mid = "gray85", high = significant,
                                     midpoint = 0.001, limits = c(0, 5), oob = scales::squish, name = "-log10(adj.p.value)") +
      ggplot2::scale_size(range = c(3, 8), name = "Protein Count") +
      ggplot2::labs(x = "", y = "Enrichment") +
      ggplot2::theme_minimal()

    if (!is.na(plot_title)){
      p <- p + ggplot2::ggtitle(plot_title)
    }
  }

  return(p)
}

################
# Fetch Models #
################
#' Fetch Models
#'
#' @description
#' This function loads precomputed KOG/GO/KEGG models for an organism and print infos on data origin
#'
#' @param strain Strain identifier used to query the database for available precomputed models (see available_models$strain).
#'
#' @return A list with model information, gtf annotation and precomputed KOG/GO/KEGG models (info, gtf, transcript2protein_id, kog, go, kegg)
#'
#' @examples
#' qm6a_models <- fetch_models("Trichoderma reesei (QM6a)")
#' qm6a_gtf <- qm6a_models$gtf
#' qm6a_transcript2protein_id <- qm6a_models$transcript2protein_id
#' qm6a_kog <- qm6a_models$kog
#' qm6a_go <- qm6a_models$go
#' qm6a_kegg <- qm6a_models$kegg
#' @export

fetch_models <- function(strain){

  data(available_models, envir = environment())

  if (strain %in% available_models$strain){
    file_name <- available_models[available_models$strain == strain,]$file
    model_env <- new.env()
    data(list = file_name, envir = model_env)
    model_data <- get(file_name, envir = model_env)
    print(model_data$info)
    return(model_data)
  }

  else {
    models <- sort(available_models$strain)
    print(sprintf("Invalid strain name. Available strain models are: '%s'", paste(models, collapse = ", ")))
  }
}

#######################
# Loads EC Annotation #
#######################
#' Load EC Annotation
#'
#' @description
#' This function loads a EC (Enzyme Commission) based KEGG annotation into a dataframe.
#'
#' @param path 	Path to EC annotation
#'
#' @return A dataframe with EC annotation information
#'
#' @examples
#' ec_annotation_df <- load_ec_annotation(path/to/ec/annotation)
#' @export

load_ec_annotation <- function(path){

  annotation_data <- read.delim(path, header = TRUE, sep = "\t")

  firstup <- function(x) {
    x <- paste(substr(x, 1, 1), tolower(substr(x, 2, nchar(x))), sep = "")
    return(x)
  }

  if (dim(annotation_data)[2] != 9) {
    stop("EC annotation file must have nine columns: proteinId,	ecNum,	definition,	catalyticActivity,	cofactors,	associatedDiseases,	pathway,	pathway_class &	pathway_type")
  }

  colnames(annotation_data) <- c("proteinId",	"ecNum",	"definition",	"catalyticActivity",	"cofactors",	"associatedDiseases",	"pathway",	"pathway_class",	"pathway_type")

  annotation_data$pathway <- trimws(annotation_data$pathway)

  for (i in 1:nrow(annotation_data)){

    annotation_data[i,]$pathway <- firstup(annotation_data[i,]$pathway)
    annotation_data[i,]$pathway_class <- firstup(annotation_data[i,]$pathway_class)

    if (annotation_data[i,]$catalyticActivity == "\\N"){
      annotation_data[i,]$catalyticActivity <- NA
    }
    if (annotation_data[i,]$cofactors == "\\N"){
      annotation_data[i,]$cofactors <- NA
    }
    if (annotation_data[i,]$associatedDiseases == "\\N"){
      annotation_data[i,]$associatedDiseases <- NA
    }
    if (annotation_data[i,]$pathway == "\\n"){
      annotation_data[i,]$pathway <- NA
    }
    if (annotation_data[i,]$pathway_class == "\\n"){
      annotation_data[i,]$pathway_class <- NA
    }
    if (annotation_data[i,]$pathway_type == "\\N"){
      annotation_data[i,]$pathway_type <- NA
    }

  }

  return(annotation_data)
}


################################################
# Create KEGG Pathway Model from EC Annotation #
################################################
#' Create KEGG Pathway Model from EC Annotation
#'
#' @description
#' This function generates a background frequency model for KEGG pathway enrichment analysis from an EC (Enzyme Commission) Annotation.
#'
#' @param ec_annotation A dataframe with EC annotation information (see load_ec_annotation)
#'
#' @return A list of dataframes with KEGG pathway background frequencies per KEGG category (pathway_type, pathway_class, pathway_name, definition, complete)
#'
#' @examples
#' kegg_model_df <- create_kegg_model_from_ec(ec_annotation_df, filter = "All")
#' @export

create_kegg_model_from_ec <- function(ec_annotation_df){

  kegg_annotation <- ec_annotation_df[, c("proteinId", "ecNum")]
  kegg_annotation <- unique(kegg_annotation)
  kegg_annotation$KEGG_ko <- NA
  filter = "All"

  data(kegg_pathways, envir = environment())

  for (i in 1:nrow(kegg_annotation)){
    ec_ko <- kegg_pathways[grepl(kegg_annotation[i,]$ecNum, kegg_pathways$ecNum),]$KEGG_ko

    if (length(ec_ko) > 0){
    kegg_annotation[i,]$KEGG_ko <- paste(ec_ko, collapse = ", ")
    }
  }

  kegg_annotation <- kegg_annotation[, c("proteinId", "KEGG_ko")]
  kegg_annotation <- as.data.frame(tidyr::separate_rows(kegg_annotation, KEGG_ko, sep = ", "))

  return(create_kegg_model(kegg_annotation))
}
