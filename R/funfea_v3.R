######################################
# Package for functional annotations #
# enrichment analysis and plotting   #
#                                    #
# Written by Julien Charest &        #
# Paul Loebenstein                   #
#                                    #
# Version 3                          #
# 25.11.2024                         #
######################################

#########
# To Do #
#########

# Rename functions
# Update roxygen comments
# Handle errors

#############
# Utilities #
#############

#################################
# Loads a GTF annotation #

#' Load GTF Annotation
#'
#' @description
#' This function loads a GTF annotation into a dataframe.
#'
#' @param path Path to GTF annotation
#' @return A dataframe with GTF annotation information
#' @examples
#' gtf_df <- load_gtf_annotation(gtf_annotation_path)
#' @export

load_gtf_annotation <- function(path){

  gtf_data <- read.delim(path, header = FALSE, sep = "\t")
  colnames(gtf_data) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
  gtf_data <- gtf_data[!grepl("^#", gtf_data$seqname),]
  rownames(gtf_data) <- NULL

  return(gtf_data)
}

#####################################################
# Generates a transcript_id to protein_id dataframe #

#' Create Transcript ID to Protein ID Dataframe
#'
#' @description
#' This function parses a GTF annotation dataframe to generate a transcript ID to protein ID conversion dataframe.
#'
#' @param gtf_df A dataframe with GTF annotation information
#'
#' @return A dataframe with gene ID, transcript ID and protein ID information
#' @examples
#' transcript2protein_id_df <- create_transcript2protein_id_df(gtf_df)
#' @export

create_transcript2protein_id_df <- function(gtf_df){

  gtf_attributes <- gtf_df[gtf_df$feature == "CDS",]

  gene_ids <- c()
  gene_names <- c()
  transcript_ids <- c()
  protein_ids <- c()

  for (i in 1:nrow(gtf_attributes)){

    attribute_i <- gtf_attributes[i,]$attribute
    attribute_i <- trimws(unlist(strsplit(attribute_i, ";")))

    gene_id <- NA
    gene_name <- NA
    transcript_id <- NA
    protein_id <- NA
    genefound <- FALSE

    for (j in attribute_i){
      if (grepl("^gene_id ", j)){
        gene_id <- substr(j,9,nchar(j))
      }
      if (grepl('^gene ', j) && !genefound){
        gene_name <- substr(j,6,nchar(j))
        genefound <- TRUE
      }
      if (grepl("^transcript_id ", j)){
        transcript_id <- substr(j,15,nchar(j))
      }
      if (grepl("^protein_id ", j)){
        protein_id <- substr(j,12,nchar(j))
      }
    }

    gene_ids<- append(gene_ids, gene_id)
    gene_names <- append(gene_names, gene_name)
    transcript_ids <- append(transcript_ids, transcript_id)
    protein_ids <- append(protein_ids, protein_id)
  }

  transcript2protein_id_df <- data.frame(gene_ids, gene_names, transcript_ids, protein_ids)
  transcript2protein_id_df <- unique(transcript2protein_id_df)
  row.names(transcript2protein_id_df) <- NULL

  return(transcript2protein_id_df)
}

####################################################
# Converts a list of transcript ids to protein ids #

#' Transcript ID to Protein ID
#'
#' @description
#' Converts a list of Transcript IDs to Protein Ids.
#'
#' @param transcript_ids A vector of transcript IDs
#' @param transcript2protein_id_df A dataframe with gene ID, transcript ID and protein ID information
#'
#' @return A vector of protein IDs
#'
#' @examples
#' protein_ids <- transcript2protein_id(transcript_ids, transcript2protein_id_df)
#' @export

transcript2protein_id <- function(transcript_ids, transcript2protein_id_df){

  protein_ids <- c()

  for (i in transcript_ids){
    protein_ids <- append(protein_ids, transcript2protein_id_df[transcript2protein_id_df$transcript_ids == i,]$protein_ids)
  }

  return(protein_ids)
}


##############################################
# Converts a list of gene ids to protein ids #

#' Gene ID to Protein ID
#'
#' @description
#' Converts a list of Gene IDs and/or Gene Names to Protein Ids.
#'
#' @param transcript_ids A vector of gene IDs
#' @param transcript2protein_id_df A dataframe with gene ID, transcript ID and protein ID information
#'
#' @return A vector of protein IDs
#'
#' @examples
#' protein_ids <- gene2protein_id(genes, transcript2protein_id_df)
#' @export

gene2protein_id <- function(genes, transcript2protein_id_df) {

  protein_ids <- c()

  for (i in toupper(genes)) {

    # Check in gene_ids
    if (i %in% toupper(transcript2protein_id_df$gene_ids)) {
      cleaned_df <- transcript2protein_id_df[!is.na(transcript2protein_id_df$gene_ids), ]
      match <- cleaned_df[toupper(cleaned_df$gene_ids) == i, ]$protein_ids

      # Take the first match or handle as required
      protein_ids <- append(protein_ids, match[1])
      next
    }

    # Check in gene_names
    else if (i %in% toupper(transcript2protein_id_df$gene_names)) {
      cleaned_df <- transcript2protein_id_df[!is.na(transcript2protein_id_df$gene_names), ]
      match <- cleaned_df[toupper(cleaned_df$gene_names) == i, ]$protein_ids

      # Take the first match or handle as required
      protein_ids <- append(protein_ids, match[1])
      next
    }

    # If no match, append NA
    else {
      protein_ids <- append(protein_ids, NA)
    }
  }

  return(protein_ids)
}


###############################
# COG/KOG enrichment analysis #
###############################

########################################
# Loads COG/KOG annotation file #

#' Load COG/KOG Annotation
#'
#' @description
#' This function loads a COG/KOG (Clusters of Orthologous Genes) annotation into a dataframe.
#'
#' @param path Path to COG/KOG annotation
#'
#' @return A dataframe with COG/KOG annotation information
#'
#' @examples
#' kog_annotation_df <- load_kog_annotation(kog_annotation_path)
#' @export

load_kog_annotation <- function(path){

  annotation_data <- read.delim(path, header = TRUE, sep = "\t")

  if (dim(annotation_data)[2] != 6) {
    stop("KOG annotation file must have six columns: transcriptId,	proteinId,	kogid,	kogdefline,	kogClass	& kogGroup")
  }

  colnames(annotation_data) <- c("transcriptId",	"proteinId",	"kogid",	"kogdefline",	"kogClass",	"kogGroup")

  return(annotation_data)
}

########################
# Create COG/KOG model #

#' Create COG/KOG Model
#'
#' @description
#' This function generates a background frequency model for COG/KOG categories enrichment analysis.
#'
#' @param kog_annotation A dataframe with COG/KOG annotation information
#'
#' @return A dataframe with COG/KOG categories background frequencies
#'
#' @examples
#' kog_model_df <- create_kog_model(kog_annotation_df)
#' @export

create_kog_model <- function(kog_annotation){

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
              "UNANNOTATED",
              "ANNOTATED")

  sum_annotated_class <- c()

  for (group in groups[1:length(groups)]){
    kogclass_sum <- sum(grepl(group, kog_annotation$kogClass))
    sum_annotated_class <- append(sum_annotated_class, kogclass_sum)
  }

  kog_model <- data.frame(kogClass = groups, kogModel = sum_annotated_class)
  kog_model[kog_model$kogClass == "ANNOTATED",]$kogModel <- length(unique(kog_annotation$proteinId))
  kog_model[kog_model$kogClass == "CELLULAR PROCESSES AND SIGNALING",]$kogModel <- sum(kog_model$kogModel[2:10])
  kog_model[kog_model$kogClass == "INFORMATION STORAGE AND PROCESSING",]$kogModel <- sum(kog_model$kogModel[12:16])
  kog_model[kog_model$kogClass == "METABOLISM",]$kogModel <- sum(kog_model$kogModel[18:26])
  kog_model[kog_model$kogClass == "POORLY CHARACTERIZED",]$kogModel <- sum(kog_model$kogModel[28:29])

  return(kog_model)
}

###################################
# Create COG/KOG enrichment table #

#' Create COG/KOG Enrichment Dataframe
#'
#' @description
#' This function generates a dataframe containing background/sample frequencies and enrichment statistics per COG/KOG category.
#'
#' @param kog_model A dataframe with COG/KOG categories background frequencies
#' @param kog_annotation_df A dataframe with COG/KOG annotation information
#' @param protein_ids A vector of protein IDs
#' @param test A statistical method for testing the null of independence of rows and columns in a contingency table (default = "fisher")
#' @param p.adjust.method A statistical method to adjust p-values for multiple comparisons (default = "BH")
#'
#' @return A dataframe with COG/KOG categories enrichment statistics
#'
#' @examples
#' kog_enrichment_df <- kog_enrichment(kog_model_df, kog_annotation_df, protein_ids)
#' kog_enrichment_df <- kog_enrichment(kog_model_df, kog_annotation_df, protein_ids, test = "fisher", p.adjust.method = "BH")
#' kog_enrichment_df <- kog_enrichment(kog_model_df, kog_annotation_df, protein_ids, test = "chisq", p.adjust.method = "bonferroni")
#' kog_enrichment_df <- kog_enrichment(kog_model_df, kog_annotation_df, protein_ids, test = "chisq", p.adjust.method = "none")
#' @export

kog_enrichment <- function(kog_model, kog_annotation_df, protein_ids, test= "fisher", p.adjust.method= "BH"){

  kog_count_df <- kog_model
  kog_count_df$proteinCount <- 0
  protein_ids <- unique(protein_ids)

  for (i in 1:length(protein_ids)){
    if (protein_ids[i] %in% kog_annotation_df$proteinId){
      kog_classes_i <- kog_annotation_df[kog_annotation_df$proteinId == protein_ids[i],]$kogClass
      kog_count_df[kog_count_df$kogClass == "ANNOTATED",]$proteinCount <- kog_count_df[kog_count_df$kogClass == "ANNOTATED",]$proteinCount + 1

      for (j in 1:length(kog_classes_i)){
        kog_count_df[kog_count_df$kogClass == trimws(kog_classes_i[j]),]$proteinCount <- kog_count_df[kog_count_df$kogClass == trimws(kog_classes_i[j]),]$proteinCount + 1
      }
    }

    else{
      kog_count_df[kog_count_df$kogClass == "UNANNOTATED",]$proteinCount <- kog_count_df[kog_count_df$kogClass == "UNANNOTATED",]$proteinCount + 1
    }
  }

  kog_count_df[kog_count_df$kogClass == "CELLULAR PROCESSES AND SIGNALING",]$proteinCount <- sum(kog_count_df$proteinCount[2:10])
  kog_count_df[kog_count_df$kogClass =="INFORMATION STORAGE AND PROCESSING",]$proteinCount <- sum(kog_count_df$proteinCount[12:16])
  kog_count_df[kog_count_df$kogClass =="METABOLISM",]$proteinCount <- sum(kog_count_df$proteinCount[18:26])
  kog_count_df[kog_count_df$kogClass =="POORLY CHARACTERIZED",]$proteinCount <- sum(kog_count_df$proteinCount[28:29])

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

  return(enrichment_df)
}

####################################
# Generate COG/KOG enrichment plot #

#' Generate COG/KOG Enrichment Plot
#'
#' @description
#' This function generates a bar plot for visualization of COG/KOG categories enrichment statistics.
#'
#' @param enrichment_df A dataframe with COG/KOG enrichment statistics
#' @param significant A color assignment for statistically significant COG/KOG categories (default = "#0B775E")
#' @param plot_title A plot title (default = NA)
#'
#' @return A ggplot2 bar plot object
#'
#' @examples
#' generate_kog_plot(kog_enrichment_df)
#' generate_kog_plot(kog_enrichment_df, plot_title = "KOG Enrichment")
#' generate_kog_plot(kog_enrichment_df, significant = "blue")
#' @export

generate_kog_plot <- function(enrichment_df, significant = "#0B775E", plot_title = NA){

  fig_tab <- enrichment_df[1:29,]
  fig_tab$color <- "grey"

  for (i in 1:nrow(fig_tab)){
    if (fig_tab[i,]$adj.p.value <= 0.05 & fig_tab[i,]$enrichment > 1){
      fig_tab[i,]$color <- significant
    }
  }

  fig_tab$kogClass <- factor(fig_tab$kogClass, levels = rev(fig_tab$kogClass))

  p <- ggplot2::ggplot(data = fig_tab, ggplot2::aes(x = kogClass, y = enrichment)) + ggplot2::geom_bar(stat = "identity", width = 0.75, fill = fig_tab$color) +
    ggplot2::coord_flip() + ggplot2::xlab("") + ggplot2::ylab("Enrichment") + ggplot2::ylim(0, (max(fig_tab$enrichment, na.rm = TRUE) + 0.03)) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "#B40F20")

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

  return(p)
}


##########################
# GO enrichment analysis #
##########################

###################################
# Loads GO annotation file #

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
#' go_annotation_df <- load_go_annotation(go_annotation_path)
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
# Create GO model #

#' Create GO Model
#'
#' @description
#' This function generates a background frequency model for GO term enrichment analysis.
#'
#' @param go_annotation A dataframe with GO annotation information
#'
#' @return A list of dataframes with GO terms background frequencies per GO term type
#'
#' @examples
#' go_model_df <- create_go_model(go_annotation_df)
#' @export

create_go_model <- function(go_annotation){

  groups <- unique(go_annotation$goName)

  sum_annotated_goname <- c()
  go_term_ids <- c()
  go_term_types <- c()
  go_acc <- c()

  for (group in groups[1:length(groups)]){
    if (group %in% go_annotation$goName){
      goname_sum <- length(go_annotation[go_annotation$goName == group,]$proteinId)
      sum_annotated_goname <- append(sum_annotated_goname, goname_sum)
      go_term_ids <- append(go_term_ids, unique(go_annotation[go_annotation$goName == group,]$gotermId))
      go_term_types <- append(go_term_types, unique(go_annotation[go_annotation$goName == group,]$gotermType))
      go_acc <- append(go_acc, unique(go_annotation[go_annotation$goName == group,]$goAcc))
    }
    else {
      sum_annotated_goname <- append(sum_annotated_goname, NA)
      go_term_ids <- append(go_term_ids, NA)
      go_term_types <- append(go_term_types, NA)
      go_acc <- append(go_acc, NA)
    }
  }

  go_model <- data.frame(gotermId = go_term_ids, goName = groups, gotermType = go_term_types, goAcc = go_acc, goModel = sum_annotated_goname)
  annotated <- c(NA,	"ANNOTATED",	NA,	NA, 0)
  unannotated <- c(NA,	"UNANNOTATED",	NA,	NA, 0)

  biological_process <- go_model[go_model$gotermType == "biological_process",]
  biological_process <- rbind(biological_process, annotated)
  biological_process[biological_process$goName == "ANNOTATED",]$goModel <- length(unique(go_annotation[go_annotation$gotermType == "biological_process",]$proteinId))
  biological_process <- rbind(biological_process, unannotated)
  rownames(biological_process) <- NULL

  cellular_component <- go_model[go_model$gotermType == "cellular_component",]
  cellular_component <- rbind(cellular_component, annotated)
  cellular_component[cellular_component$goName == "ANNOTATED",]$goModel <- length(unique(go_annotation[go_annotation$gotermType == "cellular_component",]$proteinId))
  cellular_component <- rbind(cellular_component, unannotated)
  rownames(cellular_component) <- NULL

  molecular_function <- go_model[go_model$gotermType == "molecular_function",]
  molecular_function <- rbind(molecular_function, annotated)
  molecular_function[molecular_function$goName == "ANNOTATED",]$goModel <- length(unique(go_annotation[go_annotation$gotermType == "molecular_function",]$proteinId))
  molecular_function <- rbind(molecular_function, unannotated)
  rownames(molecular_function) <- NULL

  go_model <- list(biological_process = biological_process, cellular_component = cellular_component, molecular_function = molecular_function)

  return(go_model)
}

##############################
# Create GO enrichment table #

#' Create GO Enrichment Dataframe
#'
#' @description
#' This function generates a dataframe containing background/sample frequencies and enrichment statistics per GO term.
#'
#' @param go_model 	A list of dataframes with GO terms background frequencies per GO term type
#' @param go_annotation_df A dataframe with GO annotation information
#' @param protein_ids A vector of protein IDs
#' @param test A statistical method for testing the null of independence of rows and columns in a contingency table (default = "fisher")
#' @param p.adjust.method A statistical method to adjust p-values for multiple comparisons (default = "BH")
#'
#' @return A list of dataframes with GO terms enrichment statistics per GO term type
#'
#' @examples
#' go_enrichment_df <- go_enrichment(go_model_df, go_annotation_df, protein_ids)
#' go_enrichment_df <- go_enrichment(go_model_df, go_annotation_df, protein_ids, test = "fisher", p.adjust.method = "BH")
#' go_enrichment_df <- go_enrichment(go_model_df, go_annotation_df, protein_ids, test = "chisq", p.adjust.method = "bonferroni")
#' go_enrichment_df <- go_enrichment(go_model_df, go_annotation_df, protein_ids, test = "chisq", p.adjust.method = "none")
#' @export

go_enrichment <- function(go_model, go_annotation_df, protein_ids, test= "fisher", p.adjust.method= "BH"){

  go_count_df <- go_model
  go_groups <- c("biological_process", "cellular_component", "molecular_function")
  protein_ids <- unique(protein_ids)

  for (group in go_groups){
    go_count_df[[group]]$proteinCount <- 0

    for (i in 1:length(protein_ids)){
      if (protein_ids[i] %in% go_annotation_df[go_annotation_df$gotermType == group,]$proteinId){
        go_classes_i <- go_annotation_df[go_annotation_df$proteinId == protein_ids[i],]$goName
        go_count_df[[group]][go_count_df[[group]]$goName == "ANNOTATED",]$proteinCount <- go_count_df[[group]][go_count_df[[group]]$goName == "ANNOTATED",]$proteinCount + 1

        for (j in 1:length(go_classes_i)){
          go_count_df[[group]][go_count_df[[group]]$goName == trimws(go_classes_i[j]),]$proteinCount <- go_count_df[[group]][go_count_df[[group]]$goName == trimws(go_classes_i[j]),]$proteinCount + 1
        }
      }

      else{
        go_count_df[[group]][go_count_df[[group]]$goName == "UNANNOTATED",]$proteinCount <- go_count_df[[group]][go_count_df[[group]]$goName == "UNANNOTATED",]$proteinCount + 1
      }
    }

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
  row.names(enrichment_df[[group]]) <- NULL

  enrichment_df[[group]]$adj.p.value <- p.adjust(enrichment_df[[group]]$p.value, method = p.adjust.method)

  }

  return(enrichment_df)
}

###############################
# Generate GO enrichment plot #

#' Generate GO Terms Enrichment Plot
#'
#' @description
#' This function generates a bar plot for visualization of GO terms enrichment statistics.
#'
#'
#' @param enrichment_df A list of dataframes with GO terms enrichment statistics per GO term type
#' @param gotermType A GO term type ("biological_process", "molecular_function" or "cellular_component"; default = "biological_process")
#' @param n A number of top-ranking GO terms to include in the plot (default = 10)
#' @param significant A color assignment for statistically significant GO terms (default = "#E58601")
#' @param plot_title A plot title (default = NA)
#'
#' @return A ggplot2 bar plot object
#'
#' @examples
#' generate_go_plot(go_enrichment_df)
#' generate_go_plot(go_enrichment_df, gotermType = "biological_process", n = 15, plot_title = "GO Enrichment: Biological Process")
#' generate_go_plot(go_enrichment_df, gotermType = "molecular_function")
#' generate_go_plot(go_enrichment_df, gotermType = "cellular_component")
#' generate_go_plot(go_enrichment_df, gotermType = "biological_process", n = 10, significant = "blue")
#' @export

generate_go_plot <- function(enrichment_df, gotermType = "biological_process", n = 10, significant = "#E58601", plot_title = NA){

  fig_tab <- enrichment_df[[gotermType]][1:n,]
  fig_tab$color <- "grey"

  for (i in 1:n){
    if (fig_tab[i,]$adj.p.value <= 0.05 & fig_tab[i,]$enrichment > 1){
      fig_tab[i,]$color <- significant
    }
  }

  fig_tab$goName <- factor(fig_tab$goName, levels = rev(fig_tab$goName))

  if (gotermType == "biological_process"){
    go_term_type <- "Biological Process"
  }
  else if (gotermType == "molecular_function"){
    go_term_type <- "Molecular Function"
  }
  else if (gotermType == "cellular_component"){
    go_term_type <- "Cellular Component"
  }

  p <- ggplot2::ggplot(data = fig_tab, ggplot2::aes(x = goName, y = enrichment)) + ggplot2::geom_bar(stat = "identity", width = 0.75, fill = fig_tab$color) +
    ggplot2::coord_flip() + ggplot2::xlab("") + ggplot2::ylab("Enrichment") + ggplot2::ylim(0, (max(fig_tab$enrichment, na.rm = TRUE) + 0.03)) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "#B40F20")

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

  return(p)
}

##########################
# EC enrichment analysis #
##########################

# Loads EC annotation file

#' Load EC Annotation
#'
#' @description
#' This function loads a EC (Enzyme Commission) annotation into a dataframe.
#'
#'
#' @param path 	Path to EC annotation
#'
#' @return A dataframe with EC annotation information
#'
#' @examples
#' ec_annotation_df <- load_ec_annotation(ec_annotation_path)
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

###################
# Create EC model #

#' Create EC Model
#'
#' @description
#' This function generates a background frequency model for EC pathway enrichment analysis.
#'
#' @param ec_annotation A dataframe with EC annotation information
#'
#' @return A list of dataframes with EC pathway background frequencies per EC pathway type
#'
#' @examples
#' ec_model_df <- create_ec_model(ec_annotation_df)
#' @export

create_ec_model <- function(ec_annotation){

  firstup <- function(x) {
    x <- paste(substr(x, 1, 1), tolower(substr(x, 2, nchar(x))), sep = "")
    return(x)
  }

  pathway_types_names <- names(table(factor(ec_annotation$pathway_type)))

  # For summary df of Pathway Types & Classes
  annotated_classes <- c()
  sum_annotated_class <- c()

  for (i in pathway_types_names){
    annotated_classes <- append(annotated_classes, i)
    sum_annotated_class <- append(sum_annotated_class, length(table(ec_annotation[ec_annotation$pathway_type == i,]$proteinId)))

    for (j in names(table(factor(ec_annotation[ec_annotation$pathway_type == i,]$pathway_class)))){
      annotated_classes <- append(annotated_classes, firstup(j))
      sum_annotated_class <- append(sum_annotated_class, table(factor(ec_annotation$pathway_class))[[j]])
    }
  }

  annotated_classes <- append(annotated_classes, "ANNOTATED")
  sum_annotated_class <- append(sum_annotated_class, length(unique(ec_annotation[!is.na(ec_annotation$pathway_type), ]$proteinId)))
  annotated_classes <- append(annotated_classes, "UNANNOTATED")
  sum_annotated_class <- append(sum_annotated_class, NA)

  ec_model_summary <- data.frame(ecPathway = annotated_classes, ecModel = sum_annotated_class)

  # For detailed dfs of Pathways by Class
  annotated_pathways_metabolic <- c()
  sum_annotated_pathways_metabolic <- c()

  for (i in names(table(factor(ec_annotation[ec_annotation$pathway_type == "METABOLIC",]$pathway)))){
    annotated_pathways_metabolic <- append(annotated_pathways_metabolic, i)
    sum_annotated_pathways_metabolic <- append(sum_annotated_pathways_metabolic, table(factor(ec_annotation$pathway))[[i]])
  }

  annotated_pathways_metabolic <- append(annotated_pathways_metabolic, "ANNOTATED")
  sum_annotated_pathways_metabolic <- append(sum_annotated_pathways_metabolic, length(unique(ec_annotation[ec_annotation$pathway_type == "METABOLIC", ]$proteinId)))
  annotated_pathways_metabolic <- append(annotated_pathways_metabolic, "UNANNOTATED")
  sum_annotated_pathways_metabolic <- append(sum_annotated_pathways_metabolic, NA)

  ec_model_metabolic <- data.frame(ecPathway = annotated_pathways_metabolic, ecModel = sum_annotated_pathways_metabolic)

  annotated_pathways_regulatory <- c()
  sum_annotated_pathways_regulatory <- c()

  for (i in names(table(factor(ec_annotation[ec_annotation$pathway_type == "REGULATORY",]$pathway)))){
    annotated_pathways_regulatory <- append(annotated_pathways_regulatory, i)
    sum_annotated_pathways_regulatory <- append(sum_annotated_pathways_regulatory, table(factor(ec_annotation$pathway))[[i]])
  }

  annotated_pathways_regulatory <- append(annotated_pathways_regulatory, "ANNOTATED")
  sum_annotated_pathways_regulatory <- append(sum_annotated_pathways_regulatory, length(unique(ec_annotation[ec_annotation$pathway_type == "REGULATORY", ]$proteinId)))
  annotated_pathways_regulatory <- append(annotated_pathways_regulatory, "UNANNOTATED")
  sum_annotated_pathways_regulatory <- append(sum_annotated_pathways_regulatory, NA)

  ec_model_regulatory <- data.frame(ecPathway = annotated_pathways_regulatory, ecModel = sum_annotated_pathways_regulatory)

  ec_model <- list(summary = ec_model_summary, metabolic = ec_model_metabolic, regulatory = ec_model_regulatory)

  return(ec_model)
}


###############################
# Create EC enrichment table #

#' Create EC Enrichment Dataframe
#'
#' @description
#' This function generates a dataframe containing background/sample frequencies and enrichment statistics per EC pathway type.
#'
#' @param ec_model A list of dataframes with EC pathway background frequencies per EC pathway type
#' @param ec_annotation_df A dataframe with EC annotation information
#' @param protein_ids A vector of protein IDs
#' @param test A statistical method for testing the null of independence of rows and columns in a contingency table (default = "fisher")
#' @param p.adjust.method 	A statistical method to adjust p-values for multiple comparisons (default = "BH")
#'
#' @return A list of dataframes with EC pathway enrichment statistics per EC pathway type
#'
#' @examples
#' ec_enrichment_df <- ec_enrichment(ec_model_df, ec_annotation_df, protein_ids)
#' ec_enrichment_df <- ec_enrichment(ec_model_df, ec_annotation_df, protein_ids, test= "fisher", p.adjust.method= "BH")
#' ec_enrichment_df <- ec_enrichment(ec_model_df, ec_annotation_df, protein_ids, test= "chisq", p.adjust.method= "bonferonni")
#' ec_enrichment_df <- ec_enrichment(ec_model_df, ec_annotation_df, protein_ids, test= "chisq", p.adjust.method= "none")
#' @export

ec_enrichment <- function(ec_model, ec_annotation_df, protein_ids, test= "fisher", p.adjust.method= "BH"){

  ec_count_df <- ec_model
  ec_groups <- c("metabolic", "regulatory")
  protein_ids <- unique(protein_ids)

  # Counting proteins for summary table
  ec_count_df[["summary"]]$proteinCount <- 0
  ec_pathway_classes_i <- c()
  ec_pathway_types_i <- c()

  for (i in 1:length(protein_ids)){
    if (protein_ids[i] %in% ec_annotation_df[!is.na(ec_annotation_df$pathway_type),]$proteinId){
      ec_pathway_classes_i <- append(ec_pathway_classes_i, ec_annotation_df[ec_annotation_df$proteinId == protein_ids[i],]$pathway_class)
      ec_pathway_types_i <- append(ec_pathway_types_i, names(table(ec_annotation_df[ec_annotation_df$proteinId == protein_ids[i],]$pathway_type)))
      ec_count_df[["summary"]][ec_count_df[["summary"]]$ecPathway == "ANNOTATED",]$proteinCount <- ec_count_df[["summary"]][ec_count_df[["summary"]]$ecPathway == "ANNOTATED",]$proteinCount + 1
    }
    else{
      ec_count_df[["summary"]][ec_count_df[["summary"]]$ecPathway == "UNANNOTATED",]$proteinCount <- ec_count_df[["summary"]][ec_count_df[["summary"]]$ecPathway == "UNANNOTATED",]$proteinCount + 1
    }
  }

  for (i in names(table(factor(ec_pathway_classes_i)))){
    ec_count_df[["summary"]][ec_count_df[["summary"]]$ecPathway == i,]$proteinCount <- table(factor(ec_pathway_classes_i))[i]
  }

  for (i in names(table(factor(ec_pathway_types_i)))){
    ec_count_df[["summary"]][ec_count_df[["summary"]]$ecPathway == i,]$proteinCount <- table(factor(ec_pathway_types_i))[i]
  }

  # For METABOLIC and REGULATORY groups
  for (group in ec_groups){
    ec_count_df[[group]]$proteinCount <- 0

    ec_pathway_i <- c()

    for (i in 1:length(protein_ids)){
      if (protein_ids[i] %in% ec_annotation_df[!is.na(ec_annotation_df$pathway_type),]$proteinId){
        ec_pathway_i <- append(ec_pathway_i, ec_annotation_df[ec_annotation_df$proteinId == protein_ids[i],]$pathway)
        ec_count_df[[group]][ec_count_df[[group]]$ecPathway == "ANNOTATED",]$proteinCount <- ec_count_df[[group]][ec_count_df[[group]]$ecPathway == "ANNOTATED",]$proteinCount + 1
      }
      else{
        ec_count_df[[group]][ec_count_df[[group]]$ecPathway == "UNANNOTATED",]$proteinCount <- ec_count_df[[group]][ec_count_df[[group]]$ecPathway == "UNANNOTATED",]$proteinCount + 1
      }
    }

    for (i in names(table(factor(ec_pathway_i)))){
      if (i %in% ec_count_df[[group]]$ecPathway){
      ec_count_df[[group]][ec_count_df[[group]]$ecPathway == i,]$proteinCount <- table(factor(ec_pathway_i))[[i]]
    }
    }
  }

  enrichment_df <- ec_count_df
  ec_groups <- c("summary", "metabolic", "regulatory")

  for (group in ec_groups){
    annotated_model <- as.numeric(enrichment_df[[group]][enrichment_df[[group]]$ecPathway == "ANNOTATED",]$ecModel)
    annotated_proteins <- as.numeric(enrichment_df[[group]][enrichment_df[[group]]$ecPathway == "ANNOTATED",]$proteinCount)

    enrichment_df[[group]]$enrichment <- (as.numeric(enrichment_df[[group]]$proteinCount)/annotated_proteins)/(as.numeric(enrichment_df[[group]]$ecModel)/annotated_model)

    enrichment_df[[group]]$enrichment[is.infinite(enrichment_df[[group]]$enrichment)] <- NaN
    enrichment_df[[group]]$enrichment[enrichment_df[[group]]$ecPathway == "UNANNOTATED"] <- NaN
    enrichment_df[[group]]$enrichment[enrichment_df[[group]]$ecPathway == "ANNOTATED"] <- NaN

    stat_results <- c()

    for (i in enrichment_df[[group]]$ecPathway){
      if (is.nan(enrichment_df[[group]][enrichment_df[[group]]$ecPathway == i,]$enrichment)){
        stat_results <- append(stat_results, NaN)
        next
      }

      contingency_table <- matrix(c(as.numeric(enrichment_df[[group]][enrichment_df[[group]]$ecPathway == i,]$proteinCount),
                                    annotated_proteins - as.numeric(enrichment_df[[group]][enrichment_df[[group]]$ecPathway == i,]$proteinCount),
                                    as.numeric(enrichment_df[[group]][enrichment_df[[group]]$ecPathway == i,]$ecModel),
                                    annotated_model - as.numeric(enrichment_df[[group]][enrichment_df[[group]]$ecPathway == i,]$ecModel)),2,2)

      if (test == "fisher"){
        stat_results <- append(stat_results, fisher.test(contingency_table, alternative='greater')$p.value)
      }

      else if (test == "chisq"){
        stat_results <- append(stat_results, chisq.test(contingency_table)$p.value)
      }
    }

    enrichment_df[[group]]$p.value <- stat_results

    if (group %in% c("metabolic", "regulatory")){
    enrichment_df[[group]] <- enrichment_df[[group]][order(enrichment_df[[group]]$p.value),]
    row.names(enrichment_df[[group]]) <- NULL
    }

    enrichment_df[[group]]$adj.p.value <- p.adjust(enrichment_df[[group]]$p.value, method = p.adjust.method)

  }

  return(enrichment_df)
}


###############################
# Create EC enrichment plot #

#' Generate EC Pathway Enrichment Plot
#'
#' @description
#' This function generates a bar plot for visualization of EC pathways enrichment statistics.
#'
#' @param enrichment_df A list of dataframes with EC pathway enrichment statistics per EC pathway type
#' @param pathwayType A EC pathway type ("summary", "regulatory" or "metabolic"; default = "summary")
#' @param n A number of top-ranking EC pathways to include in the plot (default = 10)
#' @param significant A color assignment for statistically significant GO terms (default = "#9986A5")
#' @param plot_title A plot title (default = NA)
#'
#' @return A ggplot2 bar plot object
#'
#' @examples
#' generate_ec_plot(ec_enrichment_df)
#' generate_ec_plot(ec_enrichment_df, pathwayType = "summary", plot_title = "EC Enrichment: Summary")
#' generate_ec_plot(ec_enrichment_df, pathwayType = "metabolic", significant = "blue")
#' generate_ec_plot(ec_enrichment_df, pathwayType = "regulatory", n = 15)
#' @export

generate_ec_plot <- function(enrichment_df, pathwayType = "summary", n = 10, significant = "#9986A5", plot_title = NA){

  if (pathwayType == "summary"){
    n <- nrow(enrichment_df[[pathwayType]]) - 2
  }
  else if (n > nrow(enrichment_df[[pathwayType]])){
    n <- nrow(enrichment_df[[pathwayType]]) - 2
  }

  fig_tab <- enrichment_df[[pathwayType]][1:n,]
  fig_tab$color <- "grey"

  for (i in 1:n){
    if (fig_tab[i,]$adj.p.value <= 0.05 & fig_tab[i,]$enrichment > 1){
      fig_tab[i,]$color <- significant
    }
  }

  fig_tab$ecPathway <- factor(fig_tab$ecPathway, levels = rev(fig_tab$ecPathway))

  if (pathwayType == "summary"){
    ec_type <- "Summary"
  }
  else if (pathwayType == "metabolic"){
    ec_type <- "Metabolic"
  }
  else if (pathwayType == "regulatory"){
    ec_type <- "Regulatory"
  }

  p <- ggplot2::ggplot(data = fig_tab, ggplot2::aes(x = ecPathway, y = enrichment)) + ggplot2::geom_bar(stat = "identity", width = 0.75, fill = fig_tab$color) +
    ggplot2::coord_flip() + ggplot2::xlab("") + ggplot2::ylab("Enrichment") + ggplot2::ylim(0, (max(fig_tab$enrichment, na.rm = TRUE) + 0.03)) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "#B40F20")

  if (!is.na(plot_title)){
    p <- p + ggplot2::ggtitle(plot_title)
  }

  for (i in fig_tab$ecPathway){
    if (fig_tab$adj.p.value[fig_tab$ecPathway == i] <= 0.05 & fig_tab$enrichment[fig_tab$ecPathway == i] > 1){
      if (fig_tab$adj.p.value[fig_tab$ecPathway == i] <= 0.001){
        p <- p + ggplot2::geom_text(x=dim(fig_tab)[1]-which(fig_tab$ecPathway == i) + 1,
                           y=fig_tab$enrichment[fig_tab$ecPathway == i] + 0.05,
                           label="***", angle = 90)
      }

      else if (fig_tab$adj.p.value[fig_tab$ecPathway == i] <= 0.01){
        p <- p + ggplot2::geom_text(x=dim(fig_tab)[1]-which(fig_tab$ecPathway == i) + 1,
                           y=fig_tab$enrichment[fig_tab$ecPathway == i] + 0.05,
                           label="**", angle = 90)
      }

      else {
        p <- p + ggplot2::geom_text(x=dim(fig_tab)[1]-which(fig_tab$ecPathway == i) + 1,
                           y=fig_tab$enrichment[fig_tab$ecPathway == i] + 0.05,
                           label="*" , , angle = 90)
      }

    }
  }

  return(p)
}


######################################
# Loading a eggNOG-mapper Annotation #
######################################

#' Load eggNOG-mapper Annotation
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
#' eggnog_annotation_df <- load_eggnog_annotation(eggnog_annotation_path)

load_eggnog_annotation <- function(path){

  annotation_data <- read.delim(path, header = FALSE, sep = "\t")

  if (dim(annotation_data)[2] != 21) {
    stop("Write Error Here")
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
#' @param eggnog_annotation A dataframe with eggNOG-mapper annotation information
#'
#' @return A dataframe with COG/KOG categories background frequencies
#' @export
#'
#' @examples
#' kog_model_df <- create_kog_model_eggnog(eggnog_annotation_df)

create_kog_model_eggnog <- function(eggnog_annotation){

  cog_classes_df <- cog_classes_df

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
              "UNANNOTATED",
              "ANNOTATED")

  kog_model <- data.frame(kogClass = groups, kogModel = 0)

  classes_model <- paste(eggnog_annotation$COG_category, collapse = '')
  n_proteins <- length(eggnog_annotation$query)
  unannotated <- length(gregexpr("-", classes_model)[[1]][gregexpr("-", classes_model)[[1]] > 0])

  for (i in cog_classes_df$code){
    kog_model[kog_model$kogClass == cog_classes_df[cog_classes_df$code == i,]$class,]$kogModel <- length(gregexpr(i, classes_model)[[1]][gregexpr(i, classes_model)[[1]] > 0])
  }

  kog_model[kog_model$kogClass == "UNANNOTATED",]$kogModel <- unannotated
  kog_model[kog_model$kogClass == "ANNOTATED",]$kogModel <- n_proteins - unannotated
  kog_model[kog_model$kogClass == "CELLULAR PROCESSES AND SIGNALING",]$kogModel <- sum(kog_model$kogModel[2:10])
  kog_model[kog_model$kogClass =="INFORMATION STORAGE AND PROCESSING",]$kogModel <- sum(kog_model$kogModel[12:16])
  kog_model[kog_model$kogClass =="METABOLISM",]$kogModel <- sum(kog_model$kogModel[18:26])
  kog_model[kog_model$kogClass =="POORLY CHARACTERIZED",]$kogModel <- sum(kog_model$kogModel[28:29])

  return(kog_model)
}


###################################
# Create COG/KOG enrichment table #

#' Create COG/KOG Enrichment Dataframe from eggNOG-mapper annotation
#'
#' @param kog_model A dataframe with COG/KOG categories background frequencies
#' @param eggnog_annotation_df A dataframe with eggNOG-mapper annotation information
#' @param protein_ids A vector of protein IDs
#' @param test A statistical method for testing the null of independence of rows and columns in a contingency table (default = "fisher")
#' @param p.adjust.method A statistical method to adjust p-values for multiple comparisons (default = "BH")
#'
#' @return A dataframe with COG/KOG categories enrichment statistics
#' @export
#'
#' @examples
#' kog_enrichment_df <- kog_enrichment_eggnog(kog_model_df, eggnog_annotation_df, protein_ids)
#' kog_enrichment_df <- kog_enrichment_eggnog(kog_model_df, eggnog_annotation_df, protein_ids, test = "fisher", p.adjust.method = "BH")
#' kog_enrichment_df <- kog_enrichment_eggnog(kog_model_df, eggnog_annotation_df, protein_ids, test = "chisq", p.adjust.method = "bonferroni")
#' kog_enrichment_df <- kog_enrichment_eggnog(kog_model_df, eggnog_annotation_df, protein_ids, test = "chisq", p.adjust.method = "none")

kog_enrichment_eggnog <- function(kog_model, eggnog_annotation_df, protein_ids, test= "fisher", p.adjust.method= "BH"){

  cog_classes_df <- cog_classes_df

  kog_count_df <- kog_model
  kog_count_df$proteinCount <- 0

  protein_ids <- unique(protein_ids)
  protein_classes <- eggnog_annotation_df[eggnog_annotation_df$query %in% protein_ids,]$COG_category

  n_proteins <- length(protein_ids)
  n_annotated <- length(protein_classes)


  for (i in names(table(factor(protein_classes)))){

    if (i == "-"){
      kog_count_df[kog_count_df$kogClass == "UNANNOTATED",]$proteinCount <- table(factor(protein_classes))[[i]] + (n_proteins - n_annotated)
    }

    else {
    kog_count_df[kog_count_df$kogClass == cog_classes_df[cog_classes_df$code == i,]$class,]$proteinCount <- table(factor(protein_classes))[[i]]
    }
  }

  kog_count_df[kog_count_df$kogClass == "CELLULAR PROCESSES AND SIGNALING",]$proteinCount <- sum(kog_count_df$proteinCount[2:10])
  kog_count_df[kog_count_df$kogClass =="INFORMATION STORAGE AND PROCESSING",]$proteinCount <- sum(kog_count_df$proteinCount[12:16])
  kog_count_df[kog_count_df$kogClass =="METABOLISM",]$proteinCount <- sum(kog_count_df$proteinCount[18:26])
  kog_count_df[kog_count_df$kogClass =="POORLY CHARACTERIZED",]$proteinCount <- sum(kog_count_df$proteinCount[28:29])
  kog_count_df[kog_count_df$kogClass =="ANNOTATED",]$proteinCount <- n_proteins - kog_count_df[kog_count_df$kogClass =="UNANNOTATED",]$proteinCount

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

  return(enrichment_df)
}

################################################################################
