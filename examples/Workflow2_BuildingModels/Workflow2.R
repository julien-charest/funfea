###########################################
# FunFEA: Example Workflow #2             #
# Building Models From Public Annotations #
###########################################

library("funfea")

##############################################################
# Building Transcript IDs To Protein ID Conversion Dataframe #
##############################################################

# Loading GTF/GFF/GFF3 Annotation (https://mycocosm.jgi.doe.gov/Trire2/Trire2.home.html)
gtf_df <- load_gtf_annotation("examples/Workflow2_BuildingModels/TreeseiV2_FilteredModelsv2.0.gtf")

# Building Transcript IDs To Protein ID Conversion Dataframe
transcript2protein_id_df <- create_transcript2protein_id_df(gtf_df)


##############################
# Building Enrichment Models #
##############################

#############
# KOG Model #

# Loading KOG Annotation (https://mycocosm.jgi.doe.gov/Trire2/Trire2.home.html)
kog_annotation_df <- load_kog_annotation("examples/Workflow2_BuildingModels/TreeseiV2_koginfo_FilteredModelsv2.0.tab")

# Building KOG Model
kog_model_df <- create_kog_model(kog_annotation_df)

############
# GO Model #

# Loading GO Annotation (https://mycocosm.jgi.doe.gov/Trire2/Trire2.home.html)
go_annotation_df <- load_go_annotation("examples/Workflow2_BuildingModels/TreeseiV2_goinfo_FilteredModelsv2.0.tab")

# Building GO Model
go_model_df <- create_go_model(go_annotation_df)

# Isolating GO Model Dataframes
go_model_biological_df <- go_model_df$biological_process
go_model_molecular_df <- go_model_df$molecular_function
go_model_cellular_df <- go_model_df$cellular_component

############################################
# KEGG Pathway Model (from JGI Annotation) #

# Loading EC Annotation (https://mycocosm.jgi.doe.gov/Trire2/Trire2.home.html)
ec_annotation_df <- load_ec_annotation("examples/Workflow2_BuildingModels/TreeseiV2_ecpathwayinfo_FilteredModelsv2.0.tab")

# Building KEGG Pathway Model from EC Annotation
kegg_model_df <- create_kegg_model_from_ec(ec_annotation_df)

#############################################
# KEGG Pathway Model (from KEGG Annotation) #

# Loading KEGG Pathway Annotation (https://rest.kegg.jp/link/ko/tre)
kegg_annotation_df  <- load_kegg_annotation("examples/Workflow2_BuildingModels/TreeseiV2_KEGG.txt")

# Formatting Protein IDs To Match Annotation (If Needed)
kegg_annotation_df$proteinId <- gsub("TRIREDRAFT_", "", kegg_annotation_df$proteinId)

# Building KEGG Pathway Model from KEGG Annotation
kegg_model_df <- create_kegg_model(kegg_annotation_df)


##################
# Preparing data #
##################

# Loading Example Gene IDs (from Derntl et al., 2017; https://doi.org/10.1073/pnas.1609348114)
gene_ids <- read.csv("examples/Workflow2_BuildingModels/gene_ids.csv")$gene_ids

# Convert Gene IDs to Protein IDs
protein_ids <- gene2protein_id(gene_ids, transcript2protein_id_df)


######################################
# Performing KOG Enrichment Analysis #
######################################

# Computing KOG Enrichment
kog_enrichment_df <- kog_enrichment(kog_model_df, protein_ids, test = "fisher", p.adjust.method = "BH")

# Generating Lollipop Plot of KOG Enrichment
generate_kog_plot(kog_enrichment_df, plot_type = "lollipop")


#####################################
# Performing GO Enrichment Analysis #
#####################################

# Computing GO Enrichment
go_enrichment_df <- go_enrichment(go_model_df, protein_ids, test = "fisher", p.adjust.method = "BH")

# Extracting GO Enrichment Table for Biological Process
go_enrichment_biological_df <- go_enrichment_df$biological_process

# Generating Lollipop Plot for GO Enrichment (Biological Process)
generate_go_plot(go_enrichment_df, gotermType = "biological_process", n = 10, plot_type = "lollipop")

# Extracting GO Enrichment Table for Molecular Function
go_enrichment_molecular_df <- go_enrichment_df$molecular_function

# Generating Lollipop Plot for GO Enrichment (Molecular Function)
generate_go_plot(go_enrichment_df, gotermType = "molecular_function", n = 10, plot_type = "lollipop")

# Extracting GO Enrichment Table for Cellular Component
go_enrichment_cellular_df <- go_enrichment_df$cellular_component

# Generating Lollipop Plot for GO Enrichment (Cellular Component)
generate_go_plot(go_enrichment_df, gotermType = "cellular_component", n = 10, plot_type = "lollipop")


###############################################
# Performing KEGG Pathway Enrichment Analysis #
###############################################

# Computing KEGG Pathway Enrichment
kegg_enrichment_df <- kegg_enrichment(kegg_model_df, protein_ids, test= "fisher", p.adjust.method= "BH")

# Extracting KEGG Pathway Enrichment Table for Pathway Type
kegg_enrichment_type_df <- kegg_enrichment_df$pathway_type

# Generating Lollipop Plot for KEGG Pathway (Pathway Type)
generate_kegg_plot(kegg_enrichment_df, model = "pathway_type", n = 10, plot_type = "lollipop")

# Extracting KEGG Pathway Enrichment Table for Pathway Class
kegg_enrichment_class_df <- kegg_enrichment_df$pathway_class

# Generating Lollipop Plot for KEGG Pathway (Pathway Class)
generate_kegg_plot(kegg_enrichment_df, model = "pathway_class", n = 10, plot_type = "lollipop")

# Extracting KEGG Pathway Enrichment Table for Pathway Name
kegg_enrichment_name_df <- kegg_enrichment_df$pathway_name

# Generating Lollipop Plot for KEGG Pathway (Pathway Name)
generate_kegg_plot(kegg_enrichment_df, model = "pathway_name", n = 10, plot_type = "lollipop")

# Extracting KEGG Pathway Enrichment Table for Enzyme Definition
kegg_enrichment_name_df <- kegg_enrichment_df$definition

# Generating Lollipop Plot for KEGG Pathway (Enzyme Definition)
generate_kegg_plot(kegg_enrichment_df, model = "definition", n = 10, plot_type = "lollipop")
