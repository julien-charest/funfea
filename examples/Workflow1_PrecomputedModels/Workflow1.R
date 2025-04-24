###################################
# FunFEA: Example Workflow #1     #
# Working with Precomputed Models #
###################################

library("funfea")

##################
# Preparing data #
##################

# Loading Example Gene IDs (from Derntl et al., 2017; https://doi.org/10.1073/pnas.1609348114)
gene_ids <- read.csv("examples/Workflow1_PrecomputedModels/gene_ids.csv")$gene_ids

# Fetch Precomputed Models for Trichoderma reesei (QM6a)
qm6a_models <- fetch_models("Trichoderma reesei (QM6a)")

# Convert Gene IDs to Protein IDs
protein_ids <- gene2protein_id(gene_ids, qm6a_models$transcript2protein_id)


######################################
# Performing KOG Enrichment Analysis #
######################################

# Computing KOG Enrichment
kog_enrichment_df <- kog_enrichment(qm6a_models$kog, protein_ids, test = "fisher", p.adjust.method = "BH")

# Generating Lollipop Plot of KOG Enrichment
generate_kog_plot(kog_enrichment_df, plot_type = "lollipop")


#####################################
# Performing GO Enrichment Analysis #
#####################################

# Computing GO Enrichment
go_enrichment_df <- go_enrichment(qm6a_models$go, protein_ids, test = "fisher", p.adjust.method = "BH")

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
kegg_enrichment_df <- kegg_enrichment(qm6a_models$kegg, protein_ids, test= "fisher", p.adjust.method= "BH")

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
