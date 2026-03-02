# NOTE: Run 01_batch_correction.R before this

# Install and load libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")

library(biomaRt)
library(dplyr)

xcell_input_dir <- "./results/xCell_input/"
dir.create(xcell_input_dir,
           showWarnings = FALSE,
           recursive = TRUE)

# TPM function
tpm <- function(counts, gene_length) {
  gene_length_kb <- gene_length / 1000
  cpk <- counts / gene_length_kb
  scaling <- colSums(cpk)
  t(t(cpk) / scaling * 1e6)
}

# Import all raw counts datasets
## GSE184316
gse184316_chp_ctrl_raw <- read.csv("./data/GSE184316/GSE184316_raw_counts_HP_Control.csv",
                         header = TRUE,
                         row.names = 1)

gse184316_chp_ipf_raw <- read.csv("./data/GSE184316/GSE184316_raw_counts_HP_IPF.csv",
                                  header = TRUE,
                                  row.names = 1)

gse184316_ipf_ctrl_raw <- read.csv("./data/GSE184316/GSE184316_raw_counts_IPF_Control.csv",
                                   header = TRUE,
                                   row.names = 1)

## GSE150910
gse150910_chp_ctrl_raw <- read.csv("./data/GSE150910/GSE150910_raw_counts_HP_Control.csv",
                                   header = TRUE,
                                   row.names = 1)

gse150910_chp_ipf_raw <- read.csv("./data/GSE150910/GSE150910_raw_counts_HP_IPF.csv",
                                  header = TRUE,
                                  row.names = 1) 

gse150910_ipf_ctrl_raw <- read.csv("./data/GSE150910/GSE150910_raw_counts_IPF_Control.csv",
                                   header = TRUE,
                                   row.names = 1)

message("All raw counts are successfully exported!")

##############################################################################
# Get transcript lengths for each gene
# Each gene can have multiple transcripts (isoforms), hence only the transcript 
# with the highest length is kept for TPM normalization

# We first start with dataset GSE184316

# Check if all row names are in same order (all must return TRUE)
all(rownames(gse184316_chp_ctrl_raw) == rownames(gse184316_chp_ipf_raw)) #Output - TRUE
all(rownames(gse184316_chp_ctrl_raw) == rownames(gse184316_ipf_ctrl_raw)) #Output - TRUE
all(rownames(gse184316_chp_ipf_raw) == rownames(gse184316_ipf_ctrl_raw)) #Output - TRUE

# Get genes
genes <- rownames(gse184316_chp_ctrl_raw)

# Connect to Ensembl
mart <- useEnsembl("genes", dataset = "hsapiens_gene_ensembl")

# Annotation using the genes vector
annot <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id", 
                 "ensembl_transcript_id", "transcript_length"),
  filters = "hgnc_symbol",
  values = genes,
  mart = mart
)

#Collapse transcripts to length per gene
gene_lengths <- annot %>%
  group_by(hgnc_symbol) %>%
  summarise(gene_length = max(transcript_length, na.rm = TRUE))

#match gene lengths back to our counts matrix
gene_length_vec <- gene_lengths$gene_length[
  match(genes, gene_lengths$hgnc_symbol)
]

#Check NA values
sum(is.na(gene_length_vec))

#Separate the symbols which has no corresponding gene length values
# eg: pseudogenes, incomplete locus transcripts, miRNA, snoRNA etc. 
missing <- setdiff(genes, gene_lengths$hgnc_symbol)
length(missing)
head(missing)

# Remove them from all the counts matrices to prevent NA
keep <- !genes %in% missing

gse184316_chp_ctrl_filt <- gse184316_chp_ctrl_raw[keep, ]
gse184316_chp_ipf_filt <- gse184316_chp_ipf_raw[keep, ]
gse184316_ipf_ctrl_filt <- gse184316_ipf_ctrl_raw[keep, ]

gene_length_vec_filtered <- gene_length_vec[keep]

#############################################################################
#perform TPM
tpm_gse184316_chp_ctrl <- tpm(gse184316_chp_ctrl_filt, gene_length_vec_filtered)
write.table(tpm_gse184316_chp_ctrl,
            file = paste0(xcell_input_dir, "TPM_GSE184316_HP_Ctrls.txt"),
            sep = '\t',
            quote = FALSE,
            col.names = NA)


tpm_gse184316_chp_ipf <- tpm(gse184316_chp_ipf_filt, gene_length_vec_filtered)
write.table(tpm_gse184316_chp_ipf,
            file = paste0(xcell_input_dir, "TPM_GSE184316_HP_IPF.txt"),
            sep = '\t',
            quote = FALSE,
            col.names = NA)

tpm_gse184316_ipf_ctrl <- tpm(gse184316_ipf_ctrl_filt, gene_length_vec_filtered)
write.table(tpm_gse184316_ipf_ctrl,
            file = paste0(xcell_input_dir, "TPM_GSE184316_IPF_Ctrls.txt"),
            sep = '\t',
            quote = FALSE,
            col.names = NA)


###########################################
###########################################

# We then move to dataset GSE150910

# Check if all row names are in same order (all must return TRUE)
all(rownames(gse150910_chp_ctrl_raw) == rownames(gse150910_chp_ipf_raw)) #Output - TRUE
all(rownames(gse150910_chp_ctrl_raw) == rownames(gse150910_ipf_ctrl_raw)) #Output - TRUE
all(rownames(gse150910_chp_ipf_raw) == rownames(gse150910_ipf_ctrl_raw)) #Output - TRUE

# Get genes
genes <- rownames(gse150910_chp_ctrl_raw)

# Connect to Ensembl
mart <- useEnsembl("genes", dataset = "hsapiens_gene_ensembl")

# Annotation using the genes vector
annot <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id", 
                 "ensembl_transcript_id", "transcript_length"),
  filters = "hgnc_symbol",
  values = genes,
  mart = mart
)

#Collapse transcripts to length per gene
gene_lengths <- annot %>%
  group_by(hgnc_symbol) %>%
  summarise(gene_length = max(transcript_length, na.rm = TRUE))

#match gene lengths back to our counts matrix
gene_length_vec <- gene_lengths$gene_length[
  match(genes, gene_lengths$hgnc_symbol)
]

#Check NA values
sum(is.na(gene_length_vec))

#Separate the symbols which has no corresponding gene length values
# eg: pseudogenes, incomplete locus transcripts, miRNA, snoRNA etc. 
missing <- setdiff(genes, gene_lengths$hgnc_symbol)
length(missing)
head(missing)

# Remove them from all the counts matrices to prevent NA
keep <- !genes %in% missing

gse150910_chp_ctrl_filt <- gse150910_chp_ctrl_raw[keep, ]
gse150910_chp_ipf_filt <- gse150910_chp_ipf_raw[keep, ]
gse150910_ipf_ctrl_filt <- gse150910_ipf_ctrl_raw[keep, ]

gene_length_vec_filtered <- gene_length_vec[keep]

#############################################################################
#perform TPM
tpm_gse150910_chp_ctrl <- tpm(gse150910_chp_ctrl_filt, gene_length_vec_filtered)
write.table(tpm_gse150910_chp_ctrl,
            file = paste0(xcell_input_dir, "TPM_GSE150910_HP_Ctrls.txt"),
            sep = '\t',
            quote = FALSE,
            col.names = NA)


tpm_gse150910_chp_ipf <- tpm(gse150910_chp_ipf_filt, gene_length_vec_filtered)
write.table(tpm_gse150910_chp_ipf,
            file = paste0(xcell_input_dir, "TPM_GSE150910_HP_IPF.txt"),
            sep = '\t',
            quote = FALSE,
            col.names = NA)

tpm_gse150910_ipf_ctrl <- tpm(gse150910_ipf_ctrl_filt, gene_length_vec_filtered)
write.table(tpm_gse150910_ipf_ctrl,
            file = paste0(xcell_input_dir, "TPM_GSE150910_IPF_Ctrls.txt"),
            sep = '\t',
            quote = FALSE,
            col.names = NA)