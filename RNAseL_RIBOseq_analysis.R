# RIBOseq analysis of WT and RNAseL-KO samples induced with 2-5A/polyI:C
# Agnes Karasik, Hernan Lorenzi, Nicholas Guydosh
# 7/17/2023

# Contact:** Nicholas Guydosh
# Contact email: nicholas.guydosh@nih.gov

# ------------- Load libraries -------------

suppressMessages(library("org.Hs.eg.db"))
suppressMessages(library("pheatmap"))
suppressMessages(library("ggpubr"))
suppressMessages(library("DESeq2"))
suppressMessages(library("ggsci"))
suppressMessages(library("RColorBrewer"))

# ------------- Required functions -------------

# Function to remove all-zero rows (by adjusting min_total_count the function can filter out rows based on total counts other than 0)
remove_all_zero_rows <- function(df, min_total_count = 0){
  df <- df[rowSums(df) > min_total_count,]
  return(df)
}

# Function to replace gene symbols with Ensembl IDs
replace_annotation_keys <- function(key_list, from, to){
  if(from %in% keytypes(org.Hs.eg.db) & to  %in% keytypes(org.Hs.eg.db)){
    if(typeof(key_list) == "S4"){
      ensemble_ids <- rownames(key_list)
      symbols <- mapIds(org.Hs.eg.db, keys = ensemble_ids, column = to, keytype = from, multiVals = "first")
      symbols <- symbols[!is.na(symbols)]
      to_name <- rownames(key_list) %in% names(symbols)
      key_list@rownames[to_name] <- as.vector(symbols)
      return(key_list)
    } else {
      to_ids <- mapIds(org.Hs.eg.db, keys = key_list, column = to, keytype = from, multiVals = "first")
      to_ids <- to_ids[!is.na(to_ids)]
      to_name <- key_list %in% names(to_ids)
      key_list[to_name] <- as.vector(to_ids)
      return(key_list)
    }
  } else {
    print("Annotation keys are not keytypes")
    stop()
  }
}

# Function to generate MA plots
make_MA_plots <- function(DESEQ_result, file_prefix = 'MAplot'){
  plotMA(DESEQ_result, ylim=c(-10,10), main = file_prefix)
  pdf(file = paste0("./Plots/",file_prefix,"_MAplot.pdf"))
  plotMA(DESEQ_result, ylim=c(-10,10))
  dev.off()
}

# Function to calculate mean DESeq2-normalized read counts per gene
get_mean_norm_counts <- function(my_dds, genotype = 'no_genotype', feature_type = 'no_feature'){
  
  mean_norm_counts.wt.polyIC <- apply(counts(my_dds[ , my_dds$Inducer %in% c("None","polyIC")]), 1, function(x) mean(x))
  write.table(x = as.data.frame(mean_norm_counts.wt.polyIC), 
              file = paste0('./Mean_norm_counts/mean_norm_counts.',
                            genotype,
                            '.polyIC.',
                            feature_type,
                            '.txt'), 
              sep = "\t" )
  
  mean_norm_counts.wt.25A <- apply(counts(my_dds[ , my_dds$Inducer %in% c("None","25A")]), 1, function(x) mean(x))
  write.table(x = as.data.frame(mean_norm_counts.wt.25A), 
              file = paste0('./Mean_norm_counts/mean_norm_counts.'
                            ,genotype,
                            '.25A.',
                            feature_type,
                            '.txt'), 
              sep = "\t" )
}

# Function for processing and saving DE result tables
sort_and_write_res_table <- function(result_table, file_name){
  # Sort genes by (padj)
  result_table_sorted <- result_table[order(result_table$padj, decreasing = FALSE),]
  # Add gene symbols
  gene_list <- rownames(result_table_sorted)
  symbol_list <- ensembl_to_symbol$gene_name[match(gene_list, ensembl_to_symbol$Ensembl_ID)]
  df <-as.data.frame(cbind(result_table_sorted, Gene_name = symbol_list))
  
  # Write sorted table to file
  write.table(df, file = paste0("./DE/",file_name,".txt"), 
              sep = "\t", col.names=NA)
  return(result_table_sorted)
}

# Function to calculate mean DESeq2-normalized read counts per gene
get_mean_norm_counts <- function(my_dds, genotype = 'no_genotype', feature_type = 'no_feature'){
  
  mean_norm_counts.wt.polyIC <- apply(counts(my_dds[ , my_dds$Inducer %in% c("None","polyIC")]), 1, function(x) mean(x))
  write.table(x = as.data.frame(mean_norm_counts.wt.polyIC), 
              file = paste0('./Mean_norm_counts/mean_norm_counts.',
                            genotype,
                            '.polyIC.',
                            feature_type,
                            '.txt'), 
              sep = "\t" )
  
  mean_norm_counts.wt.25A <- apply(counts(my_dds[ , my_dds$Inducer %in% c("None","25A")]), 1, function(x) mean(x))
  write.table(x = as.data.frame(mean_norm_counts.wt.25A), 
              file = paste0('./Mean_norm_counts/mean_norm_counts.'
                            ,genotype,
                            '.25A.',
                            feature_type,
                            '.txt'), 
              sep = "\t" )
}

# Define function for processing and saving result tables
sort_and_write_res_table <- function(result_table, file_name){
  # Sort genes by (padj)
  result_table_sorted <- result_table[order(result_table$padj, decreasing = FALSE),]
  # Add gene symbols
  gene_list <- rownames(result_table_sorted)
  symbol_list <- ensembl_to_symbol$gene_name[match(gene_list, ensembl_to_symbol$Ensembl_ID)]
  df <-as.data.frame(cbind(result_table_sorted, Gene_name = symbol_list))
  
  # Write sorted table to file
  write.table(df, file = paste0("./DE/",file_name,".txt"), 
              sep = "\t", col.names=NA)
  return(result_table_sorted)
}


# ------------- Set filtering steps -------------

# Remove genes with 0 counts across all samples
remove_genes_with_zero_counts <- TRUE

# Filter reads with less than N reads
N_reads <- 20

# ------------- Load data - RIBOseq -------------

all <- read.csv("./data/cdsfoot_round.csv", row.names = 1) 

# Keep table with Ensemble IDs and gene Symbols
gene_symbols <- replace_annotation_keys(key_list = rownames(all), from = "ENSEMBL", to = "SYMBOL")
ensembl_to_symbol <- as.data.frame(cbind("Ensembl_ID" = rownames(all), "gene_name" = gene_symbols), row.names = 1)

# Load metadata
metadata <- read.delim2("./data/Metadata_footprint.txt", sep = "\t", row.names = 1, header = T)

# Sort tables so metadata and read counts match order
metadata<-  metadata[match(colnames(all), rownames(metadata)), ]

# Add total read counts and sample id columns to metadata
metadata <- cbind(metadata, Read_counts =colSums(all), Sample_id = rownames(metadata))

if (remove_genes_with_zero_counts){
  all <- remove_all_zero_rows(all, min_total_count = 0)
}

# ------------- Analysis of expression data using DESeq2 - RIBOseq -------------

# Convert metadata to factors
for (variable in c("Read_length","Sequencing_pool","Date_1st_submitted","Species","Strain_name","Genotype","Media","Comments","Group","Colection_time","X25A","PolyIC","Inducer","Treatment","Sample_id")){
  metadata[,variable] <- as.factor(metadata[,variable])  
}

# Subset metadata and count tables by Data
meta_one_RS <- metadata
all_one_RS <- all

# I created a new column in metadata (Group_gt_ind) that concatenates the info from Genotype and Inducer columns so coefficients include genotype info.
# I also added a new column (Read_depth) to tag samples with High or Low sequencing depth so this factor can be controlled for in the design formula.
# The design formula use in DESeq2 is the following:
  
design = ~ Read_depth + Group_gt_ind

# add new factors (Group_gt_ind and Read_depth (high > 10M reads / Low < 10M reads))
meta_one_RS$Group_gt_ind <- factor(paste0(meta_one_RS$Genotype, meta_one_RS$Inducer))
meta_one_RS$Read_depth <- 'High'
meta_one_RS[meta_one_RS$Read_counts < 10e6,]$Read_depth <- 'Low'
meta_one_RS$Read_depth <- as.factor(meta_one_RS$Read_depth)

# Adding read_depth in design to control for read_depth
dds.one.RS <- DESeqDataSetFromMatrix(countData = all_one_RS, 
                                     colData = meta_one_RS,  
                                     design = ~ Read_depth + Group_gt_ind)

# ------------- Exploratory analysis with DESeq object- RIBOseq -------------

# Plot total reads per sample using barchar
p <- ggbarplot(data = meta_one_RS, 
               x = "Sample_id", 
               y = "Read_counts",
               x.text.angle = 90,
               fill = "Inducer", 
               title = "Total read counts", 
               ylab = "Read counts",
               sort.by.groups = TRUE,
               palette = "jco",
               sort.val = "asc", 
               facet.by = "Genotype")
ggsave("Plots/barplot_read_counts_RS.pdf", plot = p)
p

# Normalize counts
rlog.one <- rlog(dds.one.RS, blind=FALSE)

# Keep genes with at least 20 reads total across samples
keep <- rowSums(counts(dds.one.RS)) >= 0 # Agnes wanted to keep all genes for the analysis
dds.one.RS <- dds.one.RS[keep,]

# Calculate distances between samples
sampleDists <- dist(t(assay(rlog.one)))

# Plot inter-sample distances
old.par <- par(no.readonly=T)

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rlog.one$Sequencing_pool, rlog.one$Genotype, rlog.one$Inducer, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# PCA
my_top_genes = 500

pcaData <- plotPCA(rlog.one, intgroup=c("Genotype", "Inducer"), returnData=TRUE, ntop = my_top_genes)
percentVar <- round(100 * attr(pcaData, "percentVar"))
y.coords = c(min(pcaData$PC1, pcaData$PC2), max(pcaData$PC1, pcaData$PC2))
x.coords = y.coords
p1 <- ggplot(pcaData, aes(PC1, PC2, color=Genotype, shape=Inducer)) +
  geom_point(size=3) + scale_color_lancet() + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed(ratio = (max(pcaData$PC1)-min(pcaData$PC1))/(max(pcaData$PC2)-min(pcaData$PC2))) 

ggsave("Plots/pca_dataset_1_Induc_gt_RS.pdf", plot = p1)
p1

pcaData <- plotPCA(rlog.one, intgroup=c("Read_counts", "Inducer"), returnData=TRUE, ntop = my_top_genes)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p2 <- ggplot(pcaData, aes(PC1, PC2, color=Read_counts, shape=Inducer)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed(ratio = (max(pcaData$PC1)-min(pcaData$PC1))/(max(pcaData$PC2)-min(pcaData$PC2))) + scale_color_gradient2(high = "blue", mid = "yellow")

ggsave("Plots/pca_dataset_1_Induc_read_counts_RS.pdf", plot = p2)
p2

# PCA plots indicate that samples separated on PC1 mostly by sequencing depth, within treatments and genotypes. This implies that sequencing depth has to be controlled for by including this factor in the design formula.

# ------------- Filtering out poorly-expressed genes (less than 20 reads across all samples) - RIBOseq -------------

# Keep genes with at least 10 reads total across samples
keep <- rowSums(counts(dds.one.RS)) >= N_reads
dds.one.RS <- dds.one.RS[keep,]

# ------------- Spliting DESeq object based on genotype - RIBOseq -------------

dds.one.RS.wt <- dds.one.RS[ , dds.one.RS$Genotype == "WT"]
dds.one.RS.wt$Genotype <- droplevels( dds.one.RS.wt$Genotype)
dds.one.RS.wt$Group_gt_ind <- droplevels( dds.one.RS.wt$Group_gt_ind)
dds.one.RS.wt$Group <- droplevels( dds.one.RS.wt$Group)

dds.one.RS.ko <- dds.one.RS[ , dds.one.RS$Genotype == "RNaseL_KO"]
dds.one.RS.ko$Genotype <- droplevels( dds.one.RS.ko$Genotype)
dds.one.RS.ko$Group_gt_ind <- droplevels( dds.one.RS.ko$Group_gt_ind)
dds.one.RS.ko$Group <- droplevels( dds.one.RS.ko$Group)

# ------------- Calculate differential expression for WT - RIBOseq -------------

# Calculate DE for WT samples
dds.one.RS.wt$Group_gt_ind <- relevel(dds.one.RS.wt$Group_gt_ind, "WTNone")
dds.one.RS.wt <- DESeq(dds.one.RS.wt)
resultsNames(dds.one.RS.wt)

res_wtIC_vs_wtNone <- results(dds.one.RS.wt, list(c( "Group_gt_ind_WTpolyIC_vs_WTNone")))
res_wt25A_vs_wtNone <- results(dds.one.RS.wt, list(c( "Group_gt_ind_WT25A_vs_WTNone")))


# Generate MA plots
make_MA_plots(res_wtIC_vs_wtNone, file_prefix = "res_wtIC_vs_wtNone_RIBOseq_cds")
make_MA_plots(res_wt25A_vs_wtNone, file_prefix = "res_wt25A_vs_wtNone_RIBOseq_cds")

summary(res_wtIC_vs_wtNone, alpha = 0.05)
summary(res_wt25A_vs_wtNone, alpha = 0.05)

# ------------- Calculate differential expression for RNaseL_KO - RIBOseq -------------

dds.one.RS.ko$Group_gt_ind <- relevel(dds.one.RS.ko$Group_gt_ind, "RNaseL_KONone")
design(dds.one.RS.ko) <- ~Group_gt_ind # Changing design given that there is no KO sample with Low read counts, otherwise you get an error message 
# Error: full model matrix is less than full rank
dds.one.RS.ko <- DESeq(dds.one.RS.ko)
resultsNames(dds.one.RS.ko)

res_koIC_vs_koNone <- results(dds.one.RS.ko, list(c( "Group_gt_ind_RNaseL_KOpolyIC_vs_RNaseL_KONone")))
res_ko25A_vs_koNone <- results(dds.one.RS.ko,list(c( "Group_gt_ind_RNaseL_KO25A_vs_RNaseL_KONone")))

# Generate MA plots
make_MA_plots(res_koIC_vs_koNone, file_prefix = "res_koIC_vs_koNone_RIBOseq_cds")
make_MA_plots(res_ko25A_vs_koNone, file_prefix = "res_ko25A_vs_koNone_RIBOseq_cds")

summary(res_koIC_vs_koNone, alpha = 0.05)
summary(res_ko25A_vs_koNone, alpha = 0.05)

# ------------- Calculate mean normalized read counts per gene -------------

dir.create('./Mean_norm_counts', showWarnings = FALSE)

get_mean_norm_counts(dds.one.RS.wt, genotype = 'wt', feature_type = 'RS')
get_mean_norm_counts(dds.one.RS.ko, genotype = 'ko', feature_type = 'RS')

# ------------- Write DE tables to file - RIBOseq -------------

# Sort results by Log2FC
res_wtIC_vs_wtNone.logfc_sorted <- sort_and_write_res_table(res_wtIC_vs_wtNone, "DE_wtIC_vs_wtNone_RS")
res_koIC_vs_koNone.logfc_sorted <- sort_and_write_res_table(res_koIC_vs_koNone, "DE_koIC_vs_koNone_RS")
res_wt25A_vs_wtNone.logfc_sorted <- sort_and_write_res_table(res_wt25A_vs_wtNone, "DE_wt25A_vs_wtNone_RS")
res_ko25A_vs_koNone.logfc_sorted <- sort_and_write_res_table(res_ko25A_vs_koNone, "DE_ko25A_vs_koNone_RS")

# Save sorted files as a list
DE_results <- list()
DE_results[["wtIC_vs_wtNone_RS"]]  <- res_wtIC_vs_wtNone.logfc_sorted
DE_results[["koIC_vs_koNone_RS"]] <- res_koIC_vs_koNone.logfc_sorted
DE_results[["wt25A_vs_wtNone_RS"]] <- res_wt25A_vs_wtNone.logfc_sorted
DE_results[["ko25A_vs_koNone_RS"]] <- res_ko25A_vs_koNone.logfc_sorted


sessionInfo()
