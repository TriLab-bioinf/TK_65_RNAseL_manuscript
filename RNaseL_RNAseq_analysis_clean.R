# RNAseq analysis of WT, RNAseL-KO and RNAseA samples induced with 2-5A/polyI:C
# Agnes Karasik, Hernan Lorenzi, Nicholas Guydosh
# 4/1/2024

# Contact: Nicholas Guydosh
# Contact email: nicholas.guydosh@nih.gov


suppressMessages(library("org.Hs.eg.db"))
suppressMessages(library("pheatmap"))
suppressMessages(library("EnhancedVolcano"))
suppressMessages(library("ggpubr"))
suppressMessages(library("DESeq2"))
suppressMessages(library("stringr"))
suppressMessages(library("biomaRt"))
suppressMessages(library("tidyverse"))
suppressMessages(library("clusterProfiler"))
suppressMessages(library("ggsci"))
suppressMessages(library("viridis"))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("msigdbr"))
suppressMessages(library("vioplot"))
suppressMessages(library("cowplot"))
suppressMessages(library("broom"))

#  ------------- Control filtering steps and shrinkage -------------

# Remove genes with 0 counts across all samples
remove_genes_with_zero_counts <- TRUE

# Filter reads with less than N reads
N_reads <- 20


# ------------- Required functions -------------

# Function to remove all-zero rows (by adjusting min_total_count the function can filter out rows based on total counts other than 0)
remove_all_zero_rows <- function(df, min_total_count = 0){
  df <- df[rowSums(df) > min_total_count,]
  return(df)
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

# Function for retrieving list of significantly DE genes from a DESeq object
fetch_significant_genes <- function(deseq_res, adj_p = 0.05, log_fc, change_direction = 'both'){
  if(change_direction == 'up'){
    res <- rownames(deseq_res[!is.na(deseq_res$padj) & 
                                deseq_res$padj <= 0.05 & 
                                deseq_res$log2FoldChange > 1, ])
  } 
  else if(change_direction == 'down'){
    res <- rownames(deseq_res[!is.na(deseq_res$padj) & 
                                deseq_res$padj <= 0.05 & 
                                deseq_res$log2FoldChange < 1, ])
  }
  else {
    # both
    res <- rownames(deseq_res[!is.na(deseq_res$padj) & 
                                deseq_res$padj <= 0.05 & 
                                abs(deseq_res$log2FoldChange) > 1, ])
  }
  return(res)
}

# Function to plot correlation between 2 log2FC vectors
plot_correlation <- function(df, my_group){
  df$min=min(df$adjp_pIC,df$adjp_25A)
  df <- df[order(df$min, decreasing = T),]
  df$Significance = ifelse(df$adjp_pIC <= 0.05 | df$adjp_25A <= 0.05, "adj.p < 0.05","N.S." )
  my.p <- ggscatter(df, x = "Log2_pIC", y = "Log2_25A", 
                    color = "Significance", alpha=0.5, palette = c("darkred","lightgray"),
                    title = paste('Group =',my_group),
                    add = "reg.line",                                 
                    conf.int = TRUE,                                  
                    add.params = list(color = "blue",
                                     fill = "gray", legend.title="")
  ) + 
    stat_cor(method = "pearson", label.x = mean(df$Log2_pIC), label.y = max(df$Log2_25A)+1) +
    geom_abline(slope = 1, intercept = 0, col = "black", size=0.5, linetype="dashed") +
    geom_hline(yintercept = 0, size=0.5, linetype="dashed", color = 'gray') +
    geom_vline(xintercept = 0, size=0.5, linetype="dashed", color = 'gray') 
  
  ggsave2(filename = paste0("./Plots/corr_plot_",str_replace_all(my_group,' ','_'),".pdf"))
  print(my.p)
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

# Function to plot heat maps
plot_heat_map <- function( my_vstd, gene_list, file_name, variables){
  
  # Replace gene IDs by gene symbols
  gene_list <- replace_annotation_keys(gene_list, from = 'ENSEMBL', to = 'SYMBOL')
  rownames(my_vstd) <- replace_annotation_keys(rownames(my_vstd), from = 'ENSEMBL', to = 'SYMBOL')
  
  # Plot the heat map
  hmp <- pheatmap(assay(my_vstd)[gene_list,], cluster_rows=T, show_rownames=TRUE,
                  cluster_cols=T, annotation_col = variables, fontsize_col = 5, tile = file_name)
  
  ggsave(filename = paste0(file_name,"_heatmap.eps"), plot = hmp, width = 8.5, height = 11, units = "in")
  
  print(hmp)
}

# Function to make volcane plots
generate_volcano_plot_simple <- function(res.tmp, my_file_name, log_scale = FALSE){
  res.tmp <- replace_annotation_keys(res.tmp, from = 'ENSEMBL', to = 'SYMBOL')
  
  vp <- EnhancedVolcano(res.tmp,
                        lab = rownames(res.tmp),
                        x = 'log2FoldChange',
                        y = 'padj',
                        #ylim = c(-0.5, 30),
                        pCutoff = 0.05,
                        FCcutoff = 1,
                        pointSize = 2,
                        #shape = c(1, 1, 1, 1),
                        col = c("grey", "grey","grey","light blue"),
                        colAlpha = 4/5,
                        labSize = 6,  # Controls labels size
                        labCol = "black",
                        title = res.tmp@elementMetadata$description[2],
                        titleLabSize = 10,
                        subtitle = '', # add subtitle here
                        subtitleLabSize = 10,
                        legendPosition = 'right',
                        legendLabSize = 15,
                        legendIconSize = 6.0,
                        axisLabSize = 20,
                        drawConnectors = TRUE,
                        selectLab = c("JUN", "FOSB", "FOSL", "IL6", "CXCL8", "ATF3", "PPPR15A", "GDF15", "CXCL2", "CCL20", "DUSP1", "DUSP8", "FOXA2", "ATF5"), # vector of gene symbols to label on volcanoplot
                        boxedLabels = FALSE,
                        gridlines.major = FALSE,
                        gridlines.minor = FALSE
  )
  
  if (log_scale){
    vp <- vp + scale_x_log10()
  }
  ggsave(filename = paste0(my_file_name,"_simple.pdf"), height=10, width=10,  plot = vp )
  
}

# ------------- Create output dirs -------------

dir.create(path = "./Plots", showWarnings = F)
dir.create(path = "./DE", showWarnings = F)
dir.create('./Mean_norm_counts', showWarnings = FALSE)
dir.create('./Tables', showWarnings = FALSE)

# ------------- Load RNAseq CDS data -------------

all <- read.delim2("./data/cdsrna_round.csv", sep = ",", header = TRUE, row.names = 1, comment.char = c("#"))
all.tmp <- as.data.frame(lapply(all, as.numeric)) %>% round()
all[, 1:ncol(all)] <- all.tmp

# Keep table with Ensemble IDs and gene Symbols
gene_symbols <- replace_annotation_keys(rownames(all),from = 'ENSEMBL', to = 'SYMBOL')
ensembl_to_symbol <- as.data.frame(cbind("Ensembl_ID" = rownames(all), "gene_name" = gene_symbols), row.names = 1)

# Load metadata
metadata <- read.delim2("./data/Metadata_RNAseq.txt", sep = "\t", row.names = 1, header = T)

# keep only samples that are present in all object
metadata <- metadata[colnames(all),]

# Add total read counts and sample id columns to metadata
metadata <- cbind(metadata, Read_counts =colSums(all), Sample_id = rownames(metadata))

# add new factors (Group_gt_ind for CDS based counts)
metadata <- metadata %>% mutate(Group_gt_ind = paste0(Genotype, Inducer))

# Add Sequencing pool groups (Batch) to metadata: G69 = GL069, G72_83 =  GL072 + GL083, G136 = GL136
metadata <- metadata %>% mutate(Batch = case_when(
                                Sequencing_pool == "GL069" ~ "G69",
                                Sequencing_pool %in% c("GL072","GL083") ~ "G72_83",
                                Sequencing_pool == "GL136" ~ "G136"
                                )
                              )

#Remove all zero rows from read counts table (all)
if (remove_genes_with_zero_counts){
  all <- remove_all_zero_rows(all, min_total_count = 0)
}

# ------------- Analysis of expression data using DESeq2 - CDS -------------

# Convert variables into factors
to_factor <- c("Genotype","Inducer","Batch","Group_gt_ind","Group")
metadata[to_factor] <- lapply(metadata[to_factor], factor)

# Rename metadata and all dataframes
meta_one_cds <- metadata
all_one_cds <- all

# Adding read_depth in design to control for read_depth
dds.one.cds <- DESeqDataSetFromMatrix(countData = all_one_cds,
                                      colData = meta_one_cds,
                                      design = ~ Batch + Group_gt_ind)

# Normalize counts
vsd.one <- vst(dds.one.cds, blind=FALSE)
rlog.one <- rlog(dds.one.cds, blind=FALSE)

# Keep genes with at least N reads total in at least 3 samples
keep <- rowSums(counts(dds.one.cds) >= N_reads) >= 3
dds.one.cds <- dds.one.cds[keep,]

# Calculate distances between samples
sampleDists <- dist(t(assay(vsd.one)))

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
pcaData <- plotPCA(rlog.one, intgroup=c("Genotype", "Inducer"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
y.coords = c(min(pcaData$PC1, pcaData$PC2), max(pcaData$PC1, pcaData$PC2))
x.coords = y.coords
p1 <- ggplot(pcaData, aes(PC1, PC2, color=Genotype, shape=Inducer)) +
  geom_point(size=3) + scale_color_lancet() + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed(ratio = (max(pcaData$PC1)-min(pcaData$PC1))/(max(pcaData$PC2)-min(pcaData$PC2))) 

ggsave("Plots/pca_dataset_1_Induc_gt.pdf", plot = p1)
p1

pcaData <- plotPCA(rlog.one, intgroup=c("Batch", "Inducer"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p2 <- ggplot(pcaData, aes(PC1, PC2, color=Batch, shape=Inducer)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed(ratio = (max(pcaData$PC1)-min(pcaData$PC1))/(max(pcaData$PC2)-min(pcaData$PC2))) 
p2
ggsave("Plots/pca_dataset_1_Induc_batch.pdf", plot = p2)


# Spliting DESeq object based on genotype - CDS
dds.one.cds.wt <- dds.one.cds[ , dds.one.cds$Genotype == "WT"]
dds.one.cds.wt$Genotype <- droplevels( dds.one.cds.wt$Genotype)
dds.one.cds.wt$Group_gt_ind <- droplevels( dds.one.cds.wt$Group_gt_ind)
dds.one.cds.wt$Group <- droplevels( dds.one.cds.wt$Group)

dds.one.cds.ko <- dds.one.cds[ , dds.one.cds$Genotype == "RNaseL_KO"]
dds.one.cds.ko$Genotype <- droplevels( dds.one.cds.ko$Genotype)
dds.one.cds.ko$Group_gt_ind <- droplevels( dds.one.cds.ko$Group_gt_ind)
dds.one.cds.ko$Group <- droplevels( dds.one.cds.ko$Group)

# ------------- Calculate DE for WT samples -------------

dds.one.cds.wt$Group_gt_ind <- relevel(dds.one.cds.wt$Group_gt_ind, "WTNone")
dds.one.cds.wt <- DESeq(dds.one.cds.wt)
resultsNames(dds.one.cds.wt)

# Using results function 
res_wtIC_vs_wtNone <- results(dds.one.cds.wt, list(c( "Group_gt_ind_WTpolyIC_vs_WTNone" )))
res_wt25A_vs_wtNone <- results(dds.one.cds.wt, list(c( "Group_gt_ind_WT25A_vs_WTNone")))

# ------------- Calculate DE for KO samples -------------

# Changing design formula given that there is no KO sample with Low read counts, otherwise you get an error message 

dds.one.cds.ko$Group_gt_ind <- relevel(dds.one.cds.ko$Group_gt_ind, "RNaseL_KONone")
dds.one.cds.ko$Batch <- droplevels(dds.one.cds.ko$Batch)
dds.one.cds.ko <- DESeq(dds.one.cds.ko)
resultsNames(dds.one.cds.ko)

res_koIC_vs_koNone <- results(dds.one.cds.ko, list(c( "Group_gt_ind_RNaseL_KOpolyIC_vs_RNaseL_KONone")))
res_ko25A_vs_koNone <- results(dds.one.cds.ko, list(c( "Group_gt_ind_RNaseL_KO25A_vs_RNaseL_KONone")))

# ------------- Calculate mean normalized read counts per gene -------------

# This is for TE plot for filtering, Figure 5A (we only plot genes with mean normalized counts >50)

# Subset DESeq object for 25A treatment and control for WT
dds.one.cds.wt.25A <- dds.one.cds.wt[,dds.one.cds.wt$Inducer %in% c("None", "25A")]
mean_norm_counts.wt.25A <- apply(counts(dds.one.cds.wt.25A, normalized=TRUE), 1, function(x) mean(x))

# Subset DESeq object for polyIC treatment and control for WT
dds.one.cds.wt.pIC <- dds.one.cds.wt[,dds.one.cds.wt$Inducer %in% c("None", "polyIC")]
mean_norm_counts.wt.pIC <- apply(counts(dds.one.cds.wt.pIC, normalized=TRUE), 1, function(x) mean(x))

# Subset DESeq object for polyIC treatment and control for RNAseL KO samples
dds.one.cds.ko.pIC <- dds.one.cds.ko[,dds.one.cds.ko$Inducer %in% c("None", "polyIC")]
mean_norm_counts.ko.pIC <- apply(counts(dds.one.cds.ko.pIC, normalized=TRUE), 1, function(x) mean(x))

write.table(x = as.data.frame(mean_norm_counts.wt.25A), file = './Mean_norm_counts/mean_norm_counts.wt.cds.25A.txt', sep = "\t",  )
write.table(x = as.data.frame(mean_norm_counts.wt.pIC), file = './Mean_norm_counts/mean_norm_counts.wt.cds.IC', sep = "\t" )
write.table(x = as.data.frame(mean_norm_counts.ko.pIC), file = './Mean_norm_counts/mean_norm_counts.ko.cds.IC', sep = "\t" )

# ------------- Write DE tables to file -------------

# Sort results by Log2FC
res_wtIC_vs_wtNone.logfc_sorted <- sort_and_write_res_table(res_wtIC_vs_wtNone, "DE_wtIC_vs_wtNone_cds")
res_koIC_vs_koNone.logfc_sorted <- sort_and_write_res_table(res_koIC_vs_koNone, "DE_koIC_vs_koNone_cds")
res_wt25A_vs_wtNone.logfc_sorted <- sort_and_write_res_table(res_wt25A_vs_wtNone, "DE_wt25A_vs_wtNone_cds")
res_ko25A_vs_koNone.logfc_sorted <- sort_and_write_res_table(res_ko25A_vs_koNone, "DE_ko25A_vs_koNone_cds")

# ------------- Store sorted DE results into DE_results list -------------

DE_results = list()
DE_results[["wtIC_vs_wtNone_cds"]]  <- res_wtIC_vs_wtNone.logfc_sorted
DE_results[["koIC_vs_koNone_cds"]] <- res_koIC_vs_koNone.logfc_sorted
DE_results[["wt25A_vs_wtNone_cds"]] <- res_wt25A_vs_wtNone.logfc_sorted
DE_results[["ko25A_vs_koNone_cds"]] <- res_ko25A_vs_koNone.logfc_sorted


# ------------- Fetch genes significantly DE (upregulated, padj <= 0.05 and log2FC >1) -------------

# Genes significantly DE in WT and induced by polyIC
wt.g <- fetch_significant_genes(deseq_res = res_wtIC_vs_wtNone.logfc_sorted, adj_p = 0.05, change_direction = 'up')

# Genes significantly DE in KO and induced by polyIC
ko.g <- fetch_significant_genes(deseq_res = res_koIC_vs_koNone.logfc_sorted, adj_p = 0.05, change_direction = 'up')

# Genes significantly DE in WT and induced by 25A
wt.25A.g <- fetch_significant_genes(deseq_res = res_wt25A_vs_wtNone.logfc_sorted, adj_p = 0.05, change_direction = 'up')

# Genes significantly DE in KO and induced by 25A
ko.25A.g <- fetch_significant_genes(deseq_res = res_ko25A_vs_koNone.logfc_sorted, adj_p = 0.05, change_direction = 'up')

# save upregulated genes into a csv file
write_csv(as.data.frame(wt.g), "./DE/WT_IC_DEgenes.csv")
write_csv(as.data.frame(ko.g), "./DE/KO_IC_DEgenes.csv")
write_csv(as.data.frame(wt.25A.g), "./DE/WT_25A_DEgenes.csv")

# ------------- Plot Venn diagram to depiction overlap of significant DE genes -------------

VennDiagram::venn.diagram(x = list(wt.g, ko.g, wt.25A.g ),
                          category.names = c("WT poly I:C", "KO polyI:C", "WT 25A"),
                          filename = "./Plots/Inducer_DE_venn_0.05.png",
                          imagetype = "png",
                          output=TRUE,
                          height = 480 , 
                          width = 480 , 
                          resolution = 300,
                          compression = "lzw",
                          lwd = 1,
                          col=c("#440154ff", '#21908dff', '#fde725ff'),
                          fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3)),
                          cex = 0.5,
                          fontfamily = "sans",
                          cat.cex = 0.3,
                          cat.default.pos = "outer",
                          cat.pos = c(-27, 27, 135),
                          cat.dist = c(0.055, 0.055, 0.085),
                          cat.fontfamily = "sans",
                          cat.col = c("#440154ff", '#21908dff', '#fde725ff'),
                          rotation = 1
)

# ------------- Plot Heatmap -------------

# Fetch metadata for heatmap
my_metadata <- select(meta_one_cds, Inducer)

# Plot heatmaps
top_60 <- list('25A_WT'=wt.25A.g, 'polyIC_WT'=wt.g, '25A_RNaseL_KO'=wt.25A.g, 'polyIC_RNaseL_KO'=ko.g)
for (inducer in c('25A','polyIC')){
  for (genotype in c('WT','RNaseL_KO')){
    id <- paste0(inducer,'_',genotype)
    vsd.one.subset <- vsd.one[,vsd.one$Inducer  %in% c(inducer,'None') & 
                                vsd.one$Genotype  %in% genotype]
    plot_heat_map(my_vstd = vsd.one.subset, 
                  gene_list = top_60[[id]][1:60], 
                  file_name = paste0("Plots/top60_",inducer,"_treatment_",genotype,"_clustered"), 
                  variables = my_metadata)
  }
}


# ------------- Plot volcano plots -------------

generate_volcano_plot_simple(res.tmp = res_wtIC_vs_wtNone, my_file_name = "./Plots/res_wtIC_vs_wtNone.logfc_sorted")
generate_volcano_plot_simple(res_koIC_vs_koNone, my_file_name = "./Plots/res_koIC_vs_koNone.logfc_sorted")

# Figure 1B
generate_volcano_plot_simple(res_wt25A_vs_wtNone, my_file_name = "./Plots/res_wt25A_vs_wtNone.logfc_sorted")

# ------------- Correlation -------------

# Generate result table of DE genes in  res_wtIC_vs_wtNone but not in res_koIC_vs_koNone
# or in res_wt25A_vs_wtNone. 
# This are genes that might require both an induction via dsRNA/DNA and a functional RNAseL

# Filter 1: adj.pvalues > 0.25 in other samples
res_wt25A_vs_wtNone.sig.id <- rownames(res_wt25A_vs_wtNone)[which(res_wt25A_vs_wtNone$padj <= 0.25)]
res_wtIC_vs_wtNone.sig.id <- rownames(res_wtIC_vs_wtNone)[which(res_wtIC_vs_wtNone$padj <= 0.05)]
res_koIC_vs_koNone.sig.id <- rownames(res_koIC_vs_koNone)[which(res_koIC_vs_koNone$padj <= 0.25)]

# Genes significantly DE in polyIC-WT (adj.p <= 0.05) but not in polyIC-KO or 25A-WT (adj.p > 0.25)
filtered_ids <- res_wtIC_vs_wtNone.sig.id[!(res_wtIC_vs_wtNone.sig.id %in% c(res_wt25A_vs_wtNone.sig.id, res_koIC_vs_koNone.sig.id))]

res_wtIC_vs_wtNone.filtered <- subset(res_wtIC_vs_wtNone, rownames(res_wtIC_vs_wtNone) %in% filtered_ids)
res_koIC_vs_koNone.filtered <- subset(res_koIC_vs_koNone, rownames(res_koIC_vs_koNone) %in% filtered_ids)
res_wt25A_vs_wtNone.filtered <- subset(res_wt25A_vs_wtNone, rownames(res_wt25A_vs_wtNone) %in% filtered_ids)

# Filter 2: abs(Log2FC) res_wt25A_vs_wtNone.filtered < 0.1 
filter2 <- abs(res_wt25A_vs_wtNone.filtered$log2FoldChange) < 0.1

# Filter 3: opposite direction of Log2FC between polyIC-WT and 25A-WT
filter3 <- res_wtIC_vs_wtNone.filtered$log2FoldChange * res_wt25A_vs_wtNone.filtered$log2FoldChange < 0

# Apply filter2 and filter3 to res_wtIC_vs_wtNone.filtered
res_wtIC_vs_wtNone.filtered23 <-  res_wtIC_vs_wtNone.filtered[filter2 | filter3, ]

# Save res_wtIC_vs_wtNone.filtered object in DE_results list
DE_results[["wtIC_vs_wtNone_KOns_25Ans"]] <- res_wtIC_vs_wtNone.filtered
DE_results[["wtIC_vs_wtNone_KOns_25Ans_f23"]] <- res_wtIC_vs_wtNone.filtered23
res_wtIC_vs_wtNone.filtered.logfc_sorted <- sort_and_write_res_table(res_wtIC_vs_wtNone.filtered, "wtIC_vs_wtNone_KOns_25Ans")
res_wtIC_vs_wtNone.filtered23.logfc_sorted <- sort_and_write_res_table(res_wtIC_vs_wtNone.filtered23, "wtIC_vs_wtNone_KOns_25Ans_f23")

# -------------  Correlation analysis -------------

# Calculate DE between polyIC WT and 25A WT samples
dds.one.cds.wt$Group_gt_ind <- relevel(dds.one.cds.wt$Group_gt_ind, "WT25A")
dds.one.cds.wt <- DESeq(dds.one.cds.wt)
resultsNames(dds.one.cds.wt)

res_wtIC_vs_wt25A <- results(dds.one.cds.wt, list(c( "Group_gt_ind_WTpolyIC_vs_WT25A" )))

# GeneIDs from 25A-WT vs None-WT with |Log2FC| <= 0.1
gene_ids_wt25A_vs_wtNone_lowFC <- rownames(subset(res_wt25A_vs_wtNone,
                                                  abs(res_wt25A_vs_wtNone$log2FoldChange) <= 0.1)
                                           )

# GeneIDs from pIC-WT vs None-WT with |Log2FC| <= 0.1
gene_ids_wtIC_vs_wtNone_lowFC <- rownames(subset(res_wtIC_vs_wtNone,
                                                  abs(res_wtIC_vs_wtNone$log2FoldChange) <= 0.1)
)

# Filter wtIC_vs_wt25A results keeping genes with padj <= 0.05 
# and that do not change in  wt25A_vs_wtNone (Log2FC ~ 0)
# or are up in wtIC_vs_None and down in wt25A_vs_None (res_wtIC_vs_wt25A_sig_fc1) or that are 
# down in wtIC_vs_None and up in wt25A_vs_None (res_wtIC_vs_wt25A_sig_fc2)
# These are genes highly likely to respond differently between polyIC and 25A
res_wtIC_vs_wt25A_sig_fc1 <- subset(res_wtIC_vs_wt25A, res_wtIC_vs_wt25A$padj <= 0.05 &
                                      (rownames(res_wtIC_vs_wt25A) %in%
                                      gene_ids_wt25A_vs_wtNone_lowFC | (res_wt25A_vs_wtNone$log2FoldChange < 0 & res_wtIC_vs_wtNone$log2FoldChange > 0)))

res_wtIC_vs_wt25A_sig_fc2 <- subset(res_wtIC_vs_wt25A, res_wtIC_vs_wt25A$padj <= 0.05 &
                                      (rownames(res_wtIC_vs_wt25A) %in%
                                         gene_ids_wtIC_vs_wtNone_lowFC | (res_wt25A_vs_wtNone$log2FoldChange > 0 & res_wtIC_vs_wtNone$log2FoldChange < 0)))

# NOTE: res_wtIC_vs_wt25A_sig_fc2 == 0 genes and therefore is not analyzed further

# Prepare dataframe for correlation analysis
corr_input.df <- as.data.frame(cbind(
                        Log2_pIC = res_wtIC_vs_wtNone$log2FoldChange,
                        Log2_25A = res_wt25A_vs_wtNone$log2FoldChange,
                        adjp_pIC = res_wtIC_vs_wtNone$padj,
                        adjp_25A = res_wt25A_vs_wtNone$padj
                      )
                    )
rownames(corr_input.df) <- rownames(res_wtIC_vs_wtNone)
corr_input.df$adjp_pIC[is.na(corr_input.df$adjp_pIC)] <- 1
corr_input.df$adjp_25A[is.na(corr_input.df$adjp_25A)] <- 1

# Add gene groups to df where:
# "pIC vs 25A WT" = DE genes between pIC and 25A with Log2FC in 25A_vs_None between (-0.1,0.1) or Log2FC is up in pIC_vs_None and negative in 25A_vs_None
# "pIC:25A vs None WT" = DE genes between pIC vs None AND + between 25A vs None in WT
# "pIC vs None WT" = DE genes ONLY in pIC vs None in WT
# "25A vs None WT" = DE genes ONLY in 25A vs None in WT
# "n.s." = genes that are not DE in either condition vs None in WT
my_groups = ifelse(rownames(corr_input.df) %in% rownames(res_wtIC_vs_wt25A_sig_fc1),
                   "pIC vs 25A WT",
                   ifelse(corr_input.df$adjp_pIC <= 0.05 & corr_input.df$adjp_25A <= 0.05,
                          "pIC:25A vs None WT",
                          ifelse(corr_input.df$adjp_pIC <= 0.05,
                                 "pIC vs None WT",
                                 ifelse(corr_input.df$adjp_25A <= 0.05,
                                        "25A vs None WT",
                                        "n.s" 
                                 )
                          )
                   )
            )
corr_input.df$Groups = my_groups


## ADD CORRELATION for the whole set
# Run correlation analysis between pIC and 25A withing groups

corr.df <- data.frame()
idx = 1
for (i in levels(factor(corr_input.df$Groups))){
  df.sub <- subset(corr_input.df, corr_input.df$Groups == i)
  plot_correlation(df.sub, my_group = df.sub$Groups[1] )
  x <- cor.test(df.sub$Log2_pIC, df.sub$Log2_25A, method = "pearson")
  corr.df[idx, c("Group","Num_genes","R", "pval")] <- c(i,dim(df.sub)[1] ,round(x$estimate, digits = 2), round(x$p.value,digits = 3))
  idx <- idx + 1
}
# Add correlation across all genes
x <- cor.test(corr_input.df$Log2_pIC, corr_input.df$Log2_25A, method = "pearson")
corr.df[idx, c("Group","Num_genes","R", "pval")] <- c('All', dim(corr_input.df)[1], round(x$estimate, digits = 2), round(x$p.value,digits = 3))
corr.df$Description <- c("Genes DE only in 2:5A vs control in WT",
                  "genes that were not DE in either treatment vs control in WT",
                  "Genes DE between polyIC and 2:5A treatments having a Log2FC in 2:5A vs control between (-0.1,0.1) or the Log2FC is up in polyIC vs control and negative in 2:5A vs control",
                  "Genes DE only in polyIC vs control in WT",
                  "Genes DE between polyIC vs control and also between 2:5A vs control in WT",
                  "All genes")

tibble(corr.df)
write.table(x = corr.df, file = "./Tables/corr_analysis.csv", sep = ",", row.names = FALSE)

# Correlation across all genes ignoring Groups
plot_correlation(corr_input.df[order(corr_input.df$adjp_pIC, decreasing = TRUE),], my_group = 'pIC_vs_None vs 25A_vs_None')


# -------------  Functional analysis -------------

# import GeneRIF and hallmark datasets that are used below

# hallmark
hallmark_IFN_gamma <- read.csv("./data/hallmark_IFNgamma.csv", header=FALSE, colClasses = c('character','character','character','character'))

IFN_gamma.tbl <- tibble(ensembl_id=replace_annotation_keys(key_list = hallmark_IFN_gamma$V2,
                                                           from = "ENTREZID", to = "ENSEMBL"),
                        gene_name = hallmark_IFN_gamma$V1, entrez_id=hallmark_IFN_gamma$V2)
IFN_gamma_ens <- IFN_gamma.tbl$ensembl_id

my_new_ids <- c("ENSG00000173020", "ENSG00000118523", "ENSG00000142871", "ENSG00000107201", "ENSG00000188486", "ENSG00000136999", "ENSG00000143207", "ENSG00000162430", "ENSG00000104415", "ENSG00000091436", "ENSG00000184293", "ENSG00000204628", "ENSG00000170854", "ENSG00000184640")
names(my_new_ids) <- c("ADRBK1", "CTGF", "CYR61", "DDX58", "H2AFX", "NOV", "RFWD2", "SEPN1", "WISP1", "ZAK", "CLECL1", "GNB2L1", "MINA", "`SEPT9")

p38_GR <- read.csv("./data/p38_GR.csv", header=FALSE)
p38_GR.tbl <- tibble(ensembl_id=replace_annotation_keys(key_list = p38_GR$V1,
                                                        from = "SYMBOL", 
                                                        to = "ENSEMBL"),
                     gene_name = p38_GR$V1)

p38_GR.tbl$ensembl_id[p38_GR.tbl$ensembl_id %in% names(my_new_ids)] <- my_new_ids[p38_GR.tbl$ensembl_id[p38_GR.tbl$ensembl_id %in% names(my_new_ids)]]

JNK <- read.csv("./data/JNK_GeneRIF.csv", header=FALSE)
JNK.tbl <- tibble(ensembl_id=replace_annotation_keys(key_list = JNK$V1,
                                                        from = "SYMBOL", 
                                                        to = "ENSEMBL"),
                     gene_name = JNK$V1)

JNK.tbl$ensembl_id[JNK.tbl$ensembl_id %in% names(my_new_ids)] <- my_new_ids[JNK.tbl$ensembl_id[JNK.tbl$ensembl_id %in% names(my_new_ids)]]


JNK_p38 <- unique(p38_GR.tbl$ensembl_id, JNK.tbl$ensembl_id)
JNK_p38 <- JNK_p38[str_starts(JNK_p38, pattern = "ENSG")]

# Fig 2C: correlation plot between poly I:C and 2-5A showing JNK and IFN genes

corr_input.filtered.df <- filter(corr_input.df, adjp_pIC<=1|adjp_25A<=1) # Plotting all dots (no filtering)
my_groups2 = ifelse(rownames(corr_input.filtered.df) %in% JNK_p38,
                    "JNK/p38", 
                    ifelse(rownames(corr_input.filtered.df) %in% IFN_gamma_ens,
                           "IFN", "other"))
corr_input.filtered.df$Groups = my_groups2

corr_input.filtered.df <- corr_input.filtered.df[order(corr_input.filtered.df$Groups, decreasing = T),]

# add symbol column
corr_input.filtered.df$symbol <- replace_annotation_keys(rownames(corr_input.filtered.df),from = 'ENSEMBL', to = 'SYMBOL')

# Relabel gene PPP1R15A with its synonym GADD34
corr_input.filtered.df$symbol[corr_input.filtered.df$symbol == 'PPP1R15A'] <- 'GADD34'

p3 <- ggplot(corr_input.filtered.df, aes(x=Log2_pIC, y=Log2_25A, label=symbol)) + 
  labs(title = "Comparison of Log2FC between pIC WT and 25A WT treated cells") +
  geom_point(aes(colour=Groups), size = 3, alpha = 0.5)  + 
  xlim(-2.5,10) + ylim(-2.5,10) +
  xlab(bquote('Poly I:C WT vs control WT'* ~Log[2]*' fold change')) +
  ylab(bquote('2-5A WT vs control WT'* ~Log[2]*' fold change')) +
  geom_abline(slope = 1, intercept = 0, col = "black", size=0.5, linetype="dashed") +
  theme_minimal() + scale_color_manual(values = c("red",
                                                  "blue",
                                                  "grey")) +
  geom_text_repel(data = subset(corr_input.filtered.df, symbol %in% c('GADD34')),
                  nudge_y = 3, segment.size  = 0.6, segment.color = "grey50") +
  geom_point(data = subset(corr_input.filtered.df, symbol %in% c('GADD34')), 
             color='black', size = 4, alpha = 1)  + 
  scale_x_continuous(breaks=seq(-4,15,2)) +
  scale_y_continuous(breaks=seq(-4,15,2)) +
  theme(
    legend.direction = "vertical",  
    legend.position = c(.15, .95),
    legend.justification = c("left", "top"),
    legend.box.just = "right",
    legend.margin = margin(3, 3, 3, 3),
    legend.text = element_text(size = 16),
    axis.line = element_line(colour = "black"),
    
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color="black", size=0.5, linetype=1), 
  )
ggsave(filename = "pICwt_vs_25Awt_Log2FC_correlation_JNK_IFN_hallmark_all_filtered.pdf", plot = p3, device = "pdf", path = "./Plots/", width = 8, height = 8)


# ------------- violin plots -------------

#filtering step
padj <- 0.05
log2fc_treshold_up <- 1
log2fc_treshold_down <- -1

# creating list for p38/JNK used in violoin plot (Figure 2D and 3D)
pval_wt_25A <- row.names(res_wt25A_vs_wtNone[which(res_wt25A_vs_wtNone$padj <= padj), ])
fc_up_wt_25A <- row.names(res_wt25A_vs_wtNone[which(res_wt25A_vs_wtNone$log2FoldChange >=  log2fc_treshold_up), ])
up_wt_25A <- intersect(pval_wt_25A, fc_up_wt_25A)
JNK_p38_25A_up <- intersect(up_wt_25A, JNK_p38)

WT_25A <-as.data.frame(res_wt25A_vs_wtNone) 
WT_25A <- WT_25A %>% mutate(ID=row.names(WT_25A))
WT_IC <-as.data.frame(res_wtIC_vs_wtNone) 
WT_IC <- WT_IC %>% mutate(ID=row.names(WT_IC))
KO_IC <-as.data.frame(res_koIC_vs_koNone) 
KO_IC <- KO_IC %>% mutate(ID=row.names(KO_IC))
KO_25A <-as.data.frame(res_ko25A_vs_koNone) 
KO_25A <- KO_25A %>% mutate(ID=row.names(KO_25A))

# Fig. 2D: Combined JNK/p38

WT_25A_JNK_p38_filtered <- merge(WT_25A,as.data.frame(JNK_p38_25A_up), by.x="ID", by.y="JNK_p38_25A_up")
WT_IC_JNK_p38_filtered <- merge(WT_IC,as.data.frame(JNK_p38_25A_up), by.x="ID", by.y="JNK_p38_25A_up")
KO_IC_JNK_p38_filtered <- merge(KO_IC,as.data.frame(JNK_p38_25A_up), by.x="ID", by.y="JNK_p38_25A_up")
KO_25A_JNK_p38_filtered <- merge(KO_25A,as.data.frame(JNK_p38_25A_up), by.x="ID", by.y="JNK_p38_25A_up")
sum_JNK_p38 <- data.frame(KO_25A_JNK_p38_filtered$log2FoldChange,WT_25A_JNK_p38_filtered$log2FoldChange, KO_IC_JNK_p38_filtered$log2FoldChange, WT_IC_JNK_p38_filtered$log2FoldChange) 
write.csv(sum_JNK_p38, "./Tables/filtered_JNK_FC.csv")

### add this to the functions?
setEPS()
postscript( "./Plots/violinplot_JNK_p38_combined_filtered.eps")
vioplot(sum_JNK_p38,
        col = "white",         # Color of the area
        rectCol = "grey",      # Color of the rectangle
        lineCol = "black",     # Color of the line
        colMed = "red",        # Pch symbol color
        border = "black",      # Color of the border of the violin
        pchMed = 16,           # Pch symbol for the median
        plotCentre = "points", # If "line", plots a median line
        names=str_replace(string = names(sum_JNK_p38), 
                          pattern = '_JNK_p38_filtered.log2FoldChange',
                          replacement = ''), ylab = "Log2FC", main = "JNK/p38 pathway"
        ) 
        
dev.off()

# t-test 

ttest <- 
  rbind(
        tidy(t.test(sum_JNK_p38$KO_25A_JNK_p38_filtered.log2FoldChange, sum_JNK_p38$WT_25A_JNK_p38_filtered.log2FoldChange, paired=TRUE)) %>% mutate(comparison = "JNK/p38  2-5A_KO vs 2-5A_WT"),
        tidy(t.test(sum_JNK_p38$WT_IC_JNK_p38_filtered.log2FoldChange, sum_JNK_p38$WT_25A_JNK_p38_filtered.log2FoldChange, paired=TRUE)) %>% mutate(comparison = "JNK/p38  pIC_WT vs 2-5A_WT"),
        tidy(t.test(sum_JNK_p38$WT_IC_JNK_p38_filtered.log2FoldChange, sum_JNK_p38$KO_IC_JNK_p38_filtered.log2FoldChange, paired=TRUE)) %>% mutate(comparison = "JNK/p38  pIC_WT vs pIC_KO")
  )
colnames(ttest) <- gsub(pattern = 'parameter', replacement = 'degrees of freedom', x = colnames(ttest))
write_csv(x = ttest, file = "./Tables/t.test_JNK_p38.csv")

# Fig. 2D: IFN
# filter IFN list it for genes that are DE in WT poly I:C
pval_wt_IC <- row.names(res_wtIC_vs_wtNone[which(res_wtIC_vs_wtNone$padj <= padj), ])
fc_up_wt_IC <- row.names(res_wtIC_vs_wtNone[which(res_wtIC_vs_wtNone$log2FoldChange >=  log2fc_treshold_up), ])
up_wt_IC <- intersect(pval_wt_IC, fc_up_wt_IC)
IFN_IC_up <- intersect(up_wt_IC, IFN_gamma_ens)

WT_25A_IFN_filtered <- merge(WT_25A,as.data.frame(IFN_IC_up), by.x="ID", by.y="IFN_IC_up")
WT_IC_IFN_filtered <- merge(WT_IC,as.data.frame(IFN_IC_up), by.x="ID", by.y="IFN_IC_up")
KO_IC_IFN_filtered <- merge(KO_IC,as.data.frame(IFN_IC_up), by.x="ID", by.y="IFN_IC_up")
KO_25A_IFN_filtered <- merge(KO_25A,as.data.frame(IFN_IC_up), by.x="ID", by.y="IFN_IC_up")
sum_IFN_hallmark <- data.frame(KO_25A_IFN_filtered$log2FoldChange,WT_25A_IFN_filtered$log2FoldChange, KO_IC_IFN_filtered$log2FoldChange, WT_IC_IFN_filtered$log2FoldChange) 

setEPS()
postscript( "./Plots/violinplot_IFN_gamma_hallmark.eps")
vioplot(sum_IFN_hallmark,col = "white",               # Color of the area
        rectCol = "grey",                             # Color of the rectangle
        lineCol = "black",                            # Color of the line
        colMed = "red",                               # Pch symbol color
        border = "black",                             # Color of the border of the violin
        pchMed = 16,                                  # Pch symbol for the median
        plotCentre = "points",                        # If "line", plots a median line
        names=str_replace(string = names(sum_IFN_hallmark), 
                          pattern = '_IFN_filtered.log2FoldChange',
                          replacement = ''), ylab = "Log2FC", main = "IFN-gamma hallmark"
        ) 
dev.off()

ttest <- 
  rbind(
    tidy(t.test(sum_IFN_hallmark$KO_25A_IFN_filtered.log2FoldChange, sum_IFN_hallmark$WT_25A_IFN_filtered.log2FoldChange, paired=TRUE)) %>% mutate(comparison = "IFN  2-5A_KO vs 2-5A_WT"),
    tidy(t.test(sum_IFN_hallmark$WT_IC_IFN_filtered.log2FoldChange, sum_IFN_hallmark$WT_25A_IFN_filtered.log2FoldChange, paired=TRUE)) %>% mutate(comparison = "IFN  pIC_WT vs 2-5A_WT"),
    tidy(t.test(sum_IFN_hallmark$WT_IC_IFN_filtered.log2FoldChange, sum_IFN_hallmark$KO_IC_IFN_filtered.log2FoldChange, paired=TRUE)) %>% mutate(comparison = "IFN  pIC_WT vs pIC_KO")
  )
write_csv(x = ttest, file = "./Tables/t.test_IFN.csv", append = TRUE)

# ------------- Gene ontology enrichment analysis -------------

# GO for 2-5A only (Fig. 1C)
GO_25A  <- DE_results$wt25A_vs_wtNone
original_gene_list <- select(as.data.frame(GO_25A), log2FoldChange)
gene_list <- na.omit(original_gene_list$log2FoldChange)
names(gene_list) = rownames(original_gene_list)
gene_list <- sort(gene_list, decreasing = TRUE)

gse <- gseGO(geneList=gene_list,
             ont ="MF",
             keyType = 'ENSEMBL',
             nPerm = 10000,
             minGSSize = 3,
             maxGSSize = 800,
             pvalueCutoff = 0.05,
             verbose = TRUE,
             OrgDb = org.Hs.eg.db,
             pAdjustMethod = "BH")

write.csv(as.data.frame(gse), "./Tables/GO_25A.csv")

# ------------- RNase A analysis -------------

#import data set
all2 <- read.csv("./data/RNA_count_RNaseA.csv", header = TRUE, row.names = 1, comment.char = c("#"))
all2 <-round(all2)

#Remove all zero rows
if (remove_genes_with_zero_counts){
  all2 <- remove_all_zero_rows(all2, min_total_count = 0)
}

# Load metadata
metadata2 <- read.delim2("./data/RNAseqA_metadata.txt", sep = "\t", row.names = 1, header = T)

# Keep only samples that are present in all
metadata2 <- metadata2[colnames(all2),]

# Add total read counts and sample id columns to metadata
metadata2 <- cbind(metadata2, Read_counts =colSums(all2), Sample_id = rownames(metadata2))

# ------------- Analysis of expression data using DESeq2 - CDS RNase A -------------

# These are all high read counts
for (variable in c("Treatment")){
  metadata2[,variable] <- as.factor(metadata2[,variable])  
}


dds.RNaseA <- DESeqDataSetFromMatrix(countData = all2,
                                      colData = metadata2,
                                      design = ~Treatment)

# Normalize counts
vsd.RNaseA <- vst(dds.RNaseA, blind=FALSE)

# Keep genes with > 0 reads total across samples
keep <- rowSums(counts(dds.RNaseA)) >= N_reads
dds.RNaseA <- dds.RNaseA[keep,]

dds.RNaseA <- DESeq(dds.RNaseA)
res_RNaseA_vs_BSA <- results(dds.RNaseA)
res_RNaseA_vs_BSA.logfc_sorted <- sort_and_write_res_table(res_RNaseA_vs_BSA, "DE_RNaseA_vs_BSA_cds")

# filter for DE genes > 2-fold change upregulated and p value <0.05
rnasea.g <- rownames(res_RNaseA_vs_BSA.logfc_sorted[!is.na(res_RNaseA_vs_BSA.logfc_sorted$padj) & res_RNaseA_vs_BSA.logfc_sorted$padj <= 0.05 & res_RNaseA_vs_BSA.logfc_sorted$log2FoldChange > 1, ])
write_csv(as.data.frame(rnasea.g), "./DE/RNaseA_DEgenes.csv")

df.ns3 <- as.data.frame(colData(dds.RNaseA)[,c("Type", "Genotype")])

#  ------------- Heatmaps -------------

# For highest in RNaseA samples
plot_heat_map(my_vstd = vsd.RNaseA, gene_list = rnasea.g[1:60], file_name = "Plots/top60_RNaseA_clustered", variables = select(df.ns3, Type))

# Highest in RNase L samples
plot_heat_map(my_vstd = vsd.RNaseA, gene_list = wt.25A.g[1:60], file_name = "Plots/top60_RNaseA_top25A_clustered", variables = select(df.ns3, Type))

# ------------- Plot Venn diagram Fig 3B comparing RNAseL vs RNAseA -------------

# Venn-diagram showing the number of upregulated genes (defined as padjusted value <0.05 and log2fold changes >1) in RNase A electroporated samples (purple) and 2-5A treated cells (yellow)

# Genes significantly DE in WT and induced by 25A
wt.25A.up.g <- fetch_significant_genes(deseq_res = res_wt25A_vs_wtNone.logfc_sorted, adj_p = 0.05, change_direction = 'up')
wt.rnaseA.up.g <- fetch_significant_genes(deseq_res = res_RNaseA_vs_BSA.logfc_sorted, adj_p = 0.05, change_direction = 'up')

VennDiagram::venn.diagram(x = list(wt.rnaseA.up.g, wt.25A.up.g ),
                          category.names = c("RNase A", "RNase L"),
                          filename = "./Plots/RNaseA_vs_RNAseL_DE_venn_0.05.png",
                          imagetype = "png",
                          output=TRUE,
                          height = 480 , 
                          width = 480 , 
                          resolution = 300,
                          compression = "lzw",
                          lwd = 1,
                          col=c("#440154ff", '#fde725ff'),
                          fill = c(alpha("#440154ff",0.3), alpha('#fde725ff',0.3)),
                          cex = 0.5,
                          fontfamily = "sans",
                          cat.cex = 0.3,
                          cat.default.pos = "outer",
                          cat.pos = c(-27, 27),
                          cat.dist = c(0.055, 0.055),
                          cat.fontfamily = "sans",
                          cat.col = c("#440154ff", '#fde725ff')
                          )

#  ------------- violin plot  -------------

RNaseA <-as.data.frame(res_RNaseA_vs_BSA ) 
RNaseA <- RNaseA %>% mutate(ID=row.names(RNaseA))
RNaseA_JNK_p38_filtered <- merge(RNaseA,as.data.frame(JNK_p38_25A_up), by.x="ID", by.y="JNK_p38_25A_up")
sum_JNK_p38_RNaseA <- data.frame(KO_25A_JNK_p38_filtered$log2FoldChange,WT_25A_JNK_p38_filtered$log2FoldChange, RNaseA_JNK_p38_filtered$log2FoldChange) 

# Fig. 3D
setEPS()
postscript( "./Plots/violinplot_JNK_p38_RNaseA_filtered.eps")
vioplot(sum_JNK_p38_RNaseA,col = "white",    # Color of the area
        rectCol = "grey",                    # Color of the rectangle
        lineCol = "black",                   # Color of the line
        colMed = "red",                      # Pch symbol color
        border = "black",                    # Color of the border of the violin
        pchMed = 16,                         # Pch symbol for the median
        plotCentre = "points", ylab="Log2FC", main="JNK/p38 pathway (RNaseA-treatment)",
        names=str_replace(string = names(sum_JNK_p38_RNaseA), 
                          pattern = '_JNK_p38_filtered.log2FoldChange',
                          replacement = '')
)
dev.off()

# t-test 
RNAseqA.ttest <- 
rbind(tidy(t.test(sum_JNK_p38_RNaseA$KO_25A_JNK_p38_filtered.log2FoldChange, sum_JNK_p38_RNaseA$RNaseA_JNK_p38_filtered.log2FoldChange, paired=TRUE)) %>% mutate(comparison = "JNK/p38 RNAseA vs 2-5A_KO"),
tidy(t.test(sum_JNK_p38_RNaseA$WT_25A_JNK_p38_filtered.log2FoldChange, sum_JNK_p38_RNaseA$RNaseA_JNK_p38_filtered.log2FoldChange, paired=TRUE)) %>% mutate(comparison = "JNK/p38 RNAseA vs 2-5A_WT"))

colnames(RNAseqA.ttest) <- gsub(pattern = 'parameter', replacement = 'degrees of freedom', x = colnames(RNAseqA.ttest))
write_csv(x = RNAseqA.ttest, file = "./Tables/t.test_JNK_p38.csv", append = TRUE)

# correlation
DE_RNaseA_vs_BSA_cds <- read.delim("./DE/DE_RNaseA_vs_BSA_cds.txt")
DE_wt25A_vs_wtNone_cds <- read.delim("./DE/DE_wt25A_vs_wtNone_cds.txt")
res_RNaseA_25A <- merge(DE_RNaseA_vs_BSA_cds, DE_wt25A_vs_wtNone_cds, "X")

# Fig. 3C
RNaseA_25A.corr <- tidy(cor.test(res_RNaseA_25A$log2FoldChange.x, res_RNaseA_25A$log2FoldChange.y, method="pearson"))
write_csv(RNaseA_25A.corr, file = "./Tables/corr_RNaseA_25A.csv")

my_groups3 = ifelse(res_RNaseA_25A$X %in% JNK_p38,
                    "JNK/p38", "other")
res_RNaseA_25A$Groups = my_groups3
res_RNaseA_25A <- res_RNaseA_25A[order(res_RNaseA_25A$Groups,decreasing = T),]

p5 <- ggplot(res_RNaseA_25A, aes(x=log2FoldChange.x, y=log2FoldChange.y)) + 
  labs(title = "Comparison of Log2FC between RNAseA-vs-BSA (RNAse L KO cells) and 2-5A-vs-control (WT  cells)") +
  geom_point(aes(colour=Groups), size = 3, alpha = 0.5)  + 
  xlab(bquote('RNAse A vs BSA'* ~Log[2]*' fold change')) +
  ylab(bquote('2-5A WT vs control WT'* ~Log[2]*' fold change')) +
  geom_abline(slope = 1, intercept = 0, col = "black", size=0.5, linetype="dashed") +
  theme_minimal() + scale_color_manual(values = c(
                                                  "blue",
                                                  "grey")) +  
  scale_x_continuous(breaks=seq(-8,15,2)) +
  scale_y_continuous(breaks=seq(-8,15,2)) +
  theme(
    legend.direction = "vertical",  
    legend.position = c(.25, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(3, 3, 3, 3),
    legend.text = element_text(size = 16),
    axis.line = element_line(colour = "black"),
    
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color="black", size=0.5, linetype=1),
  )

ggsave(filename = "RNaseA_vs_25Awt_Log2FC_correlation_JNK_IFN_hallmark_all_filtered.pdf", plot = p5, device = "pdf", path = "./Plots/", width = 8, height = 8)
p5

# volcano plot
generate_volcano_plot_simple(res.tmp = res_RNaseA_vs_BSA, my_file_name = "./Plots/res_RNaseA_vs_BSA")

sessionInfo()