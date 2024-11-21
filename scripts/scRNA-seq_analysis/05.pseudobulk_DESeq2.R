### https://hbctraining.github.io/scRNA-seq/lessons/pseudobulk_DESeq2_scrnaseq.html
### https://github.com/hbctraining/scRNA-seq/blob/master/lessons/pseudobulk_DESeq2_scrnaseq.md

# Load libraries
library(scater)
library(Seurat)
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(dplyr)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)
library(ggthemes)

# Bring in Seurat object    # normal　と　TumorのうちinferCNVScoreでMalignant判定だったもの (02を使用)
seurat <- readRDS("../Harmony_RNA_round1.epi.inferCNVScore02.subcluster.rds")

hpv_expression <- FetchData(seurat, vars=c("HPV16-E1", "HPV16-E2", "HPV16-E5", "HPV16-L2", "HPV16-L1", "HPV16-E6", "HPV16-E7"))
hpv_expression$HPV <- hpv_expression$`HPV16-E1` + hpv_expression$`HPV16-E2` + hpv_expression$`HPV16-E5` + hpv_expression$`HPV16-L2` + hpv_expression$`HPV16-L1` + hpv_expression$`HPV16-E6` + hpv_expression$`HPV16-E7`
seurat$HPV_group <- ifelse(hpv_expression$HPV > 0, "positive", "negative")

seurat
##  An object of class Seurat
##  25523 features across 11608 samples within 1 assay
##  Active assay: RNA (25523 features, 0 variable features)
##   1 dimensional reduction calculated: umap

table(seurat@meta.data %>% select(Sample, epithelial_infercnvScore))
##            epithelial_infercnvScore
##  Sample     Normal Tumor
##    OP13         19    27
##    OP14        208    96
##    OP17        157   142
##    OP20        277  2053
##    OP33        476  2708
##    OP33norm     92     0
##    OP34        343  1271
##    OP34norm     99     0
##    OP35BOT      30   211
##    OP35LN        0    88
##    OP35norm     24     0
##    OP4         301  1189
##    OP5          35   281
##    OP6          49   283
##    OP9         526   623

table(seurat@meta.data %>% select(HPV_group, epithelial_infercnvScore))
##            epithelial_infercnvScore
##  HPV_group  Normal Tumor
##    negative   2486  3012
##    positive    150  5960


seurat.sub <- seurat
seurat.sub@meta.data <- seurat.sub@meta.data %>% unite(sample_new, Sample, epithelial_infercnvScore, sep="_", remove=FALSE)

seurat.sub
##  An object of class Seurat
##  25523 features across 11608 samples within 1 assay
##  Active assay: RNA (25523 features, 0 variable features)
##   1 dimensional reduction calculated: umap

### 除外
seurat.sub <- subset(x = seurat.sub, subset = subclusters %in% c("OP33norm_s1", "OP33norm_s3", "OP34norm_s2", "OP34norm_s3", "OP35LN_s4", "OP6_s4", "OP9_s1", "OP13_s1"), invert = TRUE)

seurat.sub
##  An object of class Seurat
##  25523 features across 10965 samples within 1 assay
##  Active assay: RNA (25523 features, 0 variable features)
##   1 dimensional reduction calculated: umap


### covariate用のsampleID
seurat.sub@meta.data <- seurat.sub@meta.data %>% mutate(
    cov_sample = case_when(
        Sample=="OP33norm" ~ "OP33",
        Sample=="OP34norm" ~ "OP34",
        Sample=="OP35norm" ~ "OP35",
        Sample=="OP35BOT" ~ "OP35",
        Sample=="OP35LN" ~ "OP35",
        TRUE ~ Sample))

table(seurat.sub@meta.data %>% select(Sample, cov_sample))
##            cov_sample
##  Sample     OP13 OP14 OP17 OP20 OP33 OP34 OP35  OP4  OP5  OP6  OP9
##    OP13       19    0    0    0    0    0    0    0    0    0    0
##    OP14        0  304    0    0    0    0    0    0    0    0    0
##    OP17        0    0  299    0    0    0    0    0    0    0    0
##    OP20        0    0    0 2330    0    0    0    0    0    0    0
##    OP33        0    0    0    0 3184    0    0    0    0    0    0
##    OP33norm    0    0    0    0   39    0    0    0    0    0    0
##    OP34        0    0    0    0    0 1614    0    0    0    0    0
##    OP34norm    0    0    0    0    0   40    0    0    0    0    0
##    OP35BOT     0    0    0    0    0    0  241    0    0    0    0
##    OP35LN      0    0    0    0    0    0   75    0    0    0    0
##    OP35norm    0    0    0    0    0    0   24    0    0    0    0
##    OP4         0    0    0    0    0    0    0 1490    0    0    0
##    OP5         0    0    0    0    0    0    0    0  316    0    0
##    OP6         0    0    0    0    0    0    0    0    0  283    0
##    OP9         0    0    0    0    0    0    0    0    0    0  707

table(seurat.sub$epithelial_infercnvScore)
##  Normal  Tumor
##    2475   8490

# Extract raw counts and metadata to create SingleCellExperiment object
counts <- seurat.sub@assays$RNA@counts 

metadata <- seurat.sub@meta.data


sample_celltype <- table(seurat.sub@meta.data$Sample, seurat.sub@meta.data$predicted.celltype.l1)
write.csv(sample_celltype,
    "UMI_per_Celltype_per_Sample.csv",
    quote=FALSE,
    row.names=TRUE)

# Set up metadata as desired for aggregation and DE analysis
#seurat.sub <- SetIdent(seurat.sub, value = "predicted.celltype.l1")
#metadata$cluster_id <- factor(seurat@active.ident)
metadata$cluster_id <- seurat.sub@active.ident %>% str_replace_all(pattern="_", replacement=" ") %>% factor()

metadata$sample_id <- factor(metadata$subclusters)
metadata$group_id <- seurat.sub$epithelial_infercnvScore %>% str_replace_all(pattern="-", replacement="_") %>% factor()
#metadata$institute_id <- factor(seurat$Institute)
metadata$cov_sample <- factor(seurat.sub$cov_sample)

# Create single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = counts), 
                           colData = metadata)





# Explore the raw counts for the dataset

## Check the assays present
assays(sce)

## Explore the raw counts for the dataset
dim(counts(sce))

counts(sce)[1:6, 1:6]



# Perform QC if not already performed
dim(sce)
##  [1] 25523 10965
table(sce$group_id)
##  Normal  Tumor
##    2475   8490


### sampleがDROPすることでfactorの数が減るがlevelに残ってしまっている (BCLLATLAS_06_d8kwy76j_7blj9otf)
# 現在のファクターレベルを確認
levels(sce$sample_id)
##   [1] "OP13"    "OP14"    "OP17"    "OP20"    "OP33"    "OP34"    "OP35BOT"
##   [8] "OP35LN"  "OP4"     "OP5"     "OP6"     "OP9"

# droplevels() で不要なレベルを削除
sce$sample_id <- droplevels(sce$sample_id)

# 再度ファクターレベルを確認
levels(sce$sample_id)
##   [1] "OP13"    "OP14"    "OP17"    "OP20"    "OP33"    "OP34"    "OP35BOT"
##   [8] "OP35LN"  "OP4"     "OP5"     "OP6"     "OP9"




## Explore the cellular metadata for the dataset
dim(colData(sce))
##  [1] 4879    36

head(colData(sce))[1:6,1:6]


### Acquiring necessary metrics for aggregation across cells in a sample

# Named vector of cluster names
cids <- purrr::set_names(levels(sce$cluster_id))
cids
##    Tumor
##  "Tumor"

# Total number of clusters
nc <- length(cids)
nc
##  [1] 1

# Named vector of sample names
sids <- purrr::set_names(levels(sce$sample_id))

# Total number of samples 
ns <- length(sids)
ns
##  [1] 6


### To perform sample-level differential expression analysis, we need to generate sample-level metadata. To do this, we will reorder samples in the single-cell metadata to match the order of the factor levels of the sample ID, then extract only the sample-level information from the first cell corresponding to that sample.
# Generate sample level metadata

## Determine the number of cells per sample
table(sce$sample_id)

## Turn named vector into a numeric vector of number of cells per sample
n_cells <- as.numeric(table(sce$sample_id))

## Determine how to reoder the samples (rows) of the metadata to match the order of sample names in sids vector
m <- match(sids, sce$sample_id)

## Create the sample level metadata by combining the reordered metadata with the number of cells corresponding to each sample.
ei <- data.frame(colData(sce)[m, ], 
                  n_cells, row.names = NULL) %>% 
                select(-"cluster_id")
ei[1:6,1:6]


### Count aggregation to sample level
# Aggregate the counts per sample_id and cluster_id

# Subset metadata to only include the cluster and sample IDs to aggregate across
groups <- colData(sce)[, c("cluster_id", "sample_id")]
groups

# Aggregate across cluster-sample groups
pb <- aggregate.Matrix(t(counts(sce)), 
                       groupings = groups, fun = "sum") 

class(pb)

dim(pb)
##  [1]    31 11050

pb[1:6, 1:6]

### To perform DE analysis on a per cell type basis, we need to wrangle our data in a couple ways. We need to do the following steps:
## 1. Split our data by cell type
## 2. Transform the matrix so that the genes are the row names and the samples are the column names

# Not every cluster is present in all samples; create a vector that represents how to split samples
splitf <- sapply(stringr::str_split(rownames(pb), 
                                    pattern = "_",
                                    n = 2),
                `[`, 1)
head(splitf)

pb <- split.data.frame(pb, 
                       factor(splitf)) %>%
        lapply(function(u) 
                set_colnames(t(u), 
                             sapply(stringr::str_split(rownames(u), pattern = "_", n = 2), `[`, 2)))

class(pb)

# Explore the different components of list
str(pb)
##  List of 1
##   $ TRUE:Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
##    .. ..@ i       : int [1:295169] 15 17 20 21 23 25 26 28 33 34 ...
##    .. ..@ p       : int [1:36] 0 10030 11000 12253 14837 19055 21803 26806 33558 37027 ...
##    .. ..@ Dim     : int [1:2] 31844 35
##    .. ..@ Dimnames:List of 2
##    .. .. ..$ : chr [1:31844] "MIR1302-2HG" "AL627309.1" "AL627309.3" "AL627309.4" ...
##    .. .. ..$ : chr [1:35] "BCLLATLAS_05_izi9unx1_8qdzhivu" "BCLLATLAS_06_d8kwy76j_7blj9otf" "BCLLATLAS_06_giz3qso4_783dbpu6" "BCLLATLAS_07_i5udk3x0_57gv6ncx" ...
##    .. ..@ x       : num [1:295169] 1 2 6 64 4 6 3 1 6 1 ...
##    .. ..@ factors : list()


# Print out the table of cells in each cluster-sample group
options(width = 100)
table(sce$cluster_id, sce$sample_id)[1,1:6]






### Differential gene expression with DESeq2
# DESeq2では、まずサンプル間のライブラリサイズやRNA組成の違いを考慮し、カウントデータを正規化します。次に、正規化されたカウントを使用して、遺伝子レベルおよびサンプルレベルでのQCのためのいくつかのプロットを作成します。最後のステップは、DESeq2パッケージの適切な関数を使用して、差分発現解析を実行することです。


### Sample-level metadata
# Get sample names for each of the cell type clusters

# prep. data.frame for plotting
get_sample_ids <- function(x){
        pb[[x]] %>%
                colnames()
}

de_samples <- map(1:length(cids), get_sample_ids) %>%
        unlist()
head(de_samples)


# Get cluster IDs for each of the samples

samples_list <- map(1:length(cids), get_sample_ids)

get_cluster_ids <- function(x){
        rep(names(pb)[x], 
            each = length(samples_list[[x]]))
}

de_cluster_ids <- map(1:length(cids), get_cluster_ids) %>%
        unlist()
head(de_cluster_ids)

# Create a data frame with the sample IDs, cluster IDs and condition
#factorにしておかないと下のlevels()がNULLを返してしまう
gg_df <- data.frame(cluster_id = factor(de_cluster_ids),
                    sample_id = factor(de_samples))

gg_df <- left_join(gg_df, ei[, c("sample_id", "group_id", "cov_sample")]) 


metadata <- gg_df %>%
        dplyr::select(cluster_id, sample_id, group_id, cov_sample)

head(metadata)     


write.csv(metadata, "metadata_for_DESeq2.csv", quote=FALSE, row.names=FALSE)


saveRDS(pb, "pb.rds")
saveRDS(metadata, "metadata.rds")
# pb <- readRDS("pb.rds")
# metadata <- readRDS("metadata.rds")





### Subsetting dataset to cluster(s) of interest
# Generate vector of cluster IDs
clusters <- levels(metadata$cluster_id)
clusters


#  Let’s perform the DE analysis on "B naive"
clusters[1]

# Subset the metadata to only the B cells                                       ### SG_SL_KR_00011ではB naiveなし
cluster_metadata <- metadata[which(metadata$cluster_id == clusters[1]), ]
head(cluster_metadata)
nrow(cluster_metadata)


# Assign the rownames of the metadata to be the sample IDs
rownames(cluster_metadata) <- cluster_metadata$sample_id
head(cluster_metadata)

# Subset the counts to only the B cells
counts <- pb[[clusters[1]]]
head(counts)

cluster_counts <- data.frame(counts[, which(colnames(counts) %in% rownames(cluster_metadata))])
head(cluster_counts)

# Check that all of the row names of the metadata are the same and in the same order as the column names of the counts in order to use as input to DESeq2
all(rownames(cluster_metadata) == colnames(cluster_counts))        
##  [1] TRUE


#------------------------#
### filtering by count ###
#------------------------#
summary(rowSums(cluster_counts))

filtered_counts <- cluster_counts

#   #-----------------------#
#   ### Yuruでは下記を飛ばす ###
#   #-----------------------#
#   
#   ### Tumor/Normalで2以上のカウントを持つ遺伝子が少なくとも2サンプルである場合にフィルタリング
#   table(cluster_metadata$group_id)
#   ##  Normal  Tumor
#   ##      16     39
#   # Tumor/Normalを抽出
#   normal_n=nrow(cluster_metadata[cluster_metadata$group_id=="Normal",])
#   tumor_n=nrow(cluster_metadata[cluster_metadata$group_id=="Tumor",])
#   normal_counts <- filtered_counts[, 1:normal_n]
#   tumor_counts <- filtered_counts[, (normal_n+1):(normal_n+tumor_n)]
#   
#   # Tumor/Normalで2以上のカウントを持つ遺伝子が少なくとも半分のサンプルである遺伝子
#   normal_threshold <- rowSums(normal_counts >= 2) >= (normal_n / 2)
#   tumor_threshold <- rowSums(tumor_counts >= 2) >= (tumor_n / 2)
#   total_threshold <- rowSums(filtered_counts >= 2) >= (nrow(cluster_metadata) / 2)
#   
#   table(normal_threshold)
#   ##  FALSE  TRUE
#   ##   1022 11820
#   table(tumor_threshold)
#   ##  FALSE  TRUE
#   ##   1417 11425
#   table(normal_threshold | tumor_threshold)
#   ##  FALSE  TRUE
#   ##    137 11168
#   
#   #   # TumorまたはNormalのどちらかで条件を満たす遺伝子を抽出
#   #   filtered2_counts <- filtered_counts[normal_threshold | tumor_threshold, ]
#   
#   ## 今回はnormalで半数以上で発現がある遺伝子に絞るべきだろう
#   filtered2_counts <- filtered_counts[normal_threshold, ]


### Create DESeq2 object

dds <- DESeqDataSetFromMatrix(filtered_counts, 
                              colData = cluster_metadata, 
                              design = ~ group_id + cov_sample)
##  converting counts to integer mode
##              ### scale(metadata$age)
##              ##    the design formula contains one or more numeric variables with integer values,
##              ##    specifying a model with increasing fold change for higher values.
##              ##    did you mean for this to be a factor? if so, first convert
##              ##    this variable to a factor using the factor() function
##              ##    the design formula contains one or more numeric variables that have mean or
##              ##    standard deviation larger than 5 (an arbitrary threshold to trigger this message).
##              ##    Including numeric variables with large mean can induce collinearity with the intercept.
##              ##    Users should center and scale numeric variables in the design to improve GLM convergence.   scale() or 離散化（年代別にわける）

######################
### gene filtering ###  https://support.bioconductor.org/p/65256/
######################  ATMを残してCADM1やZBTB16が落ちそう
dds <- estimateSizeFactors(dds)

normal_n=nrow(cluster_metadata[cluster_metadata$group_id=="Normal",])
tumor_n=nrow(cluster_metadata[cluster_metadata$group_id=="Tumor",])
ifelse(normal_n <= tumor_n, small_samplesize <- normal_n, small_samplesize <- tumor_n)

keep_normalized_count = 10

keep <- rowSums( counts(dds, normalized=TRUE) >= keep_normalized_count ) >= small_samplesize
dds <- dds[keep,]

table(keep)
##  keep
##  FALSE  TRUE
##  17388  8135

saveRDS(keep, "keep.rds")


### ATM
# Normal
counts(dds, normalized=TRUE) %>% data.frame() %>% subset(rownames(data)=="ATM") %>% select(cluster_metadata[cluster_metadata$group_id=="Normal",]$sample_id)
counts(dds, normalized=TRUE) %>% data.frame() %>% subset(rownames(data)=="ATM") %>% select(cluster_metadata[cluster_metadata$group_id=="Normal",]$sample_id) %>% rowMeans()

# Tumor
counts(dds, normalized=TRUE) %>% data.frame() %>% subset(rownames(data)=="ATM") %>% select(cluster_metadata[cluster_metadata$group_id=="Tumor",]$sample_id)
counts(dds, normalized=TRUE) %>% data.frame() %>% subset(rownames(data)=="ATM") %>% select(cluster_metadata[cluster_metadata$group_id=="Tumor",]$sample_id) %>% rowMeans()

### ZBTB16
# Normal
counts(dds, normalized=TRUE) %>% data.frame() %>% subset(rownames(data)=="ZBTB16") %>% select(cluster_metadata[cluster_metadata$group_id=="Normal",]$sample_id)
counts(dds, normalized=TRUE) %>% data.frame() %>% subset(rownames(data)=="ZBTB16") %>% select(cluster_metadata[cluster_metadata$group_id=="Tumor",]$sample_id) %>% rowMeans()

# Tumor
counts(dds, normalized=TRUE) %>% data.frame() %>% subset(rownames(data)=="ZBTB16") %>% select(cluster_metadata[cluster_metadata$group_id=="Tumor",]$sample_id)
counts(dds, normalized=TRUE) %>% data.frame() %>% subset(rownames(data)=="ZBTB16") %>% select(cluster_metadata[cluster_metadata$group_id=="Tumor",]$sample_id) %>% rowMeans()

### CADM1
# Normal
counts(dds, normalized=TRUE) %>% data.frame() %>% subset(rownames(data)=="CADM1") %>% select(cluster_metadata[cluster_metadata$group_id=="Normal",]$sample_id)
counts(dds, normalized=TRUE) %>% data.frame() %>% subset(rownames(data)=="CADM1") %>% select(cluster_metadata[cluster_metadata$group_id=="Normal",]$sample_id) %>% rowMeans()

# Tumor
counts(dds, normalized=TRUE) %>% data.frame() %>% subset(rownames(data)=="CADM1") %>% select(cluster_metadata[cluster_metadata$group_id=="Tumor",]$sample_id)
counts(dds, normalized=TRUE) %>% data.frame() %>% subset(rownames(data)=="CADM1") %>% select(cluster_metadata[cluster_metadata$group_id=="Tumor",]$sample_id) %>% rowMeans()


data <- counts(dds, normalized=TRUE) %>% data.frame()
tumor <- data %>% select(cluster_metadata[cluster_metadata$group_id=="Tumor",]$sample_id)
normal <- data %>% select(cluster_metadata[cluster_metadata$group_id=="Normal",]$sample_id)
tumor <- t(tumor) %>% data.frame() %>% rownames_to_column(var = "SampleID")
normal <- t(normal) %>% data.frame() %>% rownames_to_column(var = "SampleID")
tumor <- tumor %>% pivot_longer(!SampleID, names_to="gene", values_to="normalized_count")
normal <- normal %>% pivot_longer(!SampleID, names_to="gene", values_to="normalized_count")
tumor$SampleGroup <- "Tumor"
normal$SampleGroup <- "Normal"
df <- rbind(tumor, normal)

df_select <- df[df$gene %in% c("ATM", "BIRC2", "BIRC3", "GNAI2", "CYLD", "TRAF3"),]

p <- ggplot(df_select, aes(x=gene, y=normalized_count, color=SampleGroup)) +
        geom_boxplot(outlier.shape=NA) +
        geom_point(position=position_jitterdodge()) +
        scale_color_manual(values=c("#0A74B2", "#B93F2B")) +
        theme_bw()
ggsave(plot=p, file="5genes.normalized_count.pdf")





### Quality Control - sample level

## Principal component analysis
# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)
saveRDS(rld, "rld.rds")

# Plot PCA


### PC3 vs PC4
plotPCA = function(object, intgroup="condition", ntop=500, returnData=FALSE, ax1="PC1", ax2="PC2")
{
  # calculate the variance for each gene
  rv <- rowVars(assay(object))

  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]

  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))

  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )

  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }

  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])
  
  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=":"))
  } else {
    colData(object)[[intgroup]]
  }

  # assembly the data for the plot
  ax1_n <- str_replace_all(ax1, pattern="PC", replacement="") %>% as.numeric()
  ax2_n <- str_replace_all(ax2, pattern="PC", replacement="") %>% as.numeric()
  d <- data.frame(ax1=pca$x[,ax1_n], ax2=pca$x[,ax2_n], group=group, intgroup.df, name=colnames(object))

  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  
  ggplot(data=d, aes_string(x="ax1", y="ax2", color="group")) + geom_point(size=3) + 
    xlab(paste0(ax1,": ",round(percentVar[ax1_n] * 100),"% variance")) +
      ylab(paste0(ax2,": ",round(percentVar[ax2_n] * 100),"% variance")) +
        coord_fixed() +
        ggrepel::geom_text_repel(aes(label=row.names(d)))
}


p <- plotPCA(rld, intgroup = "group_id", ax1="PC1", ax2="PC2")
ggsave(file="DESeq2_PCA12.TRUE.group.pdf", plot=p)



p <- plotPCA(rld, intgroup = "group_id", ax1="PC3", ax2="PC4")
ggsave(file="DESeq2_PCA34.TRUE.group.pdf", plot=p)



## Hierarchical clustering
# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)

# Plot heatmap
p <- pheatmap(rld_cor, annotation = cluster_metadata[, c("group_id"), drop=F])
ggsave(file="DESeq2_Clustering.TRUE.group.pdf", plot=p)


# ここで、除去が必要な外れ値や、デザイン式に回帰させたい追加の変動要因があるかどうかを判断します。PCAや階層的クラスタリングで外れ値が検出されず、回帰すべき追加の変動要因もないため、差分発現解析の実行に進むことができます。




### Running DESeq2
#DESeq2は、ライブラリの深さの違いを考慮するために正規化因子（サイズ因子）を使用して、生のカウントをモデル化します。次に、遺伝子ごとの分散を推定し、その推定値を縮小して、より正確な分散推定値を生成してカウントをモデル化します。最後に、DESeq2は負の二項モデルを適合させ、Wald検定または尤度比検定を使用して仮説検定を実行します。これらのステップはすべて、追加資料で詳細に説明されています。

# Run DESeq2 differential expression analysis
dds <- DESeq(dds)

# Plot dispersion estimates
pdf("DESeq2_DispersionEstimates.TRUE.pdf")
plotDispEsts(dds)
dev.off()
# このプロットは、平均値の増加とともに分散が減少し、ベストフィットの線に従うと予想されるため、勇気づけられるものです。



### Results
#さて、差分発現解析を行ったので、特定の比較について結果を調べることができます。目的の比較を示すために、コントラストを指定し、log2 fold changeの縮小を行う必要があります。
#刺激されたグループとコントロールの相対的な比較をしてみましょう：
# Output results of Wald test for contrast for stim vs ctrl
levels(cluster_metadata$group_id)[2]
##  [1] "Malignant"
levels(cluster_metadata$group_id)[1]
##  [1] "nonMalignant"

contrast <- c("group_id", levels(cluster_metadata$group_id)[2], levels(cluster_metadata$group_id)[1])
contrast
##  [1] "group_id" "Malignant"  "nonMalignant"

resultsNames(dds)
##      [1] "Intercept"                "group_id_Malignant_vs_nonMalignant"

res <- results(dds, 
               contrast = contrast,
               alpha = 0.05)

res

resShrink <- lfcShrink(dds,
                coef = resultsNames(dds)[2],
                res=res)

resShrink



### Table of results for all genes
# Turn the results object into a tibble for use with tidyverse functions
res_tbl <- res %>%
        data.frame() %>%
        rownames_to_column(var="gene") %>%
        as_tibble()
resShrink_tbl <- resShrink %>%
        data.frame() %>%
        rownames_to_column(var="gene") %>%
        as_tibble()

resShrink_tbl[resShrink_tbl$gene=="ATM",]

# Check results output
res_tbl

resShrink_tbl

# Write all results to file
cluster_name <- str_replace_all(clusters[1], " ", "_")
write.csv(res_tbl,
          paste0(cluster_name, "_", levels(cluster_metadata$group_id)[2], "_vs_", levels(cluster_metadata$group_id)[1], "_all_genes.csv"),
          quote = FALSE, 
          row.names = FALSE)
write.csv(resShrink_tbl,
          paste0(cluster_name, "_", levels(cluster_metadata$group_id)[2], "_vs_", levels(cluster_metadata$group_id)[1], "_all_genes_lfcShrink.csv"),
          quote = FALSE, 
          row.names = FALSE)




# Set thresholds
padj_cutoff <- 0.05

# Subset the significant results
sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
        dplyr::arrange(padj)

sig_resShrink <- dplyr::filter(resShrink_tbl, padj < padj_cutoff) %>%
        dplyr::arrange(padj)

# Check significant genes output
sig_res

sig_resShrink

# Write significant results to file
write.csv(sig_res,
          paste0(cluster_name, "_", levels(cluster_metadata$group_id)[2], "_vs_", levels(cluster_metadata$group_id)[1], "_sig_genes.csv"),
          quote = FALSE, 
          row.names = FALSE)
write.csv(sig_resShrink,
          paste0(cluster_name, "_", levels(cluster_metadata$group_id)[2], "_vs_", levels(cluster_metadata$group_id)[1], "_sig_genes_lfcShrink.csv"),
          quote = FALSE, 
          row.names = FALSE)



### Scatterplot of normalized expression of top 20 most significant genes
## ggplot of top genes
normalized_counts <- counts(dds, 
                            normalized = TRUE)

## Order results by padj values
top20_sig_genes <- res_tbl %>%
        dplyr::arrange(padj) %>%
        dplyr::pull(gene) %>%
        head(n=20)


top20_sig_norm <- data.frame(normalized_counts) %>%
        rownames_to_column(var = "gene") %>%
        dplyr::filter(gene %in% top20_sig_genes) %>% dplyr::arrange(match(gene, top20_sig_genes))

gathered_top20_sig <- top20_sig_norm %>%
        gather(colnames(top20_sig_norm)[2:length(colnames(top20_sig_norm))], key = "samplename", value = "normalized_counts")
        
gathered_top20_sig <- inner_join(ei[, c("sample_id", "group_id" )], gathered_top20_sig, by = c("sample_id" = "samplename"))

gathered_top20_sig$gene <- factor(gathered_top20_sig$gene, levels=top20_sig_genes)

## plot using ggplot2
p <- ggplot(gathered_top20_sig) +
        geom_point(aes(x = gene, 
                       y = normalized_counts, 
                       color = group_id), 
                   position=position_jitter(w=0.1,h=0)) +
        scale_y_log10() +
        xlab("Genes") +
        ylab("log10 Normalized Counts") +
        ggtitle("Top 20 Significant DE Genes") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
        theme(plot.title = element_text(hjust = 0.5))
ggsave(file=paste0(cluster_name, "_", levels(cluster_metadata$group_id)[2], "_vs_", levels(cluster_metadata$group_id)[1], "_top20sigDEG.pdf"), plot=p)



### Heatmap of all significant genes
# Extract normalized counts for only the significant genes
sig_norm <- data.frame(normalized_counts) %>%
        rownames_to_column(var = "gene") %>%
        dplyr::filter(gene %in% sig_res$gene)
        
# Set a color palette
heat_colors <- brewer.pal(6, "YlOrRd")

# Run pheatmap using the metadata data frame for the annotation
p <- pheatmap(sig_norm[ , 2:length(colnames(sig_norm))], 
    color = heat_colors, 
    cluster_rows = T, 
    show_rownames = F,
    annotation = cluster_metadata[, c("group_id", "cluster_id")], 
    border_color = NA, 
    fontsize = 10, 
    scale = "row", 
    fontsize_row = 10, 
    height = 20)        
ggsave(file=paste0(cluster_name, "_", levels(cluster_metadata$group_id)[2], "_vs_", levels(cluster_metadata$group_id)[1], "_sigDEG_Heatmap.pdf"), plot=p)
# Error in hclust(d, method = method) : must have n >= 2 objects to cluster





### Volcano plot of results
## Obtain logical vector where TRUE values denote padj values < 0.05 and fold change > 1.5 in either direction
res_table_thres <- res_tbl[complete.cases(res_tbl),] %>% 
                  mutate(threshold = padj < 0.05 & abs(log2FoldChange) > 1)
res_table_thres$label <- NA
res_table_thres$label[res_table_thres$threshold==TRUE] <- res_table_thres$gene[res_table_thres$threshold==TRUE]
    
## Volcano plot
p <- ggplot(res_table_thres) +
    geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
    xlab("log2 fold change") + 
    ylab("-log10 adjusted p-value") +
    theme(legend.position = "none",
          plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.title = element_text(size = rel(1.25))) +
    ggrepel::geom_text_repel(aes(x = log2FoldChange, y = -log10(padj), label=label),
                        size=3,
                        max.overlaps = 50) +
    theme_bw() +
    scale_color_manual(values=c("#034789", "#E70808"))
ggsave(file=paste0(cluster_name, "_", levels(cluster_metadata$group_id)[2], "_vs_", levels(cluster_metadata$group_id)[1], "_Volcano.pdf"), plot=p)


resShrink_table_thres <- resShrink_tbl[complete.cases(resShrink_tbl),] %>% 
                  mutate(threshold = padj < 0.05 & abs(log2FoldChange) > 1)
resShrink_table_thres$label <- NA
resShrink_table_thres$label[resShrink_table_thres$threshold==TRUE] <- resShrink_table_thres$gene[resShrink_table_thres$threshold==TRUE]
                  
## Volcano plot
p <- ggplot(resShrink_table_thres) +
    geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
    xlab("log2 fold change") + 
    ylab("-log10 adjusted p-value") +
    theme(legend.position = "none",
          plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.title = element_text(size = rel(1.25))) +
    ggrepel::geom_text_repel(aes(x = log2FoldChange, y = -log10(padj), label=label),
                        size=3,
                        max.overlaps = 50) +
    theme_bw() +
    scale_color_manual(values=c("#034789", "#E70808"))
ggsave(file=paste0(cluster_name, "_", levels(cluster_metadata$group_id)[2], "_vs_", levels(cluster_metadata$group_id)[1], "_Volcano_lfcShrink.pdf"), plot=p)


### HI genes
data <- resShrink_tbl

df <- read.table("/mnt/d/Visium/Seurat_HNSCC/pseudobulk_DESeq2_epithelial_infercnvScore02_cd45n_countFilterYuru/HI.genes", header=FALSE)
df02 <- read.table("/mnt/d/Visium/Seurat_HNSCC/pseudobulk_DESeq2_epithelial_infercnvScore02_cd45n_countFilterYuru/PCGR.TSG_HI.tsv", header=TRUE)
df02$SYMBOL_gene <- str_split(df02$SYMBOL, '_', simplify=TRUE)[,1]
df02 <- df02 %>% select(SYMBOL_gene, SYMBOL)

df <- dplyr::inner_join(df02, df, by=c("SYMBOL_gene"="V1"))

library(tidyverse)
HI_data <- dplyr::left_join(df, data, by=c("SYMBOL_gene"="gene"))

# remove NA
HI_data <- HI_data[complete.cases(HI_data),]

# filtering 38 TSGsと85 signaling genes
filtered_HI_data <- HI_data[HI_data$log2FoldChange <= -1 & HI_data$padj<0.05,]

filtered_HI_data$sig <- (filtered_HI_data$SYMBOL %>% stringr::str_split(pattern = "_", simplify = TRUE))[,2]
table(filtered_HI_data$sig)
##      Sig
##   49  77

write.table(filtered_HI_data, "filtered_HI_49TSG_77sig_genes.tsv", sep="\t", row.names=FALSE, quote=FALSE)

write.table(filtered_HI_data %>% filter(sig!="Sig"), "filtered_HI_49TSG.tsv", sep="\t", row.names=FALSE, quote=FALSE)





# Figureで使ったTSG＆signalingに絞る
genelosslist <- c("ZBTB16", "HSPA8", "CADM1", "DRD2", "CDON", "PPP2R1B", "TIRAP", "GRIA4", "ATM", "GUCY1A2", "OPCML", "MMP3", "PDGFD", "BIRC2",
    "BIRC3", "TRPC6", "YAP1", "GRM5", "DLG2", "FZD4", "CAB39L", "NISCH", "BAP1", "CACNA1D", "GNAI2", "HTR2A", "MAGI1", "NPRL2", "PLXNB1", "PRICKLE2",
    "PTPRG", "RASSF1", "RB1", "RHOA", "ARRB1", "CACNA2D2", "CACNA2D3", "CISH", "DGKH", "FOXO1", "GRM2", "IP6K1", "IP6K2", "MAPKAPK3", "PBRM1", "PFKFB4",
    "RBM5", "SEMA3F", "TNFSF11", "TNNC1", "TUSC2", "WNT11", "WNT5A", "ZMYND10", "FAM107A", "LRIG1", "PDHB", "TRAF3", "CYLD")
genelosslist <- data.frame(SYMBOL_gene=genelosslist)

data <- resShrink_tbl
df <- read.table("/mnt/d/Visium/Seurat_HNSCC/pseudobulk_DESeq2_epithelial_infercnvScore02_cd45n_countFilterYuru/HI.genes", header=FALSE)
df02 <- read.table("/mnt/d/Visium/Seurat_HNSCC/pseudobulk_DESeq2_epithelial_infercnvScore02_cd45n_countFilterYuru/PCGR.TSG_HI.tsv", header=TRUE)
df02$SYMBOL_gene <- str_split(df02$SYMBOL, '_', simplify=TRUE)[,1]
df02 <- df02 %>% select(SYMBOL_gene, SYMBOL)

df02 <- dplyr::inner_join(df02, genelosslist, by="SYMBOL_gene")

df <- dplyr::inner_join(df02, df, by=c("SYMBOL_gene"="V1"))

library(tidyverse)
HI_data <- dplyr::left_join(df, data, by=c("SYMBOL_gene"="gene"))

# remove NA
HI_data <- HI_data[complete.cases(HI_data),]

# filtering 38 TSGsと87 signaling genes
filtered_HI_data <- HI_data[HI_data$log2FoldChange <= -1 & HI_data$padj<0.05,]
##     SYMBOL_gene    SYMBOL   baseMean log2FoldChange     lfcSE       pvalue         padj
##  2          ATM       ATM   8.563333      -1.180018 0.2451245 9.109667e-08 5.621314e-07
##  4        BIRC2 BIRC2_Sig  24.735189      -1.200553 0.1786789 1.047794e-12 1.474344e-11
##  14        DGKH  DGKH_Sig  12.855354      -1.414288 0.2621038 1.117933e-09 1.001117e-08
##  17     FAM107A   FAM107A   8.815055      -2.863104 0.7324105 4.902416e-07 2.643989e-06
##  18       FOXO1     FOXO1  10.477581      -1.506077 0.2545268 4.198725e-10 3.970724e-09
##  20       GNAI2 GNAI2_Sig  79.792953      -2.617859 0.2025130 3.803803e-42 5.624787e-40
##  30       MAGI1 MAGI1_Sig  15.898874      -1.675770 0.2274546 1.126048e-14 2.081397e-13
##  44      RASSF1    RASSF1  16.324936      -1.343102 0.1617281 1.178523e-17 3.184360e-16
##  47        RHOA      RHOA 241.106785      -1.059659 0.1217206 3.622074e-19 1.155229e-17
##  48      SEMA3F    SEMA3F  19.462862      -1.654425 0.3423423 1.382806e-07 8.308454e-07
##  56       WNT5A     WNT5A  12.121594      -1.373452 0.3777459 3.293694e-05 1.208827e-04

### 11q TSG
# ATM, HIPred 0.75, pHaplo 0.88
# CADM1, HIPred 0.66, pHaplo 0.98
# ZBTB16, HIPred 0.83, pHaplo 0.98
### 11q sig
# ARRB1, HIPred 0.73, pHaplo 0.71
# BIRC2, HIPred 0.73, pHaplo 0.67
# FZD4, HIPred 0.70, pHaplo 0.99
# GUCY1A2, HIPred 0.75, pHaplo 0.77
# PDGFD, HIPred 0.60, pHaplot 0.51
# PPP2R1B, HIPred 0.62, pHaplo 0.82

### plot Shrinkのthresholdで色
filtered_HI_data <- HI_data[abs(HI_data$log2FoldChange) > 1 & HI_data$padj<0.05,]

resShrink_table_thres <- resShrink_tbl[complete.cases(resShrink_tbl),] %>% 
                  mutate(threshold = case_when(
                    padj < 0.05 & log2FoldChange > 1 ~ "Upregulated",
                    padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
                    TRUE ~ "Not significant"
                    )) %>% select(gene, threshold)
res_table_thres <- res_tbl[complete.cases(res_tbl),]
res_table_thres <- dplyr::left_join(res_table_thres, resShrink_table_thres, by="gene")
label <- filtered_HI_data %>% select(SYMBOL_gene)
label$label <- label$SYMBOL_gene
res_table_thres <- dplyr::left_join(res_table_thres, label, by=c("gene"="SYMBOL_gene"))
res_table_thres$threshold <- factor(res_table_thres$threshold, levels=c("Upregulated", "Downregulated", "Not significant"))
  
## Volcano plot
p <- ggplot(res_table_thres) +
    geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
    xlab("log2 fold change") + 
    ylab("-log10 adjusted p-value") +
    theme(legend.position = "none",
          plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.title = element_text(size = rel(1.25))) +
    ggrepel::geom_text_repel(aes(x = log2FoldChange, y = -log10(padj), label=label),
                        size=3,
                        max.overlaps = 50,
                        force = 80, max.iter = 3000, min.segment.length = 0.1
                        ) +
    theme_bw() +
    scale_color_manual(values=c("#B93F2B", "#0A74B2", "#AFA3A0"))
ggsave(file=paste0(cluster_name, "_", levels(cluster_metadata$group_id)[2], "_vs_", levels(cluster_metadata$group_id)[1], "_Volcano.Shrink_thres.labeled.40perc.abs.pdf"), plot=p)








