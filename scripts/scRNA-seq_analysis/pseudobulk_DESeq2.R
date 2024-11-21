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
##  List of length 1
##  names(1): counts

## Explore the raw counts for the dataset
dim(counts(sce))
##  [1] 25523 10965

counts(sce)[1:6, 1:6]
##      6 x 6 sparse Matrix of class "dgCMatrix"
##                  HN12_CD45p_Tumor_AAAGTAGGTACAAGTA-1
##      MIR1302-2HG                                   .
##      AL627309.1                                    .
##      AL627309.3                                    .
##      AL627309.4                                    .
##      AL732372.1                                    .
##      AC114498.1                                    .
##                  HN12_CD45p_Tumor_ACGAGCCGTCTCTCTG-1
##      MIR1302-2HG                                   .
##      AL627309.1                                    .
##      AL627309.3                                    .
##      AL627309.4                                    .
##      AL732372.1                                    .
##      AC114498.1                                    .
##                  HN12_CD45p_Tumor_ACTGAGTTCTTGTATC-1
##      MIR1302-2HG                                   .
##      AL627309.1                                    .
##      AL627309.3                                    .
##      AL627309.4                                    .
##      AL732372.1                                    .
##      AC114498.1                                    .
##                  HN12_CD45p_Tumor_CCACGGAAGGCAGGTT-1
##      MIR1302-2HG                                   .
##      AL627309.1                                    .
##      AL627309.3                                    .
##      AL627309.4                                    .
##      AL732372.1                                    .
##      AC114498.1                                    .
##                  HN12_CD45p_Tumor_CTCGGAGAGGGATCTG-1
##      MIR1302-2HG                                   .
##      AL627309.1                                    .
##      AL627309.3                                    .
##      AL627309.4                                    .
##      AL732372.1                                    .
##      AC114498.1                                    .
##                  HN12_CD45p_Tumor_GACAGAGTCGGAGGTA-1
##      MIR1302-2HG                                   .
##      AL627309.1                                    .
##      AL627309.3                                    .
##      AL627309.4                                    .
##      AL732372.1                                    .
##      AC114498.1                                    .




# Perform QC if not already performed
dim(sce)
##  [1] 25523 10965
table(sce$group_id)
##  Normal  Tumor
##    2475   8490

#   # Calculate quality control (QC) metrics
#       #sce <- calculateQCMetrics(sce)
#       ##  Error: 'calculateQCMetrics' is defunct.
#       ##  Use 'perCellQCMetrics' instead.
#       ##  See help("Defunct")
#   per.cell <- perCellQCMetrics(sce)
#   per.cell
#   ##  DataFrame with 5966 rows and 3 columns
#   ##                                                                 sum  detected
#   ##                                                           <numeric> <integer>
#   ##  HN12_CD45n_Tumor_AAACCTGAGAACAATC-1                           7763      2609
#   ##  HN12_CD45n_Tumor_AAACCTGAGGAATGGA-1                           2830      1262
#   ##  HN12_CD45n_Tumor_AAACCTGTCGTTTGCC-1                           2359      1111
#   ##  HN12_CD45n_Tumor_AAACCTGTCTTTACAC-1                           1652       732
#   ##  HN12_CD45n_Tumor_AAACGGGAGACCTAGG-1                           2839      1390
#   ##  ...                                                            ...       ...
#   ##  BCLLATLAS_57_nqyw13tk_o7guqugw_Normal_CCGATGGCAGGACTAG-1      4785      1475
#   ##  BCLLATLAS_57_nqyw13tk_o7guqugw_Normal_GGTGGCTAGTACGAGC-1      2549       986
#   ##  BCLLATLAS_57_p4b145yq_4z8crwqq_Normal_AGGATAACAAGTATAG-1      1476       266
#   ##  BCLLATLAS_57_p4b145yq_4z8crwqq_Normal_CCAATGAGTATGATCC-1      4717       621
#   ##  BCLLATLAS_57_p4b145yq_4z8crwqq_Normal_TTGGGATCATGTGGCC-1      4565      1611
#   ##                                                               total
#   ##                                                           <numeric>
#   ##  HN12_CD45n_Tumor_AAACCTGAGAACAATC-1                           7763
#   ##  HN12_CD45n_Tumor_AAACCTGAGGAATGGA-1                           2830
#   ##  HN12_CD45n_Tumor_AAACCTGTCGTTTGCC-1                           2359
#   ##  HN12_CD45n_Tumor_AAACCTGTCTTTACAC-1                           1652
#   ##  HN12_CD45n_Tumor_AAACGGGAGACCTAGG-1                           2839
#   ##  ...                                                            ...
#   ##  BCLLATLAS_57_nqyw13tk_o7guqugw_Normal_CCGATGGCAGGACTAG-1      4785
#   ##  BCLLATLAS_57_nqyw13tk_o7guqugw_Normal_GGTGGCTAGTACGAGC-1      2549
#   ##  BCLLATLAS_57_p4b145yq_4z8crwqq_Normal_AGGATAACAAGTATAG-1      1476
#   ##  BCLLATLAS_57_p4b145yq_4z8crwqq_Normal_CCAATGAGTATGATCC-1      4717
#   ##  BCLLATLAS_57_p4b145yq_4z8crwqq_Normal_TTGGGATCATGTGGCC-1      4565
#   
#   # Get cells w/ few/many detected genes
#   sce$is_outlier <- isOutlier(
#           metric = per.cell$total,
#           nmads = 2, type = "both", log = TRUE)
#   
#   # Remove outlier cells
#   sce <- sce[, !sce$is_outlier]
#   dim(sce)
#   ##  [1] 25523 10635
#   table(sce$group_id)
#   ##  Normal  Tumor
#   ##    2248   8389
#   
#   ## Remove lowly expressed genes which have less than 10 cells with any counts
#   sce <- sce[rowSums(counts(sce) > 1) >= 10, ]
#   
#   dim(sce)
#   ##  [1] 12842 10635


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
##      DataFrame with 6 rows and 6 columns
##                                             orig.ident nCount_RNA nFeature_RNA
##                                            <character>  <numeric>    <integer>
##      HN12_CD45p_Tumor_AAAGTAGGTACAAGTA-1 SeuratProject       7099         1726
##      HN12_CD45p_Tumor_ACGAGCCGTCTCTCTG-1 SeuratProject      12567         3371
##      HN12_CD45p_Tumor_ACTGAGTTCTTGTATC-1 SeuratProject       1585          685
##      HN12_CD45p_Tumor_CCACGGAAGGCAGGTT-1 SeuratProject       7102         1534
##      HN12_CD45p_Tumor_CTCGGAGAGGGATCTG-1 SeuratProject       1162          516
##      HN12_CD45p_Tumor_GACAGAGTCGGAGGTA-1 SeuratProject        588          322
##                                               Sample            Barcode SampleGroup
##                                          <character>        <character> <character>
##      HN12_CD45p_Tumor_AAAGTAGGTACAAGTA-1  HN12_CD45p AAAGTAGGTACAAGTA-1       Tumor
##      HN12_CD45p_Tumor_ACGAGCCGTCTCTCTG-1  HN12_CD45p ACGAGCCGTCTCTCTG-1       Tumor
##      HN12_CD45p_Tumor_ACTGAGTTCTTGTATC-1  HN12_CD45p ACTGAGTTCTTGTATC-1       Tumor
##      HN12_CD45p_Tumor_CCACGGAAGGCAGGTT-1  HN12_CD45p CCACGGAAGGCAGGTT-1       Tumor
##      HN12_CD45p_Tumor_CTCGGAGAGGGATCTG-1  HN12_CD45p CTCGGAGAGGGATCTG-1       Tumor
##      HN12_CD45p_Tumor_GACAGAGTCGGAGGTA-1  HN12_CD45p GACAGAGTCGGAGGTA-1       Tumor





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
##   BCLLATLAS_05_izi9unx1_8qdzhivu  BCLLATLAS_06_d8kwy76j_7blj9otf
##                              395                             883
##   BCLLATLAS_06_giz3qso4_783dbpu6  BCLLATLAS_07_i5udk3x0_57gv6ncx
##                             1024                             554
##   BCLLATLAS_07_umt51kfr_p8ei65ms BCLLATLAS_131_rfs8oamh_0obpdn7k
##                              608                             129
##  BCLLATLAS_131_xxfne43y_x9dh95oz  BCLLATLAS_16_bw94nf57_vm85woki
##                              110                             998
##   BCLLATLAS_16_ggq3nifm_jkilwp1x  BCLLATLAS_19_dvcbn9p8_ix0j3k8b
##                             1245                             758
##   BCLLATLAS_19_ff8s19u3_7e96iusr  BCLLATLAS_21_dvdzq8et_eot75su8
##                              806                             551
##   BCLLATLAS_21_md651vbh_eymr91s7  BCLLATLAS_21_n1b3su0a_l7shyi35
##                             1602                            1040
##   BCLLATLAS_21_x739d5z1_dsamhgey  BCLLATLAS_25_bz5rpwtv_kg7w108r
##                              479                             985
##   BCLLATLAS_25_wf4su8ny_h4yj8bv7  BCLLATLAS_34_kjzv2rwx_sfomyxok
##                             1043                            4578
##   BCLLATLAS_34_v8g80gtx_ps9bamz7  BCLLATLAS_41_ejto2bae_y5mydeam
##                             4506                            2419
##   BCLLATLAS_41_z3of7uaq_mzbhy4tt  BCLLATLAS_54_altbaco5_45sf3wul
##                             2360                            1520
##   BCLLATLAS_54_c3ftguo2_xxwxs507  BCLLATLAS_57_nqyw13tk_o7guqugw
##                             1578                            3091
##   BCLLATLAS_57_p4b145yq_4z8crwqq                      HN12_CD45p
##                             3191                             242
##                       HN14_CD45p                      HN16_CD45p
##                               27                             108
##                       HN17_CD45p                      HN18_CD45p
##                              543                               1

## Turn named vector into a numeric vector of number of cells per sample
n_cells <- as.numeric(table(sce$sample_id))

## Determine how to reoder the samples (rows) of the metadata to match the order of sample names in sids vector
m <- match(sids, sce$sample_id)

## Create the sample level metadata by combining the reordered metadata with the number of cells corresponding to each sample.
ei <- data.frame(colData(sce)[m, ], 
                  n_cells, row.names = NULL) %>% 
                select(-"cluster_id")
ei[1:6,1:6]
##           orig.ident nCount_RNA nFeature_RNA                          Sample
##      1 SeuratProject        938          433  BCLLATLAS_05_izi9unx1_8qdzhivu
##      2 SeuratProject        832          500  BCLLATLAS_06_d8kwy76j_7blj9otf
##      3 SeuratProject       3069         1253  BCLLATLAS_06_giz3qso4_783dbpu6
##      4 SeuratProject        951          504  BCLLATLAS_07_i5udk3x0_57gv6ncx
##      5 SeuratProject        501          316  BCLLATLAS_07_umt51kfr_p8ei65ms
##      6 SeuratProject       3395         1319 BCLLATLAS_131_rfs8oamh_0obpdn7k
##                   Barcode SampleGroup
##      1 AAATGGAGTTCGAAGG-1      Normal
##      2 CTCAGGGAGGACGGAG-1      Normal
##      3 TCAGCAACATAAGCGG-1      Normal
##      4 AGGTTGTGTGCTCCGA-1      Normal
##      5 AAATGGAGTTCTCTCG-1      Normal
##      6 AGATCGTTCAGCAGAG-1      Normal



### QC https://bioconductor.org/packages/devel/bioc/vignettes/scuttle/inst/doc/qc.html
# SKIP


### Count aggregation to sample level
# Aggregate the counts per sample_id and cluster_id

# Subset metadata to only include the cluster and sample IDs to aggregate across
groups <- colData(sce)[, c("cluster_id", "sample_id")]
groups
##  DataFrame with 5840 rows and 2 columns
##                                                           cluster_id
##                                                             <factor>
##  HN12_CD45p_Tumor_AAAGTAGGTACAAGTA-1                            TRUE
##  HN12_CD45p_Tumor_ACGAGCCGTCTCTCTG-1                            TRUE
##  HN12_CD45p_Tumor_ACTGAGTTCTTGTATC-1                            TRUE
##  HN12_CD45p_Tumor_CACACTCGTAGGGTAC-1                            TRUE
##  HN12_CD45p_Tumor_CCACGGAAGGCAGGTT-1                            TRUE
##  ...                                                             ...
##  BCLLATLAS_57_nqyw13tk_o7guqugw_Normal_CCGATGGCAGGACTAG-1       TRUE
##  BCLLATLAS_57_nqyw13tk_o7guqugw_Normal_GGTGGCTAGTACGAGC-1       TRUE
##  BCLLATLAS_57_p4b145yq_4z8crwqq_Normal_AGGATAACAAGTATAG-1       TRUE
##  BCLLATLAS_57_p4b145yq_4z8crwqq_Normal_CCAATGAGTATGATCC-1       TRUE
##  BCLLATLAS_57_p4b145yq_4z8crwqq_Normal_TTGGGATCATGTGGCC-1       TRUE
##                                                                                sample_id
##                                                                                 <factor>
##  HN12_CD45p_Tumor_AAAGTAGGTACAAGTA-1                                          HN12_CD45p
##  HN12_CD45p_Tumor_ACGAGCCGTCTCTCTG-1                                          HN12_CD45p
##  HN12_CD45p_Tumor_ACTGAGTTCTTGTATC-1                                          HN12_CD45p
##  HN12_CD45p_Tumor_CACACTCGTAGGGTAC-1                                          HN12_CD45p
##  HN12_CD45p_Tumor_CCACGGAAGGCAGGTT-1                                          HN12_CD45p
##  ...                                                                                 ...
##  BCLLATLAS_57_nqyw13tk_o7guqugw_Normal_CCGATGGCAGGACTAG-1 BCLLATLAS_57_nqyw13tk_o7guqugw
##  BCLLATLAS_57_nqyw13tk_o7guqugw_Normal_GGTGGCTAGTACGAGC-1 BCLLATLAS_57_nqyw13tk_o7guqugw
##  BCLLATLAS_57_p4b145yq_4z8crwqq_Normal_AGGATAACAAGTATAG-1 BCLLATLAS_57_p4b145yq_4z8crwqq
##  BCLLATLAS_57_p4b145yq_4z8crwqq_Normal_CCAATGAGTATGATCC-1 BCLLATLAS_57_p4b145yq_4z8crwqq
##  BCLLATLAS_57_p4b145yq_4z8crwqq_Normal_TTGGGATCATGTGGCC-1 BCLLATLAS_57_p4b145yq_4z8crwqq

# Aggregate across cluster-sample groups
pb <- aggregate.Matrix(t(counts(sce)), 
                       groupings = groups, fun = "sum") 

class(pb)
##  [1] "dgCMatrix"
##  attr(,"package")
##  [1] "Matrix"

dim(pb)
##  [1]    31 11050

pb[1:6, 1:6]
##      6 x 6 sparse Matrix of class "dgCMatrix"
##                                                 MIR1302-2HG AL627309.1 AL627309.3
##      Epithelial_BCLLATLAS_05_izi9unx1_8qdzhivu            .          .          .
##      Epithelial_BCLLATLAS_06_d8kwy76j_7blj9otf            .          .          .
##      Epithelial_BCLLATLAS_06_giz3qso4_783dbpu6            .          .          .
##      Epithelial_BCLLATLAS_07_i5udk3x0_57gv6ncx            .          .          .
##      Epithelial_BCLLATLAS_07_umt51kfr_p8ei65ms            .          .          .
##      Epithelial_BCLLATLAS_131_rfs8oamh_0obpdn7k           .          .          .
##                                                 AL627309.4 AL732372.1 AC114498.1
##      Epithelial_BCLLATLAS_05_izi9unx1_8qdzhivu           .          .          .
##      Epithelial_BCLLATLAS_06_d8kwy76j_7blj9otf           .          .          .
##      Epithelial_BCLLATLAS_06_giz3qso4_783dbpu6           .          .          .
##      Epithelial_BCLLATLAS_07_i5udk3x0_57gv6ncx           .          .          .
##      Epithelial_BCLLATLAS_07_umt51kfr_p8ei65ms           .          .          .
##      Epithelial_BCLLATLAS_131_rfs8oamh_0obpdn7k          .          .          .


### To perform DE analysis on a per cell type basis, we need to wrangle our data in a couple ways. We need to do the following steps:
## 1. Split our data by cell type
## 2. Transform the matrix so that the genes are the row names and the samples are the column names

# Not every cluster is present in all samples; create a vector that represents how to split samples
splitf <- sapply(stringr::str_split(rownames(pb), 
                                    pattern = "_",
                                    n = 2),
                `[`, 1)
head(splitf)
##  [1] "Epithelial" "Epithelial" "Epithelial" "Epithelial" "Epithelial"
##  [6] "Epithelial"

# Turn into a list and split the list into components for each cluster and transform, so rows are genes and columns are samples and make rownames as the sample IDs
#pb <- split.data.frame(pb, 
#                       factor(splitf)) %>%
#        lapply(function(u) 
#                set_colnames(t(u), 
#                             stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))
pb <- split.data.frame(pb, 
                       factor(splitf)) %>%
        lapply(function(u) 
                set_colnames(t(u), 
                             sapply(stringr::str_split(rownames(u), pattern = "_", n = 2), `[`, 2)))

class(pb)
##  [1] "list"

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
##   BCLLATLAS_05_izi9unx1_8qdzhivu  BCLLATLAS_06_d8kwy76j_7blj9otf  BCLLATLAS_06_giz3qso4_783dbpu6
##                               77                               3                               1
##   BCLLATLAS_07_i5udk3x0_57gv6ncx  BCLLATLAS_07_umt51kfr_p8ei65ms BCLLATLAS_131_rfs8oamh_0obpdn7k
##                                7                               9                               5








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
##  [1] "SG_SL_KR_00117" "SG_SL_OI_00010" "SG_SL_OI_00046" "SG_SL_OI_00053" "SG_SL_UO_00159"
##  [6] "SG_SL_KR_00010"
# はじめのASDCは5例でしかカウントデータなし


# Get cluster IDs for each of the samples

samples_list <- map(1:length(cids), get_sample_ids)

get_cluster_ids <- function(x){
        rep(names(pb)[x], 
            each = length(samples_list[[x]]))
}

de_cluster_ids <- map(1:length(cids), get_cluster_ids) %>%
        unlist()
head(de_cluster_ids)
##  [1] "ASDC"           "ASDC"           "ASDC"           "ASDC"           "ASDC"
##  [6] "B intermediate"


# Create a data frame with the sample IDs, cluster IDs and condition
#factorにしておかないと下のlevels()がNULLを返してしまう
gg_df <- data.frame(cluster_id = factor(de_cluster_ids),
                    sample_id = factor(de_samples))

gg_df <- left_join(gg_df, ei[, c("sample_id", "group_id", "cov_sample")]) 


metadata <- gg_df %>%
        dplyr::select(cluster_id, sample_id, group_id, cov_sample)

head(metadata)     
##    cluster_id                       sample_id group_id
##  1       TRUE  BCLLATLAS_05_izi9unx1_8qdzhivu   Normal
##  2       TRUE  BCLLATLAS_06_d8kwy76j_7blj9otf   Normal
##  3       TRUE  BCLLATLAS_06_giz3qso4_783dbpu6   Normal
##  4       TRUE  BCLLATLAS_07_i5udk3x0_57gv6ncx   Normal
##  5       TRUE  BCLLATLAS_07_umt51kfr_p8ei65ms   Normal
##  6       TRUE BCLLATLAS_131_rfs8oamh_0obpdn7k   Normal


#       metadata02 <- read.table("metadata02.tsv", sep="\t", header=T)
#       
#       
#       
#       metadata <- left_join(metadata, metadata02)
#       metadata$sex <- factor(metadata$sex, levels=c("M", "F"))
#       metadata$scaled_age <- scale(metadata$age)


write.csv(metadata, "metadata_for_DESeq2.csv", quote=FALSE, row.names=FALSE)


saveRDS(pb, "pb.rds")
saveRDS(metadata, "metadata.rds")
# pb <- readRDS("pb.rds")
# metadata <- readRDS("metadata.rds")





### Subsetting dataset to cluster(s) of interest
# Generate vector of cluster IDs
clusters <- levels(metadata$cluster_id)
clusters
##      [1] "Epithelial"


#  Let’s perform the DE analysis on "B naive"
clusters[1]
##      [1] "Epithelial"

# Subset the metadata to only the B cells                                       ### SG_SL_KR_00011ではB naiveなし
cluster_metadata <- metadata[which(metadata$cluster_id == clusters[1]), ]
head(cluster_metadata)
##    cluster_id                       sample_id group_id
##  1       TRUE  BCLLATLAS_05_izi9unx1_8qdzhivu   Normal
##  2       TRUE  BCLLATLAS_06_d8kwy76j_7blj9otf   Normal
##  3       TRUE  BCLLATLAS_06_giz3qso4_783dbpu6   Normal
##  4       TRUE  BCLLATLAS_07_i5udk3x0_57gv6ncx   Normal
##  5       TRUE  BCLLATLAS_07_umt51kfr_p8ei65ms   Normal
##  6       TRUE BCLLATLAS_131_rfs8oamh_0obpdn7k   Normal
nrow(cluster_metadata)
##  [1] 30


# Assign the rownames of the metadata to be the sample IDs
rownames(cluster_metadata) <- cluster_metadata$sample_id
head(cluster_metadata)
##                                  cluster_id                       sample_id group_id
##  BCLLATLAS_05_izi9unx1_8qdzhivu        TRUE  BCLLATLAS_05_izi9unx1_8qdzhivu   Normal
##  BCLLATLAS_06_d8kwy76j_7blj9otf        TRUE  BCLLATLAS_06_d8kwy76j_7blj9otf   Normal
##  BCLLATLAS_06_giz3qso4_783dbpu6        TRUE  BCLLATLAS_06_giz3qso4_783dbpu6   Normal
##  BCLLATLAS_07_i5udk3x0_57gv6ncx        TRUE  BCLLATLAS_07_i5udk3x0_57gv6ncx   Normal
##  BCLLATLAS_07_umt51kfr_p8ei65ms        TRUE  BCLLATLAS_07_umt51kfr_p8ei65ms   Normal
##  BCLLATLAS_131_rfs8oamh_0obpdn7k       TRUE BCLLATLAS_131_rfs8oamh_0obpdn7k   Normal

# Subset the counts to only the B cells
counts <- pb[[clusters[1]]]
head(counts)
##      6 x 34 sparse Matrix of class "dgCMatrix"
##        [[ suppressing 34 column names ‘BCLLATLAS_05_izi9unx1_8qdzhivu’, ‘BCLLATLAS_06_d8kwy76j_7blj9otf’, ‘BCLLATLAS_06_giz3qso4_783dbpu6’ ... ]]
##      
##      MIR1302-2HG . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
##      AL627309.1  . . . . . . . . . . . . . . . 1 . . . . . . . . . 8 . 2 . . 1 1 1 .
##      AL627309.3  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
##      AL627309.4  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
##      AL732372.1  . . . . . . . . . . . . . . . 1 . . . . . . . . . . . . . . 1 1 . .
##      AC114498.1  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 1 1 . .

cluster_counts <- data.frame(counts[, which(colnames(counts) %in% rownames(cluster_metadata))])
head(cluster_counts)
##                  BCLLATLAS_05_izi9unx1_8qdzhivu BCLLATLAS_06_d8kwy76j_7blj9otf
##      MIR1302-2HG                              0                              0
##      AL627309.1                               0                              0
##      AL627309.3                               0                              0
##      AL627309.4                               0                              0
##      AL732372.1                               0                              0
##      AC114498.1                               0                              0
##                  BCLLATLAS_06_giz3qso4_783dbpu6 BCLLATLAS_07_i5udk3x0_57gv6ncx
##      MIR1302-2HG                              0                              0
##      AL627309.1                               0                              0
##      AL627309.3                               0                              0
##      AL627309.4                               0                              0
##      AL732372.1                               0                              0
##      AC114498.1                               0                              0
##                  BCLLATLAS_07_umt51kfr_p8ei65ms BCLLATLAS_131_rfs8oamh_0obpdn7k
##      MIR1302-2HG                              0                               0
##      AL627309.1                               0                               0
##      AL627309.3                               0                               0
##      AL627309.4                               0                               0
##      AL732372.1                               0                               0
##      AC114498.1                               0                               0
##                  BCLLATLAS_131_xxfne43y_x9dh95oz BCLLATLAS_16_bw94nf57_vm85woki
##      MIR1302-2HG                               0                              0
##      AL627309.1                                0                              0
##      AL627309.3                                0                              0
##      AL627309.4                                0                              0
##      AL732372.1                                0                              0
##      AC114498.1                                0                              0
##                  BCLLATLAS_16_ggq3nifm_jkilwp1x BCLLATLAS_19_dvcbn9p8_ix0j3k8b
##      MIR1302-2HG                              0                              0
##      AL627309.1                               0                              0
##      AL627309.3                               0                              0
##      AL627309.4                               0                              0
##      AL732372.1                               0                              0
##      AC114498.1                               0                              0
##                  BCLLATLAS_19_ff8s19u3_7e96iusr BCLLATLAS_21_dvdzq8et_eot75su8
##      MIR1302-2HG                              0                              0
##      AL627309.1                               0                              0
##      AL627309.3                               0                              0
##      AL627309.4                               0                              0
##      AL732372.1                               0                              0
##      AC114498.1                               0                              0
##                  BCLLATLAS_21_md651vbh_eymr91s7 BCLLATLAS_21_n1b3su0a_l7shyi35
##      MIR1302-2HG                              0                              0
##      AL627309.1                               0                              0
##      AL627309.3                               0                              0
##      AL627309.4                               0                              0
##      AL732372.1                               0                              0
##      AC114498.1                               0                              0
##                  BCLLATLAS_21_x739d5z1_dsamhgey BCLLATLAS_25_bz5rpwtv_kg7w108r
##      MIR1302-2HG                              0                              0
##      AL627309.1                               0                              1
##      AL627309.3                               0                              0
##      AL627309.4                               0                              0
##      AL732372.1                               0                              1
##      AC114498.1                               0                              0
##                  BCLLATLAS_25_wf4su8ny_h4yj8bv7 BCLLATLAS_34_kjzv2rwx_sfomyxok
##      MIR1302-2HG                              0                              0
##      AL627309.1                               0                              0
##      AL627309.3                               0                              0
##      AL627309.4                               0                              0
##      AL732372.1                               0                              0
##      AC114498.1                               0                              0
##                  BCLLATLAS_34_v8g80gtx_ps9bamz7 BCLLATLAS_41_ejto2bae_y5mydeam
##      MIR1302-2HG                              0                              0
##      AL627309.1                               0                              0
##      AL627309.3                               0                              0
##      AL627309.4                               0                              0
##      AL732372.1                               0                              0
##      AC114498.1                               0                              0
##                  BCLLATLAS_41_z3of7uaq_mzbhy4tt BCLLATLAS_54_altbaco5_45sf3wul
##      MIR1302-2HG                              0                              0
##      AL627309.1                               0                              0
##      AL627309.3                               0                              0
##      AL627309.4                               0                              0
##      AL732372.1                               0                              0
##      AC114498.1                               0                              0
##                  BCLLATLAS_54_c3ftguo2_xxwxs507 BCLLATLAS_57_nqyw13tk_o7guqugw
##      MIR1302-2HG                              0                              0
##      AL627309.1                               0                              0
##      AL627309.3                               0                              0
##      AL627309.4                               0                              0
##      AL732372.1                               0                              0
##      AC114498.1                               0                              0
##                  BCLLATLAS_57_p4b145yq_4z8crwqq HN12_CD45n HN12_CD45p HN13_CD45n HN14_CD45n HN16_CD45n
##      MIR1302-2HG                              0          0          0          0          0          0
##      AL627309.1                               0          8          0          2          0          0
##      AL627309.3                               0          0          0          0          0          0
##      AL627309.4                               0          0          0          0          0          0
##      AL732372.1                               0          0          0          0          0          0
##      AC114498.1                               0          0          0          0          0          0
##                  HN17_CD45n HN17_CD45p HN18_CD45n HN18_CD45p
##      MIR1302-2HG          0          0          0          0
##      AL627309.1           1          1          1          0
##      AL627309.3           0          0          0          0
##      AL627309.4           0          0          0          0
##      AL732372.1           1          1          0          0
##      AC114498.1           1          1          0          0

# Check that all of the row names of the metadata are the same and in the same order as the column names of the counts in order to use as input to DESeq2
all(rownames(cluster_metadata) == colnames(cluster_counts))        
##  [1] TRUE


#------------------------#
### filtering by count ###
#------------------------#
summary(rowSums(cluster_counts))
##      Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
##        63      994     2328    13489     5375 10630954

### カウント合計が10以上の遺伝子だけを残す
cutoff <- 10
table(rowSums(cluster_counts) >= cutoff)
##  FALSE  TRUE
##   5875 19648

#filtered_counts <- cluster_counts[rowSums(cluster_counts) >= cutoff, ]
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

#   > keep["ATM"]
#    ATM
#   TRUE
#   > keep["CADM1"]
#   CADM1
#   FALSE
#   > keep["ZBTB16"]
#   ZBTB16
#    FALSE
#   > keep["ARRB1"]
#   ARRB1
#   FALSE
#   > keep["BIRC2"]
#   BIRC2
#    TRUE

saveRDS(keep, "keep.rds")


### ATM
# Normal
counts(dds, normalized=TRUE) %>% data.frame() %>% subset(rownames(data)=="ATM") %>% select(cluster_metadata[cluster_metadata$group_id=="Normal",]$sample_id)
##      OP13_s2  OP14_s1  OP14_s2  OP17_s1  OP20_s4  OP33_s4 OP33norm_s2 OP33norm_s4  OP34_s3
##  ATM       0 19.63438 7.798658 10.76861 11.25542 11.34045    32.16547           0 13.32537
##      OP34norm_s1 OP35BOT_s4 OP35norm   OP4_s2  OP5_s5   OP9_s2   OP9_s3   OP9_s5
##  ATM     21.7032          0 22.09937 14.91688 8.70393 21.48134 18.02215 20.64444
counts(dds, normalized=TRUE) %>% data.frame() %>% subset(rownames(data)=="ATM") %>% select(cluster_metadata[cluster_metadata$group_id=="Normal",]$sample_id) %>% rowMeans()
##       ATM
##  13.75645

# Tumor
counts(dds, normalized=TRUE) %>% data.frame() %>% subset(rownames(data)=="ATM") %>% select(cluster_metadata[cluster_metadata$group_id=="Tumor",]$sample_id)
##      OP14_s3  OP14_s4  OP17_s2 OP17_s3  OP20_s1  OP20_s2  OP20_s3  OP20_s5  OP20_s6  OP20_s7
##  ATM       0 1.926994 4.071093 4.45987 13.01274 8.382454 8.827991 4.939107 6.037842 2.231336
##       OP33_s1  OP33_s2  OP33_s3 OP33_s5  OP33_s6  OP34_s1  OP34_s2  OP34_s4 OP35BOT_s1 OP35BOT_s2
##  ATM 16.91709 15.39012 8.381166 7.43386 7.339431 5.675074 5.103013 3.881482          0          0
##      OP35BOT_s3 OP35LN_s1 OP35LN_s2 OP35LN_s3   OP4_s1   OP4_s3  OP4_s4  OP4_s5   OP4_s6   OP4_s7
##  ATM   4.796842  10.19156         0  22.62211 5.057271 6.848189 8.11186 5.50279 4.201155 3.556312
##        OP4_s8   OP4_s9   OP5_s1   OP5_s2   OP5_s3   OP5_s4   OP6_s1   OP6_s2   OP6_s3   OP9_s4
##  ATM 7.049878 3.833412 11.19881 7.790135 6.049358 3.821937 3.541533 6.998752 2.260343 6.807374
counts(dds, normalized=TRUE) %>% data.frame() %>% subset(rownames(data)=="ATM") %>% select(cluster_metadata[cluster_metadata$group_id=="Tumor",]$sample_id) %>% rowMeans()
##       ATM
##  6.356257

### ZBTB16
# Normal
counts(dds, normalized=TRUE) %>% data.frame() %>% subset(rownames(data)=="ZBTB16") %>% select(cluster_metadata[cluster_metadata$group_id=="Normal",]$sample_id)
##         OP13_s2  OP14_s1 OP14_s2  OP17_s1  OP20_s4  OP33_s4 OP33norm_s2 OP33norm_s4  OP34_s3
##  ZBTB16       0 1.682947       0 1.076861 5.627708 3.240128           0    36.07137 1.615197
##         OP34norm_s1 OP35BOT_s4 OP35norm   OP4_s2   OP5_s5   OP9_s2 OP9_s3 OP9_s5
##  ZBTB16    24.80366          0 29.46582 10.26746 17.40786 8.162909      0      0
counts(dds, normalized=TRUE) %>% data.frame() %>% subset(rownames(data)=="ZBTB16") %>% select(cluster_metadata[cluster_metadata$group_id=="Tumor",]$sample_id) %>% rowMeans()
##   ZBTB16
##  8.20129

# Tumor
counts(dds, normalized=TRUE) %>% data.frame() %>% subset(rownames(data)=="ZBTB16") %>% select(cluster_metadata[cluster_metadata$group_id=="Tumor",]$sample_id)
##         OP14_s3 OP14_s4 OP17_s2 OP17_s3  OP20_s1  OP20_s2  OP20_s3 OP20_s5 OP20_s6 OP20_s7   OP33_s1
##  ZBTB16       0       0       0       0 3.203136 1.635601 3.433108       0       0       0 0.4126119
##           OP33_s2   OP33_s3 OP33_s5 OP33_s6 OP34_s1 OP34_s2 OP34_s4 OP35BOT_s1 OP35BOT_s2 OP35BOT_s3
##  ZBTB16 0.1399102 0.1611763       0       0       0       0       0          0          0          0
##         OP35LN_s1 OP35LN_s2 OP35LN_s3    OP4_s1 OP4_s3 OP4_s4 OP4_s5 OP4_s6 OP4_s7 OP4_s8    OP4_s9
##  ZBTB16         0         0         0 0.2593473      0      0      0      0      0      0 0.7666825
##         OP5_s1 OP5_s2 OP5_s3 OP5_s4    OP6_s1 OP6_s2 OP6_s3    OP9_s4
##  ZBTB16      0      0      0      0 0.5902555      0      0 0.7563749
counts(dds, normalized=TRUE) %>% data.frame() %>% subset(rownames(data)=="ZBTB16") %>% select(cluster_metadata[cluster_metadata$group_id=="Tumor",]$sample_id) %>% rowMeans()
##     ZBTB16
##  0.2839551

### CADM1
# Normal
counts(dds, normalized=TRUE) %>% data.frame() %>% subset(rownames(data)=="CADM1") %>% select(cluster_metadata[cluster_metadata$group_id=="Normal",]$sample_id)
##         OP13_s2 OP14_s1  OP14_s2  OP17_s1  OP20_s4  OP33_s4 OP33norm_s2 OP33norm_s4  OP34_s3
##  CADM1 3.261301 5.04884 31.19463 12.92233 11.25542 6.300249    16.08274           0 2.018996
##        OP34norm_s1 OP35BOT_s4 OP35norm   OP4_s2 OP5_s5    OP9_s2   OP9_s3   OP9_s5
##  CADM1    6.200914   18.32164        0 9.492561      0 0.4296268 13.05052 1.423755
counts(dds, normalized=TRUE) %>% data.frame() %>% subset(rownames(data)=="CADM1") %>% select(cluster_metadata[cluster_metadata$group_id=="Normal",]$sample_id) %>% rowMeans()
##    CADM1
##  8.05903

# Tumor
counts(dds, normalized=TRUE) %>% data.frame() %>% subset(rownames(data)=="CADM1") %>% select(cluster_metadata[cluster_metadata$group_id=="Tumor",]$sample_id)
##        OP14_s3 OP14_s4 OP17_s2 OP17_s3  OP20_s1  OP20_s2  OP20_s3  OP20_s5  OP20_s6  OP20_s7
##  CADM1       0       0       0       0 18.81843 14.72041 17.65598 5.556496 4.830273 1.115668
##          OP33_s1 OP33_s2   OP33_s3  OP33_s5 OP33_s6   OP34_s1  OP34_s2   OP34_s4 OP35BOT_s1
##  CADM1 0.1650447       0 0.4835288 2.262479       0 0.5870766 1.342898 0.7762965          0
##        OP35BOT_s2 OP35BOT_s3 OP35LN_s1 OP35LN_s2 OP35LN_s3    OP4_s1 OP4_s3 OP4_s4 OP4_s5   OP4_s6
##  CADM1          0   4.796842         0         0         0 0.7780418      0      0      0 1.260347
##        OP4_s7 OP4_s8   OP4_s9   OP5_s1 OP5_s2 OP5_s3 OP5_s4 OP6_s1 OP6_s2 OP6_s3   OP9_s4
##  CADM1      0      0 1.533365 2.036148      0      0      0      0      0      0 23.44762
counts(dds, normalized=TRUE) %>% data.frame() %>% subset(rownames(data)=="CADM1") %>% select(cluster_metadata[cluster_metadata$group_id=="Tumor",]$sample_id) %>% rowMeans()
##     CADM1
##  2.554173


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
#RNA-seq解析の最初のステップとして、サンプル間の全体的な類似性を評価することが有効です：
#どのサンプルが似ていて、どのサンプルが違うか？
#これは実験計画からの予想と合っているか？
#データセットの主な変動要因は何か？
#サンプルの類似性を調べるために、主成分分析（PCA）と階層的クラスタリング法を用いてサンプルレベルのQCを実施する予定です。サンプルレベルのQCでは、レプリカがどの程度まとまっているかを確認し、実験条件がデータの主な変動要因であるかどうかを観察することができます。また、サンプルレベルのQCを行うことで、サンプルの外れ値を特定することができます。この外れ値は、DE分析の前に除去する必要があるかどうかを判断するために、さらに調査する必要があるかもしれません。

#これらの教師なしクラスタリング手法を使用する場合、カウントの正規化とlog2変換を行うことで、距離やクラスタリングが改善され可視化されます。DESeq2では、カウントの正規化にmedian of ratiosメソッドを使用し、サンプルレベルのQCに正規化カウントの正規化ログ変換（rlog）を使用することで、平均値に対する分散を緩和し、クラスタリングを改善します。


## Principal component analysis
# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)
saveRDS(rld, "rld.rds")

# Plot PCA

#p <- DESeq2::plotPCA(rld, intgroup = "group_id")
#ggsave(file="DESeq2_PCA.Bnaive.group.pdf", plot=p)
#
#p <- DESeq2::plotPCA(rld, intgroup = "institute_id")
#ggsave(file="DESeq2_PCA.Bnaive.institute.pdf", plot=p)
#
#p <- DESeq2::plotPCA(rld, intgroup = "sex")
#ggsave(file="DESeq2_PCA.Bnaive.sex.pdf", plot=p)
#
#p <- DESeq2::plotPCA(rld, intgroup = "scaled_age")
#ggsave(file="DESeq2_PCA.Bnaive.scaled_age.pdf", plot=p)


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

#p <- plotPCA(rld, intgroup = "institute_id", ax1="PC1", ax2="PC2")
#ggsave(file="DESeq2_PCA12.Bnaive.institute.pdf", plot=p)
#
#p <- plotPCA(rld, intgroup = "sex", ax1="PC1", ax2="PC2")
#ggsave(file="DESeq2_PCA12.Bnaive.sex.pdf", plot=p)
#
#p <- plotPCA(rld, intgroup = "scaled_age", ax1="PC1", ax2="PC2")
#ggsave(file="DESeq2_PCA12.Bnaive.scaled_age.pdf", plot=p)




p <- plotPCA(rld, intgroup = "group_id", ax1="PC3", ax2="PC4")
ggsave(file="DESeq2_PCA34.TRUE.group.pdf", plot=p)

#p <- plotPCA(rld, intgroup = "institute_id", ax1="PC3", ax2="PC4")
#ggsave(file="DESeq2_PCA34.Bnaive.institute.pdf", plot=p)
#
#p <- plotPCA(rld, intgroup = "sex", ax1="PC3", ax2="PC4")
#ggsave(file="DESeq2_PCA34.Bnaive.sex.pdf", plot=p)
#
#p <- plotPCA(rld, intgroup = "scaled_age", ax1="PC3", ax2="PC4")
#ggsave(file="DESeq2_PCA34.Bnaive.scaled_age.pdf", plot=p)





## Hierarchical clustering
# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)

# Plot heatmap
p <- pheatmap(rld_cor, annotation = cluster_metadata[, c("group_id"), drop=F])
ggsave(file="DESeq2_Clustering.TRUE.group.pdf", plot=p)

#p <- pheatmap(rld_cor, annotation = cluster_metadata[, c("institute_id"), drop=F])
#ggsave(file="DESeq2_Clustering.Bnaive.institute.pdf", plot=p)
#
#p <- pheatmap(rld_cor, annotation = cluster_metadata[, c("sex"), drop=F])
#ggsave(file="DESeq2_Clustering.Bnaive.sex.pdf", plot=p)
#
#p <- pheatmap(rld_cor, annotation = cluster_metadata[, c("scaled_age"), drop=F])
#ggsave(file="DESeq2_Clustering.Bnaive.scaled_age.pdf", plot=p)

# -> SG_SL_OI_00046（Ctrl, OI）が全然違うが、、、

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
##      log2 fold change (MLE): group_id Tumor vs Normal
##      Wald test p-value: group id Tumor vs Normal
##      DataFrame with 31844 rows and 6 columns
##                    baseMean log2FoldChange     lfcSE      stat    pvalue      padj
##                   <numeric>      <numeric> <numeric> <numeric> <numeric> <numeric>
##      MIR1302-2HG 0.00000000             NA        NA        NA        NA        NA
##      AL627309.1  0.01299569       -2.50459   2.55112 -0.981759  0.326219        NA
##      AL627309.3  0.00000000             NA        NA        NA        NA        NA
##      AL627309.4  0.00000000             NA        NA        NA        NA        NA
##      AL732372.1  0.00489816       -3.83322   3.33086 -1.150818  0.249807        NA
##      ...                ...            ...       ...       ...       ...       ...
##      AC011751.1    0.000000             NA        NA        NA        NA        NA
##      AC007244.1    0.000000             NA        NA        NA        NA        NA
##      AC010889.2    0.000000             NA        NA        NA        NA        NA
##      AC009494.2    0.000000             NA        NA        NA        NA        NA
##      AC004556.3    0.107952       -4.42858   3.33126   -1.3294  0.183716   0.33431

resShrink <- lfcShrink(dds,
                coef = resultsNames(dds)[2],
                res=res)

resShrink
##      log2 fold change (MAP): group_id Tumor vs Normal
##      Wald test p-value: group id Tumor vs Normal
##      DataFrame with 31844 rows and 5 columns
##                    baseMean log2FoldChange      lfcSE    pvalue      padj
##                   <numeric>      <numeric>  <numeric> <numeric> <numeric>
##      MIR1302-2HG 0.00000000             NA         NA        NA        NA
##      AL627309.1  0.01299569    9.49884e-08 0.00144269  0.326219        NA
##      AL627309.3  0.00000000             NA         NA        NA        NA
##      AL627309.4  0.00000000             NA         NA        NA        NA
##      AL732372.1  0.00489816   -1.01278e-07 0.00144270  0.249807        NA
##      ...                ...            ...        ...       ...       ...
##      AC011751.1    0.000000             NA         NA        NA        NA
##      AC007244.1    0.000000             NA         NA        NA        NA
##      AC010889.2    0.000000             NA         NA        NA        NA
##      AC009494.2    0.000000             NA         NA        NA        NA
##      AC004556.3    0.107952   -1.72624e-07  0.0014427  0.183716   0.33431



#   We will output our significant genes and perform a few different visualization techniques to explore our results:
#   
#   Table of results for all genes
#   Table of results for significant genes (padj < 0.05)
#   Scatterplot of normalized expression of top 20 most significant genes
#   Heatmap of all significant genes
#   Volcano plot of results



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
##  # A tibble: 1 × 6
##    gene  baseMean log2FoldChange lfcSE       pvalue        padj
##    <chr>    <dbl>          <dbl> <dbl>        <dbl>       <dbl>
##  1 ATM       8.56          -1.18 0.245 0.0000000911 0.000000562

# Check results output
res_tbl
##      # A tibble: 31,844 × 7
##         gene        baseMean log2FoldChange lfcSE    stat pvalue   padj
##         <chr>          <dbl>          <dbl> <dbl>   <dbl>  <dbl>  <dbl>
##       1 MIR1302-2HG 0               NA      NA    NA      NA     NA
##       2 AL627309.1  0.0130          -2.50    2.55 -0.982   0.326 NA
##       3 AL627309.3  0               NA      NA    NA      NA     NA
##       4 AL627309.4  0               NA      NA    NA      NA     NA
##       5 AL732372.1  0.00490         -3.83    3.33 -1.15    0.250 NA
##       6 AC114498.1  0.000613        -3.81    3.33 -1.14    0.253 NA
##       7 AL669831.2  0.000487        -3.85    3.33 -1.15    0.248 NA
##       8 AL669831.5  0.0935          -0.0491  1.42 -0.0346  0.972 NA
##       9 FAM87B      0               NA      NA    NA      NA     NA
##      10 LINC00115   0.101           -0.318   1.22 -0.260   0.795  0.880
##      # … with 31,834 more rows
##      # ℹ Use `print(n = ...)` to see more rows

resShrink_tbl
##      # A tibble: 31,844 × 6
##         gene        baseMean log2FoldChange    lfcSE pvalue   padj
##         <chr>          <dbl>          <dbl>    <dbl>  <dbl>  <dbl>
##       1 MIR1302-2HG 0         NA            NA       NA     NA
##       2 AL627309.1  0.0130     0.0000000950  0.00144  0.326 NA
##       3 AL627309.3  0         NA            NA       NA     NA
##       4 AL627309.4  0         NA            NA       NA     NA
##       5 AL732372.1  0.00490   -0.000000101   0.00144  0.250 NA
##       6 AC114498.1  0.000613   0.0000000259  0.00144  0.253 NA
##       7 AL669831.2  0.000487   0.0000000175  0.00144  0.248 NA
##       8 AL669831.5  0.0935     0.00000105    0.00144  0.972 NA
##       9 FAM87B      0         NA            NA       NA     NA
##      10 LINC00115   0.101      0.00000108    0.00144  0.795  0.880
##      # … with 31,834 more rows
##      # ℹ Use `print(n = ...)` to see more rows

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
##      # A tibble: 4,081 × 7
##         gene    baseMean log2FoldChange lfcSE  stat    pvalue      padj
##         <chr>      <dbl>          <dbl> <dbl> <dbl>     <dbl>     <dbl>
##       1 PTPRC      35.3           -9.40 0.420 -22.4 5.53e-111 8.04e-107
##       2 UBE2V1     12.4           -4.65 0.214 -21.8 3.89e-105 2.83e-101
##       3 EEF1G      70.4           -5.86 0.270 -21.7 9.81e-105 4.76e-101
##       4 SNHG29     34.2          -11.7  0.580 -20.2 6.91e- 91 2.51e- 87
##       5 ATP6V0C    15.9           -9.20 0.484 -19.0 1.71e- 80 4.97e- 77
##       6 SEPTIN7    22.4          -11.0  0.587 -18.7 4.67e- 78 1.13e- 74
##       7 GAS5       25.4          -11.1  0.598 -18.6 2.97e- 77 6.17e- 74
##       8 SNHG6      18.6          -10.7  0.585 -18.3 5.67e- 75 1.03e- 71
##       9 EVI2B       8.66          -7.16 0.397 -18.0 1.17e- 72 1.89e- 69
##      10 TOMM6      14.3          -10.4  0.580 -17.9 2.15e- 71 3.12e- 68
##      # … with 4,071 more rows
##      # ℹ Use `print(n = ...)` to see more rows

sig_resShrink
##      # A tibble: 4,081 × 6
##         gene    baseMean log2FoldChange lfcSE    pvalue      padj
##         <chr>      <dbl>          <dbl> <dbl>     <dbl>     <dbl>
##       1 PTPRC      35.3           -9.64 0.441 5.53e-111 8.04e-107
##       2 UBE2V1     12.4           -4.65 0.220 3.89e-105 2.83e-101
##       3 EEF1G      70.4           -5.84 0.279 9.81e-105 4.76e-101
##       4 SNHG29     34.2          -16.5  3.62  6.91e- 91 2.51e- 87
##       5 ATP6V0C    15.9           -9.71 0.577 1.71e- 80 4.97e- 77
##       6 SEPTIN7    22.4          -15.7  3.54  4.67e- 78 1.13e- 74
##       7 GAS5       25.4          -15.9  3.57  2.97e- 77 6.17e- 74
##       8 SNHG6      18.6          -15.5  3.52  5.67e- 75 1.03e- 71
##       9 EVI2B       8.66          -7.34 0.420 1.17e- 72 1.89e- 69
##      10 TOMM6      14.3          -15.0  3.45  2.15e- 71 3.12e- 68
##      # … with 4,071 more rows
##      # ℹ Use `print(n = ...)` to see more rows

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

### PaperとしてはDEG解析から HI gene (HIPred or DECIPHER) かつ adj.P<0.05かつlogFC<0　を抽出してOncoprintをfilteringするという流れにすべき
# /work23/home/nsasa/data/HPVrOPC_WGS_Analysis_forRevision_EGA/Fig8_GS18001T_noCloneSig_01Tmut3/nf-kb_genes/PCGR.groupBy_OncoPrint_HI/oncoprint_OnlyPrimary_OncoTSGs_amp/PCGR.TSG_HI.tsv
# /work23/home/nsasa/data/HPVrOPC_WGS_Analysis_forRevision_EGA/Fig8_GS18001T_noCloneSig_01Tmut3/nf-kb_genes/notPCGR_genes_Data/tmp/HI.genes
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


### plot Shrinkのthresholdで色
resShrink_table_thres <- resShrink_tbl[complete.cases(resShrink_tbl),] %>% 
                  mutate(threshold = padj < 0.05 & log2FoldChange <= -1) %>% select(gene, threshold)
res_table_thres <- res_tbl[complete.cases(res_tbl),]
res_table_thres <- dplyr::left_join(res_table_thres, resShrink_table_thres, by="gene")
label <- filtered_HI_data %>% filter(sig!="Sig") %>% select(SYMBOL_gene)
label$label <- label$SYMBOL_gene
res_table_thres <- dplyr::left_join(res_table_thres, label, by=c("gene"="SYMBOL_gene"))
    
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
ggsave(file=paste0(cluster_name, "_", levels(cluster_metadata$group_id)[2], "_vs_", levels(cluster_metadata$group_id)[1], "_Volcano.Shrink_thres.labeled.pdf"), plot=p)

p <- ggplot(res_table_thres) +
    geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
    xlab("log2 fold change") + 
    ylab("-log10 adjusted p-value") +
    theme(legend.position = "none",
          plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.title = element_text(size = rel(1.25))) +
    theme_bw() +
    scale_color_manual(values=c("#034789", "#E70808"))
ggsave(file=paste0(cluster_name, "_", levels(cluster_metadata$group_id)[2], "_vs_", levels(cluster_metadata$group_id)[1], "_Volcano.Shrink_thres.no_labeled.pdf"), plot=p)





# Figureで使ったTSG＆signalingに絞る    -> HIで絞らなくても下記の5つのみ！
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


#   # ちなみに消えたTRAF3は
#   resShrink_tbl[resShrink_tbl$gene=="TRAF3",]
#   ##  # A tibble: 1 × 6
#   ##    gene  baseMean log2FoldChange lfcSE  pvalue   padj
#   ##    <chr>    <dbl>          <dbl> <dbl>   <dbl>  <dbl>
#   ##  1 TRAF3     1.43         -0.950 0.469 0.00906 0.0301
#   
#   filtered_HI_data$sig <- (filtered_HI_data$SYMBOL %>% stringr::str_split(pattern = "_", simplify = TRUE))[,2]
#   table(filtered_HI_data$sig)
#   ##      Sig
#   ##    2   2
#   
#   write.table(filtered_HI_data, "filtered_HI_2TSG_2sig_genes.tsv", sep="\t", row.names=FALSE, quote=FALSE)
#   
#   write.table(filtered_HI_data %>% filter(sig!="Sig"), "filtered_HI_2TSG.tsv", sep="\t", row.names=FALSE, quote=FALSE)

### plot Shrinkのthresholdで色
resShrink_table_thres <- resShrink_tbl[complete.cases(resShrink_tbl),] %>% 
                  mutate(threshold = padj < 0.05 & log2FoldChange <= -1) %>% select(gene, threshold)
res_table_thres <- res_tbl[complete.cases(res_tbl),]
res_table_thres <- dplyr::left_join(res_table_thres, resShrink_table_thres, by="gene")
label <- filtered_HI_data %>% select(SYMBOL_gene)
label$label <- label$SYMBOL_gene
res_table_thres <- dplyr::left_join(res_table_thres, label, by=c("gene"="SYMBOL_gene"))
    
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
    scale_color_manual(values=c("#034789", "#E70808"))
ggsave(file=paste0(cluster_name, "_", levels(cluster_metadata$group_id)[2], "_vs_", levels(cluster_metadata$group_id)[1], "_Volcano.Shrink_thres.labeled.40perc.pdf"), plot=p)








### 逆に唯一あがっているのがBIRC3でありこれも面白い
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


















### ATM Vlnplog
vlns <- VlnPlot(seurat.sub, features=c("ATM", "BIRC2", "BIRC3", "CYLD", "GNAI2", "TRAF3"), group.by="epithelial_infercnvScore", pt.size=0, combine=FALSE, slot = "counts")
p <- CombinePlots(plots=vlns, ncol=1)
ggsave(file = "VlnPlot_5genes.pdf", plot = p , width = 10, height = 20)

seurat.sub <- SetIdent(seurat.sub, value="SampleGroup")
seurat.sub.sub <- subset(seurat.sub, idents="Tumor")

vlns <- VlnPlot(seurat.sub.sub, features=c("ATM", "BIRC2", "BIRC3", "CYLD", "GNAI2", "TRAF3"), group.by="hpv_samplegroup", pt.size=0, combine=FALSE, slot = "counts")
p <- CombinePlots(plots=vlns, ncol=1)
ggsave(file = "VlnPlot_5genes.tumor.pdf", plot = p , width = 10, height = 20)





















### CNA Gain + MUT (oncogene) ###
# Figureで使ったTSG＆signalingに絞る
oncogenelist <- c("PIK3CA", "BCL6", "MECOM", "ECT2", "SOX2", "PRKCI", "CIP2A", "CBLB", "TFG", "IQANK1", "BCL2L1", "BCL3", "MYC", "TP53INP2",
    "CD274", "FGFR3", "PIGU", "PLAGL2", "SKP2", "AKT2", "FASN", "NFATC2", "PABPC1", "ZNF217", "ZNF703", "ALK", "AXL", "CCND1", "CCNE1", "DDX5", "EGF",
    "EPCAM", "GOLPH3", "JAK2", "MTDH", "NCOA3", "NOTCH3", "NSD3", "PRAM1", "PRMT5", "SALL4", "TERT", "URI1", "YTHDF1", "ABL1", "AURKA", "BIRC5",
    "CBL", "CDH17", "CRP", "EEF1A2", "ERBB3", "FGFR1", "FOXK1", "FOXM1", "GNAS", "KAT6A", "KDM5A", "KITLG", "LSM1", "MEIS1", "MOS", "MST1R", "MYBL2",
    "NFIB", "NFKB1", "PBK", "PBX1", "PCNA", "PI3", "PRDM16", "PSIP1", "PTPN1", "PWWP3A", "RAB1A", "REL", "RET", "RUNX2", "SET", "SRC", "SRSF6", "TPD52",
    "UBE2C", "USP6", "VWA5A")
oncogenelist <- data.frame(gene=oncogenelist)

data <- resShrink_tbl

data <- dplyr::inner_join(data, oncogenelist, by="gene")

# remove NA
data <- data[complete.cases(data),]

# filtering 38 TSGsと87 signaling genes
filtered_data <- data[data$log2FoldChange >= 1 & data$padj<0.05,]
##      gene  baseMean log2FoldChange     lfcSE       pvalue         padj
##  2  EPCAM 3.4729875       2.345687 0.8982852 3.893917e-03 1.487110e-02
##  4  RAB1A 7.3722242       1.600812 0.2201281 7.715496e-14 2.213788e-12
##  11 PRKCI 1.7887359       1.654685 0.5472935 1.712944e-03 7.476379e-03
##  14  SOX2 5.4937741       2.930074 1.0830926 2.640330e-03 1.081814e-02
##  16 FGFR3 1.3733610       2.738562 1.0773714 1.062756e-02 3.403500e-02
##  23  LSM1 4.1056074       1.592771 0.2582640 4.468941e-10 8.794841e-09
##  30   MYC 7.1388191       3.282481 0.4575454 1.566769e-13 4.416677e-12
##  34  NFIB 2.9730001       2.146619 0.7254275 1.482617e-03 6.599110e-03
##  45 PRMT5 2.0039185       1.735349 0.5108395 5.040503e-04 2.638601e-03
##  49  PCNA 9.4712421       2.157958 0.3049562 3.569548e-13 9.680809e-12
##  52  PIGU 0.7737526       1.867710 0.6591902 8.304540e-03 2.799117e-02
##  55 SRSF6 5.2911572       1.184416 0.3645967 6.184795e-04 3.146109e-03
##  58 UBE2C 3.9620964       1.219284 0.5014186 8.875626e-03 2.948316e-02


write.table(filtered_data, "filtered_13oncogenes.tsv", sep="\t", row.names=FALSE, quote=FALSE)

### plot Shrinkのthresholdで色
resShrink_table_thres <- resShrink_tbl[complete.cases(resShrink_tbl),] %>% 
                  mutate(threshold = padj < 0.05 & log2FoldChange >= 1) %>% select(gene, threshold)
res_table_thres <- res_tbl[complete.cases(res_tbl),]
res_table_thres <- dplyr::left_join(res_table_thres, resShrink_table_thres, by="gene")
label <- filtered_data %>% select(gene)
label$label <- label$gene
res_table_thres <- dplyr::left_join(res_table_thres, label, by=c("gene"))
    
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
    scale_color_manual(values=c("#034789", "#E70808"))
ggsave(file="Oncogene_Volcano.Shrink_thres.labeled.pdf", plot=p)



















table(seurat@meta.data %>% select(predicted.celltype.l1, SampleGroup))
##                         SampleGroup
##  predicted.celltype.l1   Normal  Tumor
##    aDC                      570     57
##    B activated             9593     98
##    B memory               54760   3044
##    B naive               114607   2373
##    CD4 naive              40696   2109
##    CD4 Non-TFH             4135    342
##    CD4 TCM                21818   3577
##    CD4 TFH                36963   1242
##    CD4 TFH Mem             2468      6
##    CD4 TREG                8231   4262
##    CD8 naive               6659    244
##    CD8 T                   9909   6044
##    CD8 TCM                 2274    526
##    Cycling DZ GCB         40529    552
##    Cycling FDC              559     27
##    Cycling myeloid           64     38
##    Cycling T               1025    229
##    DC                      1488    425
##    dnT                      995     27
##    DZ GCB                  8172     19
##    DZtoLZ GCB transition  29908    195
##    Epithelial               953   6750
##    FCRL4/5+ B memory       9889     89
##    FDC                     2322   1860
##    Granulocytes             667    554
##    ILC                      740     25
##    LZ GCB                  7348     74
##    LZtoDZ GCB transition   4424     36
##    MAIT/TRDV2+ gdT          304     55
##    Mast                     122     79
##    Mono/Macro              6135   1507
##    NK                       251    420
##    NK_CD56bright            428     16
##    non-TRDV2+ gdT           524      7
##    PB                       296     16
##    PC                     13422   2056
##    PC/doublet               260      1
##    PDC                      622    172
##    preB/T                    96      0
##    preGCB                 10037    584
##    preMBC/doublet          1910     56
##    prePB                    250      0


for (celltype in c("B memory", "B naive", "CD4 TCM", "CD8 T")) {
    celltype02 = trimws(celltype)
    new_col = paste0("col_", celltype02)
    seurat@meta.data <- seurat@meta.data %>% mutate(
        celltype_trimmed = trimws(predicted.celltype.l1))
    seurat@meta.data <- seurat@meta.data %>% mutate(
                        !!new_col := case_when(
                            celltype_trimmed==celltype02 & SampleGroup=="Tumor" & HPV_group=="negative" & Sample!="HN12_CD45n" & Sample!="HN13_CD45n" & Sample!="HN14_CD45n" & Sample!="HN16_CD45n" & Sample!="HN17_CD45n" & Sample!="HN18_CD45n"  ~ celltype02,            ### TumorをHPV発現関係なくEpithelialかつCD45nで絞った
                            celltype_trimmed==celltype02 & SampleGroup=="Normal" ~ celltype02,
                            TRUE ~ "remove"
                            )
                        ) 
    
    seurat <- SetIdent(seurat, value=new_col)
    seurat.sub <- subset(seurat, idents=celltype02)
    
    
    
    # Extract raw counts and metadata to create SingleCellExperiment object
    counts <- seurat.sub@assays$RNA@counts 
    
    metadata <- seurat.sub@meta.data
    
    
    sample_celltype <- table(seurat.sub@meta.data$Sample, seurat.sub@meta.data$predicted.celltype.l1)
    write.csv(sample_celltype,
        paste0(celltype02, ".UMI_per_Celltype_per_Sample.csv"),
        quote=FALSE,
        row.names=TRUE)
    
    # Set up metadata as desired for aggregation and DE analysis
    #seurat.sub <- SetIdent(seurat.sub, value = "predicted.celltype.l1")
    #metadata$cluster_id <- factor(seurat@active.ident)
    metadata$cluster_id <- seurat.sub@active.ident %>% str_replace_all(pattern="_", replacement=" ") %>% factor()
    
    metadata$sample_id <- factor(metadata$Sample)
    metadata$group_id <- seurat.sub$SampleGroup %>% str_replace_all(pattern="-", replacement="_") %>% factor()
    #metadata$institute_id <- factor(seurat$Institute)
    
    
    # Create single cell experiment object
    sce <- SingleCellExperiment(assays = list(counts = counts), 
                               colData = metadata)
    
    per.cell <- perCellQCMetrics(sce)
    
    # Get cells w/ few/many detected genes
    sce$is_outlier <- isOutlier(
            metric = per.cell$total,
            nmads = 2, type = "both", log = TRUE)
    
    # Remove outlier cells
    sce <- sce[, !sce$is_outlier]
    
    ## Remove lowly expressed genes which have less than 10 cells with any counts
    sce <- sce[rowSums(counts(sce) > 1) >= 10, ]
    
    
    
    ### sampleがDROPすることでfactorの数が減るがlevelに残ってしまっている (BCLLATLAS_06_d8kwy76j_7blj9otf)
    # droplevels() で不要なレベルを削除
    sce$sample_id <- droplevels(sce$sample_id)
    
    ### Acquiring necessary metrics for aggregation across cells in a sample
    
    # Named vector of cluster names
    cids <- purrr::set_names(levels(sce$cluster_id))
    
    # Total number of clusters
    nc <- length(cids)
    
    # Named vector of sample names
    sids <- purrr::set_names(levels(sce$sample_id))
    
    # Total number of samples 
    ns <- length(sids)
    
    
    ### To perform sample-level differential expression analysis, we need to generate sample-level metadata. To do this, we will reorder samples in the single-cell metadata to match the order of the factor levels of the sample ID, then extract only the sample-level information from the first cell corresponding to that sample.
    # Generate sample level metadata
        
    ## Turn named vector into a numeric vector of number of cells per sample
    n_cells <- as.numeric(table(sce$sample_id))
    
    ## Determine how to reoder the samples (rows) of the metadata to match the order of sample names in sids vector
    m <- match(sids, sce$sample_id)
    
    ## Create the sample level metadata by combining the reordered metadata with the number of cells corresponding to each sample.
    ei <- data.frame(colData(sce)[m, ], 
                      n_cells, row.names = NULL) %>% 
                    select(-"cluster_id")    
    
    
    ### QC https://bioconductor.org/packages/devel/bioc/vignettes/scuttle/inst/doc/qc.html
    # SKIP
    
    
    ### Count aggregation to sample level
    # Aggregate the counts per sample_id and cluster_id
    
    # Subset metadata to only include the cluster and sample IDs to aggregate across
    groups <- colData(sce)[, c("cluster_id", "sample_id")]
    
    # Aggregate across cluster-sample groups
    pb <- aggregate.Matrix(t(counts(sce)), 
                           groupings = groups, fun = "sum") 
    
    ### To perform DE analysis on a per cell type basis, we need to wrangle our data in a couple ways. We need to do the following steps:
    ## 1. Split our data by cell type
    ## 2. Transform the matrix so that the genes are the row names and the samples are the column names
    
    # Not every cluster is present in all samples; create a vector that represents how to split samples
    splitf <- sapply(stringr::str_split(rownames(pb), 
                                        pattern = "_",
                                        n = 2),
                    `[`, 1)

    pb <- split.data.frame(pb, 
                           factor(splitf)) %>%
            lapply(function(u) 
                    set_colnames(t(u), 
                                 sapply(stringr::str_split(rownames(u), pattern = "_", n = 2), `[`, 2)))
    
    
    # Print out the table of cells in each cluster-sample group
    options(width = 100)
    
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
    
    # Get cluster IDs for each of the samples
    
    samples_list <- map(1:length(cids), get_sample_ids)
    
    get_cluster_ids <- function(x){
            rep(names(pb)[x], 
                each = length(samples_list[[x]]))
    }
    
    de_cluster_ids <- map(1:length(cids), get_cluster_ids) %>%
            unlist()
    
    
    # Create a data frame with the sample IDs, cluster IDs and condition
    #factorにしておかないと下のlevels()がNULLを返してしまう
    gg_df <- data.frame(cluster_id = factor(de_cluster_ids),
                        sample_id = factor(de_samples))
    
    gg_df <- left_join(gg_df, ei[, c("sample_id", "group_id")]) 
    
    
    metadata <- gg_df %>%
            dplyr::select(cluster_id, sample_id, group_id)
    
    
    
    write.csv(metadata, paste0(celltype02, ".metadata_for_DESeq2.csv"), quote=FALSE, row.names=FALSE)
    
    #saveRDS(pb, "pb.rds")
    #saveRDS(metadata, "metadata.rds")
    # pb <- readRDS("pb.rds")
    # metadata <- readRDS("metadata.rds")
    
    
    ### Subsetting dataset to cluster(s) of interest
    # Generate vector of cluster IDs
    clusters <- levels(metadata$cluster_id)
    
    # Subset the metadata to only the B cells                                       ### SG_SL_KR_00011ではB naiveなし
    cluster_metadata <- metadata[which(metadata$cluster_id == clusters[1]), ]
    
    # Assign the rownames of the metadata to be the sample IDs
    rownames(cluster_metadata) <- cluster_metadata$sample_id
    
    # Subset the counts to only the B cells
    counts <- pb[[clusters[1]]]
    
    cluster_counts <- data.frame(counts[, which(colnames(counts) %in% rownames(cluster_metadata))])
    
    # Check that all of the row names of the metadata are the same and in the same order as the column names of the counts in order to use as input to DESeq2
    all(rownames(cluster_metadata) == colnames(cluster_counts))        
    ##  [1] TRUE
    
    
    #------------------------#
    ### filtering by count ###
    #------------------------#
    summary(rowSums(cluster_counts))
    ##      Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
    ##        63      994     2328    13489     5375 10630954
    
    ### カウント合計が10以上の遺伝子だけを残す
    cutoff <- 10
    table(rowSums(cluster_counts) >= cutoff)
    ##   TRUE
    ##  11047
    filtered_counts <- cluster_counts[rowSums(cluster_counts) >= cutoff, ]
    
    
    ### Create DESeq2 object
    
    dds <- DESeqDataSetFromMatrix(filtered_counts, 
                                  colData = cluster_metadata, 
                                  design = ~ group_id)
    
    ### Quality Control - sample level
    #RNA-seq解析の最初のステップとして、サンプル間の全体的な類似性を評価することが有効です：
    #どのサンプルが似ていて、どのサンプルが違うか？
    #これは実験計画からの予想と合っているか？
    #データセットの主な変動要因は何か？
    #サンプルの類似性を調べるために、主成分分析（PCA）と階層的クラスタリング法を用いてサンプルレベルのQCを実施する予定です。サンプルレベルのQCでは、レプリカがどの程度まとまっているかを確認し、実験条件がデータの主な変動要因であるかどうかを観察することができます。また、サンプルレベルのQCを行うことで、サンプルの外れ値を特定することができます。この外れ値は、DE分析の前に除去する必要があるかどうかを判断するために、さらに調査する必要があるかもしれません。
    
    #これらの教師なしクラスタリング手法を使用する場合、カウントの正規化とlog2変換を行うことで、距離やクラスタリングが改善され可視化されます。DESeq2では、カウントの正規化にmedian of ratiosメソッドを使用し、サンプルレベルのQCに正規化カウントの正規化ログ変換（rlog）を使用することで、平均値に対する分散を緩和し、クラスタリングを改善します。
    
    
    ## Principal component analysis
    # Transform counts for data visualization
    rld <- rlog(dds, blind=TRUE)
    
    
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
    ggsave(file=paste0(celltype02, ".DESeq2_PCA12.TRUE.group.pdf"), plot=p)
    
    p <- plotPCA(rld, intgroup = "group_id", ax1="PC3", ax2="PC4")
    ggsave(file=paste0(celltype02, ".DESeq2_PCA34.TRUE.group.pdf"), plot=p)
    
    
    
    
    ## Hierarchical clustering
    # Extract the rlog matrix from the object and compute pairwise correlation values
    rld_mat <- assay(rld)
    rld_cor <- cor(rld_mat)
    
    # Plot heatmap
    p <- pheatmap(rld_cor, annotation = cluster_metadata[, c("group_id"), drop=F])
    ggsave(file=paste0(celltype02, ".DESeq2_Clustering.TRUE.group.pdf"), plot=p)
    
    # ここで、除去が必要な外れ値や、デザイン式に回帰させたい追加の変動要因があるかどうかを判断します。PCAや階層的クラスタリングで外れ値が検出されず、回帰すべき追加の変動要因もないため、差分発現解析の実行に進むことができます。
    
    
    
    ### Running DESeq2
    #DESeq2は、ライブラリの深さの違いを考慮するために正規化因子（サイズ因子）を使用して、生のカウントをモデル化します。次に、遺伝子ごとの分散を推定し、その推定値を縮小して、より正確な分散推定値を生成してカウントをモデル化します。最後に、DESeq2は負の二項モデルを適合させ、Wald検定または尤度比検定を使用して仮説検定を実行します。これらのステップはすべて、追加資料で詳細に説明されています。
    
    # Run DESeq2 differential expression analysis
    dds <- DESeq(dds)
    
    # Plot dispersion estimates
    pdf(paste0(celltype02, ".DESeq2_DispersionEstimates.TRUE.pdf"))
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
    ##      log2 fold change (MLE): group_id Tumor vs Normal
    ##      Wald test p-value: group id Tumor vs Normal
    ##      DataFrame with 31844 rows and 6 columns
    ##                    baseMean log2FoldChange     lfcSE      stat    pvalue      padj
    ##                   <numeric>      <numeric> <numeric> <numeric> <numeric> <numeric>
    ##      MIR1302-2HG 0.00000000             NA        NA        NA        NA        NA
    ##      AL627309.1  0.01299569       -2.50459   2.55112 -0.981759  0.326219        NA
    ##      AL627309.3  0.00000000             NA        NA        NA        NA        NA
    ##      AL627309.4  0.00000000             NA        NA        NA        NA        NA
    ##      AL732372.1  0.00489816       -3.83322   3.33086 -1.150818  0.249807        NA
    ##      ...                ...            ...       ...       ...       ...       ...
    ##      AC011751.1    0.000000             NA        NA        NA        NA        NA
    ##      AC007244.1    0.000000             NA        NA        NA        NA        NA
    ##      AC010889.2    0.000000             NA        NA        NA        NA        NA
    ##      AC009494.2    0.000000             NA        NA        NA        NA        NA
    ##      AC004556.3    0.107952       -4.42858   3.33126   -1.3294  0.183716   0.33431
    
    resShrink <- lfcShrink(dds,
                    coef = resultsNames(dds)[2],
                    res=res)
    
    resShrink
    ##      log2 fold change (MAP): group_id Tumor vs Normal
    ##      Wald test p-value: group id Tumor vs Normal
    ##      DataFrame with 31844 rows and 5 columns
    ##                    baseMean log2FoldChange      lfcSE    pvalue      padj
    ##                   <numeric>      <numeric>  <numeric> <numeric> <numeric>
    ##      MIR1302-2HG 0.00000000             NA         NA        NA        NA
    ##      AL627309.1  0.01299569    9.49884e-08 0.00144269  0.326219        NA
    ##      AL627309.3  0.00000000             NA         NA        NA        NA
    ##      AL627309.4  0.00000000             NA         NA        NA        NA
    ##      AL732372.1  0.00489816   -1.01278e-07 0.00144270  0.249807        NA
    ##      ...                ...            ...        ...       ...       ...
    ##      AC011751.1    0.000000             NA         NA        NA        NA
    ##      AC007244.1    0.000000             NA         NA        NA        NA
    ##      AC010889.2    0.000000             NA         NA        NA        NA
    ##      AC009494.2    0.000000             NA         NA        NA        NA
    ##      AC004556.3    0.107952   -1.72624e-07  0.0014427  0.183716   0.33431
    
    
    
    #   We will output our significant genes and perform a few different visualization techniques to explore our results:
    #   
    #   Table of results for all genes
    #   Table of results for significant genes (padj < 0.05)
    #   Scatterplot of normalized expression of top 20 most significant genes
    #   Heatmap of all significant genes
    #   Volcano plot of results
    
    
    
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
    
    # Check results output
    res_tbl
    ##      # A tibble: 31,844 × 7
    ##         gene        baseMean log2FoldChange lfcSE    stat pvalue   padj
    ##         <chr>          <dbl>          <dbl> <dbl>   <dbl>  <dbl>  <dbl>
    ##       1 MIR1302-2HG 0               NA      NA    NA      NA     NA
    ##       2 AL627309.1  0.0130          -2.50    2.55 -0.982   0.326 NA
    ##       3 AL627309.3  0               NA      NA    NA      NA     NA
    ##       4 AL627309.4  0               NA      NA    NA      NA     NA
    ##       5 AL732372.1  0.00490         -3.83    3.33 -1.15    0.250 NA
    ##       6 AC114498.1  0.000613        -3.81    3.33 -1.14    0.253 NA
    ##       7 AL669831.2  0.000487        -3.85    3.33 -1.15    0.248 NA
    ##       8 AL669831.5  0.0935          -0.0491  1.42 -0.0346  0.972 NA
    ##       9 FAM87B      0               NA      NA    NA      NA     NA
    ##      10 LINC00115   0.101           -0.318   1.22 -0.260   0.795  0.880
    ##      # … with 31,834 more rows
    ##      # ℹ Use `print(n = ...)` to see more rows
    
    resShrink_tbl
    ##      # A tibble: 31,844 × 6
    ##         gene        baseMean log2FoldChange    lfcSE pvalue   padj
    ##         <chr>          <dbl>          <dbl>    <dbl>  <dbl>  <dbl>
    ##       1 MIR1302-2HG 0         NA            NA       NA     NA
    ##       2 AL627309.1  0.0130     0.0000000950  0.00144  0.326 NA
    ##       3 AL627309.3  0         NA            NA       NA     NA
    ##       4 AL627309.4  0         NA            NA       NA     NA
    ##       5 AL732372.1  0.00490   -0.000000101   0.00144  0.250 NA
    ##       6 AC114498.1  0.000613   0.0000000259  0.00144  0.253 NA
    ##       7 AL669831.2  0.000487   0.0000000175  0.00144  0.248 NA
    ##       8 AL669831.5  0.0935     0.00000105    0.00144  0.972 NA
    ##       9 FAM87B      0         NA            NA       NA     NA
    ##      10 LINC00115   0.101      0.00000108    0.00144  0.795  0.880
    ##      # … with 31,834 more rows
    ##      # ℹ Use `print(n = ...)` to see more rows
    
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
    ##      # A tibble: 4,081 × 7
    ##         gene    baseMean log2FoldChange lfcSE  stat    pvalue      padj
    ##         <chr>      <dbl>          <dbl> <dbl> <dbl>     <dbl>     <dbl>
    ##       1 PTPRC      35.3           -9.40 0.420 -22.4 5.53e-111 8.04e-107
    ##       2 UBE2V1     12.4           -4.65 0.214 -21.8 3.89e-105 2.83e-101
    ##       3 EEF1G      70.4           -5.86 0.270 -21.7 9.81e-105 4.76e-101
    ##       4 SNHG29     34.2          -11.7  0.580 -20.2 6.91e- 91 2.51e- 87
    ##       5 ATP6V0C    15.9           -9.20 0.484 -19.0 1.71e- 80 4.97e- 77
    ##       6 SEPTIN7    22.4          -11.0  0.587 -18.7 4.67e- 78 1.13e- 74
    ##       7 GAS5       25.4          -11.1  0.598 -18.6 2.97e- 77 6.17e- 74
    ##       8 SNHG6      18.6          -10.7  0.585 -18.3 5.67e- 75 1.03e- 71
    ##       9 EVI2B       8.66          -7.16 0.397 -18.0 1.17e- 72 1.89e- 69
    ##      10 TOMM6      14.3          -10.4  0.580 -17.9 2.15e- 71 3.12e- 68
    ##      # … with 4,071 more rows
    ##      # ℹ Use `print(n = ...)` to see more rows
    
    sig_resShrink
    ##      # A tibble: 4,081 × 6
    ##         gene    baseMean log2FoldChange lfcSE    pvalue      padj
    ##         <chr>      <dbl>          <dbl> <dbl>     <dbl>     <dbl>
    ##       1 PTPRC      35.3           -9.64 0.441 5.53e-111 8.04e-107
    ##       2 UBE2V1     12.4           -4.65 0.220 3.89e-105 2.83e-101
    ##       3 EEF1G      70.4           -5.84 0.279 9.81e-105 4.76e-101
    ##       4 SNHG29     34.2          -16.5  3.62  6.91e- 91 2.51e- 87
    ##       5 ATP6V0C    15.9           -9.71 0.577 1.71e- 80 4.97e- 77
    ##       6 SEPTIN7    22.4          -15.7  3.54  4.67e- 78 1.13e- 74
    ##       7 GAS5       25.4          -15.9  3.57  2.97e- 77 6.17e- 74
    ##       8 SNHG6      18.6          -15.5  3.52  5.67e- 75 1.03e- 71
    ##       9 EVI2B       8.66          -7.34 0.420 1.17e- 72 1.89e- 69
    ##      10 TOMM6      14.3          -15.0  3.45  2.15e- 71 3.12e- 68
    ##      # … with 4,071 more rows
    ##      # ℹ Use `print(n = ...)` to see more rows
    
    # Write significant results to file
    write.csv(sig_res,
              paste0(cluster_name, "_", levels(cluster_metadata$group_id)[2], "_vs_", levels(cluster_metadata$group_id)[1], "_sig_genes.csv"),
              quote = FALSE, 
              row.names = FALSE)
    write.csv(sig_resShrink,
              paste0(cluster_name, "_", levels(cluster_metadata$group_id)[2], "_vs_", levels(cluster_metadata$group_id)[1], "_sig_genes_lfcShrink.csv"),
              quote = FALSE, 
              row.names = FALSE)
}


