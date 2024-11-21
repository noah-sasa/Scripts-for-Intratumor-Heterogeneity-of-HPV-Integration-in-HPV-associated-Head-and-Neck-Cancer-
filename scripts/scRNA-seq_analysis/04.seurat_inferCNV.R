library(Seurat)

# Bring in Seurat object
seurat <- readRDS("../harmony_seurat02/Harmony_RNA_round1.rds")

seurat
##  An object of class Seurat
##  25523 features across 31286 samples within 1 assay
##  Active assay: RNA (25523 features, 0 variable features)
##   1 dimensional reduction calculated: umap

table(seurat$SampleGroup)
##  Normal  Tumor
##    3813  27473

### HPV positive or negative
hpv_expression <- FetchData(seurat, vars=c("HPV16-E1", "HPV16-E2", "HPV16-E5", "HPV16-L2", "HPV16-L1", "HPV16-E6", "HPV16-E7"))
hpv_expression$HPV <- hpv_expression$`HPV16-E1` + hpv_expression$`HPV16-E2` + hpv_expression$`HPV16-E5` + hpv_expression$`HPV16-L2` + hpv_expression$`HPV16-L1` + hpv_expression$`HPV16-E6` + hpv_expression$`HPV16-E7`
seurat$HPV_group <- ifelse(hpv_expression$HPV > 0, "positive", "negative")

df <- seurat@meta.data
table(df[df$SampleGroup=="Tumor" & df$HPV_group=="positive",]$predicted.celltype.l1)


table(seurat.sub$Sample)


### HPV-positiveでEpithelial以外は除外


table(df[df$SampleGroup=="Tumor" & df$predicted.celltype.l1=="Epithelial",]$HPV_group)
##  negative positive
##      5294     6099

### Tumor Epithelialの半分はHPV-negativeであり、Malignantかどうかを知りたい。

table(df[df$SampleGroup=="Normal" & df$predicted.celltype.l1=="Epithelial",]$HPV_group)
##  negative positive
##       204       11


library(tidyverse)



### predicted.celltype.l1でEpithelial に絞って解析
seurat <- SetIdent(seurat, value="predicted.celltype.l1")
seurat.sub <- subset(seurat, idents="Epithelial")

seurat.sub
##  An object of class Seurat
##  25523 features across 11608 samples within 1 assay
##  Active assay: RNA (25523 features, 0 variable features)
##   1 dimensional reduction calculated: umap

df <- seurat.sub@meta.data

table(seurat.sub$SampleGroup)


table(df %>% select(SampleGroup, HPV_group))
##             HPV_group
##  SampleGroup negative positive
##       Normal      204       11
##       Tumor      5294     6099

counts_matrix = GetAssayData(seurat.sub, slot="counts")
annot_df <- df[,"Sample", drop=FALSE]
write.table(annot_df, "annot_sample.txt", sep="\t", quote=FALSE, col.names=FALSE)

### 所属するCellが1個だけのannotationがある場合にERROR
##  Error in base::rowMeans(x, na.rm = na.rm, dims = dims, ...) :
##    'x' must be an array of at least two dimensions

annot_df <- df[,"SampleCase", drop=FALSE]
write.table(annot_df, "annot_samplecase.txt", sep="\t", quote=FALSE, col.names=FALSE)


### inferCNA
library(infercnv)


ref_names <- c("OP33norm", "OP34norm", "OP35norm")

infercnv_obj = CreateInfercnvObject(
    raw_counts_matrix = counts_matrix,
    annotations_file = "annot_sample.txt",
    gene_order_file = "/mnt/d/Visium/Seurat_HNSCC/inferCNV/hg38_gencode_v27.txt",
    ref_group_names = ref_names
    )
##  INFO [2024-10-01 18:17:52] Parsing gene order file: /mnt/d/Visium/Seurat_HNSCC/inferCNV/hg38_gencode_v27.txt
##  INFO [2024-10-01 18:17:52] Parsing cell annotations file: annot_sample.txt
##  INFO [2024-10-01 18:17:52] ::order_reduce:Start.
##  INFO [2024-10-01 18:17:52] .order_reduce(): expr and order match.
##  INFO [2024-10-01 18:17:53] ::process_data:order_reduce:Reduction from positional data, new dimensions (r,c) = 23204,5330 Total=39930422 Min=0 Max=14337.
##  INFO [2024-10-01 18:17:53] num genes removed taking into account provided gene ordering list: 5545 = 23.8967419410446% removed.
##  INFO [2024-10-01 18:17:53] -filtering out cells < 100 or > Inf, removing 0 % of cells
##  INFO [2024-10-01 18:17:55] validating infercnv_obj

out_dir = "output_dir"
infercnv_obj = infercnv::run(
    infercnv_obj,
    cutoff = 0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
    out_dir = out_dir,
    cluster_by_groups = TRUE,
    denoise = TRUE,
    HMM = TRUE
    )


### Tumor heterogeneity
### https://github.com/broadinstitute/inferCNV/wiki/infercnv-tumor-subclusters
infercnv_obj_subclusters = CreateInfercnvObject(
    raw_counts_matrix = counts_matrix,
    annotations_file = "annot_sample.txt",
    gene_order_file = "/mnt/d/Visium/Seurat_HNSCC/inferCNV/hg38_gencode_v27.txt",
    ref_group_names = ref_names
    )
out_dir = "output_dir_subclusters"
infercnv_obj_subclusters = infercnv::run(
    infercnv_obj_subclusters,
    cutoff = 0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
    out_dir = out_dir,
    cluster_by_groups = TRUE,
    denoise = TRUE,
    HMM = TRUE,
    analysis_mode='subclusters',
    output_format='pdf'
    )

saveRDS(infercnv_obj_subclusters, "infercnv_obj_subclusters.rds")


### discover subgroups
# infercnv.20_HMM_predHMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations_dendrogram.txt　の各行からnwk作成

# cell id確認
colnames(infercnv_obj_subclusters@expr.data)

# subcluster確認
# infercnv_obj_subclusters@tumor_subclusters$subclusters$

# HN12のnormal疑い
infercnv_obj_subclusters@tumor_subclusters$subclusters$HN12$HN12_s6
# HN13のnormal疑い NASI
# HN14のnormal疑い NASI
# HN16のnormal疑い
infercnv_obj_subclusters@tumor_subclusters$subclusters$HN16$HN16_s3
infercnv_obj_subclusters@tumor_subclusters$subclusters$HN16$HN16_s4
# HN17のnormal疑い HN17はCD45pとCD45nが重複している気がしないでもない、、、
infercnv_obj_subclusters@tumor_subclusters$subclusters$HN17$HN17_s7
# HN18のnormal疑い NASI






#--------------------#
### x1以外の遺伝子数 ###
#--------------------#
###  based on the HMM predictions that give you specific CNV fold change levels 
ref <- read.table("output_dir_subclusters/infercnv.20_HMM_predHMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.references.txt", sep=" ", header=TRUE, row.names=1, check.names=FALSE)
obs <- read.table("output_dir_subclusters/infercnv.20_HMM_predHMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations.txt", sep=" ", header=TRUE, row.names=1, check.names=FALSE)

ref_sum <-  data.frame(sapply(ref, function(x) sum(x != 1)))
library(tidyverse)
ref_tbl <- ref_sum %>%
        data.frame() %>%
        rownames_to_column(var="cell") %>%
        as_tibble()
ref_tbl$sample <- sapply(ref_tbl$cell, function(text) {
    split_text <- str_split(text, "_", simplify=TRUE)
    str_c(split_text[-length(split_text)], collapse="_")
    })
ref_tbl$sampleCase <- sapply(ref_tbl$cell, function(text) {
    split_text <- str_split(text, "_", simplify=TRUE)
    str_c(split_text[1:2], collapse="_")
    })
ref_tbl$samplegroup <- "Normal"
ref_tbl <- ref_tbl %>% mutate(
    subclusters = case_when(
        ref_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP33norm$OP33norm_s1] ~ "OP33norm_s1",
        ref_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP33norm$OP33norm_s2] ~ "OP33norm_s2",
        ref_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP33norm$OP33norm_s3] ~ "OP33norm_s3",
        ref_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP33norm$OP33norm_s4] ~ "OP33norm_s4",
        ref_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP34norm$OP34norm_s1] ~ "OP34norm_s1",
        ref_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP34norm$OP34norm_s2] ~ "OP34norm_s2",
        ref_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP34norm$OP34norm_s3] ~ "OP34norm_s3",
        ref_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP35norm$OP35norm] ~ "OP35norm"
        ))
colnames(ref_tbl) <- c("cell", "score", "sample", "sampleCase", "samplegroup", "subclusters")


obs_sum <- data.frame(sapply(obs, function(x) sum(x != 1)))
library(tidyverse)
obs_tbl <- obs_sum %>%
        data.frame() %>%
        rownames_to_column(var="cell") %>%
        as_tibble()
obs_tbl$sample <- sapply(obs_tbl$cell, function(text) {
    split_text <- str_split(text, "_", simplify=TRUE)
    str_c(split_text[-length(split_text)], collapse="_")
    })
obs_tbl$sampleCase <- sapply(obs_tbl$cell, function(text) {
    split_text <- str_split(text, "_", simplify=TRUE)
    str_c(split_text[1], collapse="_")
    })
obs_tbl$samplegroup <- "Tumor"
obs_tbl <- obs_tbl %>% mutate(
    subclusters = case_when(
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP33$OP33_s1] ~ "OP33_s1",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP33$OP33_s2] ~ "OP33_s2",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP33$OP33_s3] ~ "OP33_s3",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP33$OP33_s4] ~ "OP33_s4",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP33$OP33_s5] ~ "OP33_s5",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP33$OP33_s6] ~ "OP33_s6",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP34$OP34_s1] ~ "OP34_s1",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP34$OP34_s2] ~ "OP34_s2",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP34$OP34_s3] ~ "OP34_s3",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP34$OP34_s4] ~ "OP34_s4",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP35BOT$OP35BOT_s1] ~ "OP35BOT_s1",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP35BOT$OP35BOT_s2] ~ "OP35BOT_s2",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP35BOT$OP35BOT_s3] ~ "OP35BOT_s3",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP35BOT$OP35BOT_s4] ~ "OP35BOT_s4",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP35LN$OP35LN_s1] ~ "OP35LN_s1",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP35LN$OP35LN_s2] ~ "OP35LN_s2",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP35LN$OP35LN_s3] ~ "OP35LN_s3",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP35LN$OP35LN_s4] ~ "OP35LN_s4",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP4$OP4_s1] ~ "OP4_s1",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP4$OP4_s2] ~ "OP4_s2",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP4$OP4_s3] ~ "OP4_s3",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP4$OP4_s4] ~ "OP4_s4",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP4$OP4_s5] ~ "OP4_s5",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP4$OP4_s6] ~ "OP4_s6",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP4$OP4_s7] ~ "OP4_s7",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP4$OP4_s8] ~ "OP4_s8",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP4$OP4_s9] ~ "OP4_s9",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP5$OP5_s1] ~ "OP5_s1",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP5$OP5_s2] ~ "OP5_s2",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP5$OP5_s3] ~ "OP5_s3",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP5$OP5_s4] ~ "OP5_s4",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP5$OP5_s5] ~ "OP5_s5",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP6$OP6_s1] ~ "OP6_s1",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP6$OP6_s2] ~ "OP6_s2",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP6$OP6_s3] ~ "OP6_s3",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP6$OP6_s4] ~ "OP6_s4",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP9$OP9_s1] ~ "OP9_s1",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP9$OP9_s2] ~ "OP9_s2",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP9$OP9_s3] ~ "OP9_s3",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP9$OP9_s4] ~ "OP9_s4",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP9$OP9_s5] ~ "OP9_s5",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP13$OP13_s1] ~ "OP13_s1",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP13$OP13_s2] ~ "OP13_s2",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP14$OP14_s1] ~ "OP14_s1",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP14$OP14_s2] ~ "OP14_s2",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP14$OP14_s3] ~ "OP14_s3",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP14$OP14_s4] ~ "OP14_s4",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP17$OP17_s1] ~ "OP17_s1",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP17$OP17_s2] ~ "OP17_s2",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP17$OP17_s3] ~ "OP17_s3",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP20$OP20_s1] ~ "OP20_s1",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP20$OP20_s2] ~ "OP20_s2",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP20$OP20_s3] ~ "OP20_s3",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP20$OP20_s4] ~ "OP20_s4",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP20$OP20_s5] ~ "OP20_s5",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP20$OP20_s6] ~ "OP20_s6",
        obs_tbl$cell %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP20$OP20_s7] ~ "OP20_s7"
        ))
colnames(obs_tbl) <- c("cell", "score", "sample", "sampleCase", "samplegroup", "subclusters")

df <- rbind(ref_tbl, obs_tbl)
summary(ref_tbl$score)
##     Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
##    0.000   0.000   0.000   8.781  32.000  32.000


summary(obs_tbl$score)
##     Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
##        0    1026    1531    1359    1755    2826

summary(df$score)
##     Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
##        0     915    1531    1334    1755    2826

df$subclusters <- factor(df$subclusters, levels=c(
    "OP4_s1", "OP4_s2", "OP4_s3", "OP4_s4", "OP4_s5", "OP4_s6", "OP4_s7", "OP4_s8", "OP4_s9",
    "OP5_s1", "OP5_s2", "OP5_s3", "OP5_s4", "OP5_s5",
    "OP6_s1", "OP6_s2", "OP6_s3", "OP6_s4",
    "OP9_s1", "OP9_s2", "OP9_s3", "OP9_s4", "OP9_s5",
    "OP13_s1", "OP13_s2",
    "OP14_s1", "OP14_s2", "OP14_s3", "OP14_s4",
    "OP17_s1", "OP17_s2", "OP17_s3",
    "OP20_s1", "OP20_s2", "OP20_s3", "OP20_s4", "OP20_s5", "OP20_s6", "OP20_s7",
    "OP33norm_s1", "OP33norm_s2", "OP33norm_s3", "OP33norm_s4",
    "OP33_s1", "OP33_s2", "OP33_s3", "OP33_s4", "OP33_s5", "OP33_s6",
    "OP34norm_s1", "OP34norm_s2", "OP34norm_s3",
    "OP34_s1", "OP34_s2", "OP34_s3", "OP34_s4",
    "OP35norm", "OP35BOT_s1", "OP35BOT_s2", "OP35BOT_s3", "OP35BOT_s4",
    "OP35LN_s1", "OP35LN_s2", "OP35LN_s3", "OP35LN_s4"
    ))

library(ggplot2)
p <- ggplot(data=df, aes(x=subclusters, y=score, fill = samplegroup)) +
        stat_summary(fun = mean, geom = "bar", color = 'black', width = 0.5) +
        theme_classic() +
        scale_fill_manual(values=c("#034789", "#E70808")) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        geom_hline(yintercept=397.00*3, linetype='dotted', col = 'black')
ggsave(file="output_dir_subclusters/HMM_score.v2.pdf", plot=p, width=14)

library(ggplot2)
p <- ggplot(data=df, aes(x=sampleCase, y=score, fill = samplegroup)) +
        geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0, color = 'black') +
        geom_point() +
        theme_classic() +
        scale_fill_manual(values=c("#034789", "#E70808")) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        geom_hline(yintercept=397.00*3, linetype='dotted', col = 'black')
ggsave(file="output_dir_subclusters/HMM_score.v2.sampleCase.pdf", plot=p)


# HPV 確認
df <- dplyr::left_join(df, seurat.sub@meta.data %>% select(HPV_group) %>% rownames_to_column(), by=c("cell"="rowname")) 
table(df %>% select(subclusters, HPV_group))

# HPV-positive率
hpv_posrate <- df %>%
  group_by(subclusters) %>%
  summarise(
    positive_rate = mean(HPV_group == "positive")  # positiveの割合を計算
  )


df <- dplyr::left_join(df, hpv_posrate, by="subclusters")


p <- ggplot(data=df, aes(x=subclusters, y=score, fill = positive_rate)) +
        stat_summary(fun = mean, geom = "bar", color = 'black', width = 0.5) +
        theme_classic() +
        scale_fill_gradient(low="#F2E5F4", high="#4A138A") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(file="output_dir_subclusters/HMM_score.v2.gradient.pdf", plot=p, width=14)


# Subcluster内のcell数
cell_n <- table(df %>% select(subclusters)) %>% data.frame
colnames(cell_n) <- c("subclusters", "cell_n")

df <- dplyr::left_join(df, cell_n, by="subclusters")
p <- ggplot(data=df) +
        stat_summary(aes(x=subclusters, y=score, fill = positive_rate), fun = mean, geom = "bar", color = 'black', width = 0.5) +
        stat_summary(aes(x=subclusters, y=-cell_n, fill = positive_rate), fun = mean, geom = "bar", color = 'black', width = 0.5) +
        theme_classic() +
        scale_fill_gradient(low="#F2E5F4", high="#4A138A") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(file="output_dir_subclusters/HMM_score.v2.gradient.cell_n.pdf", plot=p, width=14)

p <- ggplot(data=df) +
        stat_summary(aes(x=subclusters, y=score, fill = cell_n), fun = mean, geom = "bar", color = 'black', width = 0.5) +
        stat_summary(aes(x=subclusters, y=- positive_rate*1000, fill = cell_n), fun = mean, geom = "bar", color = 'black', width = 0.5) +
        theme_classic() +
        scale_fill_gradient(low="#FAE9E6", high="#BF350C") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(file="output_dir_subclusters/HMM_score.v2.gradient.positive_rate.pdf", plot=p, width=14)






### subcluster情報をseuratに
library(tidyverse)
seurat.sub@meta.data <- seurat.sub@meta.data %>% mutate(
    subclusters = case_when(
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP33$OP33_s1] ~ "OP33_s1",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP33$OP33_s2] ~ "OP33_s2",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP33$OP33_s3] ~ "OP33_s3",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP33$OP33_s4] ~ "OP33_s4",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP33$OP33_s5] ~ "OP33_s5",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP33$OP33_s6] ~ "OP33_s6",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP34$OP34_s1] ~ "OP34_s1",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP34$OP34_s2] ~ "OP34_s2",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP34$OP34_s3] ~ "OP34_s3",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP34$OP34_s4] ~ "OP34_s4",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP35BOT$OP35BOT_s1] ~ "OP35BOT_s1",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP35BOT$OP35BOT_s2] ~ "OP35BOT_s2",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP35BOT$OP35BOT_s3] ~ "OP35BOT_s3",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP35BOT$OP35BOT_s4] ~ "OP35BOT_s4",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP35LN$OP35LN_s1] ~ "OP35LN_s1",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP35LN$OP35LN_s2] ~ "OP35LN_s2",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP35LN$OP35LN_s3] ~ "OP35LN_s3",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP35LN$OP35LN_s4] ~ "OP35LN_s4",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP4$OP4_s1] ~ "OP4_s1",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP4$OP4_s2] ~ "OP4_s2",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP4$OP4_s3] ~ "OP4_s3",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP4$OP4_s4] ~ "OP4_s4",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP4$OP4_s5] ~ "OP4_s5",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP4$OP4_s6] ~ "OP4_s6",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP4$OP4_s7] ~ "OP4_s7",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP4$OP4_s8] ~ "OP4_s8",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP4$OP4_s9] ~ "OP4_s9",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP5$OP5_s1] ~ "OP5_s1",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP5$OP5_s2] ~ "OP5_s2",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP5$OP5_s3] ~ "OP5_s3",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP5$OP5_s4] ~ "OP5_s4",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP5$OP5_s5] ~ "OP5_s5",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP6$OP6_s1] ~ "OP6_s1",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP6$OP6_s2] ~ "OP6_s2",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP6$OP6_s3] ~ "OP6_s3",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP6$OP6_s4] ~ "OP6_s4",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP9$OP9_s1] ~ "OP9_s1",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP9$OP9_s2] ~ "OP9_s2",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP9$OP9_s3] ~ "OP9_s3",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP9$OP9_s4] ~ "OP9_s4",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP9$OP9_s5] ~ "OP9_s5",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP13$OP13_s1] ~ "OP13_s1",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP13$OP13_s2] ~ "OP13_s2",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP14$OP14_s1] ~ "OP14_s1",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP14$OP14_s2] ~ "OP14_s2",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP14$OP14_s3] ~ "OP14_s3",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP14$OP14_s4] ~ "OP14_s4",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP17$OP17_s1] ~ "OP17_s1",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP17$OP17_s2] ~ "OP17_s2",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP17$OP17_s3] ~ "OP17_s3",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP20$OP20_s1] ~ "OP20_s1",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP20$OP20_s2] ~ "OP20_s2",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP20$OP20_s3] ~ "OP20_s3",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP20$OP20_s4] ~ "OP20_s4",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP20$OP20_s5] ~ "OP20_s5",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP20$OP20_s6] ~ "OP20_s6",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP20$OP20_s7] ~ "OP20_s7",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP33norm$OP33norm_s1] ~ "OP33norm_s1",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP33norm$OP33norm_s2] ~ "OP33norm_s2",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP33norm$OP33norm_s3] ~ "OP33norm_s3",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP33norm$OP33norm_s4] ~ "OP33norm_s4",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP34norm$OP34norm_s1] ~ "OP34norm_s1",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP34norm$OP34norm_s2] ~ "OP34norm_s2",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP34norm$OP34norm_s3] ~ "OP34norm_s3",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP35norm$OP35norm] ~ "OP35norm",
        TRUE ~ "no_subcluster"
        ))

seurat.sub@meta.data <- seurat.sub@meta.data %>% mutate(
    epithelial_infercnvScore = case_when(
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP13$OP13_s2] ~ "Normal",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP14$OP14_s1] ~ "Normal",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP14$OP14_s2] ~ "Normal",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP17$OP17_s1] ~ "Normal",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP20$OP20_s4] ~ "Normal",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP33$OP33_s4] ~ "Normal",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP33norm$OP33norm_s1] ~ "Normal",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP33norm$OP33norm_s2] ~ "Normal",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP33norm$OP33norm_s3] ~ "Normal",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP33norm$OP33norm_s4] ~ "Normal",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP34$OP34_s3] ~ "Normal",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP34norm$OP34norm_s1] ~ "Normal",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP34norm$OP34norm_s2] ~ "Normal",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP34norm$OP34norm_s3] ~ "Normal",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP35BOT$OP35BOT_s4] ~ "Normal",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP35norm$OP35norm] ~ "Normal",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP4$OP4_s2] ~ "Normal",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP5$OP5_s5] ~ "Normal",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP6$OP6_s4] ~ "Normal",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP9$OP9_s2] ~ "Normal",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP9$OP9_s3] ~ "Normal",
        rownames(seurat.sub@meta.data) %in% colnames(infercnv_obj_subclusters@expr.data)[infercnv_obj_subclusters@tumor_subclusters$subclusters$OP9$OP9_s5] ~ "Normal",
        TRUE ~ "Tumor"
        ))


saveRDS(seurat.sub, "Harmony_RNA_round1.epi.inferCNVScore02.subcluster.rds")












### CD45n と　predicted.celltype.l1でEpithelial に絞って解析
seurat <- SetIdent(seurat, value="epithelial")
seurat.sub <- subset(seurat, idents="Epithelial")

seurat.sub
##  An object of class Seurat
##  31844 features across 6225 samples within 1 assay
##  Active assay: RNA (31844 features, 0 variable features)
##   1 dimensional reduction calculated: umap

df <- seurat.sub@meta.data

counts_matrix = GetAssayData(seurat.sub, slot="counts")
annot_df <- df[,"Sample", drop=FALSE]
write.table(annot_df, "annot_sample.CD45n.txt", sep="\t", quote=FALSE, col.names=FALSE)

### 所属するCellが1個だけのannotationがある場合にERROR
##  Error in base::rowMeans(x, na.rm = na.rm, dims = dims, ...) :
##    'x' must be an array of at least two dimensions

annot_df <- df[,"SampleCase", drop=FALSE]
write.table(annot_df, "annot_samplecase.CD45n.txt", sep="\t", quote=FALSE, col.names=FALSE)


### inferCNA
library(infercnv)

ref_names <- c("BCLLATLAS_05", "BCLLATLAS_06", "BCLLATLAS_07", "BCLLATLAS_131", "BCLLATLAS_16", "BCLLATLAS_19",
    "BCLLATLAS_21", "BCLLATLAS_25", "BCLLATLAS_34", "BCLLATLAS_41", "BCLLATLAS_54", "BCLLATLAS_57")


### Tumor heterogeneity
### https://github.com/broadinstitute/inferCNV/wiki/infercnv-tumor-subclusters
infercnv_obj_cd45n_subclusters = CreateInfercnvObject(
    raw_counts_matrix = counts_matrix,
    annotations_file = "annot_samplecase.CD45n.txt",
    gene_order_file = "hg38_gencode_v27.txt",
    ref_group_names = ref_names
    )
out_dir = "output_dir_CD45n_subclusters"
infercnv_obj_cd45n_subclusters = infercnv::run(
    infercnv_obj_cd45n_subclusters,
    cutoff = 0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
    out_dir = out_dir,
    cluster_by_groups = TRUE,
    denoise = TRUE,
    HMM = TRUE,
    analysis_mode='subclusters'
    )

saveRDS(infercnv_obj_cd45n_subclusters, "infercnv_obj_cd45n_subclusters.rds")
