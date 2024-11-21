seurat.combined.Scdsed.filt <- readRDS("../seurat_qc02/seurat.combined.Scdsed.filt")




# harmony + seurat https://github.com/satijalab/seurat-wrappers/blob/master/docs/harmony.md
# seurat https://kero.hgc.jp/open_tutorials/learning/contents/single2_article.html
#       https://github.com/NBISweden/workshop-scRNAseq/blob/master/labs/compiled/seurat/seurat_02_dim_reduction.md
#       https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

library(harmony)
library(Seurat)

# Harmony + seuratから
#   pbmcsca <- NormalizeData(pbmcsca) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)
#   pbmcsca <- RunHarmony(pbmcsca, group.by.vars = "Method")
#   pbmcsca <- RunUMAP(pbmcsca, reduction = "harmony", dims = 1:30)
#   pbmcsca <- FindNeighbors(pbmcsca, reduction = "harmony", dims = 1:30) %>% FindClusters()
#   DimPlot(pbmcsca, group.by = c("Method", "ident", "CellType"), ncol = 3)

### SCTを用いて正規化  https://hbctraining.github.io/scRNA-seq_online/lessons/06a_integration_harmony.html



#   ### remove SG_SL_UO_00170
#   seurat.combined.miQCed.Scdsed.filt <- subset(seurat.combined.miQCed.Scdsed.filt, subset = Sample != "SG_SL_UO_00170")



### Harmonyは先にmerge推奨
#     # Merge raw samples
#     merged_seurat <- merge(x = raw_seurat_list[[1]],
#                    y = raw_seurat_list[2:length(raw_seurat_list)],
#                    merge.data = TRUE)
#     
#     # Perform log-normalization and feature selection, as well as SCT normalization on global object
#     merged_seurat <- merged_seurat %>%
#         NormalizeData() %>%
#         FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
#         ScaleData() %>%
#         SCTransform(vars.to.regress = c("mitoRatio"))
#     
#     # Calculate PCs using variable features determined by SCTransform (3000 by default)
#     merged_seurat <- RunPCA(merged_seurat, assay = "SCT", npcs = 50)


### SCTransformは個別に正規化推奨
##      bm40k.list <- SplitObject(hcabm40k, split.by = "orig.ident")
##      bm40k.list <- lapply(X = bm40k.list, FUN = function(x) {
##          x <- SCTransform(x, verbose = FALSE, do.scale = FALSE, do.center = FALSE)
##      })
##      features <- SelectIntegrationFeatures(object.list = bm40k.list)
##      bm40k.list <- lapply(X = bm40k.list, FUN = function(x) {
##          x <- ScaleData(x, features = features, verbose = FALSE)
##          x <- RunPCA(x, features = features, verbose = FALSE)
##      })





### 順を追って処理 https://github.com/immunogenomics/harmony/issues/41#issuecomment-642862186

### SplitObject
seurat.list <- SplitObject(object = seurat.combined.Scdsed.filt, split.by = "Sample")
rm(seurat.combined.Scdsed.filt)

### SCTransform normalization separately
for (i in 1:length(seurat.list)) {
    seurat.list[[i]] <- SCTransform(seurat.list[[i]], method = "glmGamPoi", vars.to.regress = "percent.mt", vst.flavor = "v2", verbose = TRUE)
}

### Integration Features without V/D/J genes & no_named_gene       # 予定のn_featuresより1000多くしている(V/D/J遺伝子などを除去するため)
for (i in 1:length(seurat.list)) {
    DefaultAssay(seurat.list[[i]]) <- "SCT"
}
seurat.features <- SelectIntegrationFeatures(object.list = seurat.list, nfeatures = 3000)

library(tidyverse)
seurat.features_tibble <- seurat.features %>% as_tibble
out_gene <- c("^TRAV","^TRAJ","^TRBV","^TRDV","^TRBJ","^IGH","^IGL","^IGK")                                              # B cell immunoglobulin data: IGH, IGK, IGL
seurat.features_tibble %>% filter(nchar(value) == 10) %>% filter(substr(.$value,9,9) == ".") %>% .$value -> noname_gene   # no_named_genes 
if (length(noname_gene) > 0) {
    seurat.features_tibble %>% filter(!str_detect(value, paste(out_gene, collapse = "|"))) %>% filter(!str_detect(value, paste(noname_gene, collapse = "|"))) %>% .$value -> gene_filtered
} else {
    seurat.features_tibble %>% filter(!str_detect(value, paste(out_gene, collapse = "|"))) %>% .$value -> gene_filtered
}

seurat.features <- gene_filtered[1:2000]

write(seurat.features, file = "2000_features_fromSCT.txt") 




saveRDS(seurat.list, "seurat.list_SCTed")
# seurat.list <- readRDS("seurat.list_SCTed")






### Azimuth annotation per Sample
library(Azimuth)

# Load the reference
# Change the file path based on where the reference is located on your system.
reference <- LoadReference(path = "/mnt/e/scRNA_SLE/azimuth_references/human_tonsil/v2.0.0")

# Preprocess with SCTransform
for (i in 1:length(seurat.list)) {
    seurat.list[[i]] <- SCTransform(
        object = seurat.list[[i]],
        assay = "RNA",
        new.assay.name = "refAssay",
        residual.features = rownames(x = reference$map),
        reference.SCT.model = reference$map[["refAssay"]]@SCTModel.list$refmodel,
        method = 'glmGamPoi',
        ncells = 2000,
        n_genes = 2000,
        do.correct.umi = FALSE,
        do.scale = FALSE,
        do.center = TRUE
        )    
}

# Find anchors between query and reference
anchors <- list()
for (i in 1:length(seurat.list)) {
    anchors[[i]] <- FindTransferAnchors(
        reference = reference$map,
        query = seurat.list[[i]],
        k.filter = NA,
        reference.neighbors = "refdr.annoy.neighbors",
        reference.assay = "refAssay",
        query.assay = "refAssay",
        reference.reduction = "refDR",
        normalization.method = "SCT",
        features = intersect(rownames(x = reference$map), VariableFeatures(object = seurat.list[[i]])),
        dims = 1:50,
        n.trees = 20,
        mapping.score.k = 100
        )
}

# Transfer cell type labels and impute protein expression
#
# Transferred labels are in metadata columns named "predicted.*"
# The maximum prediction score is in a metadata column named "predicted.*.score"
# The prediction scores for each class are in an assay named "prediction.score.*"
# The imputed assay is named "impADT" if computed

refdata <- lapply(X = c("celltype.l2", "celltype.l1"), function(x) {
  reference$map[[x, drop = TRUE]]
})
names(refdata) <- c("celltype.l2", "celltype.l1")
#   if (TRUE) {
#     refdata[["impADT"]] <- GetAssayData(
#       object = reference$map[['ADT']],
#       slot = 'data'
#     )
#   }

for (i in 1:length(seurat.list)) {
    seurat.list[[i]] <- TransferData(
        reference = reference$map,
        query = seurat.list[[i]],
        dims = 1:50,
        anchorset = anchors[[i]],
        refdata = refdata,
        n.trees = 20,
        store.weights = TRUE
        )
}
rm(refdata)

# Calculate the embeddings of the query data on the reference SPCA
for (i in 1:length(seurat.list)) {
    seurat.list[[i]] <- IntegrateEmbeddings(
        anchorset = anchors[[i]],
        reference = reference$map,
        query = seurat.list[[i]],
        reductions = "pcaproject",
        reuse.weights.matrix = TRUE
        )
}

# Calculate the query neighbors in the reference
# with respect to the integrated embeddings
for (i in 1:length(seurat.list)) {
    seurat.list[[i]][["query_ref.nn"]] <- FindNeighbors(
        object = Embeddings(reference$map[["refDR"]]),
        query = Embeddings(seurat.list[[i]][["integrated_dr"]]),
        return.neighbor = TRUE,
        l2.norm = TRUE
        )
}

# The reference used in the app is downsampled compared to the reference on which
# the UMAP model was computed. This step, using the helper function NNTransform,
# corrects the Neighbors to account for the downsampling.

## なぜかNNTransform functionが呼び出せないのでhelpers.Rからコピペ
#' Transform an NN index
#'
#' @param object Seurat object
#' @param meta.data Metadata
#' @param neighbor.slot Name of Neighbor slot
#' @param key Column of metadata to use
#'
#' @return \code{object} with transfomed neighbor.slot
#'
#' @importFrom SeuratObject Indices
#'
#' @keywords internal
#'
NNTransform <- function(
  object,
  meta.data,
  neighbor.slot = "query_ref.nn",
  key = 'ori.index'
) {
  on.exit(expr = gc(verbose = FALSE))
  ind <- SeuratObject::Indices(object[[neighbor.slot]])
  ori.index <- t(x = sapply(
    X = 1:nrow(x = ind),
    FUN = function(i) {
      return(meta.data[ind[i, ], key])
    }
  ))
  rownames(x = ori.index) <- rownames(x = ind)
  slot(object = object[[neighbor.slot]], name = "nn.idx") <- ori.index
  return(object)
}

for (i in 1:length(seurat.list)) {
    seurat.list[[i]] <- NNTransform(
        object = seurat.list[[i]],
        meta.data = reference$map[[]]
        )
}

# Project the query to the reference UMAP.
for (i in 1:length(seurat.list)) {
    seurat.list[[i]][["proj.umap"]] <- RunUMAP(
        object = seurat.list[[i]][["query_ref.nn"]],
        reduction.model = reference$map[["refUMAP"]],
        reduction.key = 'UMAP_'
        )
}

# Calculate mapping score and add to metadata
for (i in 1:length(seurat.list)) {
    seurat.list[[i]] <- AddMetaData(
        object = seurat.list[[i]],
        metadata = MappingScore(anchors = anchors[[i]]),
        col.name = "mapping.score"
        )
}
rm(anchors)
saveRDS(seurat.list,"seurat.list_Azimuthed")


# VISUALIZATIONS

# First predicted metadata field, change to visualize other predicted metadata
celltypes <- c("celltype.l2", "celltype.l1")

dir.create("Azimuth_Plot", recursive = T)
dir.create("Azimuth_Plot/Prediction_Score", recursive = T)
dir.create("Azimuth_Plot/Mapping_Score", recursive = T)
dir.create("Azimuth_Plot/DimPlot", recursive = T)
dir.create("Azimuth_Plot/RNA", recursive = T)


for (n in 1:length(celltypes)) {
    id <- celltypes[[n]]
    predicted.id <- paste0("predicted.", id)

    # DimPlot of the reference
    png(paste0("Azimuth_Plot/DimPlot/reference.",id,".png"), width = 3200, height = 1600)
    p <- DimPlot(object = reference$plot, reduction = "refUMAP", group.by = id, label = TRUE, repel=TRUE)
    p
    dev.off()

    for (i in 1:length(seurat.list)) {
        sampleid <- names(seurat.list[i])

        # DimPlot of the query, colored by predicted cell type
        png(paste0("Azimuth_Plot/DimPlot/",sampleid,".",predicted.id,".png"), width = 1600, height = 500)
        p <- DimPlot(object = seurat.list[[i]], reduction = "proj.umap", group.by = predicted.id, label = TRUE, repel=TRUE)
        plot(p)
        dev.off()
        
        # Plot the score for the predicted cell type of the query
        png(paste0("Azimuth_Plot/Prediction_Score/",sampleid,".",predicted.id,".png"), width = 2000, height = 500)
        p <- FeaturePlot(object = seurat.list[[i]], features = paste0(predicted.id, ".score"), reduction = "proj.umap")
        p <- p + VlnPlot(object = seurat.list[[i]], features = paste0(predicted.id, ".score"), group.by = predicted.id) + NoLegend()
        plot(p)
        dev.off()
        
        # Plot the mapping score
        png(paste0("Azimuth_Plot/Mapping_Score/",sampleid,".",predicted.id,".png"), width = 1600, height = 500)
        p <- FeaturePlot(object = seurat.list[[i]], features = "mapping.score", reduction = "proj.umap")
        p <- p + VlnPlot(object = seurat.list[[i]], features = "mapping.score", group.by = predicted.id) + NoLegend()
        plot(p)
        dev.off()
        
        ## Plot the prediction score for the class CD16 Mono
        #dir.create(paste0("Azimuth_Plot/Prediction_Score/","CD16_Mono"), recursive = T)
        #png(paste0("Azimuth_Plot/Prediction_Score/","CD16_Mono","/",sampleid,".",predicted.id,".png"), width = 800, height = 500)
        #FeaturePlot(object = seurat.list[[i]], features = "CD16 Mono", reduction = "proj.umap")
        #VlnPlot(object = seurat.list[[i]], features = "CD16 Mono", group.by = predicted.id) + NoLegend()
        #dev.off()
        
        # Plot an RNA feature
        dir.create(paste0("Azimuth_Plot/RNA/","PIK3CA"), recursive = T)
        png(paste0("Azimuth_Plot/RNA/","PIK3CA","/",sampleid,".",predicted.id,".png"), width = 1600, height = 500)
        p <- FeaturePlot(object = seurat.list[[i]], features = "PIK3CA", reduction = "proj.umap")
        p <- p + VlnPlot(object = seurat.list[[i]], features = "PIK3CA", group.by = predicted.id) + NoLegend()
        dev.off()
        
        # Plot an RNA feature
        dir.create(paste0("Azimuth_Plot/RNA/","ATM"), recursive = T)
        png(paste0("Azimuth_Plot/RNA/","ATM","/",sampleid,".",predicted.id,".png"), width = 1600, height = 500)
        p <- FeaturePlot(object = seurat.list[[i]], features = "ATM", reduction = "proj.umap")
        p <- p + VlnPlot(object = seurat.list[[i]], features = "ATM", group.by = predicted.id) + NoLegend()
        dev.off()
    }
}



rm(reference)

for (i in 1:length(seurat.list)) {
    DefaultAssay(seurat.list[[i]]) <- "SCT"
    seurat.list[[i]][["refAssay"]] <- NULL
    seurat.list[[i]][["prediction.score.celltype.l2"]] <- NULL
    seurat.list[[i]][["prediction.score.celltype.l1"]] <- NULL
}
saveRDS(seurat.list,"seurat.list_Azimuthed_slim")



### merge
seurat_merged <- merge(x = seurat.list[[1]], y = seurat.list[2:length(seurat.list)], project = "natgenet", merge.data = TRUE)
DefaultAssay(seurat_merged) <- "SCT"
saveRDS(seurat_merged,"seurat_merged")


### Manually set variable features of merged Seurat object        highly variable gene (HVG)
VariableFeatures(seurat_merged) <- seurat.features

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_merged), 10)

vdj_gene <- seurat.features_tibble %>% filter(str_detect(value, paste(out_gene, collapse = "|"))) %>% .$value
repel_gene = c(top10, vdj_gene, noname_gene)

# plot variable features with and without labels
for (i in 1:length(seurat.list)) {
    #png(paste0("SCT.VariableFeatures.",i,".png"), width=500, height=500)
    plot1 <- VariableFeaturePlot(seurat.list[[i]], assay = "SCT", selection.method="sct")
    #plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
    #print(plot2)
    #dev.off()
    plot_final <- paste0("plot_", i)
    assign(plot_final, LabelPoints(plot = plot1, points = repel_gene, repel = TRUE, xnudge = 0, ynudge = 0))
}
png("SCT.VariableFeatures.png", width=4000, height=1000)
plot_1 + plot_2 + plot_3 + plot_4 + plot_5 + plot_6 + plot_7
dev.off()

rm(seurat.list)


### PCA
seurat_merged <- RunPCA(object = seurat_merged, assay = "SCT", npcs = 50)
##  PC_ 1
##  Positive:  KRT19, KRT15, KRT5, CSTA, S100A2, FXYD3, KRT17, SFN, S100A8, TACSTD2
##             HSPB1, S100A9, KRT6A, S100A14, PERP, CCL20, CALML3, CDKN2A, KRT13, KRT8
##             SPINT2, ELF3, CD24, GSTP1, SLPI, CSTB, KRT18, CLDN7, PTTG1, MAL2
##  Negative:  IGFBP7, VIM, B2M, TMSB4X, SPARC, HLA-B, CD74, CST3, HLA-DRB1, HLA-A
##             TIMP1, HLA-DRA, COL1A1, GNG11, S100A4, CALD1, LGALS1, SERPING1, RAMP2, HLA-DPA1
##             HLA-DPB1, HLA-C, IGFBP4, IFITM2, RNASE1, SRGN, COL4A2, CRIP2, RGCC, TIMP3
##  PC_ 2
##  Positive:  CXCR4, HLA-DRA, CD52, SRGN, CD74, RGS1, CORO1A, HLA-DPB1, HLA-DPA1, HLA-DRB1
##             HLA-DQA1, LTB, TMSB4X, ALOX5AP, CD3D, HLA-DQB1, CREM, TRBC2, LAPTM5, TRAC
##             CCL5, HCST, GPR183, NKG7, DUSP2, CD7, BTG1, TYROBP, CD2, RGS2
##  Negative:  IGFBP7, SPARC, IFITM3, COL1A1, CALD1, IFI27, LGALS1, TIMP1, IGFBP4, SERPING1
##             S100A6, GNG11, VIM, CTHRC1, CST3, RAMP2, ID3, CRIP2, ADIRF, PRKCDBP
##             IGFBP3, IGFBP5, TIMP3, CAV1, CCDC80, TPM1, HSPB1, COL4A2, HTRA1, ID1
##  PC_ 3
##  Positive:  COL1A1, LGALS1, CALD1, S100A4, SPARC, CTHRC1, HMGB2, CPE, CCDC80, HIST1H4C
##             TPM1, FBLN1, MT2A, TOP2A, UBE2C, PRRX1, TRAC, CENPF, PTTG1, RPLP1
##             NUSAP1, RPS27, CCL5, MFGE8, CD3D, NBL1, CD7, CXCR4, TRBC2, SERPING1
##  Negative:  CD74, HLA-DRA, HLA-DRB1, HLA-DPA1, S100A9, HLA-DPB1, S100A8, HLA-DQA1, HLA-DQB1, LYZ
##             C15orf48, TYROBP, AIF1, RNASE1, SAT1, RAMP2, FCER1G, CST3, KRT6A, TACSTD2
##             CSTA, CSTB, SLPI, TMSB4X, TMSB10, LST1, HLA-DMA, CD14, CD68, TM4SF1
##  PC_ 4
##  Positive:  RAMP2, GNG11, RNASE1, CRIP2, HYAL2, HSPG2, TM4SF1, IGFBP7, SLC9A3R2, PRCP
##             GIMAP7, FKBP1A, TSPAN7, INSR, BCAM, FLT1, COL4A2, GIMAP4, HMGB2, EPAS1
##             STMN1, TGFBR2, PDLIM1, PTMA, ADAM15, HIST1H4C, SPRY1, ITGA6, CAV1, IFI27
##  Negative:  COL1A1, FTH1, LGALS1, S100A4, FTL, SERPING1, SERPINF1, S100A9, CST3, LGALS3BP
##             S100A8, C15orf48, SOD2, TYROBP, CTHRC1, CALD1, TIMP1, FCER1G, LYZ, PLAUR
##             APOE, CCDC80, AIF1, MT2A, CD68, CTSB, CD63, NBL1, CXCL14, FBLN1
##  PC_ 5
##  Positive:  HLA-DRA, CD74, HMGB2, HLA-DRB1, UBE2C, HLA-DPA1, TOP2A, HLA-DQA1, HLA-DPB1, PTTG1
##             LYZ, TUBA1B, CENPF, CST3, AIF1, NUSAP1, TYROBP, HIST1H4C, HLA-DQB1, LGALS1
##             FTL, BIRC5, CCNB1, FCER1G, HMGN2, STMN1, CDC20, KRT15, CDK1, AURKB
##  Negative:  S100A8, S100A9, KRT6A, ANXA1, SLPI, KRT13, TACSTD2, CXCR4, CD3D, C19orf33
##             CSTA, CCL5, TRBC2, CD7, LYPD3, CLDN7, CD2, TRAC, PDZK1IP1, CSTB
##             RGCC, SPINK5, CXCL17, CD52, SDCBP2, NEAT1, ISG20, TSC22D3, MALL, NKG7

#print(seurat_merged[["pca"]], dims = 1:5, nfeatures = 5)


png("SCT.VizDimLoadings.png", width = 800, height = 500)
VizDimLoadings(seurat_merged, dims = 1:2, reduction = "pca")
dev.off()

png("SCT.pca.Sample.png", width = 1200, height = 500)
DimPlot(seurat_merged, reduction = "pca", group.by="Sample")
dev.off()

png("SCT.pca.Group.png", width = 1200, height = 500)
DimPlot(seurat_merged, reduction = "pca", group.by="SampleGroup")
dev.off()

png("SCT.pca.png", width = 1200, height = 500)
DimPlot(seurat_merged, reduction = "pca", group.by="Sample", split.by="SampleGroup")
dev.off()

png("SCT.heatmap_dim_1_15.png", width = 800, height = 1500)
DimHeatmap(seurat_merged, dims = 1:15, cells = 5000, balanced = TRUE)
dev.off()

### SCTでは各遺伝子の分散が等しいという仮定が上手くいかないのでJackStrawは動作しないようになっている。　https://github.com/satijalab/seurat/issues/1871

### ElbowPlot
png("SCT.ElbowPlot.png", width=800, height=500)
ElbowPlot(seurat_merged, ndims = 50, reduction = "pca")
dev.off()


### Harmony         ElbowPlotをみてとりあえずPC50まで全部にした
seurat_merged <- RunHarmony(object = seurat_merged,
                                    assay.use = "SCT",
                                    reduction = "pca",
                                    dims.use = 1:50,
                                    group.by.vars = c("SampleGroup", "Sample", "SampleCase"),
                                    plot_convergence = TRUE)
##  Harmony converged after 6 iterations
##  Warning: Invalid name supplied, making object name syntactically valid. New object name is Seurat..ProjectDim.SCT.harmony; see ?make.names for more details on syntax validity

### UMAP with harmony       
seurat_merged <- RunUMAP(object = seurat_merged, assay = "SCT", reduction = "harmony", dims = 1:30)
seurat_merged <- FindNeighbors(object = seurat_merged, assay = "SCT", reduction = "harmony", dims = 1:30)
seurat_merged <- FindClusters(object = seurat_merged, resolution = c(0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2))


### Change default assay to "RNA"; normalize then generate FeaturePlots and perform differential expression analysis    https://github.com/satijalab/seurat/issues/1836
DefaultAssay(seurat_merged) <-"RNA"
seurat_merged <- NormalizeData(seurat_merged, verbose = TRUE) 

saveRDS(seurat_merged,"Harmony_round1.rds")

seurat_merged
##  An object of class Seurat
##  47342 features across 31286 samples within 2 assays
##  Active assay: RNA (25523 features, 0 variable features)
##   1 other assay present: SCT
##   3 dimensional reductions calculated: pca, harmony, umap

#   seurat_merged[["refAssay"]] <- NULL
#   seurat_merged[["prediction.score.celltype.l2"]] <- NULL
#   seurat_merged[["prediction.score.celltype.l1"]] <- NULL
#   seurat_merged[["prediction.score.celltype.l3"]] <- NULL
#   seurat_merged[["impADT"]] <- NULL
#   
#   seurat_merged
#   ##  An object of class Seurat
#   ##  49185 features across 66915 samples within 2 assays
#   ##  Active assay: RNA (27204 features, 0 variable features)
#   ##   1 other assay present: SCT
#   ##   3 dimensional reductions calculated: pca, harmony, umap

#   library(SeuratDisk)
#   SaveH5Seurat(seurat_merged, filename = "Harmony_round1.h5Seurat")
#   ##  Creating h5Seurat file for version 3.1.5.9900
#   ##  Adding counts for RNA
#   ##  Adding data for RNA
#   ##  No variable features found for RNA
#   ##  No feature-level metadata found for RNA
#   ##  Adding counts for SCT
#   ##  Adding data for SCT
#   ##  Adding scale.data for SCT
#   ##  Adding variable features for SCT
#   ##  No feature-level metadata found for SCT
#   ##  Writing out SCTModel.list for SCT
#   ##  Adding cell embeddings for pca
#   ##  Adding loadings for pca
#   ##  No projected loadings for pca
#   ##  Adding standard deviations for pca
#   ##  No JackStraw data for pca
#   ##  Adding cell embeddings for harmony
#   ##  Adding loadings for harmony
#   ##  Adding projected loadings for harmony
#   ##  No standard deviations for harmony
#   ##  No JackStraw data for harmony
#   ##  Adding cell embeddings for umap
#   ##  No loadings for umap
#   ##  No projected loadings for umap
#   ##  No standard deviations for umap
#   ##  No JackStraw data for umap
#   
#   # hfile <- Connect("Harmony_round1.h5Seurat")
#   # hfile$index()
#   # ##    Data for assay RNA★ (default assay)
#   # ##       counts      data    scale.data
#   # ##         ✔          ✔          ✖
#   # ##    Data for assay SCT
#   # ##       counts      data    scale.data
#   # ##         ✔          ✔          ✔
#   # ##    Dimensional reductions:
#   # ##               Embeddings  Loadings  Projected  JackStraw
#   # ##     harmony:      ✔          ✔          ✔          ✖
#   # ##     pca:          ✔          ✔          ✖          ✖
#   # ##     umap:         ✔          ✖          ✖          ✖
#   # ##    Graphs:
#   # ##     ─ SCT_nn
#   # ##     ─ SCT_snn
#   
#   Convert("Harmony_round1.h5Seurat", dest = "h5ad")
#   ##  Validating h5Seurat file
#   ##  Adding data from RNA as X
#   ##  Adding counts from RNA as raw
#   ##  Transfering meta.data to obs
#   ##  Adding dimensional reduction information for umap (global)

seurat_merged[["SCT"]] <- NULL

seurat_merged
##  An object of class Seurat
##  25523 features across 31286 samples within 1 assay
##  Active assay: RNA (25523 features, 0 variable features)
##   1 dimensional reduction calculated: umap

saveRDS(seurat_merged,"Harmony_RNA_round1.rds")

#   library(SeuratDisk)
#   SaveH5Seurat(seurat_merged, filename = "Harmony_RNA_round1.h5Seurat")
#   ##  Creating h5Seurat file for version 3.1.5.9900
#   ##  Adding counts for RNA
#   ##  Adding data for RNA
#   ##  No variable features found for RNA
#   ##  No feature-level metadata found for RNA
#   ##  Adding cell embeddings for umap
#   ##  No loadings for umap
#   ##  No projected loadings for umap
#   ##  No standard deviations for umap
#   ##  No JackStraw data for umap
#   
#   
#   # hfile <- Connect("Harmony_round1.h5Seurat")
#   # hfile$index()
#   # hfile$index()
#   # ##    Data for assay RNA★ (default assay)
#   # ##       counts      data    scale.data
#   # ##         ✔          ✔          ✖
#   
#   Convert("Harmony_RNA_round1.h5Seurat", dest = "h5ad")
#   ##  Validating h5Seurat file
#   ##  Warning: Cannot find assay SCT in the H5Seurat file
#   ##  Adding data from RNA as X
#   ##  Adding counts from RNA as raw
#   ##  Transfering meta.data to obs

### Azimuth
## DL results

#UMAP (Seurat Reduction RDS)
#projected.umap <- readRDS('azimuth_umap.Rds')
#object <- object[, Cells(projected.umap)]
#object[['umap.proj']] <- projected.umap

#Imputed protein (Seurat Assay RDS)
#imputed.assay <- readRDS('azimuth_impADT.Rds')
#object <- object[, Cells(imputed.assay)]
#object[['impADT']] <- imputed.assay

#Predicted cell types and scores (TSV)
#predictions <- read.delim('azimuth_pred.tsv', row.names = 1)
#object <- AddMetaData(
#    object = object,
#    metadata = predictions)

#-------------------------------#
#             graph             #
#-------------------------------#            
p_list <- list()
for(meta in c("Sample","SampleGroup")){
    seurat_merged <- SetIdent(seurat_merged, value = meta)
    p <- DimPlot(seurat_merged, reduction = "umap", group.by = meta) +
            theme(axis.text=element_blank(),axis.ticks.length = unit(0, "cm"))
    ggsave(file = paste0("UMAP_", meta,".png"), plot = p, width =12, height = 10, dpi=220)
    
    p_list[[meta]] <- p

    }

for(meta in c("SCT_snn_res.0.2","SCT_snn_res.0.4","SCT_snn_res.0.1","SCT_snn_res.0.6","SCT_snn_res.0.8","SCT_snn_res.1.0","SCT_snn_res.1.2")){
    seurat_merged <- SetIdent(seurat_merged, value = meta)
    p <- DimPlot(seurat_merged, reduction = "umap", group.by = meta) + NoLegend() +
             theme(axis.text=element_blank(),axis.ticks.length = unit(0, "cm"))
    p <- LabelClusters(plot = p, id = meta, size=4.5)
    ggsave(file = paste0("UMAP_", meta,".png"), plot = p, width =10, height = 10, dpi=220)
    
    p_list[[meta]] <- p

    }

meta="predicted.celltype.l1"
seurat_merged <- SetIdent(seurat_merged, value = meta)
p <- DimPlot(seurat_merged, reduction = "umap", group.by = meta) +
         theme(axis.text=element_blank(),axis.ticks.length = unit(0, "cm"))
p <- LabelClusters(plot = p, id = meta, size=4.5)
ggsave(file = paste0("UMAP_", meta,".pdf"), plot = p, width =14, height = 10)

#p_seurat_merged <- (p_list[["predicted.celltype.l1"]] + p_list[["predicted.celltype.l2"]]) /
#           (p_list[["SCT_snn_res.0.2"]] + p_list[["SCT_snn_res.0.4"]])
#
#ggsave(file = "UMAP_clustering_seurat_merged.png", plot = p_seurat_merged, width =15, height = 12.5, dpi=220)


# FeaturePlot
fp.plt <- FeaturePlot(seurat_merged, reduction = "umap", c("PPBP", "GP9"))
ggsave(file = "FeaturePlot_plt.png", plot = fp.plt , width = 10, height = 4.5, dpi = 220)

fp.Innate <- FeaturePlot(seurat_merged, reduction = "umap", c("CD14", "FCGR3A","CD68", "LYZ"))
ggsave(file = "FeaturePlot_Innate.png", plot = fp.Innate , width = 10, height = 9, dpi = 220)

fp.DC <- FeaturePlot(seurat_merged, reduction = "umap", c("CD1C","CLEC10A","LILRA4","IL3RA"))
ggsave(file = "FeaturePlot_DC.png", plot = fp.DC , width = 10, height = 9, dpi = 220)

fp.T_NK <- FeaturePlot(seurat_merged, reduction = "umap", c("CD3E", "CD8A", "CD4", "KLRF1"))
ggsave(file = "FeaturePlot_T_NK.png", plot = fp.T_NK , width = 10, height = 9, dpi = 220)

fp.T2 <- FeaturePlot(seurat_merged, reduction = "umap", c("FOXP3","SLC4A10","TRAV1-2", "MKI67"))
ggsave(file = "FeaturePlot_T_minor.png", plot = fp.T2 , width = 10, height = 9, dpi = 220)

fp.B <- FeaturePlot(seurat_merged, reduction = "umap", c("CD79A", "MS4A1", "MZB1", "S100A9"))           # S100A9: Doublets (Plasmablast & Monocytes) 
ggsave(file = "FeaturePlot_B.png", plot = fp.B , width = 10, height = 9, dpi = 220)

fp.HPV <- FeaturePlot(seurat_merged, reduction = "umap", c("HPV16-E1", "HPV16-E2", "HPV16-E5", "HPV16-L2", "HPV16-L1", "HPV16-E6", "HPV16-E7"))
ggsave(file = "FeaturePlot_HPV.png", plot = fp.HPV , width = 18, height = 18, dpi = 220)

fp.key <- FeaturePlot(seurat_merged, reduction = "umap", c("PIK3CA", "ATM"))
ggsave(file = "FeaturePlot_PIK3CA_ATM.png", plot = fp.key , width = 10, height = 4.5, dpi = 220)

fp.epi <- FeaturePlot(seurat_merged, reduction = "umap", c("EPCAM", "KRT8", "KRT14", "CDH1"))
ggsave(file = "FeaturePlot_epi.png", plot = fp.epi , width = 10, height = 9, dpi = 220)

fp.epiTumor <- FeaturePlot(seurat_merged, reduction = "umap", c("EPCAM", "MDK", "SOX4"))
ggsave(file = "FeaturePlot_epiTumor.png", plot = fp.epi , width = 10, height = 9, dpi = 220)

fp.epiTumor <- FeaturePlot(seurat_merged, reduction = "umap", c("EPCAM", "MDK", "SOX4"))
ggsave(file = "FeaturePlot_epiTumor.png", plot = fp.epi , width = 10, height = 9, dpi = 220)

fp.HNC_epi <- FeaturePlot(seurat_merged, reduction = "umap", c("b-actenin", "E-cadherin", "ZO1"))
ggsave(file = "FeaturePlot_HNC_epi.png", plot = fp.epi , width = 10, height = 9, dpi = 220)

vlns <- VlnPlot(seurat_merged, features=c("PIK3CA", "ATM", "HPV16-E7"), group.by="HPV_or_epithelial", pt.size=0, combine=FALSE)
p <- CombinePlots(plots=vlns, ncol=1)
ggsave(file = "VlnPlot_PIK3CA_ATM_E7.png", plot = p , width = 10, height = 20, dpi = 220)



Idents(object = seurat_merged) <- seurat_merged@meta.data$'SCT_snn_res.0.2'
deg.SCT_snn_res.0.2 <- FindAllMarkers(object = seurat_merged, only.pos=T)
write.table(deg.SCT_snn_res.0.2, "deg.SCT_snn_res.0.2.tsv", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE, )

for (i in 0:13) {
    tmp.deg <- deg.SCT_snn_res.0.2[deg.SCT_snn_res.0.2$cluster==i & deg.SCT_snn_res.0.2$p_val_adj<0.05,]
    dir.create("deg.SCT_snn_res.0.2", recursive=TRUE)
    write.table(tmp.deg, paste0("deg.SCT_snn_res.0.2/cluster", i, "_deg.tsv"), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
}


seurat_merged <- SetIdent(seurat_merged, value = 'SCT_snn_res.0.1')
deg.SCT_snn_res.0.1 <- FindAllMarkers(object = seurat_merged, only.pos=T)
write.table(deg.SCT_snn_res.0.1, "deg.SCT_snn_res.0.1.tsv", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE, )
for (i in 0:9) {
    tmp.deg <- deg.SCT_snn_res.0.1[deg.SCT_snn_res.0.1$cluster==i & deg.SCT_snn_res.0.1$p_val_adj<0.05,]
    dir.create("deg.SCT_snn_res.0.1", recursive=TRUE)
    write.table(tmp.deg, paste0("deg.SCT_snn_res.0.1/cluster", i, "_deg.tsv"), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
}

for (meta in c("SCT_snn_res.0.2","SCT_snn_res.0.4","SCT_snn_res.0.1","SCT_snn_res.0.6","SCT_snn_res.0.8","SCT_snn_res.1.2")) {
    seurat_merged <- SetIdent(seurat_merged, value = meta)
    p <- VlnPlot(seurat_merged, features = c("E1", "E2", "E6", "E7", "PIK3CA", "ATM"))
    ggsave(plot=p, file=paste0("VlnPlot_PIK3CA_ATM.",meta,".pdf"), width=14)
}


