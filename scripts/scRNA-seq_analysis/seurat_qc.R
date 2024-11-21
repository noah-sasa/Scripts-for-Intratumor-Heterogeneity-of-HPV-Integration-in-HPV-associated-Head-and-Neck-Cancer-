library(data.table)
library(Seurat)
library(miQC)
library(SeuratWrappers)
library(flexmix)
library(SingleCellExperiment)
library(Matrix)
library(stringr)


for (id in c("OP33", "OP33norm", "OP34", "OP34norm", "OP35BOT", "OP35LN", "OP35norm")) {
    sce = readRDS(paste0("../scds/", id, "/sce"))
    seurat <- paste0("seurat_", id)
    assign(seurat, CreateSeuratObject(counts=counts(sce)))
    # add metadata to seurat (変数名がiの影響を受けるときは下記のようにする必要がある)
    seurat.tmp <- get(seurat)
    seurat.tmp$Sample <- sce$Sample
    seurat.tmp$Barcode <- sce$Barcode
    #seurat.tmp$SampleGroup <- sce$SampleGroup      # sceのOP34がミスってNormal
    if (id %in% c("OP33", "OP34", "OP35BOT", "OP35LN")) {
        seurat.tmp$SampleGroup <- "Tumor"
    } else {
        seurat.tmp$SampleGroup <- "Normal"
    }
    if (id %in% c("OP33", "OP33norm")) {
        seurat.tmp$SampleCase <- "OP33"
    } else if (id %in% c("OP34", "OP34norm")) {
        seurat.tmp$SampleCase <- "OP34"
    } else {
        seurat.tmp$SampleCase <- "OP35"
    }
    seurat.tmp$Scds_call <- sce$hybrid_call
    assign(seurat, seurat.tmp)
}

for (id in c("OP4", "OP5", "OP6", "OP9", "OP13", "OP14", "OP17", "OP20")) {
    sce = readRDS(paste0("../scds/", id, "/sce"))
    seurat <- paste0("seurat_", id)
    assign(seurat, CreateSeuratObject(counts=counts(sce)))
    # add metadata to seurat (変数名がiの影響を受けるときは下記のようにする必要がある)
    seurat.tmp <- get(seurat)
    seurat.tmp$Sample <- sce$Sample
    seurat.tmp$Barcode <- sce$Barcode
    seurat.tmp$SampleGroup <- sce$SampleGroup
    seurat.tmp$SampleCase <- sce$SampleGroup
    seurat.tmp$Scds_call <- sce$hybrid_call
    assign(seurat, seurat.tmp)
}


head(seurat_OP33)
##                          orig.ident nCount_RNA nFeature_RNA Sample
##  OP33_AAACCTGAGATGGCGT-1       OP33       3800         1397   OP33
##  OP33_AAACCTGAGCTACCTA-1       OP33       6141         1542   OP33
##  OP33_AAACCTGAGTTTGCGT-1       OP33       9033         2877   OP33
##  OP33_AAACCTGCAAACTGCT-1       OP33       9677         3121   OP33
##  OP33_AAACCTGCACCGCTAG-1       OP33       5213         1617   OP33
##  OP33_AAACCTGCACGCGAAA-1       OP33       3895         1365   OP33
##  OP33_AAACCTGCAGTGACAG-1       OP33       8777         2485   OP33
##  OP33_AAACCTGCATGGGACA-1       OP33       6316         1926   OP33
##  OP33_AAACCTGGTCTGGAGA-1       OP33       5341         1950   OP33
##  OP33_AAACCTGGTGCACTTA-1       OP33       4669         1477   OP33
##                                          Barcode SampleGroup SampleCase
##  OP33_AAACCTGAGATGGCGT-1 OP33_AAACCTGAGATGGCGT-1       Tumor       OP33
##  OP33_AAACCTGAGCTACCTA-1 OP33_AAACCTGAGCTACCTA-1       Tumor       OP33
##  OP33_AAACCTGAGTTTGCGT-1 OP33_AAACCTGAGTTTGCGT-1       Tumor       OP33
##  OP33_AAACCTGCAAACTGCT-1 OP33_AAACCTGCAAACTGCT-1       Tumor       OP33
##  OP33_AAACCTGCACCGCTAG-1 OP33_AAACCTGCACCGCTAG-1       Tumor       OP33
##  OP33_AAACCTGCACGCGAAA-1 OP33_AAACCTGCACGCGAAA-1       Tumor       OP33
##  OP33_AAACCTGCAGTGACAG-1 OP33_AAACCTGCAGTGACAG-1       Tumor       OP33
##  OP33_AAACCTGCATGGGACA-1 OP33_AAACCTGCATGGGACA-1       Tumor       OP33
##  OP33_AAACCTGGTCTGGAGA-1 OP33_AAACCTGGTCTGGAGA-1       Tumor       OP33
##  OP33_AAACCTGGTGCACTTA-1 OP33_AAACCTGGTGCACTTA-1       Tumor       OP33
##                          Scds_call
##  OP33_AAACCTGAGATGGCGT-1     FALSE
##  OP33_AAACCTGAGCTACCTA-1     FALSE
##  OP33_AAACCTGAGTTTGCGT-1     FALSE
##  OP33_AAACCTGCAAACTGCT-1     FALSE
##  OP33_AAACCTGCACCGCTAG-1     FALSE
##  OP33_AAACCTGCACGCGAAA-1     FALSE
##  OP33_AAACCTGCAGTGACAG-1     FALSE
##  OP33_AAACCTGCATGGGACA-1      TRUE
##  OP33_AAACCTGGTCTGGAGA-1     FALSE
##  OP33_AAACCTGGTGCACTTA-1     FALSE



seurat.combined <- merge(seurat_OP33,
    y = c(seurat_OP33norm, seurat_OP34, seurat_OP34norm, seurat_OP35BOT, seurat_OP35LN, seurat_OP35norm,
        seurat_OP4, seurat_OP5, seurat_OP6, seurat_OP9, seurat_OP13, seurat_OP14, seurat_OP17, seurat_OP20
        ),
    add.cell.ids = c(seurat_OP33$Sample[1], seurat_OP33norm$Sample[1], seurat_OP34$Sample[1], seurat_OP34norm$Sample[1],
                        seurat_OP35BOT$Sample[1], seurat_OP35LN$Sample[1], seurat_OP35norm$Sample[1],
                        seurat_OP4$Sample[1], seurat_OP5$Sample[1], seurat_OP6$Sample[1], seurat_OP9$Sample[1],
                        seurat_OP13$Sample[1], seurat_OP14$Sample[1], seurat_OP17$Sample[1], seurat_OP20$Sample[1]
                        ),
    project = "natgenet")

### HPV positive or negative
hpv_expression <- FetchData(seurat.combined, vars=c("HPV16-E1", "HPV16-E2", "HPV16-E5", "HPV16-L2", "HPV16-L1", "HPV16-E6", "HPV16-E7"))
hpv_expression$HPV <- hpv_expression$`HPV16-E1` + hpv_expression$`HPV16-E2` + hpv_expression$`HPV16-E5` + hpv_expression$`HPV16-L2` + hpv_expression$`HPV16-L1` + hpv_expression$`HPV16-E6` + hpv_expression$`HPV16-E7`
seurat.combined$HPV_group <- ifelse(hpv_expression$HPV > 0, "positive", "negative")

library(tidyverse)
table(seurat.combined@meta.data %>% select(Sample, HPV_group))
##            HPV_group
##  Sample     negative positive
##    OP13          711       33
##    OP14         1217       60
##    OP17         1046      117
##    OP20         3882     2305
##    OP33         3774     3026
##    OP33norm      412        0
##    OP34         3227      779
##    OP34norm     3213       13
##    OP35BOT       265      141
##    OP35LN         45       81
##    OP35norm      635        0
##    OP4          1332     1211
##    OP5          2357      307
##    OP6          1138      336
##    OP9          3381      648

seurat.combined
##  An object of class Seurat
##  28074 features across 35692 samples within 1 assay
##  Active assay: RNA (28074 features, 0 variable features)



# QC https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_01_qc.html
#   https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
# MT
seurat.combined <- PercentageFeatureSet(object = seurat.combined, pattern = "^MT-", col.name = "percent.mt")
# RIBO
seurat.combined <- PercentageFeatureSet(object = seurat.combined, pattern = "^RP[SL]", col.name = "percent.ribo")
# Hb
seurat.combined <- PercentageFeatureSet(object = seurat.combined, pattern = "^HB[^(P)]", col.name = "percent.hb")
# Plat
seurat.combined <- PercentageFeatureSet(object = seurat.combined, pattern = "PECAM1|PF4", col.name = "percent.plat")


dim(seurat.combined)
##  [1]   28074 35692

rowSums(seurat.combined)
## -> HPV由来リードそれなりにあり

png("QC_metrics.violin.perSample.png", width=3000, height=1000)
VlnPlot(seurat.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb", "percent.plat"), group.by="Sample", log = FALSE, pt.size=0)
dev.off()

png("QC_metrics.violin.perSample.log.png", width=3000, height=1000)
VlnPlot(seurat.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb", "percent.plat"), group.by="Sample", log = TRUE, pt.size=0)
dev.off()

png("QC_metrics.violin.Group.png", width=1000, height=1000)
VlnPlot(seurat.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb", "percent.plat"), ncol=3, group.by="Sample", split.by="SampleGroup", log = TRUE, pt.size=0)
dev.off()

png("QC_metrics.violin.Doublet.png", width=1000, height=1000)
VlnPlot(seurat.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb", "percent.plat"), ncol=3, group.by="Sample", split.by="Scds_call", log = TRUE, pt.size=0)
dev.off()

png("QC_metrics.violin.Doublet.notlog.png", width=1000, height=1000)
VlnPlot(seurat.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb", "percent.plat"), ncol=3, group.by="Sample", split.by="Scds_call", log = FALSE, pt.size=0)
dev.off()


library(patchwork)
png("QC_metrics.scatter.nCount_RNA.perSample.5.png", width=2000, height=500)
plot1 <- FeatureScatter(seurat.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by="Sample", raster=FALSE)
plot2 <- FeatureScatter(seurat.combined, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by="Sample", raster=FALSE)
plot3 <- FeatureScatter(seurat.combined, feature1 = "nCount_RNA", feature2 = "percent.ribo", group.by="Sample", raster=FALSE)
plot4 <- FeatureScatter(seurat.combined, feature1 = "nCount_RNA", feature2 = "percent.hb", group.by="Sample", raster=FALSE)
plot5 <- FeatureScatter(seurat.combined, feature1 = "nCount_RNA", feature2 = "percent.plat", group.by="Sample", raster=FALSE)
plot1 + plot2 + plot3 + plot4 + plot5 + plot_layout( guides = "collect" )
dev.off()

png("QC_metrics.scatter.nFeature_RNA.perSample.5.png", width=2000, height=500)
plot1 <- FeatureScatter(seurat.combined, feature1 = "nFeature_RNA", feature2 = "nCount_RNA", group.by="Sample", raster=FALSE)
plot2 <- FeatureScatter(seurat.combined, feature1 = "nFeature_RNA", feature2 = "percent.mt", group.by="Sample", raster=FALSE)
plot3 <- FeatureScatter(seurat.combined, feature1 = "nFeature_RNA", feature2 = "percent.ribo", group.by="Sample", raster=FALSE)
plot4 <- FeatureScatter(seurat.combined, feature1 = "nFeature_RNA", feature2 = "percent.hb", group.by="Sample", raster=FALSE)
plot5 <- FeatureScatter(seurat.combined, feature1 = "nFeature_RNA", feature2 = "percent.plat", group.by="Sample", raster=FALSE)
plot1 + plot2 + plot3 + plot4 + plot5 + plot_layout( guides = "collect" )
dev.off()

png("QC_metrics.scatter.percent.mt.perSample.5.png", width=2000, height=500)
plot1 <- FeatureScatter(seurat.combined, feature1 = "percent.mt", feature2 = "nFeature_RNA", group.by="Sample", raster=FALSE)
plot2 <- FeatureScatter(seurat.combined, feature1 = "percent.mt", feature2 = "nCount_RNA", group.by="Sample", raster=FALSE)
plot3 <- FeatureScatter(seurat.combined, feature1 = "percent.mt", feature2 = "percent.ribo", group.by="Sample", raster=FALSE)
plot4 <- FeatureScatter(seurat.combined, feature1 = "percent.mt", feature2 = "percent.hb", group.by="Sample", raster=FALSE)
plot5 <- FeatureScatter(seurat.combined, feature1 = "percent.mt", feature2 = "percent.plat", group.by="Sample", raster=FALSE)
plot1 + plot2 + plot3 + plot4 + plot5 + plot_layout( guides = "collect" )
dev.off()

png("QC_metrics.scatter.percent.ribo.perSample.5.png", width=2000, height=500)
plot1 <- FeatureScatter(seurat.combined, feature1 = "percent.ribo", feature2 = "nFeature_RNA", group.by="Sample", raster=FALSE)
plot2 <- FeatureScatter(seurat.combined, feature1 = "percent.ribo", feature2 = "nCount_RNA", group.by="Sample", raster=FALSE)
plot3 <- FeatureScatter(seurat.combined, feature1 = "percent.ribo", feature2 = "percent.mt", group.by="Sample", raster=FALSE)
plot4 <- FeatureScatter(seurat.combined, feature1 = "percent.ribo", feature2 = "percent.hb", group.by="Sample", raster=FALSE)
plot5 <- FeatureScatter(seurat.combined, feature1 = "percent.ribo", feature2 = "percent.plat", group.by="Sample", raster=FALSE)
plot1 + plot2 + plot3 + plot4 + plot5 + plot_layout( guides = "collect" )
dev.off()

png("QC_metrics.scatter.percent.hb.perSample.5.png", width=2000, height=500)
plot1 <- FeatureScatter(seurat.combined, feature1 = "percent.hb", feature2 = "nFeature_RNA", group.by="Sample", raster=FALSE)
plot2 <- FeatureScatter(seurat.combined, feature1 = "percent.hb", feature2 = "nCount_RNA", group.by="Sample", raster=FALSE)
plot3 <- FeatureScatter(seurat.combined, feature1 = "percent.hb", feature2 = "percent.mt", group.by="Sample", raster=FALSE)
plot4 <- FeatureScatter(seurat.combined, feature1 = "percent.hb", feature2 = "percent.ribo", group.by="Sample", raster=FALSE)
plot5 <- FeatureScatter(seurat.combined, feature1 = "percent.hb", feature2 = "percent.plat", group.by="Sample", raster=FALSE)
plot1 + plot2 + plot3 + plot4 + plot5 + plot_layout( guides = "collect" )
dev.off()

png("QC_metrics.scatter.percent.plat.perSample.5.png", width=2000, height=500)
plot1 <- FeatureScatter(seurat.combined, feature1 = "percent.plat", feature2 = "nFeature_RNA", group.by="Sample", raster=FALSE)
plot2 <- FeatureScatter(seurat.combined, feature1 = "percent.plat", feature2 = "nCount_RNA", group.by="Sample", raster=FALSE)
plot3 <- FeatureScatter(seurat.combined, feature1 = "percent.plat", feature2 = "percent.mt", group.by="Sample", raster=FALSE)
plot4 <- FeatureScatter(seurat.combined, feature1 = "percent.plat", feature2 = "percent.ribo", group.by="Sample", raster=FALSE)
plot5 <- FeatureScatter(seurat.combined, feature1 = "percent.plat", feature2 = "percent.hb", group.by="Sample", raster=FALSE)
plot1 + plot2 + plot3 + plot4 + plot5 + plot_layout( guides = "collect" )
dev.off()

#   # filter out low quality cells using miQC
#   seurat.combined <- RunMiQC(seurat.combined, percent.mt = "percent.mt", nFeature_RNA = "nFeature_RNA", posterior.cutoff = 0.75, model.slot = "flexmix_model")
#   #   Warning: Adding a command log without an assay associated with it
#   
#   pdf("PlotMiQC.pdf")
#   p <- PlotMiQC(seurat.combined, percent.mt="percent.mt", nFeature_RNA="nFeature_RNA", model.slot="flexmix_model", color.by="miQC.probability")
#   p
#   dev.off()
#   
#   dim(seurat.combined)
#   ##  [1]  36792 100368
#   
#   seurat.combined.miQCed <- subset(seurat.combined, miQC.keep == "keep")
#   
#   dim(seurat.combined.miQCed)
#   ##  [1] 36792 81615
#   
#   
#   
#   png("miQCed.QC_metrics.violin.perSample.png", width=3000, height=1000)
#   VlnPlot(seurat.combined.miQCed, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb", "percent.plat"), group.by="Sample", log = FALSE, pt.size=0.1)
#   dev.off()
#   
#   png("miQCed.QC_metrics.violin.perSample.log.png", width=3000, height=1000)
#   VlnPlot(seurat.combined.miQCed, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb", "percent.plat"), group.by="Sample", log = TRUE, pt.size=0.1)
#   dev.off()
#   
#   png("miQCed.QC_metrics.violin.Institute.png", width=1000, height=1000)
#   VlnPlot(seurat.combined.miQCed, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb", "percent.plat"), ncol=3, group.by="Institute", split.by="SampleGroup", log = TRUE)
#   dev.off()
#   
#   png("miQCed.QC_metrics.violin.Institute.Doublet.png", width=1000, height=1000)
#   VlnPlot(seurat.combined.miQCed, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb", "percent.plat"), ncol=3, group.by="Institute", split.by="Scds_call", log = TRUE)
#   dev.off()
#   
#   png("miQCed.QC_metrics.violin.Institute.Doublet.notlog.png", width=1000, height=1000)
#   VlnPlot(seurat.combined.miQCed, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb", "percent.plat"), ncol=3, group.by="Institute", split.by="Scds_call", log = FALSE)
#   dev.off()



seurat.combined.Scdsed <- subset(seurat.combined, Scds_call == "FALSE")

dim(seurat.combined)
##  [1]  28074 35692
dim(seurat.combined.Scdsed)
##  [1]  28074 32331

png("Scdsed.QC_metrics.violin.perSample.png", width=3000, height=1000)
VlnPlot(seurat.combined.Scdsed, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb", "percent.plat"), group.by="Sample", log = FALSE, pt.size=0)
dev.off()

png("Scdsed.QC_metrics.violin.perSample.log.png", width=3000, height=1000)
VlnPlot(seurat.combined.Scdsed, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb", "percent.plat"), group.by="Sample", log = TRUE, pt.size=0)
dev.off()

png("Scdsed.QC_metrics.violin.Group.png", width=1000, height=1000)
VlnPlot(seurat.combined.Scdsed, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb", "percent.plat"), ncol=3, group.by="Sample", split.by="SampleGroup", log = TRUE, pt.size=0)
dev.off()



selected_gene <- rownames(seurat.combined.Scdsed)[Matrix::rowSums(seurat.combined.Scdsed) >= 3]

seurat.combined.Scdsed.filt <- subset(seurat.combined.Scdsed, features=selected_gene, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA < 60000 & percent.ribo > 5 & percent.hb < 10)

dim(seurat.combined.Scdsed.filt)
##  [1]  25523 31286

png("Scdsed.filt.QC_metrics.violin.perSample.png", width=3000, height=1000)
VlnPlot(seurat.combined.Scdsed.filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb", "percent.plat"), group.by="Sample", log = FALSE, pt.size=0)
dev.off()

png("Scdsed.filt.QC_metrics.violin.perSample.log.png", width=3000, height=1000)
VlnPlot(seurat.combined.Scdsed.filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb", "percent.plat"), group.by="Sample", log = TRUE, pt.size=0)
dev.off()

png("Scdsed.filt.QC_metrics.violin.Group.png", width=1000, height=1000)
VlnPlot(seurat.combined.Scdsed.filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb", "percent.plat"), ncol=3, group.by="Sample", split.by="SampleGroup", log = TRUE, pt.size=0)
dev.off()

png("Scdsed.filt.QC_metrics.violin.Group.notlog.png", width=1000, height=1000)
VlnPlot(seurat.combined.Scdsed.filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb", "percent.plat"), ncol=3, group.by="Sample", split.by="SampleGroup", log = FALSE, pt.size=0)
dev.off()


saveRDS(seurat.combined, file="seurat.combined")
saveRDS(seurat.combined.Scdsed.filt, file="seurat.combined.Scdsed.filt")
