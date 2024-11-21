## https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/v0.0.4/Scds.html
## https://auroratummy.hatenablog.com/entry/2020/10/28/072113

library(dplyr)
library(tidyr)
library(tidyverse)
library(scds)
library(Seurat)
library(SingleCellExperiment)
set.seed(114)


## Read in data
counts <- Read10X(".", gene.column = 1)

# extract samples
seurat <- CreateSeuratObject(counts=counts, project="natgenet", min.cells=3, min.features=200)

table(seurat@meta.data$orig.ident)
##  
##       OP10      OP12      OP13      OP14 OP14CD45P      OP16      OP17      OP19
##       7078      6911       744      1277      1961      3330      1163      2274
##       OP20      OP33  OP33norm      OP34  OP34norm   OP35BOT    OP35LN  OP35norm
##       6187      6800       412      4006      3226       406       126       635
##        OP4  OP4CD45P       OP5       OP6  OP6CD45P       OP8       OP9    OP9GEX
##       2543      2463      2664      1474      1795      6964      4029      2502

for (id in c("OP33", "OP33norm", "OP34", "OP34norm", "OP35BOT", "OP35LN", "OP35norm")) {
    seurat <- SetIdent(seurat, value="orig.ident")
    seurat.sub <- subset(seurat, idents=id)

    if (id %in% c("OP33", "OP34", "OP35BOT", "OP35LN")) {
        seurat.sub$group <- "Tumor"
    } else {
        seurat.sub$group <- "Normal"
    }

    sce <- as.SingleCellExperiment(seurat.sub)
    
    ## Annotate doublet using binary classification based doublet scoring:
    sce = bcds(sce, retRes = TRUE, estNdbl=TRUE)
    
    ## Annotate doublet using co-expression based doublet scoring:
    try({
        sce = cxds(sce, retRes = TRUE, estNdbl=TRUE)
    })
    
    ### If cxds worked, run hybrid, otherwise use bcds annotations
    if ("cxds_score" %in% colnames(colData(sce))) {
        ## Combine both annotations into a hybrid annotation
        sce = cxds_bcds_hybrid(sce, estNdbl=TRUE)
        Doublets <- as.data.frame(cbind(rownames(colData(sce)), colData(sce)$hybrid_score, colData(sce)$hybrid_call))
    } else {
        print("this pool failed cxds so results are just the bcds calls")
        Doublets <- as.data.frame(cbind(rownames(colData(sce)), colData(sce)$bcds_score, colData(sce)$bcds_call))
    }
    
    ## Doublet scores are now available via colData:
    colnames(Doublets) <- c("Barcode","scds_score","scds_DropletType")
    Doublets$scds_DropletType <- gsub("FALSE","singlet",Doublets$scds_DropletType)
    Doublets$scds_DropletType <- gsub("TRUE","doublet",Doublets$scds_DropletType)
    
    message("writing output")
    dir.create("scds", recursive=TRUE)
    dir.create(paste0("scds/", id), recursive=TRUE)
    write_delim(Doublets, paste0("scds/", id, "/scds_doublets_singlets.tsv"), "\t")
    
    
    summary <- as.data.frame(table(Doublets$scds_DropletType))
    colnames(summary) <- c("Classification", "Droplet N")
    write_delim(summary, paste0("scds/", id, "/scds_doublet_summary.tsv"), "\t")

    # add sample and barcode to colData
    colData(sce) <- colData(sce) %>%
        as.data.frame() %>%
        mutate(Sample=orig.ident) %>%
        mutate(Barcode=rownames(colData(sce))) %>%
        mutate(SampleGroup=group) %>%
        DataFrame()
    
    
    saveRDS(sce, file = paste0("scds/", id, "/sce"))
}



for (id in c("OP4", "OP5", "OP6", "OP9", "OP13", "OP14", "OP17", "OP20")) {
    seurat <- SetIdent(seurat, value="orig.ident")
    seurat.sub <- subset(seurat, idents=id)

    if (id %in% c("OP4", "OP5", "OP6", "OP9", "OP13", "OP14", "OP17", "OP20")) {
        seurat.sub$group <- "Tumor"
    } else {
        seurat.sub$group <- "Normal"
    }

    sce <- as.SingleCellExperiment(seurat.sub)
    
    ## Annotate doublet using binary classification based doublet scoring:
    sce = bcds(sce, retRes = TRUE, estNdbl=TRUE)
    
    ## Annotate doublet using co-expression based doublet scoring:
    try({
        sce = cxds(sce, retRes = TRUE, estNdbl=TRUE)
    })
    
    ### If cxds worked, run hybrid, otherwise use bcds annotations
    if ("cxds_score" %in% colnames(colData(sce))) {
        ## Combine both annotations into a hybrid annotation
        sce = cxds_bcds_hybrid(sce, estNdbl=TRUE)
        Doublets <- as.data.frame(cbind(rownames(colData(sce)), colData(sce)$hybrid_score, colData(sce)$hybrid_call))
    } else {
        print("this pool failed cxds so results are just the bcds calls")
        Doublets <- as.data.frame(cbind(rownames(colData(sce)), colData(sce)$bcds_score, colData(sce)$bcds_call))
    }
    
    ## Doublet scores are now available via colData:
    colnames(Doublets) <- c("Barcode","scds_score","scds_DropletType")
    Doublets$scds_DropletType <- gsub("FALSE","singlet",Doublets$scds_DropletType)
    Doublets$scds_DropletType <- gsub("TRUE","doublet",Doublets$scds_DropletType)
    
    message("writing output")
    dir.create("scds", recursive=TRUE)
    dir.create(paste0("scds/", id), recursive=TRUE)
    write_delim(Doublets, paste0("scds/", id, "/scds_doublets_singlets.tsv"), "\t")
    
    
    summary <- as.data.frame(table(Doublets$scds_DropletType))
    colnames(summary) <- c("Classification", "Droplet N")
    write_delim(summary, paste0("scds/", id, "/scds_doublet_summary.tsv"), "\t")

    # add sample and barcode to colData
    colData(sce) <- colData(sce) %>%
        as.data.frame() %>%
        mutate(Sample=orig.ident) %>%
        mutate(Barcode=rownames(colData(sce))) %>%
        mutate(SampleGroup=group) %>%
        DataFrame()
    
    
    saveRDS(sce, file = paste0("scds/", id, "/sce"))
}


