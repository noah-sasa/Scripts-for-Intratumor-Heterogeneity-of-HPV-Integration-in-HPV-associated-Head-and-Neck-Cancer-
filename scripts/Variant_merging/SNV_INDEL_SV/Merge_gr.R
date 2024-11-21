suppressPackageStartupMessages(require(StructuralVariantAnnotation))
suppressPackageStartupMessages(require(VariantAnnotation))

#SNVも含めても結局indelのみのgrになるっぽい。
m2_vcf <- VariantAnnotation::readVcf("../M2/M2.PASS.Norm.noMNP.indel.vcf", "hg38")
# Some SV callers don't report QUAL so wee need to use a proxy
VariantAnnotation::fixed(m2_vcf)$QUAL <- info(m2_vcf)$GERMQ
m2_gr <- breakpointRanges(m2_vcf)
m2_gr$caller <- "Mutect2"

s2_vcf <- VariantAnnotation::readVcf("../S2/S2.INDEL.PASS.Norm.name-changed.vcf", "hg38")
# Some SV callers don't report QUAL so wee need to use a proxy
VariantAnnotation::fixed(s2_vcf)$QUAL <- info(s2_vcf)$QSI
s2_gr <- breakpointRanges(s2_vcf)
s2_gr$caller <- "Strelka2"

gridss_vcf <- VariantAnnotation::readVcf("../SV_VCF/GRIDSS.ohb.PASS.vcf", "hg38")
gridss_bpgr <- breakpointRanges(gridss_vcf)
### single-breakend
gridss_begr <- breakendRanges(gridss_vcf)
### Pseudo breakend ... partner=itself
gridss_begr$partner <- names(gridss_begr)
### Merge
gridss_gr <- c(gridss_bpgr, gridss_begr)
gridss_gr$caller <- "GRIDSS"

manta_vcf <- VariantAnnotation::readVcf("../SV_VCF/Manta.PASS.vcf", "hg38")
# Some SV callers don't report QUAL so wee need to use a proxy
VariantAnnotation::fixed(manta_vcf)$QUAL <- info(manta_vcf)$SOMATICSCORE
manta_bpgr <- breakpointRanges(manta_vcf)
### single-breakend
manta_begr <- breakendRanges(manta_vcf)
### Pseudo breakend ... partner=itself
manta_begr$partner <- names(manta_begr)
### Merge
manta_gr <- c(manta_bpgr, manta_begr)
manta_gr$caller <- "Manta"

delly_vcf <- VariantAnnotation::readVcf("../SV_VCF/Delly.PASS.vcf", "hg38")
delly_bpgr <- breakpointRanges(delly_vcf)
### single-breakend
delly_begr <- breakendRanges(delly_vcf)
### Pseudo breakend ... partner=itself
delly_begr$partner <- names(delly_begr)
### Merge
delly_gr <- c(delly_bpgr, delly_begr)
delly_gr$caller <- "delly"


#INDEL match条件1
m2_gr$s2.indel1_matches <- countBreakpointOverlaps(m2_gr, s2_gr,
  maxgap=0,
  sizemargin=NULL,
  restrictMarginToSizeMultiple=NULL,
  countOnlyBest=TRUE)

m2_gr$gridss.indel1_matches <- countBreakpointOverlaps(m2_gr, gridss_gr,
  maxgap=0,
  sizemargin=NULL,
  restrictMarginToSizeMultiple=NULL,
  countOnlyBest=TRUE)

m2_gr$manta.indel1_matches <- countBreakpointOverlaps(m2_gr, manta_gr,
  maxgap=0,
  sizemargin=NULL,
  restrictMarginToSizeMultiple=NULL,
  countOnlyBest=TRUE)

m2_gr$delly.indel1_matches <- countBreakpointOverlaps(m2_gr, delly_gr,
  maxgap=0,
  sizemargin=NULL,
  restrictMarginToSizeMultiple=NULL,
  countOnlyBest=TRUE)


s2_gr$m2.indel1_matches <- countBreakpointOverlaps(s2_gr, m2_gr,
  maxgap=0,
  sizemargin=NULL,
  restrictMarginToSizeMultiple=NULL,
  countOnlyBest=TRUE)

s2_gr$gridss.indel1_matches <- countBreakpointOverlaps(s2_gr, gridss_gr,
  maxgap=0,
  sizemargin=NULL,
  restrictMarginToSizeMultiple=NULL,
  countOnlyBest=TRUE)

s2_gr$manta.indel1_matches <- countBreakpointOverlaps(s2_gr, manta_gr,
  maxgap=0,
  sizemargin=NULL,
  restrictMarginToSizeMultiple=NULL,
  countOnlyBest=TRUE)

s2_gr$delly.indel1_matches <- countBreakpointOverlaps(s2_gr, delly_gr,
  maxgap=0,
  sizemargin=NULL,
  restrictMarginToSizeMultiple=NULL,
  countOnlyBest=TRUE)


gridss_gr$m2.indel1_matches <- countBreakpointOverlaps(gridss_gr, m2_gr,
  maxgap=0,
  sizemargin=NULL,
  restrictMarginToSizeMultiple=NULL,
  countOnlyBest=TRUE)

gridss_gr$s2.indel1_matches <- countBreakpointOverlaps(gridss_gr, s2_gr,
  maxgap=0,
  sizemargin=NULL,
  restrictMarginToSizeMultiple=NULL,
  countOnlyBest=TRUE)

gridss_gr$manta.indel1_matches <- countBreakpointOverlaps(gridss_gr, manta_gr,
  maxgap=0,
  sizemargin=NULL,
  restrictMarginToSizeMultiple=NULL,
  countOnlyBest=TRUE)

gridss_gr$delly.indel1_matches <- countBreakpointOverlaps(gridss_gr, delly_gr,
  maxgap=0,
  sizemargin=NULL,
  restrictMarginToSizeMultiple=NULL,
  countOnlyBest=TRUE)


manta_gr$m2.indel1_matches <- countBreakpointOverlaps(manta_gr, m2_gr,
  maxgap=0,
  sizemargin=NULL,
  restrictMarginToSizeMultiple=NULL,
  countOnlyBest=TRUE)

manta_gr$s2.indel1_matches <- countBreakpointOverlaps(manta_gr, s2_gr,
  maxgap=0,
  sizemargin=NULL,
  restrictMarginToSizeMultiple=NULL,
  countOnlyBest=TRUE)

manta_gr$gridss.indel1_matches <- countBreakpointOverlaps(manta_gr, gridss_gr,
  maxgap=0,
  sizemargin=NULL,
  restrictMarginToSizeMultiple=NULL,
  countOnlyBest=TRUE)

manta_gr$delly.indel1_matches <- countBreakpointOverlaps(manta_gr, delly_gr,
  maxgap=0,
  sizemargin=NULL,
  restrictMarginToSizeMultiple=NULL,
  countOnlyBest=TRUE)


delly_gr$m2.indel1_matches <- countBreakpointOverlaps(delly_gr, m2_gr,
  maxgap=0,
  sizemargin=NULL,
  restrictMarginToSizeMultiple=NULL,
  countOnlyBest=TRUE)

delly_gr$s2.indel1_matches <- countBreakpointOverlaps(delly_gr, s2_gr,
  maxgap=0,
  sizemargin=NULL,
  restrictMarginToSizeMultiple=NULL,
  countOnlyBest=TRUE)

delly_gr$gridss.indel1_matches <- countBreakpointOverlaps(delly_gr, gridss_gr,
  maxgap=0,
  sizemargin=NULL,
  restrictMarginToSizeMultiple=NULL,
  countOnlyBest=TRUE)

delly_gr$manta.indel1_matches <- countBreakpointOverlaps(delly_gr, manta_gr,
  maxgap=0,
  sizemargin=NULL,
  restrictMarginToSizeMultiple=NULL,
  countOnlyBest=TRUE)



#INDEL match条件2
m2_gr$s2.indel2_matches <- countBreakpointOverlaps(m2_gr, s2_gr,
  maxgap=5,
  sizemargin=0.25,
  restrictMarginToSizeMultiple=0.05,
  countOnlyBest=TRUE)

m2_gr$gridss.indel2_matches <- countBreakpointOverlaps(m2_gr, gridss_gr,
  maxgap=5,
  sizemargin=0.25,
  restrictMarginToSizeMultiple=0.05,
  countOnlyBest=TRUE)

m2_gr$manta.indel2_matches <- countBreakpointOverlaps(m2_gr, manta_gr,
  maxgap=5,
  sizemargin=0.25,
  restrictMarginToSizeMultiple=0.05,
  countOnlyBest=TRUE)

m2_gr$delly.indel2_matches <- countBreakpointOverlaps(m2_gr, delly_gr,
  maxgap=5,
  sizemargin=0.25,
  restrictMarginToSizeMultiple=0.05,
  countOnlyBest=TRUE)


s2_gr$m2.indel2_matches <- countBreakpointOverlaps(s2_gr, m2_gr,
  maxgap=5,
  sizemargin=0.25,
  restrictMarginToSizeMultiple=0.05,
  countOnlyBest=TRUE)

s2_gr$gridss.indel2_matches <- countBreakpointOverlaps(s2_gr, gridss_gr,
  maxgap=5,
  sizemargin=0.25,
  restrictMarginToSizeMultiple=0.05,
  countOnlyBest=TRUE)

s2_gr$manta.indel2_matches <- countBreakpointOverlaps(s2_gr, manta_gr,
  maxgap=5,
  sizemargin=0.25,
  restrictMarginToSizeMultiple=0.05,
  countOnlyBest=TRUE)

s2_gr$delly.indel2_matches <- countBreakpointOverlaps(s2_gr, delly_gr,
  maxgap=5,
  sizemargin=0.25,
  restrictMarginToSizeMultiple=0.05,
  countOnlyBest=TRUE)


gridss_gr$m2.indel2_matches <- countBreakpointOverlaps(gridss_gr, m2_gr,
  maxgap=5,
  sizemargin=0.25,
  restrictMarginToSizeMultiple=0.05,
  countOnlyBest=TRUE)

gridss_gr$s2.indel2_matches <- countBreakpointOverlaps(gridss_gr, s2_gr,
  maxgap=5,
  sizemargin=0.25,
  restrictMarginToSizeMultiple=0.05,
  countOnlyBest=TRUE)

gridss_gr$manta.indel2_matches <- countBreakpointOverlaps(gridss_gr, manta_gr,
  maxgap=5,
  sizemargin=0.25,
  restrictMarginToSizeMultiple=0.05,
  countOnlyBest=TRUE)

gridss_gr$delly.indel2_matches <- countBreakpointOverlaps(gridss_gr, delly_gr,
  maxgap=5,
  sizemargin=0.25,
  restrictMarginToSizeMultiple=0.05,
  countOnlyBest=TRUE)


manta_gr$m2.indel2_matches <- countBreakpointOverlaps(manta_gr, m2_gr,
  maxgap=5,
  sizemargin=0.25,
  restrictMarginToSizeMultiple=0.05,
  countOnlyBest=TRUE)

manta_gr$s2.indel2_matches <- countBreakpointOverlaps(manta_gr, s2_gr,
  maxgap=5,
  sizemargin=0.25,
  restrictMarginToSizeMultiple=0.05,
  countOnlyBest=TRUE)

manta_gr$gridss.indel2_matches <- countBreakpointOverlaps(manta_gr, gridss_gr,
  maxgap=5,
  sizemargin=0.25,
  restrictMarginToSizeMultiple=0.05,
  countOnlyBest=TRUE)

manta_gr$delly.indel2_matches <- countBreakpointOverlaps(manta_gr, delly_gr,
  maxgap=5,
  sizemargin=0.25,
  restrictMarginToSizeMultiple=0.05,
  countOnlyBest=TRUE)


delly_gr$m2.indel2_matches <- countBreakpointOverlaps(delly_gr, m2_gr,
  maxgap=5,
  sizemargin=0.25,
  restrictMarginToSizeMultiple=0.05,
  countOnlyBest=TRUE)

delly_gr$s2.indel2_matches <- countBreakpointOverlaps(delly_gr, s2_gr,
  maxgap=5,
  sizemargin=0.25,
  restrictMarginToSizeMultiple=0.05,
  countOnlyBest=TRUE)

delly_gr$gridss.indel2_matches <- countBreakpointOverlaps(delly_gr, gridss_gr,
  maxgap=5,
  sizemargin=0.25,
  restrictMarginToSizeMultiple=0.05,
  countOnlyBest=TRUE)

delly_gr$manta.indel2_matches <- countBreakpointOverlaps(delly_gr, manta_gr,
  maxgap=5,
  sizemargin=0.25,
  restrictMarginToSizeMultiple=0.05,
  countOnlyBest=TRUE)



# SV match条件
m2_gr$s2_matches <- countBreakpointOverlaps(m2_gr, s2_gr,
  maxgap=100,
  sizemargin=0.25,
  restrictMarginToSizeMultiple=0.5,
  countOnlyBest=TRUE)

m2_gr$gridss_matches <- countBreakpointOverlaps(m2_gr, gridss_gr,
  maxgap=100,
  sizemargin=0.25,
  restrictMarginToSizeMultiple=0.5,
  countOnlyBest=TRUE)

m2_gr$manta_matches <- countBreakpointOverlaps(m2_gr, manta_gr,
  maxgap=100,
  sizemargin=0.25,
  restrictMarginToSizeMultiple=0.5,
  countOnlyBest=TRUE)

m2_gr$delly_matches <- countBreakpointOverlaps(m2_gr, delly_gr,
  maxgap=100,
  sizemargin=0.25,
  restrictMarginToSizeMultiple=0.5,
  countOnlyBest=TRUE)


s2_gr$m2_matches <- countBreakpointOverlaps(s2_gr, m2_gr,
  maxgap=100,
  sizemargin=0.25,
  restrictMarginToSizeMultiple=0.5,
  countOnlyBest=TRUE)

s2_gr$gridss_matches <- countBreakpointOverlaps(s2_gr, gridss_gr,
  maxgap=100,
  sizemargin=0.25,
  restrictMarginToSizeMultiple=0.5,
  countOnlyBest=TRUE)

s2_gr$manta_matches <- countBreakpointOverlaps(s2_gr, manta_gr,
  maxgap=100,
  sizemargin=0.25,
  restrictMarginToSizeMultiple=0.5,
  countOnlyBest=TRUE)

s2_gr$delly_matches <- countBreakpointOverlaps(s2_gr, delly_gr,
  maxgap=100,
  sizemargin=0.25,
  restrictMarginToSizeMultiple=0.5,
  countOnlyBest=TRUE)


gridss_gr$m2_matches <- countBreakpointOverlaps(gridss_gr, m2_gr,
  maxgap=100,
  sizemargin=0.25,
  restrictMarginToSizeMultiple=0.5,
  countOnlyBest=TRUE)

gridss_gr$s2_matches <- countBreakpointOverlaps(gridss_gr, s2_gr,
  maxgap=100,
  sizemargin=0.25,
  restrictMarginToSizeMultiple=0.5,
  countOnlyBest=TRUE)

gridss_gr$manta_matches <- countBreakpointOverlaps(gridss_gr, manta_gr,
  maxgap=100,
  sizemargin=0.25,
  restrictMarginToSizeMultiple=0.5,
  countOnlyBest=TRUE)

gridss_gr$delly_matches <- countBreakpointOverlaps(gridss_gr, delly_gr,
  maxgap=100,
  sizemargin=0.25,
  restrictMarginToSizeMultiple=0.5,
  countOnlyBest=TRUE)


manta_gr$m2_matches <- countBreakpointOverlaps(manta_gr, m2_gr,
  maxgap=100,
  sizemargin=0.25,
  restrictMarginToSizeMultiple=0.5,
  countOnlyBest=TRUE)

manta_gr$s2_matches <- countBreakpointOverlaps(manta_gr, s2_gr,
  maxgap=100,
  sizemargin=0.25,
  restrictMarginToSizeMultiple=0.5,
  countOnlyBest=TRUE)

manta_gr$gridss_matches <- countBreakpointOverlaps(manta_gr, gridss_gr,
  maxgap=100,
  sizemargin=0.25,
  restrictMarginToSizeMultiple=0.5,
  countOnlyBest=TRUE)

manta_gr$delly_matches <- countBreakpointOverlaps(manta_gr, delly_gr,
  maxgap=100,
  sizemargin=0.25,
  restrictMarginToSizeMultiple=0.5,
  countOnlyBest=TRUE)


delly_gr$m2_matches <- countBreakpointOverlaps(delly_gr, m2_gr,
  maxgap=100,
  sizemargin=0.25,
  restrictMarginToSizeMultiple=0.5,
  countOnlyBest=TRUE)

delly_gr$s2_matches <- countBreakpointOverlaps(delly_gr, s2_gr,
  maxgap=100,
  sizemargin=0.25,
  restrictMarginToSizeMultiple=0.5,
  countOnlyBest=TRUE)

delly_gr$gridss_matches <- countBreakpointOverlaps(delly_gr, gridss_gr,
  maxgap=100,
  sizemargin=0.25,
  restrictMarginToSizeMultiple=0.5,
  countOnlyBest=TRUE)

delly_gr$manta_matches <- countBreakpointOverlaps(delly_gr, manta_gr,
  maxgap=100,
  sizemargin=0.25,
  restrictMarginToSizeMultiple=0.5,
  countOnlyBest=TRUE)






### single-breakendはassembly精度によってサイズがcaller間でかなり異なる場合がある。size制限なし。代わりにmaxgap=20
gridss_begr$manta_matches <- countBreakpointOverlaps(gridss_begr, manta_begr,
  maxgap=20,
  sizemargin=NULL,
  restrictMarginToSizeMultiple=NULL,
  ignore.strand = FALSE,
  countOnlyBest=TRUE)
gridss_begr$delly_matches <- countBreakpointOverlaps(gridss_begr, delly_begr,
  maxgap=20,
  sizemargin=NULL,
  restrictMarginToSizeMultiple=NULL,
  ignore.strand = FALSE,
  countOnlyBest=TRUE)

manta_begr$gridss_matches <- countBreakpointOverlaps(manta_begr, gridss_begr,
  maxgap=20,
  sizemargin=NULL,
  restrictMarginToSizeMultiple=NULL,
  ignore.strand = FALSE,
  countOnlyBest=TRUE)
manta_begr$delly_matches <- countBreakpointOverlaps(manta_begr, delly_begr,
  maxgap=20,
  sizemargin=NULL,
  restrictMarginToSizeMultiple=NULL,
  ignore.strand = FALSE,
  countOnlyBest=TRUE)

delly_begr$gridss_matches <- countBreakpointOverlaps(delly_begr, gridss_begr,
  maxgap=20,
  sizemargin=NULL,
  restrictMarginToSizeMultiple=NULL,
  ignore.strand = FALSE,
  countOnlyBest=TRUE)
delly_begr$manta_matches <- countBreakpointOverlaps(delly_begr, manta_begr,
  maxgap=20,
  sizemargin=NULL,
  restrictMarginToSizeMultiple=NULL,
  ignore.strand = FALSE,
  countOnlyBest=TRUE)


options(showTailLines=Inf)

write.table(m2_gr, file="m2_gr.tsv", sep="\t", quote=F)
write.table(s2_gr, file="s2_gr.tsv", sep="\t", quote=F)
write.table(gridss_gr, file="gridss_gr.tsv", sep="\t", quote=F)
write.table(gridss_begr, file="gridss_single-brealend_gr.tsv", sep="\t", quote=F)
write.table(manta_gr, file="manta_gr.tsv", sep="\t", quote=F)
write.table(manta_begr, file="manta_single-breakend_gr.tsv", sep="\t", quote=F)
write.table(delly_gr, file="delly_gr.tsv", sep="\t", quote=F)
write.table(delly_begr, file="delly_single-breakend_gr.tsv", sep="\t", quote=F)
