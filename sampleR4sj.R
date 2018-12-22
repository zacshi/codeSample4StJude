## mostly modified from an online tutorial http://www.sthda.com/english/wiki/ggbio-visualize-genomic-data
library(BayesPeak)
library(ggbio)
library(GenomicRanges)
library(Rsamtools)
library(biovizBase)
library(Homo.sapiens)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ggpubr)
library(data.table)
library(Gviz)

#my test showing covearge over all chromosomes in a bam file
bam1 <- BamFile(file="./Bioo1Sorted.bam", index="./Bioo1Sorted.bam.bai")
pdf("demoFig1.pdf", width = 15, height = 7)
autoplot(bam1, method = "estimate", main = "This is a demo to Scott")
dev.off()

#load gene symbol : GRanges, one gene/row
data(genesymbol, package = "biovizBase")

#my test showing a particular gene of interest
what <- genesymbol[c("BRAF")]
what <- range(what, ignore.strand = TRUE)
what <- keepSeqlevels(what, "chr7")
autoplot(bam1, which = what, color = "red", main = "BRAF gene")

# another gene
what <- genesymbol[c("PIK3CA")]
what <- range(what, ignore.strand = TRUE)
what <- keepSeqlevels(what, "chr3")
autoplot(bam1, which = what, color = "green", main = "PIK3CA gene")


#put probes on the graph
probesKLv1 <- rtracklayer::import.bed("./KLv1.bed")
# have to add chrX or X to the hg19sub or remove chrX from the probe bed files
WithoutX <- probesKLv1[seqnames(probesKLv1) != "chrX"]
seqlevels(WithoutX)
seqnames(WithoutX)
#newseqnames <- as.character(c(1:22))
seqlevels(WithoutX, pruning.mode="coarse") <- as.character(unique(seqnames(WithoutX)))
WithoutX
#end(WithoutX) <- start(WithoutX)
seqlengths(WithoutX) <- seqlengths(hg19sub)
KLv1probes <- ggbio() + 
    circle(WithoutX, geom = "point", color = "steelblue") + #somatic mutation
    circle(hg19sub, geom = "ideo", fill = "gray70") +#Ideogram
    circle(hg19sub, geom = "scale", size = 2) +#Scale
    circle(hg19sub, geom = "text", aes(label = seqnames), vjust = 0, size = 3)#label
KLv1probes
#anthor set of probes, my baits in ActionSeq
ActionSeqProbes <- import.bed("SHI-hg19-baits-filtration1-targetcovg-171006je.bedplus")
mybaits <- ggbio() + 
    circle(ActionSeqProbes, geom = "point", color = "steelblue") + #somatic mutation
    circle(hg19all, geom = "ideo", fill = "gray70") +#Ideogram
    circle(hg19all, geom = "scale", size = 2) +#Scale
    circle(hg19all, geom = "text", aes(label = seqnames), vjust = 0, size = 3)#label
mybaits


#Use Gviz package and a bedgraph file to visualize coverage and location
bedgraph_dt <- fread('./bio1seq.bedGraph', col.names = c('chromosome', 'start', 'end', 'value'))
# Specifiy the range to plot, "PIK3CA" gene and nearby
techr <- "chr3"
st <- 1788e5
en <-179e6
bedgraph_dt_one_chr <- bedgraph_dt[chromosome == thechr]
dtrack <- DataTrack(
    range = bedgraph_dt_one_chr,
    type = "a",
    genome = 'hg19',
    name = "Seq.Depth"
)
plotTracks(
    list(dtrack),
    from = st, to = en
)
itrack <- IdeogramTrack(
    genome = "hg19", chromosome = thechr
)
gtrack <- GenomeAxisTrack()

plotTracks(
    list(itrack, gtrack, dtrack),
    from = st, to = en
)
#add annotation
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
grtrack <- GeneRegionTrack(
    txdb,
    chromosome = thechr, start = st, end = en,
    showId = TRUE,
    name = "Gene Annotation"
)
plotTracks(
    list(itrack, gtrack, dtrack, grtrack),
    from = st, to = en
)
#graph completed. 



#last, a barplot based the function annotation of the genomic regions. 
my_GenomeReport <- read.csv("./query.output.genome_summary.csv")
str(my_GenomeReport)
myData <- as.data.frame(table(my_GenomeReport$Func.refgene))
ggbarplot(myData, x= "Var1", y = "Freq", color = "#00AFBB", fill = "#00AFBB", 
          xlab = "Locations", ylab = "Frequencies", x.text.angle=60)

