#!/usr/bin/Rscript
#=======================================================================
#
#   Analysis script for chromosomal translocations with respect to Hi-C data
#
#=======================================================================

require(RColorBrewer)   # for nice colors
require(plyr)           # count() function
require(ggplot2)        # for nice plots
require(BiocParallel)   # for parallel computing
require(GenomicRanges)  # for genomic intervals and overlap calculation
require(InteractionSet) # for HI-C data structures
require(TxDb.Hsapiens.UCSC.hg19.knownGene) # for seqinfo object
require(chromint) 	# to parse Hi-C from Rao et al. see: https://github.com/ibn-salem/chromint
require(stringr) # for str_split_fixed()
require(rtracklayer)  # for export() function to parse .bigWig files



#-----------------------------------------------------------------------
# Rao et al. 2014 Hi-C sample
#-----------------------------------------------------------------------

TICDB_FILES=c(
	gT3="data/TICDdb/gT3.txt",
	mRNAfus3="data/TICDdb/mRNAfus3.txt",
	unknown3="data/TICDdb/unknown3.txt"
)

# CTCF motifs:
CTCF_motif_file <- "../loop_prediction/data/factorbook/factorbookMotifPos.bed.CTCF"
CTCF_motif_file_ms <- "data/hg19_MA0139.1.mrkv1.ms"


CELL <- "GM12878"
HIC_RESOLUTION <- 50*10^3 # 50kb
#~ HIC_RESOLUTION <- 5*10^3 # 5kb
HIC_DATA_DIR <- "data/Rao2014"

#~ HIC_DATA_DIR_MURO="/project/jgu-cbdm/andradeLab/download/flat/databases/uncompressed/muro/ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GM12878_combined"
#~ HIC_DATA_DIR_MURO="/project/jgu-cbdm/andradeLab/download/flat/databases/uncompressed/muro/ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl"

# number of bins for distance density estimation as parameter to adjust sampling resolution
RANDOM_SEED=13521
USE_LOCAL_HIC=TRUE


VERSION="v02"

outPrefix = paste0("results/TICdb_analysis_", VERSION, "/transloc.", VERSION) 
dataPrefix = "results/TICdb_analysis_data/transloc."
# create directory, if not exist
dir.create(dirname(outPrefix), recursive=TRUE, showWarnings = FALSE)
dir.create(dirname(dataPrefix), recursive=TRUE, showWarnings = FALSE)

#-----------------------------------------------------------------------
# Options for parallel computation
#-----------------------------------------------------------------------

# use all available cores but generate random number streams on each worker
multicorParam <- MulticoreParam(RNGseed=RANDOM_SEED)   
# set options
register(multicorParam)  
# bpparam() # to print current options

#-----------------------------------------------------------------------
# load some custom functions
#-----------------------------------------------------------------------
#~ source("R/parseHiC.R")
source("R/functions.translocation.R")

txdb_hg19 <- TxDb.Hsapiens.UCSC.hg19.knownGene
seqInfo <- seqinfo(txdb_hg19)

#-----------------------------------------------------------------------
# Load Hi-C data from Rao et al. 2014
#-----------------------------------------------------------------------
if ( !USE_LOCAL_HIC) {

	# parse normalized Hi-C map from Rao et al. 
#~ 	gi = parseRaoHiCtoGI(CELL, HIC_RESOLUTION, HIC_DATA_DIR, seqInfo)
	gi = parseRaoHiCtoGI(CELL, HIC_RESOLUTION, HIC_DATA_DIR, seqInfo)
	
	# annotate GI with genomic distance:
	gi$dist <- pairdist(gi)
	gi$cis <- intrachr(gi)
	

#~ 	HiClistRaw = parseRaoHiC(CELL, HIC_RESOLUTION, HIC_DATA_DIR, seqInfo, normStr=NULL)	
	# parse expected contacts by distance for each chromosome
#~ 	expectedHiCList <- parseRaoExpected(CELL, HIC_RESOLUTION, HIC_DATA_DIR)
	
	# save data for faster loading next time
	save(gi, file=paste0(dataPrefix, ".Hi-C.", CELL, ".", HIC_RESOLUTION, ".gi.RData"))
			
}else{
	load(paste0(dataPrefix, ".Hi-C.", CELL, ".", HIC_RESOLUTION, ".gi.RData"))
}

#-----------------------------------------------------------------------
# Parse TICdb data from Enrique
#-----------------------------------------------------------------------
INFILE <- TICDB_FILES[1]

df <- read.table(INFILE, header=TRUE, sep="\t", comment.char="")

g1 <- str_split_fixed(df$gene1, ":", 4)
g2 <- str_split_fixed(df$gene2, ":", 4)

loc1 <- as.numeric(g1[,2])
loc2 <- as.numeric(g2[,2])

gr1 <- GRanges(g1[,1], IRanges(loc1, loc1), strand=g1[,3], gene=g1[,4], seqinfo=seqInfo)
gr2 <- GRanges(g2[,1], IRanges(loc2, loc2), strand=g2[,3], gene=g2[,4], seqinfo=seqInfo)

translocGI <- GInteractions(gr1, gr2)


#-----------------------------------------------------------------------
# parse CTCF motifs
#-----------------------------------------------------------------------
# parse CTCF sites as anchor GR
motifENCODE <- import(CTCF_motif_file, format="BED", seqinfo=seqInfo)

# flip strand (to match motif of Rao et al. 2014)
strand(motifENCODE) <- ifelse(strand(motifENCODE) == "+", "-", "+")
motifENCODE <- sort(motifENCODE)

# parse matrix-scan motif:
motifDF <- read.table(CTCF_motif_file_ms, 
	sep="\t", header=TRUE, comment.char=";")

motifGR <- GRanges(
	motifDF[,1],
	IRanges(motifDF$start, motifDF$end),
	strand=ifelse(motifDF$strand == "D", "+", "-"),
	motifDF[,c("sequence", "weight", "Pval", "ln_Pval", "sig")], 
	seqinfo=seqInfo
	)
#require(BSgenome.Hsapiens.UCSC.hg19)
#getSeq(Hsapiens, motifGR[1:10])
motifGR <- sort(motifGR)

#~ motifGR <- motifGR[motifGR$Pval <= 10^-4]

#~ ENCODE.uniq <- setdiff(motifENCODE, motifGR)
#~ matrix.uniq <- setdiff(motifGR, motifENCODE)
#~ ovl <- subsetByOverlaps(motifENCODE, motifGR)	

#~ c(
#~ length(ENCODE.uniq),
#~ length(ovl),
#~ length(matrix.uniq)
#~ )
#-----------------------------------------------------------------------
# Compare motifs from ENCODE with motifs from RSAT matrix-scan
#-----------------------------------------------------------------------



#-----------------------------------------------------------------------
# Analyse motif distribution
#-----------------------------------------------------------------------

#WIN_SIZE <- 10^4
WIN_SIZES <- c(10^5, 10^4, 10^3)
N_BINS <- 20

bp <- regions(translocGI)

randBp <- randReg(bp, 100, keepStrand=TRUE)

allBp <- c(bp, randBp)
allBp$group <- rep(c("transloc", "random"), c(length(bp), length(randBp)))
n <- length(allBp)


for( WIN_SIZE in WIN_SIZES){
	
	bpWin <- resize(allBp, WIN_SIZE, fix="center")
	bpBins <- unlist(tile(bpWin, N_BINS))
	
	countsPos <- countOverlaps(bpBins, motifGR[strand(motifGR) == "+"], ignore.strand=TRUE)
	countsNeg <- countOverlaps(bpBins, motifGR[strand(motifGR) == "-"], ignore.strand=TRUE)
	
	#countMat <- matrix(counts, ncol=20, byrow=TRUE)
	
	
	countDF <- data.frame(
		count=c(countsPos, countsNeg),
		bin=rep(1:N_BINS, n*2),
		bpID=rep(rep(1:n, each=N_BINS), 2),
		group=rep(rep(allBp$group, each=N_BINS), 2),
		geneStrand=rep(rep(strand(bpWin), each=N_BINS), 2),
		motifDirection=rep(c("+", "-"), each=n*N_BINS)
	)
	
	countDF$distBin <- factor(
		countDF$bin, 
		1:N_BINS, 
		seq(-WIN_SIZE/2, WIN_SIZE/2-WIN_SIZE/N_BINS, length.out=N_BINS)
		)
	
	#geneDF <- ddply(gpDF, .(group, chr, gene, pairID), summarize, nInt=mean(intPairs), total=sum(raw), mean=mean(raw))
	plotDF <- ddply(countDF, .(group, geneStrand, motifDirection, distBin), summarize, 
		mean=mean(count),
		sd=sd(count)
	)

	#-------------------------------------------------------------------
	# barplot
	#-------------------------------------------------------------------
	p <- ggplot(plotDF, aes(x=distBin, y=mean, fill=motifDirection)) + 
		geom_bar(stat="identity", position="dodge",  color="black") +
		facet_grid(geneStrand~group) + 
	    theme_bw() + 
	    theme(text = element_text(size=20), 
			axis.text.x=element_text(angle = 45, hjust = 1), 
			legend.position="bottom") +
	    ylab("CTCF motifs per bin") + xlab("Distance to breakpoint [bp]")
		
	ggsave(p, file=paste0(outPrefix, ".win_size_", WIN_SIZE, ".CTCF_motifs_at_breakpoints.byGeneStrand.barplot.pdf"), w=7,h=7)

	#-------------------------------------------------------------------
	# line
	#-------------------------------------------------------------------
	p <- ggplot(plotDF, aes(x=as.numeric(as.character(distBin)), y=mean, color=motifDirection)) + 
		geom_line(size=2) + geom_point(size=3) +
		facet_grid(geneStrand~group) + 
	    theme_bw() + 
	    theme(text = element_text(size=20), 
			axis.text.x=element_text(angle = 45, hjust = 1), 
			legend.position="bottom") +
	    ylab("CTCF motifs per bin") + xlab("Distance to breakpoint [bp]")
		
	ggsave(p, file=paste0(outPrefix, ".win_size_", WIN_SIZE, ".CTCF_motifs_at_breakpoints.byGeneStrand.lines.pdf"), w=7,h=7)


}

#-----------------------------------------------------------------------
# FINISH HERE!
stop("FINISH HERE!")
#-----------------------------------------------------------------------


#-----------------------------------------------------------------------
# Get Ensembl genes
#-----------------------------------------------------------------------

ensemblGRCh37 = useMart(host="grch37.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")

#-------------------------------------------------------------------
# get all genes with annotation:
#-------------------------------------------------------------------
geneAttributes = c("entrezgene", "external_gene_name", "chromosome_name", "transcript_start", "start_position", "end_position", "strand", "status", "gene_biotype")
geneFilters=c("chromosome_name", "biotype", "status")

# read "normal" human chromosome names (without fixes and patches) and take only KNOWN protein coding genes into account
geneValues=list(chromosome_name=c(1:22, "X", "Y"), biotype="protein_coding", status="KNOWN")

allGenes = getBM(attributes=geneAttributes, mart=ensemblGRCh37, filters=geneFilters, values=geneValues)


# unique gene entry gene id:
genes = allGenes[!duplicated(allGenes$entrezgene) & !is.na(allGenes$entrezgene),]
#~ rownames(genes) <- genes$entrezgene

# make GRanges object for gene 
geneGR = GRanges(
	paste0("chr", genes$chromosome_name),
	IRanges(genes$start_position, genes$end_position),
	strand = ifelse(genes$strand == 1, '+', '-'), 
	genes[,c("entrezgene", "external_gene_name", "transcript_start")],
	seqinfo=seqInfo
	)	
geneGR <- sort(geneGR)

plainGR <- geneGR
mcols(plainGR) <- NULL


qID <- 4297
tIDs <- c(4300, 4301, 4299, 4298)

subGenes <- geneGR[match(c(qID, tIDs), geneGR$entrezgene)]

genePairGI <- GInteractions(
	rep(match(qID, geneGR$entrezgene), length(tIDs)), 
	match(tIDs, geneGR$entrezgene), 
	plainGR
	)


#-----------------------------------------------------------------------
# Analyse Hi-C interaction between gene pairs
#-----------------------------------------------------------------------

# get random gene pairs from same chromosome
randPairs <- sampPairsFromSameChrom(genePairGI, plainGR, 1000)
randRegPairs <- sampRegionsFromSameChrom(genePairGI, plainGR, 1000)

#~ df <- DataFrame()
#~ dim(df) <- c(length(randRegPairs), 0)
#~ mcols(randRegPairs) <- DataFrame()
#~ mcols(randRegPairs) <- NULL

# drop metadata of regions:

#~ mcols(regions(genePairGI)) <- NULL
#~ mcols(regions(randPairs)) <- NULL
#~ mcols(regions(randRegPairs)) <- NULL

# make one query gene pair GI
allPairs <- c(genePairGI, randPairs, randRegPairs)
allPairs$group <- rep(c("cases", "sampled_genes", "random_regions"), 
	c(length(genePairGI), length(randPairs), length(randRegPairs)))

# annotate partners with gene names
allPairs$gene <- geneGR$external_gene_name[nearest(anchors(allPairs, type="second"), geneGR)]
allPairs$chr <- seqnames(anchors(allPairs, type="second"))
allPairs$pairID <- paste0("pID", 1:length(allPairs))


# get interactions between gene pairs
hit <- findOverlaps(gi, allPairs, ignore.strand=TRUE)

# get number of unique interaction pairs between regions
cOv <- data.frame(table(factor(subjectHits(hit), 1:length(allPairs))))$Freq

allPairs$intPairs <- cOv

# expand data.frame to have a row for all interaction pairs from gi
gp <- allPairs[c(subjectHits(hit), which(cOv == 0))]

# annotate each interaction with the raw interaction count
gp$raw <-  c(gi$raw[queryHits(hit)], rep(0, sum(cOv == 0)))

gpDF <- data.frame(mcols(gp))

#~ # add column indicating case vs. control
#~ gpDF$testGroup <- factor(gpDF$group == "cases", c(TRUE, FALSE), c("Transloc partners", "Random Control"))

#=======================================================================
# plot
#=======================================================================
COL = brewer.pal(8, "Set1")[c(3,5)]
#~ COL_GROUP = brewer.pal(8, "Dark2")[c(4,6,7)]
COL_GROUP = brewer.pal(8, "Set3")[c(4,5,1)]
COL_CHR = brewer.pal(8, "Set3")

#-----------------------------------------------------------------------
# interaction strength
#-----------------------------------------------------------------------
p <- ggplot(gpDF, aes(x=gene, fill=chr, y=raw)) + 
	geom_bar(stat="identity", color="black") +
	facet_grid(.~group, scales="free_x") + 
    scale_fill_manual(values=COL_CHR, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=10), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") +
    ylab("Hi-C interaction frequency")
ggsave(p, file=paste0(outPrefix, ".", CELL, ".", HIC_RESOLUTION, ".MLL_partners_and_rand.intFreq_byGroup.barplot.pdf"), w=14,h=7)

p <- ggplot(gpDF, aes(x=gene, fill=group, y=raw)) + 
	geom_bar(stat="identity", color="black") +
	facet_grid(.~chr, scales="free_x") +
    theme_bw() + theme(text = element_text(size=10), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom")+
    ylab("Hi-C interaction frequency")
ggsave(p, file=paste0(outPrefix, ".", CELL, ".", HIC_RESOLUTION, ".MLL_partners_and_rand.intFreq_byChr.barplot.pdf"), w=7,h=7)
  
#~ ggplot(gpDF, aes(x=group, color=group, y=interactions)) + 
#~ 	geom_jitter(position = position_jitter(w = .25, h = 0))
	

# boxplot
#~ pvalDF <- ddply(gpDF, .(chr), summarize, p=wilcox.test(raw ~ group)$p.value)
#~ pSamp <- wilcox.test(raw ~ group, data=subset(gpDF, group %in% c("cases", "sampled_genes")))$p.value
#~ pRand <- wilcox.test(raw ~ group, data=subset(gpDF, group %in% c("cases", "random_regions")))$p.value

pvalDF <- ddply(subset(gpDF, group!="cases"), .(group), summarize, p=wilcox.test(subset(gpDF, group=="cases")$raw, raw)$p.value)


p <- ggplot(gpDF, aes(x=group, color=group, y=raw)) + 
	geom_boxplot() +
#~ 	geom_jitter(position = position_jitter(w = .25, h = 0)) +
	scale_fill_manual(values=COL_CHR, guide_legend(title = "")) + 
    theme_bw() + theme(text = element_text(size=10), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom")+
    ylab("Hi-C interaction frequency") + 
	geom_text(aes(label=paste0("p=", signif(p,2)), x=c(1.5, 2.5), y=max(gpDF$raw)), data=pvalDF, size=5) 
	
ggsave(p, file=paste0(outPrefix, ".", CELL, ".", HIC_RESOLUTION, ".MLL_partners_and_rand.interaction_byGroup.boxplot.pdf"), w=3.5,h=7)


p <- ggplot(gpDF, aes(x=group, color=group, y=raw)) + 
	geom_boxplot(color="gray") +
	geom_jitter(position = position_jitter(w = .25, h = 0)) +
	facet_grid(.~chr) + 
    scale_fill_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=10), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom")+
    ylab("Hi-C interaction frequency")
ggsave(p, file=paste0(outPrefix, ".", CELL, ".", HIC_RESOLUTION, ".MLL_partners_and_rand.interaction_byChr.boxplot.pdf"), w=7,h=7)
  
#=======================================================================
# summarize by gene
geneDF <- ddply(gpDF, .(group, chr, gene, pairID), summarize, nInt=mean(intPairs), total=sum(raw), mean=mean(raw))


# add column indicating case vs. control
#~ geneDF$testGroup <- factor(geneDF$group == "cases", c(TRUE, FALSE), c("Transloc partners", "Random Control"))
#~ geneDFctl <- subset(geneDF, group!="cases")
#~ pvalDF <- ddply(geneDF, .(chr, group), summarize, p=wilcox.test(total ~ testGroup)$p.value)


#~ pvalDF <- ddply(subset(geneDF, group!="cases"), .(group), summarize, p=wilcox.test(subset(geneDF, group=="cases")$total, total)$p.value)

subSamp <- subset(geneDF, group %in% c("cases", "sampled_genes"))
subRand <- subset(geneDF, group %in% c("cases", "random_regions"))

#~ pvalDFSamp <- ddply(subSamp, .(chr), summarize, p=wilcox.test(total ~ group)$p.value)
#~ pvalDFRand <- ddply(subRand, .(chr), summarize, p=wilcox.test(total ~ group)$p.value)

empiricalP <- function(x, exp){
	mean(p < exp)
}

pvalDFSamp <- ddply(subSamp, .(chr), summarize, p=wilcox.test(total ~ group)$p.value)
pvalDFRand <- ddply(subRand, .(chr), summarize, p=wilcox.test(total ~ group)$p.value)


pvalDF <- cbind(
	rbind(pvalDFSamp, pvalDFRand),
	group = rep(c("sampled_genes", "random_regions"), c(nrow(pvalDFSamp), nrow(pvalDFRand)))
	)

nDF <- ddply(geneDF, .(chr, group), summarize, 
	n=length(pairID)
)
	
#~ # compute p-values per chromosome
#~ pvalDF <- ddply(subset(geneDF, group!="cases"), .(chr, group), summarize, p=wilcox.test(total ~ subset(geneDF, group=="cases")$total)$p.value)


p <- ggplot(geneDF, aes(x=group, color=group, y=total)) + 
	geom_boxplot() +
#~ 	geom_jitter(position = position_jitter(w = .25, h = 0)) +
#~ 	geom_violin(aes(fill=group)) +
	facet_grid(.~chr) + 
	scale_fill_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") +
	geom_text(data=pvalDF, aes(label=paste0("p=", signif(p,2)), y=1.1*max(geneDF$total)), color="black", size=3) +
	geom_text(data=nDF, aes(label=paste0("n=", n), y=-0.5), color="black", size=3) +
	geom_text(data=subset(geneDF, group=="cases"), aes(label=gene), size=3, vjust = -1) +
    ylab("Hi-C interactions with MLL") + xlab("")
    
    
ggsave(p, file=paste0(outPrefix, ".", CELL, ".", HIC_RESOLUTION, ".MLL_partners_and_rand.perGene.total_byChromAndGroup.boxplot.pdf"), w=7,h=7)

#~ p <- p + geom_text_repel(aes(label=gene))
#~ ggsave(p, file=paste0(outPrefix, ".", CELL, ".", HIC_RESOLUTION, ".MLL_partners_and_rand.perGene.total_byChromAndGroup.boxplot.withGenes.pdf"), w=7,h=7)



#=======================================================================
# OLD:
#=======================================================================

# make one query gene pair GI
allPairs <- c(genePairGI, randPairs)
allPairs$group <- rep(c("cases", "sampled_genes"), 
	c(length(genePairGI), length(randPairs)))
	
allPairs$gene <- geneGR$external_gene_name[anchors(allPairs, type="second",  id=TRUE)]
allPairs$chr <- seqnames(anchors(allPairs, type="second"))


# get interactions between gene pairs
hit <- findOverlaps(gi, allPairs, ignore.strand=TRUE)
#~ hitDF <- data.frame(hit)
#~ hitDF$raw <- gi$raw[queryHits(hit)]

#~ cOv <- countOverlaps(allPairs, gi)
cOv <- data.frame(table(factor(subjectHits(hit), 1:length(allPairs))))$Freq

allPairs$interactions <- cOv

# expand data.frame
gp <- allPairs[c(subjectHits(hit), which(cOv == 0))]

gp$raw <-  c(gi$raw[queryHits(hit)], rep(0, sum(cOv == 0)))


gpDF <- data.frame(mcols(gp))

