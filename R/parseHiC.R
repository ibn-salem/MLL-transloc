#!/usr/bin/Rscript
#=======================================================================
#
#   Provides functionality to parse Hi-C data 
#   from Rao et al 2014 and Capture Hi-C data from 
#   Mifsud et al. 2015
#
#=======================================================================


#-----------------------------------------------------------------------
# Creates a GRanges object with bins at a given resolution on given chrom
#-----------------------------------------------------------------------
getBinGR <- function(chr, resolution, seqInfo){
    options("scipen"=999)
    chrLen = seqlengths(seqInfo)[chr]
    starts = seq(1, chrLen, resolution)
    ends = starts+resolution-1
    ends = ifelse(ends < chrLen, ends, chrLen)
    n = length(starts)
    gr = GRanges(rep(chr, n), 
        IRanges(starts, ends),
        names=as.character(starts-1),
        seqinfo=seqInfo
    )
    options("scipen"=0)
    return(gr)
}
# getBinGR(chr, resolution, seqInfo)


#-----------------------------------------------------------------------
# Creates a GRanges object with bins at a given resolution on given chrom
#-----------------------------------------------------------------------
getChrToBinOffset <- function(chromosomes, resolution, seqInfo){
	
	# get length of all chromosomes
	chrLen = seqlengths(seqInfo)[chromosomes]
	
	# get number of bins for each chromosomes
	nBins <- ceiling(chrLen / resolution)
	
	# get cummulative sum for offset
	binOffset <- cumsum(c(0, nBins))[1:length(chromosomes)]
	
	# fix names
	names(binOffset) <- chromosomes

    return(binOffset)
}
	

#-----------------------------------------------------------------------
# parse Hi-C map for all chromosomes (only intra-chrom) at given resolution 
# TODO: implement KR normalization
#-----------------------------------------------------------------------
parseRaoHiCtoGI <- function(cell, resolution, dirPrefix, seqInfo, normalizeByExpected=FALSE, mapQualStr="MAPQGE30", normStr="KR"){
    
    # get all the proper file paths
    # An example path is: 
    # K562/100kb_resolution_intrachromosomal/chr1/MAPQGE30/chr1_100kb.RAWobserved
    # K562_interchromosomal/100kb_resolution_interchromosomal/chr1_chr2/MAPQGE30/chr1_2_100kb.RAWobserved
    
    # map resolution in bp to kb/mb substring:
    res2str = c(
        "1000"="1kb",
        "5000"="5kb",
        "10000"="10kb",
        "25000"="25kb",
        "50000"="50kb",
        "100000"="100kb",
        "250000"="250kb",
        "500000"="500kb",
        "1000000"="1mb"
        )

    scipenDefault <- options()$scipen
    options("scipen"=999)
    resStr = res2str[as.character(resolution)]
    options("scipen"=scipenDefault)
        
    # build path to directory with chromosome subdirectories
    chromDir = file.path(dirPrefix, cell, paste0(resStr, "_resolution_intrachromosomal"))

    # get all available chromosome names:
    chromosomes = list.dirs(path =chromDir , full.names = FALSE, recursive = FALSE)
	# sort 
	chromOrder <- match(seqnames(seqInfo), chromosomes)
	chromosomes <- chromosomes[chromOrder[!is.na(chromOrder)]]

	# create bins for given resolution and chromosome size as GRanges object:
    chromBinGR <- lapply(chromosomes, getBinGR, resolution, seqInfo)
    names(chromBinGR) <- chromosomes
	
	# combine to binGR for entire genome
	binGR <-unlist(GRangesList(chromBinGR))
	
	# get bin offsets
	binOffset <- getChrToBinOffset(chromosomes, resolution, seqInfo)

	# help function to get from bin position to index in GRanges
	posToIdx <- function(pos, resolution){
		(pos / resolution) + 1
	}

	#--------------------------------------------------------------------
	# parse intra-chromosomal interactions
	#--------------------------------------------------------------------
	
	cisDFlist <- bplapply(chromosomes, function(chr){
        
        message(paste("INFO: Begin to parse data for chromosome", chr, "..."))
		
		# get file path
        rawInteractionFile = file.path(chromDir, chr, mapQualStr, paste0(chr, "_", resStr, ".RAWobserved"))


		# parse interaction file
		intData <- data.frame(fread(rawInteractionFile))
		
		# get anchor bins as GRanges object
		anchorIdx1 <- binOffset[chr] + (intData[,1] / resolution) +1
		anchorIdx2 <- binOffset[chr] + (intData[,2] / resolution) +1
		
		# build data.frame with indexes and score
		chrDF <- data.frame(idx1=anchorIdx1, idx2=anchorIdx2, raw=intData[,3])
		
		return(chrDF)
        
    })
    Sys.sleep(1) # hack to fix problems with bplapply on MOGON
    
	#--------------------------------------------------------------------
	# parse intra-chromosomal interactions
	#--------------------------------------------------------------------
	
	
	# get all chromosome pair combination (ordered)
	chromPairs <- t(combn(chromosomes, 2))
	
	system.time(
	transDFlist <- bpmapply(function(chr1, chr2){
		message(paste("INFO: Begin to parse interactions between chromosome", chr1, "and", chr2 , "..."))

		# get file path
		# Example file:	#K562_interchromosomal/100kb_resolution_interchromosomal/chr1_chr2/MAPQGE30/chr1_2_100kb.RAWobserved

        rawInteractionFile = file.path(dirPrefix, paste0(cell, "_interchromosomal"), paste0(resStr, "_resolution_interchromosomal"), paste0(chr1, "_", chr2), mapQualStr, paste0(chr1, "_", gsub("chr", "", chr2), "_", resStr, ".RAWobserved"))

		# parse interaction file
		intData <- data.frame(fread(rawInteractionFile))
		
		# get anchor bins as GRanges object
		anchorIdx1 <- binOffset[chr1] + (intData[,1] / resolution) +1
		anchorIdx2 <- binOffset[chr2] + (intData[,2] / resolution) +1
		
		# build data.frame with indexes and score
		chrPairDF <- data.frame(idx1=anchorIdx1, idx2=anchorIdx2, raw=intData[,3])
		
		return(chrPairDF)
		
	}, chromPairs[,1], chromPairs[,2], SIMPLIFY=FALSE)
    )
    Sys.sleep(1) # hack to fix problems with bplapply on MOGON
	
	
    
    # combine all data.frames into a single one
    intDF <- rbindlist(c(cisDFlist, transDFlist))

	# build GI object
	GI <- GInteractions(intDF[,idx1], intDF[,idx2], binGR, raw=intDF[,raw], mode="strict")
	
	# get size of object
	#format(object.size(GI), units = "auto")
    
	message(paste("INFO: Finished parsing of interactions. GI object has size:", format(object.size(GI), units = "auto")))

    return(GI)
    
}
#~ rawMatrixFile
#~ cell=CELL
#~ resolution=HIC_RESOLUTION
#~ dirPrefix=HIC_DATA_DIR
#~ normalizeByExpected=FALSE
#~ mapQualStr="MAPQGE30"
#~ normStr="KR"
#~ expectedSuffix=".KRexpected"


