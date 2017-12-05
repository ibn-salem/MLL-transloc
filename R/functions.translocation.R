#!/usr/bin/Rscript
#=======================================================================
#
#   Some fucntions for analysis script for chromosomal translocations with respect to Hi-C data
#
#=======================================================================




#' returns sampling weights to sample from samp with probabilities observed in obs
#'
#' @param obs  vector of observed values according to which the sampling
#' 	should be done (e.g distances observed for paralog gene pairs).
#' @aram population vector of all values in the total population from 
#' 	which one wants to sample (e.g distances of all gene pairs) 
#' @param breaks breaks used for sampling resolution (see breaks argument in hist() function).
#' @return sampling weights to sample from samp with probabilities observed in obs
#'
weightsByBin <- function(obs, population, breaks=50){

    breaksAll <- hist(c(obs, population), breaks=breaks, plot=FALSE)$breaks
        
    # calculate the number of observation for nBin equal sized bins
    hObs <- hist(obs, breaks=breaksAll, plot=FALSE)
    
    # get for each individual in the population the bin index
    binPop <- .bincode(population, breaksAll, include.lowest=TRUE)
    
    # get counts per bin in the population
    hPop <- hist(population, breaks=breaksAll, plot=FALSE)    
    
    # get the number of observed counts normalized to one as weight
    # normalize the number of observed counts by the bias observed in the population 
    weight <- hObs$counts[binPop] / hPop$counts[binPop] 

    # remove NA's, e,g, bis not observed in obs but in population. Set their probability to zero 
    weight[is.na(weight)] = 0
        
    # normalize the weights to sum up to 1
    weightNormed <- weight / sum(weight)
    
    return(weightNormed)
}

#-----------------------------------------------------------------------
# function to place breakpoints to random regions in the genome
#-----------------------------------------------------------------------
randReg <- function(inGR, N, keepStrand=FALSE, ...){
	
	require(regioneR)	# for ranodmization of genomic regions (function randomizeRegions())
	
	randRowID <- rep(1:length(inGR), each=N)
	randGR <- randomizeRegions(inGR[randRowID], genome="hg19", ...)
	
	if ( keepStrand ){
		strand(randGR) <- rep(strand(inGR), each=N)
	}

	return(randGR)
}


#-----------------------------------------------------------------------
# Sample for each partner a n times a random regions of the same size from
# the same chromosome.
#-----------------------------------------------------------------------
sampRegionsFromSameChrom <- function(pairGI, genesGR, n=10){

	# select interaction partners
	firstGR <-  anchors(pairGI, type="first")
	partnerGR <- anchors(pairGI, type="second")
	
	# sample random regions with same size from same chrom
	randGR <- randReg(partnerGR, n, allow.overlaps=TRUE, masked=FALSE, per.chromosome=TRUE)	
	
	randRegGI <- GInteractions(rep(firstGR, each=n), randGR)
	
	return(randRegGI)

}

#-----------------------------------------------------------------------
# sample random gene pairs from same chromosomes
#-----------------------------------------------------------------------
sampPairsFromSameChrom <- function(pairGI, geneGR, n=10){
	
	# select interaction partners
	firstIDs <-  anchors(pairGI, type="first", id=TRUE)
	partners <- anchors(pairGI, type="second", id=TRUE)
	
	
	chrs <- as.vector(seqnames(geneGR))

	# iterate over all partner genes
	randPartnerIDx <- sapply(partners, function(gID){
		
		# take chr and size
		chr <- chrs[gID]
#~ 		s <- size[gID]

		# filter candidates to be on same chrom and have size in the range
		candidateIDs <- which(
			chrs == chr 
		)
		
#~ 		message("INFO: Sample for gene ", geneNames[gID], " with size ", s, " from ", length(candidateIDs), " candidates...")
		
		# sample indexes of all genes on same chrom
		sample(candidateIDs, n, replace=TRUE)
	
	}, simplify=TRUE)
	
	randPairsGI <- GInteractions(rep(firstIDs, each=n), as.vector(randPartnerIDx), geneGR)
	
	return(randPairsGI)
	
}

