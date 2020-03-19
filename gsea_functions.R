gseaScoresbpown  <- function(geneList, geneNames.perm, 
	collectionOfGeneSets, exponent = 1, nPermutations = 1000,no.cores.gsea) {
	##check arguments
	## paraCheck("genelist",geneList)
	## paraCheck("gsc",collectionOfGeneSets)
	## paraCheck("exponent",exponent)
	## paraCheck("nPermutations",nPermutations)
	if(!is.matrix(geneNames.perm))
		stop("'geneNames.perm' should be a matrix!\n")
	if(ncol(geneNames.perm)!=(nPermutations+1))
		stop("The No of columns of 'geneNames.perm' should be equal to 'nPermutations'!\n")	
	##local function for computation of gsea scores with a single core
	gseaScoresBatchLocal <- function(geneList, geneNames.perm, geneSet, 
		exponent, nPermutations) {	
		geneList.names <- names(geneList)

		##The geneSet should be a subset of the gene universe, i.e. we 
		##keep only those element of the gene set that appear in the 
		##geneList		
		geneSet <- intersect(geneList.names, geneSet)
        ##Compute the size of the gene set and of the genelist
        nh <- length(geneSet)
        N <- length(geneList)
		
        ES <- rep(0, nPermutations+1)
		Phit <- matrix(0, nrow = N, ncol = nPermutations+1)
		Pmiss <- Phit
		runningES <- NULL
		
		if(nh > N)
			stop("Gene Set is larger than Gene List")

		hits <- matrix(FALSE, nrow = N, ncol = nPermutations+1) 	
		hits[which(!is.na(match(geneNames.perm, geneSet)))] <- TRUE	
		hits <- matrix(hits, ncol = nPermutations+1, byrow = FALSE)		
		if(sum(hits[,1]) > 0) {
			junk <- sapply(1:(nPermutations+1), function(i) 
				Phit[which(hits[, i]), i] <<- 
					abs(geneList[which(hits[, i])])^exponent)	
			NR <- colSums(Phit)		
			Pmiss[which(!hits)] <- 1/(N-nh)		
			Pmiss <- sapply(1:(nPermutations+1), function(i) 
				cumsum(Pmiss[, i]))
			Phit <- sapply(1:(nPermutations+1), function(i) 
				cumsum(Phit[, i])/NR[i])		
			runningES <- Phit-Pmiss		
			ESrange <- sapply(1:(nPermutations+1), function(i) 
				range(runningES[, i]))
			ES <- sapply(1:(nPermutations+1), function(i) 
				ESrange[which.max(abs(ESrange[,i])),i])	
			if(is.list(ES)) ES <- unlist(ES)
		}	
		##Return the relevant information according to mode		
		ES <- list(scoresObserved = ES[1], scoresperm = ES[2:(nPermutations+1)])
		return(ES)	
	}
	#parallel computing
	## scores <- parSapply(getOption("cluster"), 1:length(collectionOfGeneSets), 
	## 		function(i) {
	## 			gseaScoresBatchLocal(geneList, geneNames.perm = geneNames.perm, 
	## 			geneSet = as.integer(collectionOfGeneSets[[i]]), exponent = exponent, 
	## 			nPermutations = nPermutations)
	## 		}
        ##                 )
    scores <- do.call("cbind",(mclapply(1:length(collectionOfGeneSets),function(i) {	gseaScoresBatchLocal(geneList, geneNames.perm = geneNames.perm, geneSet = as.integer(collectionOfGeneSets[[i]]), exponent = exponent, nPermutations = nPermutations)},mc.cores=no.cores.gsea)))
	return(scores)
}

    ##     


owncollectionGsea <- function(collectionOfGeneSets, geneList, exponent=1, 
	nPermutations=1000, minGeneSetSize=15, verbose=TRUE,no.cores.gsea) {
	##check input arguments
	## paraCheck("gsc", collectionOfGeneSets)
	## paraCheck("genelist", geneList)
	## paraCheck("exponent", exponent)
	## paraCheck("minGeneSetSize", minGeneSetSize)
	geneList.names <- names(geneList)
##	paraCheck("nPermutations", nPermutations)	
	##tag the gene sets that can be used in the analysis, i.e. those 
	##that are smaller than the size of the gene list and that have more 
	##than 'minGeneSetSize' elements that can be found in the geneList	
	nGeneSets <- length(collectionOfGeneSets)
	tagGeneSets <- rep(FALSE, nGeneSets)
	tagGeneSets[which(unlist(lapply(collectionOfGeneSets, length)) < 
		length(geneList))] <- TRUE
	tagGeneSets[which(unlist(lapply(lapply(collectionOfGeneSets, 
		intersect, y=geneList.names), length)) < minGeneSetSize)] <- FALSE
	##check that there are actually some gene sets that pass the max 
	##and min cutoffs
	n.tagGeneSets <- sum(tagGeneSets)
	if(n.tagGeneSets == 0) 
		warning(paste("There are no gene sets in your collection",
			" that pass the cutoffs on size", sep=""))
	if(n.tagGeneSets > 0) {
		##Generate a matrix to store the permutation-based scores, with 
		##one row for each gene set (that has been tagged) and one column 
		##for each permutation	
		scoresperm <- matrix(rep(0, (nPermutations * n.tagGeneSets)), 
			nrow=n.tagGeneSets)
		rownames(scoresperm) <- names(collectionOfGeneSets)[which(tagGeneSets)]
		##Generate a vector to store the experimental scores
		##one entry for each gene set (that has been tagged)
		scoresObserved <- rep(0, n.tagGeneSets)
		names(scoresObserved) <- names(collectionOfGeneSets)[which(tagGeneSets)]
		##Compute the scores	
		##create permutation gene list
		perm.gL <- sapply(1:nPermutations, function(n) names(geneList)[
			sample(1:length(geneList), length(geneList),replace=FALSE)])
		perm.gL<-cbind(names(geneList),perm.gL)
		##check if package snow has been loaded and a cluster object 
		##has been created for HTSanalyzeR	
			scores <- gseaScoresbpown(geneList, geneNames.perm = perm.gL,
				collectionOfGeneSets=collectionOfGeneSets[which(tagGeneSets)],
				exponent=exponent,nPermutations=nPermutations,no.cores.gsea)
			sapply(1:n.tagGeneSets, function(i) {
					scoresperm[i,]<<-unlist(scores["scoresperm",i])
					scoresObserved[i]<<-unlist(scores["scoresObserved",i])
                        }
			)
	} else {
		scoresObserved <- NULL
		scoresperm <- NULL
	}
	return(list("Observed.scores" = scoresObserved , "Permutation.scores" = scoresperm))	
}


gsea.to.frame <- function(gsea.result){
    GSCscores <- gsea.result[["scores"]]
    GSCfdrs <- gsea.result[["fdrs"]]
    GSCpvals <- gsea.result[["pvals"]]
    first.hits.go <- tibble(fdr=GSCfdrs,pvals=GSCpvals,score=GSCscores$Observed.scores,term=names(gsea.result[["pvals"]]))
    first.hits.go
}



own.gsea <- function(data.vector,ListGSC,min.set.size=10,nperm=1000,no.cores.gsea){
    hits <- names(data.vector)[data.vector >quantile(data.vector,0.5)]
    
    gsca <- new("GSCA", listOfGeneSetCollections=ListGSC,geneList=data.vector, hits=hits)
    gsca <- suppressMessages(preprocess(gsca, species="Hs", initialIDs="Ensembl.gene",keepMultipleMappings=TRUE, duplicateRemoverMethod="max",orderAbsValue=FALSE,verbose=FALSE))
    GSCscores <- owncollectionGsea(collectionOfGeneSets=ListGSC[[1]], geneList=gsca@geneList, exponent=1, nPermutations=nperm, minGeneSetSize=min.set.size,verbose=FALSE,no.cores.gsea=no.cores.gsea)
    GSCfdrs <- FDRcollectionGsea(permScores=GSCscores$Permutation.scores,dataScores=GSCscores$Observed.scores)
    GSCpvals = permutationPvalueCollectionGsea(permScores=GSCscores$Permutation.scores, dataScores=GSCscores$Observed.scores)
    list(fdrs=GSCfdrs,pvals=GSCpvals,scores=GSCscores["Observed.scores"])
}
