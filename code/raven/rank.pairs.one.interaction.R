#This function takes output from triage.one.interaction
#and ranks epistatic gene pairs.

rank.pairs.one.interaction <- function(locus1.genes, locus2.genes, tissue.adj, svm.results, 
top.n.pairs = 100){
	
	# Make a matrix describing the interactions between 
	# all gene pairs across the two interacting loci
	locus1.idx <- match(locus1.genes, rownames(tissue.adj))
	locus1.idx <- locus1.idx[which(!is.na(locus1.idx))]
	locus2.idx <- match(locus2.genes, colnames(tissue.adj))
	locus2.idx <- locus2.idx[which(!is.na(locus2.idx))]
	two.locus <- as.matrix(tissue.adj[locus1.idx,locus2.idx])
	rownames(two.locus) <- colnames(tissue.adj)[locus1.idx]
	colnames(two.locus) <- colnames(tissue.adj)[locus2.idx]

	#genes without any connections are removed from the matrix.
	#we need to remove these from our locus genes lists
	missing.locus1 <- setdiff(locus1.genes, rownames(two.locus))
	missing.locus2 <- setdiff(locus2.genes, colnames(two.locus))
	locus1.genes <- setdiff(locus1.genes, missing.locus1)
	locus2.genes <- setdiff(locus2.genes, missing.locus2)
	
	all.locus.genes <- c(locus1.genes, locus2.genes)

	#find the phenotype feature vectors for 
	pheno.feat <- svm.results$pheno.feat
	locus.pheno.feat = pheno.feat[rownames(pheno.feat) %in% all.locus.genes, ]

	locus.gene.feat = build.pairwise.features.two.locus(A=two.locus, X=locus.pheno.feat)
	
	svm.out <- svm.results$svm.out
	
	locus.predict <- predict(svm.out$opt.model, locus.gene.feat, decision.values = TRUE)

	
	dec.vals = attributes(locus.predict)$decision.values
	
	sort.dec.vals.idx = sort(dec.vals, index.return = TRUE, decreasing = TRUE)$ix
	
	max.idx = which(dec.vals == max(dec.vals))
	max.label = rownames(locus.gene.feat)[max.idx]
	
	top.pairs = rownames(locus.gene.feat)[sort.dec.vals.idx[1:top.n.pairs]]
	# top.genes = unique(as.numeric(unlist(strsplit(top.pairs, "-"))))
	pair.mat <- Reduce("rbind", strsplit(top.pairs, "-"))
	pair.mat <- apply(pair.mat, 2, as.numeric)
	final.mat <- cbind(pair.mat, dec.vals[sort.dec.vals.idx[1:top.n.pairs]])
	
	return(final.mat)
	
}