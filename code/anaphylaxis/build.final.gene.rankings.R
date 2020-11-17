
build.final.gene.rankings <- function(snp.genes, snps, neg.log.fp, results.dir, scale.scores = TRUE){
	
	
	common.snps <- intersect(snp.genes[,"SNP"], snps[,"SNP"])
	common.gene.locale <- match(common.snps, snp.genes[,"SNP"])
	# head(cbind(common.snps, snp.genes[common.gene.locale,]))
	common.snp.locale <- match(common.snps, snps[,"SNP"])
	# head(cbind(common.snps, snps[common.snp.locale,]))

	neg.log.P <- -log10(as.numeric(snps[common.snp.locale,c("p.value")]))
	snp.table <- cbind(snp.genes[common.gene.locale,c("SNP", "Nearest.Gene", "start_position", 
	"end_position")], neg.log.P)
	
	
	#find the genes that have SNPs associated with them from the EMMA analysis.
	genes.with.SNPs <- unique(snp.table[,2])
	
	#remove genes without position information
	genes.with.SNPs <- genes.with.SNPs[which(!is.na(genes.with.SNPs))]
	#for each of these genes, find the maximum -log EMMA p value and the maximum
	#functional score across modules.
	all.gene.info <- t(sapply(genes.with.SNPs, function(x) combine.snp.fp(snp.info = snp.table, 
	fp.table = neg.log.fp, x)))
	
	#remove genes without functional information
	all.gene.info <- all.gene.info[which(!is.na(all.gene.info[,4])),]
	
	if(scale.scores){
		max.emma <- max(as.numeric(all.gene.info[,"EMMA.p"]), na.rm = TRUE)
		max.fp <- max(as.numeric(all.gene.info[,"FP"]), na.rm = TRUE)
		all.gene.info[,"EMMA.p"] <- as.numeric(all.gene.info[,"EMMA.p"])/max.emma
		all.gene.info[,"FP"] <- as.numeric(all.gene.info[,"FP"])/max.fp
	}

	gene.final.score <- as.numeric(all.gene.info[,"EMMA.p"]) + as.numeric(all.gene.info[,"FP"])
	gene.order <- order(gene.final.score, decreasing = TRUE)
	all.gene.info <- cbind(all.gene.info, gene.final.score)
	all.gene.info <- all.gene.info[gene.order,]
	gene.names <- all.gene.info[,1]	
	all.gene.info <- all.gene.info[,-1] #remove the gene names from the table
	all.gene.info <- apply(all.gene.info, 2, as.numeric) #make the table numeric
	rownames(all.gene.info) <- gene.names #put the rownames back


	write.table(all.gene.info, file.path(results.dir, "Final.Gene.Rankings.txt"), sep = "\t", quote = FALSE)
	return(all.gene.info)
	}