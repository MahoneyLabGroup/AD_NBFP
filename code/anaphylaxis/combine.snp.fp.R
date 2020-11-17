#This function takes in a gene name and finds
#the most significant associated SNP as well 
#as its false positive rates across all modules.
#This function operates on specifically formatted
#tables, and is not a general function

combine.snp.fp <- function(snp.info, fp.table, gene.name){
	
	gene.snp <- gene.fp <- NA
	gene.locale <- which(snp.info[,2] == gene.name)
	if(length(gene.locale) > 0){
		gene.snps <- snp.info[gene.locale,,drop=FALSE]
		gene.pos <- mean(as.numeric(gene.snps[1,c("start_position", "end_position")]))
		gene.max.log.p <- max(as.numeric(gene.snps[,"neg.log.P"]), na.rm = TRUE)
		}
		
	gene.locale <- which(rownames(fp.table) == gene.name)
	if(length(gene.locale) > 0){
		gene.fp <- max(fp.table[gene.locale,], na.rm = TRUE)
		}
		
	final.result <- c(gene.name, gene.pos, gene.max.log.p, gene.fp)
	names(final.result) <- c("gene.name", "gene.position", "EMMA.p", "FP")
	return(final.result)
	}
