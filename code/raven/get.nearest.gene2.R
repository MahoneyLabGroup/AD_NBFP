#This function finds the nearest gene to a SNP
#if just.snp.position is TRUE, the function does
#not find the nearest gene but just returns the SNP
#positions

get.nearest.gene2 <- function(snp.names, organism = c("human", "mouse"), max.distance = 1e6, write.gene.table = TRUE, filename = "Genes.Near.SNPs.txt", restrict.to.entrez = TRUE, just.snp.position = FALSE){

	require("biomaRt")	
	require("stringr")
	organism <- organism[1]

	if(organism[1] == "mouse"){		
	# lib <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")
		lib <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host = "apr2018.archive.ensembl.org")
		snp.db = useMart(biomart="ENSEMBL_MART_SNP", dataset="mmusculus_snp", host = "apr2018.archive.ensembl.org")	
	}else{
	# lib <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")	
	lib <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "apr2018.archive.ensembl.org")
	snp.db = useMart(biomart="ENSEMBL_MART_SNP", dataset="hsapiens_snp", host = "apr2018.archive.ensembl.org")
	}

	cat("Getting SNP locations...\n")
	snp.table <- getBM(c("refsnp_id", "chr_name", "chrom_start"), filters = "snp_filter", values = snp.names, mart = snp.db)
	num.chr <- which(!is.na(as.numeric(snp.table[,2])))
	snp.table <- snp.table[num.chr,]
	
	if(just.snp.position){return(snp.table)}

	all.chr.region <- apply(snp.table, 1, function(x) paste0(x[2], ":", as.numeric(x[3])-max.distance, ":", as.numeric(x[3])+max.distance))
	
	cat("Finding genes in regions near SNPs...\n")
	region.genes <- lapply(all.chr.region, function(x) getBM(c("external_gene_name", "entrezgene","chromosome_name","start_position", "end_position"), "chromosomal_region", x, lib))


	if(restrict.to.entrez){	
		region.genes <- lapply(region.genes, function(x) x[which(!is.na(x[,"entrezgene"])),])
		}
	
	snp.gene.table <- matrix(NA, nrow = nrow(snp.table),  ncol = (ncol(snp.table)+6))
	
	for(i in 1:nrow(snp.table)){	
		snp.chr <- as.numeric(snp.table[i,"chr_name"])
		snp.locale <- as.numeric(snp.table[i,"chrom_start"])
		genes.near.snp <- region.genes[[i]]		
		
		if(nrow(genes.near.snp) > 0){
			nearest.start.idx <- get.nearest.pt(genes.near.snp[,"start_position"], snp.locale)
			start.dist <- abs(genes.near.snp[nearest.start.idx,"start_position"] - snp.locale)
			nearest.end.idx <- get.nearest.pt(genes.near.snp[,"end_position"], snp.locale)
			end.dist <- abs(genes.near.snp[nearest.end.idx,"end_position"] - snp.locale)
			
			if(start.dist < end.dist){
				gene.row <- unlist(c(snp.table[i,], genes.near.snp[nearest.start.idx,], start.dist))
				}else{
				gene.row <- unlist(c(snp.table[i,], genes.near.snp[nearest.end.idx,], end.dist))	
				}
			}else{
			gene.row <- unlist(c(snp.table[i,], rep(NA, (ncol(genes.near.snp)+1))))
			}
	snp.gene.table[i,] <- gene.row
	}
	colnames(snp.gene.table) <- c("SNP", "Chr", "Position", "Nearest.Gene", "entrez_id", "Chr", "start_position", "end_position", "Distance_to_gene")

	write.table(snp.gene.table, filename, quote = FALSE, sep = "\t", row.names = FALSE)

	invisible(snp.gene.table)
	
}