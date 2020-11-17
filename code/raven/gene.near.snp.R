#This functon finds the genes up and downstream of 
#a given SNP
#if return.nearest = TRUE, only the nearest gene
#will be returned. Otherwise all genes in the region
#will be returned
#if entrez.only is TRUE only genes with entrez ids
#will be considered, otherwise, any gene (including psuedogenes
#and predicted genes) will be returned

gene.near.snp <- function(snp.name, mart, snp.db, upstream.buffer = 5000, downstream.buffer = 5000, return.nearest = FALSE, entrez.only = FALSE){
	
	require(biomaRt)

	# organism = organism[1]
	
	# if(organism == "human"){
		# mart = useEnsembl(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
		# snp.db = useEnsembl(biomart="ENSEMBL_MART_SNP", dataset="hsapiens_snp")
		# }else{
		# mart = useEnsembl(biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")
		# snp.db = useEnsembl(biomart="ENSEMBL_MART_SNP", dataset="mmusculus_snp")
		# }

	snp.info <- as.matrix(getBM(c("refsnp_id","allele","chr_name","chrom_start"), filters = "snp_filter", values = snp.name, mart = snp.db))
	
	if(nrow(snp.info) == 0){return("cannot find SNP")}
	
	snp.pos <- as.numeric(snp.info[1,"chrom_start"])
	snp.chr <- as.numeric(snp.info[1,"chr_name"])

	chr.region <- paste0(snp.chr, ":", snp.pos-upstream.buffer, ":", snp.pos+downstream.buffer)	
		
	att = c('entrezgene','external_gene_name', 'chromosome_name','start_position','end_position', 'description', 'phenotype_description')


	genes <- getBM(att, filters = "chromosomal_region", values = chr.region, mart = mart)
	
	if(nrow(genes) == 0){
		return(list(snp.info, "No nearby genes"))
		}
	
	if(entrez.only){
		not.na.locale <- which(!is.na(genes[,1]))
		if(length(not.na.locale) == 0){
			return(list(snp.info,"No nearby entrezgenes"))
			}
		genes <- genes[not.na.locale,]
		}
	
	
	if(return.nearest){
		near.start <- get.nearest.pt(genes[,"start_position"], snp.pos)
		start.dist <- min(abs(genes[near.start,4:5] - snp.pos))
	
		near.end <- get.nearest.pt(genes[,"end_position"], snp.pos)
		end.dist <- min(abs(genes[near.end,4:5] - snp.pos))
		
		if(start.dist < end.dist){
			genes <- genes[near.start,]
			}else{
			genes <- genes[near.end,]
			}
		
		}

	result <- list(snp.info, genes)



	return(result)

}