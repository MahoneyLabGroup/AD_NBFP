#This function gets genes for linkage blocks
#in a data.obj

# setwd("~/Documents/Data/Little_Cross/little_cross_data_cape/Body_Comp_IGF_New_P")
# data.obj <- readRDS("cross.RData"); geno.obj <- NULL
# rownames(data.obj$pheno) <- rownames(data.obj$geno) <- 1:nrow(data.obj$pheno)
# data.obj <- marker2covar(data.obj, markers = c("sex", "Final.IGF.1"))
# data.obj <- get.network(data.obj, p.or.q = 0.0005, standardized = FALSE, collapse.linked.markers = FALSE, verbose = TRUE)

#fix for DO


block2gene <- function(data.obj, organism = c("mouse", "human"), collapsed.net = FALSE, remove.covar = FALSE, min.effect.size = 3, upstream.buffer = 0, downstream.buffer = 0, lookup.marker.position = FALSE){
	
	require(biomaRt)

	organism <- organism[1]

	if(organism[1] == "mouse"){		
	# lib <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")
		lib <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host = "may2017.archive.ensembl.org")
		snp.db = useMart(biomart="ENSEMBL_MART_SNP", dataset="mmusculus_snp", host = "may2017.archive.ensembl.org")	
	}else{
	# lib <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")	
	lib <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "may2017.archive.ensembl.org")
	snp.db = useMart(biomart="ENSEMBL_MART_SNP", dataset="hsapiens_snp", host = "may2017.archive.ensembl.org")
	}
	

	atts <- c("external_gene_name", "entrezgene", "ensembl_gene_id", "chromosome_name", "start_position", "end_position")
	fils <- "chromosomal_region"


	if(!collapsed.net){
		blocks <- data.obj$linkage.blocks.full
		net <- data.obj$full.net
		}else{
		blocks <- data.obj$linkage.blocks.collapsed
		net <- data.obj$collapsed.net			
		}
	
	just.int <- as.matrix(net[,1:nrow(net)])
	# hist(just.int[which(just.int != 0)], breaks = 100)
	# abline(v = min.effect.size, col = "red");abline(v = min.effect.size*-1, col = "red")
	
	if(remove.covar){
		covar.locale <- which(rownames(just.int) %in% data.obj$p.covar)
		if(length(covar.locale) > 0){
			just.int <- just.int[-covar.locale,-covar.locale,drop=FALSE]
			}
		}
	
	int.markers <- which(abs(just.int) > min.effect.size, arr.ind = TRUE)
	cat("Interactions:", nrow(int.markers), "\n")
	
	int.blocks <- apply(int.markers, 2, function(x) rownames(just.int)[x])
	
	u_blocks <- unique(c(int.blocks[,1], int.blocks[,2]))
	cat("Linkage blocks:", length(u_blocks), "\n")
	
	
	if(length(data.obj$geno.names) == 3){ #then we have a DO object
		block.alleles <- sapply(u_blocks, function(x) get.block.allele(data.obj, x, collapsed.net))
		block.markers <- sapply(u_blocks, function(x) get.block.marker(data.obj, x, collapsed.net))
		block.markers <- lapply(block.markers, function(x) unlist(lapply(strsplit(x, "_"), function(x) x[1])))
		block.chr <- lapply(block.markers, function(x) data.obj$chromosome[which(data.obj$geno.names[[3]] %in% x)])

		if(lookup.marker.position){
			cat("looking up SNP positions...\n")
			snp.info <- lapply(block.markers, function(x) as.matrix(getBM(c("refsnp_id","allele","chr_name","chrom_start"), filters = "snp_filter", values = x, mart = snp.db)))
			block.bp <- lapply(snp.info, function(x) as.numeric(x[,4]))
			}else{
			block.bp <- lapply(block.markers, function(x) data.obj$marker.location[which(data.obj$geno.names[[3]] %in% x)])
			}

		block.start <- round(unlist(lapply(block.bp, function(x) if(length(x) > 0){min(x)}))) - upstream.buffer
		block.end <- round(unlist(lapply(block.bp, function(x) if(length(x) > 0){max(x)}))) + downstream.buffer
		block.chr <- unlist(lapply(block.chr, unique))
		}else{	
		block.markers <- lapply(u_blocks, function(x) blocks[which(names(blocks) == x)])
		cat("Encompassing", length(unlist(block.markers)), "Markers\n")

		block.chr <- lapply(block.markers, function(x) data.obj$chromosome[which(data.obj$marker.num %in% x[[1]])])
		block.bp <- lapply(block.markers, function(x) data.obj$marker.location[which(data.obj$marker.num %in% x[[1]])])
		block.start <- round(unlist(lapply(block.bp, function(x) if(length(x) > 0){min(x)}))) - upstream.buffer
		block.end <- round(unlist(lapply(block.bp, function(x) if(length(x) > 0){max(x)}))) + downstream.buffer
		block.chr <- unlist(lapply(block.chr, unique))
		}

	region.table <- cbind(block.chr, block.start, block.end)
	query <- apply(region.table, 1, function(x) paste(x, collapse = ":"))
	# mm.query <- apply(region.table, 1, function(x) paste(paste(x[1:2], collapse = ":"), x[3], sep = ".."))
	cat("Finding genes...\n")
	region.genes <- getBM(atts, fils, values = query, mart = lib)
	cat("Found the following unique elements:\n")
	print(apply(region.genes, 2, function(x) length(unique(x))))
	return(region.genes)
	}