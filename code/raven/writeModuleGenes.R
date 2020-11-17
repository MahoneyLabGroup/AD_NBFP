writeModuleGenes <- function(results.dir, mart){
	
	module.dir <- get.module.dir(results.dir)
	for(i in 1:length(module.dir)){
			
		results.file <- paste0(module.dir[i], "/Module.Gene.Info.csv")	
		decomp.mat.file <- list.files(module.dir[i], pattern = "decomp", full.names = TRUE)
		if(length(decomp.mat.file) == 0){
			stop(paste0("I could not find a non.decomp.mat.RData a decomp.mat.RData file in ", dir.table[i,1], " ", dir.table[i,2], ". Please make sure generate.triage.models() was run."))
			}
		decomp.mat <- readRDS(decomp.mat.file)
		
		gene.id <- rownames(decomp.mat)
		
		gene.info <- getBM(c("external_gene_name", "entrezgene", "chromosome_name", "start_position", "end_position"), "entrezgene", values = as.numeric(gene.id), mart = mart)
		write.table(gene.info, results.file, quote = FALSE, sep = ",", row.name = FALSE)
		}
	
	
}