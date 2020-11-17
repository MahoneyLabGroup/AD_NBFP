#This function gets candidate genes above a give voting score
#alt.cutoff sets the number of top-ranked genes if there are no genes that
#pass the cutoff


get.candidate.genes <- function(vote.cutoff = 90, alt.cutoff = 1){
	
	all.files <- get.files(want = "Ranked.Genes", dont.want = c("RData", "Edgelist","Network"))
	snps <- unlist(lapply(strsplit(all.files, "\\."), function(x) x[3]))
	gene.list <- vector(mode = "list", length = length(all.files))
	names(gene.list) <- snps

	for(i in 1:length(all.files)){
		gene.table <- read.table(all.files[i], stringsAsFactors = FALSE, header = TRUE, sep = "\t")
		candidate.locale <- which(gene.table[,"gene.vote"] >= vote.cutoff)
		if(length(candidate.locale) == 0){
			candidate.locale <- 1
			}
		gene.list[[i]] <- gene.table[candidate.locale, c("external_gene_name", "entrezgene","gene.vote")]
		}
	return(gene.list)
}