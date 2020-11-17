#This function returns a table of QTL
#by chromosomal region in the Bloom data


get.QTL <- function(data.dir = "~/Documents/Data/yeast/Forsberg_et_al/", exp.type = c("full", "marginal"), buffer.bp = 2000){
	
	
	if(exp.type == "full"){
		filename = "Bloom_etal_interactions_full.csv"
		}else{
		filename = "Bloom_etal_interactions_marginal.csv"
		}

	qtl.margins <- read.csv(paste0(data.dir, "Bloom_etal_main_effects.csv"), stringsAsFactors = FALSE)
		
	int.table <- as.matrix(read.csv(paste0(data.dir, filename)))
	int.table[,4] <- gsub(pattern = "chr", "", int.table[,4])
	int.table[,6] <- gsub(pattern = "chr", "", int.table[,6])
	
		
	qtl1 <- apply(int.table, 1, function(x) paste(x[4], max(c(as.numeric(x[5])-buffer.bp), 1), as.numeric(x[5])+buffer.bp, sep = ":"))
	qtl2 <- apply(int.table, 1, function(x) paste(x[6], max(c(as.numeric(x[7])-buffer.bp), 1), as.numeric(x[7])+buffer.bp, sep = ":"))
	
	
	qtl.table <- cbind(qtl1, qtl2, int.table[,c(1,8:10)])
	
	trimmed.qtl <- trim.qtl.ends(qtl.table)
	
	return(trimmed.qtl)

}