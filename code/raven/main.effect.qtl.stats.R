
main.effect.qtl.stats <- function(){
	
	main.effects = as.matrix(read.csv("~/Documents/Data/yeast/Forsberg_et_al/Bloom_etal_main_effects.csv", stringsAsFactors = FALSE))
	
	chr.table <- as.matrix(read.table("~/Documents/Data/yeast/general_data/yeast_chromosomes.txt", sep = "\t", stringsAsFactors = FALSE))
	rownames(chr.table) <- chr.table[,1]
	chr.table <- chr.table[,-1,drop=FALSE]
	
	
	qtl.traits <- unique(main.effects[,2])
	all.qtl.size <- vector(mode = "list", length = length(qtl.traits))
	names(all.qtl.size) <- qtl.traits
	
	qtl.count.table <- matrix(NA, nrow = length(qtl.traits), ncol = nrow(chr.table))
	rownames(qtl.count.table) <- qtl.traits
	colnames(qtl.count.table) <- rownames(chr.table)
		
	for(tr in 1:length(qtl.traits)){
		trait.table <- main.effects[which(main.effects[,2] == qtl.traits[tr]),]
		trait.chr <- gsub("chr", "", trait.table[,1])
		qtl.count.table[tr,] <- unlist(lapply(rownames(chr.table), function(x) length(which(trait.chr == x))))	
		qtl.size <- as.numeric(trait.table[,10]) - as.numeric(trait.table[,8])		
		all.qtl.size[[tr]] <- log10(qtl.size)
		}

	median.qtl.size <- median(unlist(all.qtl.size))
	# pdf("Main.Effect.Stats.pdf", width = 10, height = 8)
	# par(mar = c(4,7,2,2))
	imageWithText(qtl.count.table, row.names = rownames(qtl.count.table), col.names = colnames(qtl.count.table), col.text.shift = -1, row.text.shift = -1, cex = 1, main = "Main-effect QTL counts by phenotype and chromosome")
	par(mar = c(9,4,2,2))
	boxplot(all.qtl.size, las = 2, main = "Main-effect QTL size by phenotype", ylab = "QTL size (log10 bp)")
	par(xpd = FALSE)
	abline(h = median.qtl.size , col = "red")
	# dev.off()

	return(10^median.qtl.size)
	
}