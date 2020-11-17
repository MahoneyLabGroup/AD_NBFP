#This function checks a qtl table to make sure
#none of the qtl are going off the ends of the
#chromosomes after the buffer is applied

trim.qtl.ends <- function(qtl.table){
	
	yeast.chr <- read.table("~/Documents/Data/yeast/general_data/yeast_chromosomes.txt", stringsAsFactors = FALSE, sep = "\t")
	
	
	trim.loci <- function(loci){
		split.loci <- strsplit(loci, ":")
		chr <- unlist(lapply(split.loci, function(x) x[1]))
		locus.start <- as.numeric(unlist(lapply(split.loci, function(x) x[2])))
		locus.start[which(locus.start < 1)] <- 1
		locus.end <- as.numeric(unlist(lapply(split.loci, function(x) x[3])))
		
		chr.size <- as.numeric(yeast.chr[,2])
		for(ch in 1:nrow(yeast.chr)){
			chr.locale <- which(chr == yeast.chr[ch,1])
			if(length(chr.locale) > 0){
				chr.ends <- locus.end[chr.locale]
				chr.ends[which(chr.ends > chr.size[ch])] <- chr.size[ch]
				locus.end[chr.locale] <- chr.ends
				}
			}
		
		#put the trimmed loci back together
		final.loci <- apply(cbind(chr, locus.start, locus.end), 1, function(x) paste(x[1], x[2], x[3], sep = ":"))
		return(final.loci)
	}

	qtl.table[,1] <- trim.loci(qtl.table[,1])
	qtl.table[,2] <- trim.loci(qtl.table[,2])
	return(qtl.table)
}