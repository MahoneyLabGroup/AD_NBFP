#This function looks at stats from a qtl table

qtl.stats <- function(qtl.table, qtl.genes, filename = "QTL.by.Chromosome.pdf"){
	
	yeast.chr <- read.table("~/Documents/Data/yeast/general_data/yeast_chromosomes.txt", stringsAsFactors = FALSE, sep = "\t")
	
	#============================================================
	#internal functions
	#============================================================
	get.chr <- function(loci){
		split.loci <- strsplit(unique(loci), ":")
		chr <- unlist(lapply(split.loci, function(x) x[1]))
		qtl.counts <- lapply(yeast.chr[,1], function(x) length(which(chr == x)))
		names(qtl.counts) <- yeast.chr[,1]
		return(unlist(qtl.counts))
		}
		
	get.size <- function(loci, locus.genes, plot.label = ""){
		u_loci <- unique(loci)
		split.loci <- strsplit(u_loci, ":")
		chr <- unlist(lapply(split.loci, function(x) x[1]))
		locus.start <- as.numeric(unlist(lapply(split.loci, function(x) x[2])))
		locus.end <- as.numeric(unlist(lapply(split.loci, function(x) x[3])))
		locus.size <- locus.end - locus.start
		chr.size <- yeast.chr[,2]
		max.len <- max(as.numeric(yeast.chr[,2]))
		num.chr <- nrow(yeast.chr)
		plot.new()
		plot.window(xlim = c(0, max.len), ylim = c(0,(num.chr+1)))
		for(i in 1:nrow(yeast.chr)){
			yval <- num.chr - i + 1
			draw.rectangle(1, as.numeric(yeast.chr[i,2]), yval-0.1, yval+0.1, border = "#386cb0")
			chr.locale <- which(chr == yeast.chr[i,1])
			par(xpd = TRUE)
			text(x = (0-(max.len*0.05)), y = yval, labels = yeast.chr[i,1], cex = 0.8)
			par(xpd = FALSE)
			if(length(chr.locale) > 0){
				for(j in chr.locale){
					locus.name <- paste(yeast.chr[i,1], locus.start[j], locus.end[j], sep = ":")
					locus.pos <- which(names(locus.genes) == locus.name)
					draw.rectangle(locus.start[j], locus.end[j], yval-0.1, yval+0.1, border = "#6a3d9a", fill = "#7fc97f", lwd = 0.5)
					text(x = mean(c(locus.start[j], locus.end[j])), y = yval+0.3, labels = nrow(locus.genes[[locus.pos]]), cex = 0.6, col = "#737373")
					}
				}
			}
		mtext(plot.label)
		
		names(locus.size) <- u_loci
		return(locus.size)
		}
		
	get.locus.genes <- function(loci, qtl.genes){
		u_loci <- unique(loci)
		split.loci <- strsplit(u_loci, ":")
		locus.genes <- vector(mode = "list", length = length(split.loci))
		names(locus.genes) <- u_loci
		chr <- unlist(lapply(split.loci, function(x) x[1]))
		locus.start <- as.numeric(unlist(lapply(split.loci, function(x) x[2])))
		locus.end <- as.numeric(unlist(lapply(split.loci, function(x) x[3])))
		for(i in 1:length(split.loci)){
			chr.locale <- which(qtl.genes[,"chromosome_name"] == chr[i])
			chr.genes <- qtl.genes[chr.locale,]
			overlap.locale <- apply(chr.genes, 1, function(x) segments.overlap(locus.start[i], locus.end[i], as.numeric(x[5]), as.numeric(x[6])))
			if(length(which(overlap.locale)) > 0){
				locus.genes[[i]] <- chr.genes[which(overlap.locale),]
				}
			}
		return(locus.genes)
		}
	#============================================================	
	
	locus.by.chr <- get.chr(c(qtl.table[,1], qtl.table[,2]))
	
	# pdf("Locus.Counts.by.Chromosome.pdf", width = 10, height = 4)
	# quartz(width = 8, height = 4)
	barplot(locus.by.chr, main = "Loci by Chromosome")
	# dev.off()
	
	
	loc.genes <- get.locus.genes(loci = c(qtl.table[,2], qtl.table[,1]), qtl.genes)
	num.genes <- unlist(lapply(loc.genes, nrow))
	
	# pdf(filename)
	locus.size <- get.size(loci = c(qtl.table[,2], qtl.table[,1]), loc.genes, "QTL locations with number of genes")
	# dev.off()
	
	qtl.stat.table <- cbind(names(locus.size), locus.size, num.genes)
	invisible(qtl.stat.table)

}