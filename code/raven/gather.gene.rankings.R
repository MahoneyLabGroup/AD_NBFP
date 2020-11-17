#This function takes in a directory and gathers the
#gene rankings for all subdirectories
#The parent directory contains all gene rankings
#for a single locus
#it expects two layers of directories: one for each
#gene set, and another for increasing numbers of genes
#used as the seed gene set
#This requires that there is a Low.FP table output from
#the SVM script

gather.gene.rankings <- function(parent.dir){
	require("e1071")

	setwd(parent.dir)
	all.sub.dir <- list.files()
	all.results.list <- list()
	list.names <- NULL
	
	list.count <- 1
	for(i in 1:length(all.sub.dir)){
		setwd(paste(parent.dir, all.sub.dir[i], sep = "/"))
		sub.sub <- list.files()
		for(j in 1:length(sub.sub)){
			setwd(paste(parent.dir, all.sub.dir[i], sub.sub[j], sep = "/"))
			list.names <- c(list.names, paste(all.sub.dir[i], sub.sub[j], sep = "_"))
			low.fp <- list.files(pattern = "Low.FP")
			low.fp.data <- read.table(low.fp, stringsAsFactors = FALSE, sep = "\t", header = TRUE)
			all.results.list[[list.count]] <- low.fp.data[,c("external_gene_name", "neg.log.fp")]
			list.count <- list.count + 1
			}
		}

	u_genes <- unique(unlist(lapply(all.results.list, function(x) x[,1])))
	
	get.gene.fp <- function(gene.name){
		fp <- rep(NA, length(all.results.list))
		for(i in 1:length(all.results.list)){
			gene.locale <- which(all.results.list[[i]][,1] == gene.name)
			if(length(gene.locale) > 0){
				fp[i] <- all.results.list[[i]][gene.locale[1],2]
				}
			}
		return(fp)
		}
	
	all.fp <- lapply(u_genes, get.gene.fp)
	all.fp.mat <- Reduce("rbind", all.fp)
	rownames(all.fp.mat) <- u_genes
	colnames(all.fp.mat) <- list.names
	
	na.locale <- which(is.na(all.fp.mat))
	all.fp.mat[na.locale] <- 0
	
	row.order <- hclust(dist(all.fp.mat))$order
	col.order <- hclust(dist(t(all.fp.mat)))$order
	
	all.fp.mat[na.locale] <- NA
	
	setwd(parent.dir)
	fp.means <- rowMeans(all.fp.mat)
	mean.order <- order(fp.means, decreasing = TRUE)

	pdf("All.Gene.FP.Values.pdf", width = 8, height = 40)
	imageWithText(all.fp.mat[mean.order, col.order], show.text = FALSE, row.names = rownames(all.fp.mat)[mean.order], col.names = colnames(all.fp.mat)[col.order], row.text.cex = 0.5, col.text.shift = -50)
	imageWithText(all.fp.mat[row.order, col.order], show.text = FALSE, row.names = rownames(all.fp.mat)[row.order], col.names = colnames(all.fp.mat)[col.order], row.text.cex = 0.5, col.text.shift = -50)	
	dev.off()
	
	write.table(all.fp.mat, "All.Gene.FP.Values.csv", sep = ",", quote = FALSE)
	
	gene.info <- getBM(c("external_gene_name", "start_position", "end_position"), "external_gene_name", u_genes, mus)
	gene.idx <- match(u_genes, gene.info[,1])
	gene.pos <- rowMeans(gene.info[gene.idx,c(2,3)])
	
	plot.new()
	plot.window(xlim = c(min(gene.pos), max(gene.pos)), ylim = c(min(fp.means), max(fp.means)))
	text(x = gene.pos, y = fp.means, labels = u_genes)
	axis(1);axis(2)	
	mtext("Genomic Position", side = 1, line = 2.5)
	mtext("Mean -log FP", side = 2, line = 2.5)	
	
	
	
}