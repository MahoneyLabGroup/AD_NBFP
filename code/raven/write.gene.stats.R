#This is a function to write and plot out the tables 
#from get.gene.stats in the TRiAGE pipeline


write.gene.stats <- function(gene.stats){
	library(RColorBrewer)	
	cols <- brewer.pal(8, "Dark2")

		write.table(t(gene.stats[[1]]), paste0("gene.stats.", names(gene.stats)[1], ".csv"), sep = ",", quote = FALSE)
		bar.cols <- c("#7fc97f", "#beaed4", "#fdc086")

		pdf("Genes.Selected.by.each.Algorithm.pdf", width = 10, height = 5)
		par(mfrow = c(1,2))
		barplot(gene.stats[[1]], beside = TRUE, col = bar.cols, main = "Number of Genes")
		legend("topright", fill = bar.cols, legend = c("total", "SA", "Gr"))
		perc.SA <- gene.stats[[1]][2,]/gene.stats[[1]][1,]
		perc.Gr <- gene.stats[[1]][3,]/gene.stats[[1]][1,]
		barplot(rbind(perc.SA, perc.Gr), beside = TRUE, col = bar.cols[2:3], main = "Proportion of total genes selected by each algorithm")
		legend("topright", fill = bar.cols[2:3], legend = c("SA", "Gr"))
		dev.off()
	

		sa.gene.table <- gene.stats[[4]]
		for(i in 1:nrow(sa.gene.table)){
			sa.gene.table[i,] <- paste0(sa.gene.table[i,], " (", gene.stats[[2]][i,], ")")
			}
		write.table(sa.gene.table, "gene.stats.SA.genes.csv", sep = "," ,quote = FALSE, col.names = FALSE)


		pdf("gene.stats.SA.gene.counts.pdf", width = nrow(sa.gene.table)/2, height = 5)
		for(i in 1:ncol(gene.stats[[2]])){
			barplot(gene.stats[[2]][,i], las = 2, ylim = c(0,100), main = paste("Selection", i))
			}
		dev.off()
				
		
		gr.gene.table <- gene.stats[[5]]
		for(i in 1:nrow(gr.gene.table)){
			gr.gene.table[i,] <- paste0(gr.gene.table[i,], " (", gene.stats[[3]][i,], ")")
			}
		write.table(gr.gene.table, "gene.stats.Gr.genes.csv", sep = "," ,quote = FALSE, col.names = FALSE)
		
		pdf("gene.stats.Gr.gene.counts.pdf", width = nrow(gr.gene.table)/2, height = 5)
		for(i in 1:ncol(gene.stats[[3]])){
			barplot(gene.stats[[3]][,i], las = 2, ylim = c(0,100), main = paste("Selection", i))
			}
		dev.off()
	
		total.trials <- sum(gene.stats[[2]][1,], na.rm = TRUE)
	
		pdf("gene.stats.Gr.SA.comparison.pdf")
		for(i in 1:nrow(gene.stats[[2]])){
			plot(gene.stats[[2]][i,], type = "l", ylim = c(1,total.trials), ylab = "Times Selected", xlab = "Gene Rank", col = cols[1])
			points(x = 1:ncol(gene.stats[[3]]), y = gene.stats[[3]][i,], type = "l", col = cols[2])
			legend("topright", lty = 1, col = cols[1:2], legend = c("SA", "Gr"))
			}	
		dev.off()
	}