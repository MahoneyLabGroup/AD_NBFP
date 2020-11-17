#This function plots the results of compare.gene.lists


plot.gene.comparison <- function(gene.comparison.table, net1.label = "Network1", net2.label = "Network2"){
	
	cols <- c(
	rgb(35/256,139/256,69/256, alpha = 0.5),
	rgb(106/256,81/256,163/256, alpha = 0.5))
	
	locus.names <- rownames(gene.comparison.table)[1:(nrow(gene.comparison.table)-1)]
	
	total.genes.per.locus <- gene.comparison.table[,"unique.genes"]
	num1 <- gene.comparison.table[,1]
	num2 <- gene.comparison.table[,2]
	overlap.with1 <- apply(gene.comparison.table, 1, function(x) x["num.genes1"]-x["Set1.not.in.Set2"])
	start2 <- num1 - overlap.with1
	end2 <- start2 + num2
	
	max.genes <- max(total.genes.per.locus[-length(total.genes.per.locus)])
	
	plot.new()
	plot.window(xlim = c(0,max.genes), ylim = c(0.5,nrow(gene.comparison.table)))
	
	for(i in 1:(nrow(gene.comparison.table)-1)){		
		draw.rectangle(0, num1[i], i-0.4, i+0.4, border.col = NA, fill = cols[1])
		draw.rectangle(start2[i], start2[i]+num2[i] , i-0.4, i+0.4, border.col = NA, fill = cols[2])	
		}	
	par(xpd = TRUE)
	text(x = rep(-0.1, length(locus.names)), y = 1:length(locus.names), labels = locus.names, adj = 1)
	mtext("Number of Genes", side = 1, line = 2)
		legend(x = 0, y = nrow(gene.comparison.table), fill = c("#78AC86", "#A596C6", "#6C7791"), legend = c(net1.label, net2.label, "Overlap"), yjust = 0)
	par(xpd = FALSE)
	axis(1)

}