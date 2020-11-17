plot.degree.by.rank <- function(rank.table, id.table, gene.degree, plot.label){
	
	
	all.gene.deg <- NULL
	all.gene.guesses <- NULL
	for(i in 1:nrow(rank.table)){
		gene.names <- id.table[i,]
		all.gene.deg <- c(all.gene.deg, gene.degree[match(gene.names, names(gene.degree))])
		all.gene.guesses <- c(all.gene.guesses, rank.table[i,])
		}
	model <- lm(all.gene.guesses~ all.gene.deg)
	plot(all.gene.deg, all.gene.guesses, main = plot.label)
	abline(model)

	
	
	
}