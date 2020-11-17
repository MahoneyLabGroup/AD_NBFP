#This function plots a heatmap of a matrix
#that contains NAs. It converts the NAs to 0's
#clusters the matrix, and then plots the matrix
#in the clustered order, but with the NAs returned.

heatmap.with.nas <- function(mat, cluster_cols = TRUE, cluster_rows = TRUE){
	
	new.mat <- mat
	na.locale <- which(is.na(new.mat))
	new.mat[na.locale] <- 0
	
	if(cluster_cols){
		d.col <- dist(t(new.mat))
		col.order <- hclust(d.col)$order
		}else{
		col.order <- 1:ncol(new.mat)
		}
	
	if(cluster_rows){
		d.row <- dist(new.mat)
		row.order <- hclust(d.row)$order
		}else{
		row.order <- 1:nrow(new.mat)	
		}
	
	pheatmap(mat[row.order, col.order], cluster_rows = FALSE, cluster_cols = FALSE)	
	
	invisible(mat[row.order, col.order])
	}