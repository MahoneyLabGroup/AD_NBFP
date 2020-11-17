#This function takes a list of vectors
#and plots a Venn diagram of their overlap

plotVenn <- function(Vlist, cat.names = NULL, lwd = 3, cex = 2, label.cex = 2){
	require(VennDiagram)
	all.cols <- c("#7fc97f", "#beaed4", "#fdc086", "#ffff99", "#386cb0")


	if(is.null(cat.names)){cat.names <- names(Vlist)}

	if(length(Vlist) == 2){
		shared <- intersect(Vlist[[1]], Vlist[[2]])
		list12 <- length(shared)
		
		draw.pairwise.venn(length(Vlist[[1]]), length(Vlist[[2]]), list12, category = cat.names, col = all.cols[1:2], lwd = lwd, cex = cex, cat.cex = label.cex)
	
		list1.only <- setdiff(Vlist[[1]], Vlist[[2]])
		list2.only <- setdiff(Vlist[[2]], Vlist[[1]])
	
		results <- list(list1.only, list2.only, shared)
		names(results) <- c(cat.names, "shared")
		invisible(results)
		}
		
	if(length(Vlist) == 3){
		sig12 <- length(intersect(Vlist[[1]], Vlist[[2]]))
		sig13 <- length(intersect(Vlist[[1]], Vlist[[3]]))
		sig23 <- length(intersect(Vlist[[2]], Vlist[[3]]))
		sig123 <- length(Reduce(intersect, list(Vlist[[1]], Vlist[[2]], Vlist[[3]])))
		
		draw.triple.venn(length(Vlist[[1]]), length(Vlist[[2]]), length(Vlist[[3]]), sig12, sig23, sig13, sig123, category = cat.names, col = all.cols[1:3], lwd = lwd, cex = cex, cat.cex = label.cex)

		}
		
	if(length(Vlist) == 4){
		all.area <- lapply(Vlist, length)
		n12 <- length(intersect(Vlist[[1]], Vlist[[2]]))
		n13 <- length(intersect(Vlist[[1]], Vlist[[3]]))
		n14 <- length(intersect(Vlist[[1]], Vlist[[4]]))
		n23 <- length(intersect(Vlist[[2]], Vlist[[3]]))		
		n24 <- length(intersect(Vlist[[2]], Vlist[[4]]))
		n34 <- length(intersect(Vlist[[3]], Vlist[[4]]))
		n123 <- length(Reduce(intersect, Vlist[1:3]))
		n124 <- length(Reduce(intersect, Vlist[c(1,2,4)]))
		n134 <- length(Reduce(intersect, Vlist[c(1,3,4)]))
		n234 <- length(Reduce(intersect, Vlist[c(2,3,4)]))
		n1234 <- length(Reduce(intersect, Vlist[c(1:4)]))
		
		draw.quad.venn(all.area[[1]], all.area[[2]], all.area[[3]], all.area[[4]], n12, n13, n14, n23, n24, n34, n123, n124, n134, n234, n1234, category = cat.names, col = all.cols[1:4], lwd = lwd, cex = cex, cat.cex = label.cex)

		}
		
	if(length(Vlist) == 5){
		all.area <- lapply(Vlist, length)
		n12 <- length(intersect(Vlist[[1]], Vlist[[2]]))
		n13 <- length(intersect(Vlist[[1]], Vlist[[3]]))
		n14 <- length(intersect(Vlist[[1]], Vlist[[4]]))
		n15 <- length(intersect(Vlist[[1]], Vlist[[5]]))
		n23 <- length(intersect(Vlist[[2]], Vlist[[3]]))		
		n24 <- length(intersect(Vlist[[2]], Vlist[[4]]))
		n25 <- length(intersect(Vlist[[2]], Vlist[[5]]))		
		n34 <- length(intersect(Vlist[[3]], Vlist[[4]]))
		n35 <- length(intersect(Vlist[[3]], Vlist[[5]]))
		n45 <- length(intersect(Vlist[[4]], Vlist[[5]]))
		n123 <- length(Reduce(intersect, Vlist[1:3]))
		n124 <- length(Reduce(intersect, Vlist[c(1,2,4)]))
		n125 <- length(Reduce(intersect, Vlist[c(1,2,5)]))
		n134 <- length(Reduce(intersect, Vlist[c(1,3,4)]))
		n135 <- length(Reduce(intersect, Vlist[c(1,3,5)]))
		n145 <- length(Reduce(intersect, Vlist[c(1,4,5)]))
		n234 <- length(Reduce(intersect, Vlist[c(2,3,4)]))
		n235 <- length(Reduce(intersect, Vlist[c(2,3,5)]))
		n245 <- length(Reduce(intersect, Vlist[c(2,4,5)]))
		n345 <- length(Reduce(intersect, Vlist[c(3,4,5)]))
		n1234 <- length(Reduce(intersect, Vlist[c(1:4)]))
		n1235 <- length(Reduce(intersect, Vlist[c(1,2,3,5)]))
		n1245 <- length(Reduce(intersect, Vlist[c(1,2,4,5)]))
		n1345 <- length(Reduce(intersect, Vlist[c(1,3,4,5)]))
		n2345 <- length(Reduce(intersect, Vlist[c(2,3,4,5)]))
		n12345 <- length(Reduce(intersect, Vlist[c(1:5)]))
		
		draw.quintuple.venn(all.area[[1]], all.area[[2]], all.area[[3]], all.area[[4]], all.area[[5]], n12, n13, n14, n15, n23, n24, n25, n34, n35, n45, n123, n124, n125, n134, n135, n145, n234, n235, n245, n345, n1234, n1235, n1245, n1345, n2345, n12345, category = cat.names, col = all.cols[1:5], lwd = lwd, cex = cex, cat.cex = label.cex)
	}	
		
		

	}