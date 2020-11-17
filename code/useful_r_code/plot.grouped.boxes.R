plot.grouped.boxes <- function(group.list, group.labels = names(group.list), 
group.cols = c("#a6cee3", "#b2df8a"), main = "", type = c("list", "matrix"), 
print.vals = c("mean", "median"), text.cex = 0.7, label.srt = 0){

	print.vals <- print.vals[1]
	type <- type[1]

	ymin <- min(unlist(group.list), na.rm = TRUE)*0.9
	ymax <- max(unlist(group.list), na.rm = TRUE)*1.1
	xmin = 0

	if(type == "list"){
		xmax <- sum(unlist(lapply(group.list, length))) + 1
		plot.window(xlim = c(0,xmax), ylim = c(ymin,ymax))
		max.iter <- length(group.list[[1]]) #number of elements in each list
		}else{
		xmax <- sum(unlist(lapply(group.list, function(x) dim(x)[2])))
		plot.window(xlim = c(0,xmax), ylim = c(ymin,ymax))	
		max.iter <- dim(group.list[[1]])[2]
		}
	box.pos <- 1
	plot.height <- ymax - ymin

	plot.new()
	plot.window(xlim = c(xmin, xmax), ylim = c(ymin, ymax))

	for(i in 1:max.iter){ #for each element of an individual group
		group.pos <- rep(NA, length(group.list))
		for(l in 1:length(group.list)){ #plot the ith element from group l
			
			if(type == "list"){
				data.vals <- group.list[[l]][[i]];label <- names(group.list[[l]])[i]
				}else{
				data.vals <- group.list[[l]][,i]; label <- colnames(group.list[[l]])[i]
				}
						
			boxplot(data.vals, at = box.pos, add = TRUE, col = group.cols[l], axes = FALSE, main = main)
			
			if(!is.na(print.vals)){
			par(xpd = TRUE)
			if(length(data.vals) > 0){
				if(print.vals == "mean"){
					print.val1 <- signif(mean(data.vals, na.rm = TRUE), 2)
					}else{
					print.val1 <- signif(median(data.vals, na.rm = TRUE), 2)
					}
				print.val2 <- signif(sd(data.vals, na.rm = TRUE), 2)
				text(x = box.pos, y = (ymin - (plot.height*0.05)), labels = print.val1, cex = text.cex)
				text(x = box.pos, y = (ymin - (plot.height*0.075)), labels = "+/-", cex = text.cex)
				text(x = box.pos, y = (ymin - (plot.height*0.1)), labels = print.val2, cex = text.cex)
				par(xpd = FALSE)
				}

			}
			group.pos[l] <- box.pos
			if(l < length(group.list)){
				box.pos = box.pos + 0.7
				}else{
				box.pos = box.pos + 1.3	
				}			
			} #end looping through groups
		#write a label for the group
		par(xpd = TRUE)
		text(x = mean(group.pos), y = (ymin - (plot.height*0.15)), labels = label, srt = label.srt)
		par(xpd = FALSE)
		} #end looping through group elements
	axis(2)
	legend("topleft", fill = group.cols, legend = group.labels)
	}

