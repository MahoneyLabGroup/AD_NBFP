#This function selects a random QTL from a starting point and
#a region size. The starting point is included in the region,
#but can be anywhere in the region.
#All genomic positions and sizes should be given in bp
#The default QTL size is 5Mb
#If you know the length of the chromosome in bp, put it in max.position.
#This will guarantee that no QTL go beyond the end of the chromosome.


rnd.QTL  <- function(include.pt, qtl.size = 5e6, max.position = NULL, plot.results = FALSE, plot.label = NULL){
	
	#select a point some random distance from the point we need to include
	#but within the qtl.size.
	rnd.pt <- sample((include.pt-qtl.size):(include.pt+qtl.size), 1)
	
	#if this point is upstream of the point we need to include
	#fill out the rest of the region downstream
	if(rnd.pt < include.pt){
		end.pt <- rnd.pt + qtl.size
		start.pt <- rnd.pt
		}else{
		start.pt <- rnd.pt - qtl.size
		end.pt <- rnd.pt	
		}
		
	if(start.pt < 0){
		shift.by <- abs(start.pt)
		start.pt <- start.pt + shift.by
		end.pt <- end.pt + shift.by
		}

	if(end.pt > max.position){
		shift.by <- end.pt - max.position
		end.pt <- end.pt - shift.by
		start.pt <- start.pt - shift.by
		}
	
	if(plot.results){
		plot.new()
		plot.window(xlim = c(start.pt, end.pt), ylim = c(0,1))
		points(x = c(start.pt, include.pt, end.pt), y = c(1,1,1), type = "h", col = c("red", "blue", "red"), lwd = 3)
		axis(1);axis(2)
		mtext(plot.label)
		}

	
	return(c(start.pt, end.pt))
}