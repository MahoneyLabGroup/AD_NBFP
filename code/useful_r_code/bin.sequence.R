
bin.sequence <- function(v, num.sections){

	v.chunks <- vector(mode = "list", length = num.sections)
	chunk.borders <- round(segment.region(1, length(v), num.sections+1, "ends"))

	for(i in 1:(length(chunk.borders)-1)){
		if(i == length(chunk.borders)-1){
			v.chunks[[i]] <- v[chunk.borders[i]:(chunk.borders[i+1])]	
			}else{
			v.chunks[[i]] <- v[chunk.borders[i]:(chunk.borders[i+1]-1)]
			}
		}
	return(v.chunks)

	}