source.all <-
function(){

	files <- matrix(list.files(pattern = "\\.R"), ncol = 1)
	
	for(i in 1:length(files)){
		source(files[i])
		}
	
	}
