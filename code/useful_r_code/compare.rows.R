#This function compares the rows in a matrix
#it returns the shared rows between the matrices
#the rows only in matrix1 and the rows only in matrix2
#It was adapted from an eeg script to compare probe
#sets, so that's why the word probes appears

compare.rows <- function(mat1, mat2){
		
	preserved.probes <- NULL
	preserved.index1 <- NULL
	preserved.index2 <- NULL
	probes.only.in1 <- NULL
	probes.only.in2 <- NULL

	#for each row in the first probe matrix
	for(i in 1:length(mat1[,1])){
		#look for the same probe pair in the second matrix
		pair.locale <- intersect(which(mat2[,1] == mat1[i,1]), which(mat2[,2] == mat1[i,2]))

		#if we find the same pair
		if(length(pair.locale) > 0){
			preserved.probes <- rbind(preserved.probes, mat1[i,])
			preserved.index1 <- c(preserved.index1, i)
			preserved.index2 <- c(preserved.index2, pair.locale)
			}
		}
	
	#if we the number of probe pairs that were preserved were fewer than the
	#total probe pairs in matrix 1
	#take out the pairs that were only in matrix 1
	if(length(preserved.index1) > 0 && length(preserved.index1) < length(mat1[,1])){
		probes.only.in1 <- mat1[-preserved.index1,]
		}

	if(length(preserved.index2) > 0 && length(preserved.index2) < length(mat2[,1])){
		probes.only.in2 <- mat2[-preserved.index2,]
		}

	if(length(preserved.index2) == 0 && length(preserved.index1) == 0){
		probes.only.in1 <- mat1
		probes.only.in2 <- mat2
		}
	
	results <- list(preserved.probes, probes.only.in1, probes.only.in2)
	names(results) <- c("common.pairs", "pairs.only.in.set1", "pairs.only.in.set2")
	return(results)
	
	}
