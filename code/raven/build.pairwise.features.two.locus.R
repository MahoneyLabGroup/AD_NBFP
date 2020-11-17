# Feture vectors for all gene pairs derived from connection strengths between candidates (A) and connection strengths to ontology term (X). For rectangular adjacency matrices.

build.pairwise.features.two.locus <- function(A, X){

	names.X <- rownames(X)
	# Compute pairwise connections to ontology term
	nr = dim(A)[1]
	nc = dim(A)[2]
	fX = NULL
	pair.label = NULL
	for(i in 1:nr){
		# report.progress(i, nr, percent.text = 5, percent.dot = 1)
		
		chunk = cbind(A[i, ], X[names.X %in% colnames(A), ] %*% diag(X[names.X %in% rownames(A)[i], ]))
		
		rownames(chunk) = paste(rownames(A)[i], "-", names.X[names.X %in% colnames(A)])
		
		fX = rbind(fX, chunk)
	}
	
	# Return
	return(fX)
}

# build.pairwise.features.two.locus <- function(A, X, names.X){

	# # Compute pairwise connections to ontology term
	# nr = dim(A)[1]
	# nc = dim(A)[2]
	# fX = matrix(0, nrow = nr * nc, ncol = dim(X)[2] + 1)
	# pair.label = vector(mode = "character", length = nr * nc)
	# curr.idx = 1
	# for(i in 1:nr){
		# for(j in 1:nc){
			# report.progress((i - 1) * nc + j, nr * nc, percent.text = 5, percent.dot = 1)
# #			fX = rbind(fX, c(X[names.X %in% rownames(A)[i], ], X[names.X %in% colnames(A)[j], ]))
			
			# fX[curr.idx, ] = c(A[i, j], X[names.X %in% rownames(A)[i], ] * X[names.X %in% colnames(A)[j], ])
			
			# pair.label[curr.idx] = paste(rownames(A)[i], "-", colnames(A)[j], sep = "")
			# curr.idx = curr.idx + 1
			
		# }
	# }
	
	# # Linearize adjacency between candidates and bind to other features
	# #fX = cbind(as.matrix(A[upper.tri(A)]), fX)
	
	# # Return
	# output = NULL
	# output$features = fX
	# output$pair.label = pair.label
	# return(output)
# }