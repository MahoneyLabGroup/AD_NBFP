#This function builds pairwise gene features for genes based on their
#connections to the 
#from the dimension-reduced
#adjacency matrix. If interaction.effects is TRUE, the pairwise feature 
#is the product of the gene vectors, i.e. AND operation, both genes need 
#high weights to generate a high interaction weight. If main.effects is 
#TRUE, the pairwise feature is the sum of the gene vectors. In this case, 
#only one of the genes needs a high weight for the weights to be high 
#(This will give us weights higher than 1. Should they be normalized?). 
#Also, aren't the individual genes already represented in A?
#These features are concatenated onto the A matrix, which is the TP-TP or 
#TN-TN adjacency matrix. 

build.pairwise.features <- function(A, X, interaction.effects = TRUE, main.effects = FALSE){

	# Compute pairwise connections to ontology term
	nX = dim(X)[1]	
	fX = NULL
	mX <- NULL
	for(i in 1:(nX - 1)){
		for(j in (i+1):nX){
			# print(c(i, j))
			#fX = rbind(fX, c(X[i, ], X[j, ]))
			if(interaction.effects){
				fX = rbind(fX, X[i, ] * X[j, ])
				}
			if(main.effects){
				mX  <- rbind(mX, X[i,]+X[j,])
				}
			}
		}

	# Linearize adjacency between candidates and bind to other features
	fX = cbind(as.matrix(A[upper.tri(A)]), fX)
	
	# Return
	return(fX)
}