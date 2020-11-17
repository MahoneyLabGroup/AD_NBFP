#This function builds TRiAGE models for a single epistatic
#interaction
#If dimension.reduce is TRUE the adjacency matrix is transformed
#by SVD and the first n.PC principal components are used. 
#max.pairs sets the maximum number of pairs to use, so we don't
#get overwhelmed by millions of gene pairs.
#gene pairs to train on. C.list is a vector of costs for the SVM.
#There are multiple places to potentially  have multiple trials.
#1. in the selection of TN genes
#2. in the selection of gene pairs for training (not necessary if there aren't many TP genes)
# Right now we are only allowing multiple trials for sampling gene pairs.

triage.one.interaction <- function(tissue.adj, TP.genes, dimension.reduce = TRUE, n.PC = 10, max.pairs = 1e4, interaction.effects = TRUE, main.effects = FALSE, C.list = 4^seq(-4, -1, 1), verbose = TRUE){
		
	#===============================================================		
	# Generate a matrix that describes the connections between
	# genes and the phenotype of interest.
	#===============================================================		

	#Connections between phenotype genes and all other genes
	#This matrix describes the connection of all genes to the 
	#phenotype
	TP.idx <- match(TP.genes, colnames(tissue.adj))
	TP.weights <- tissue.adj[,TP.idx]
	rownames(TP.weights) <- colnames(tissue.adj)
	#look for genes that are not in the FNTM network
	#and remove these from our list
	not.found <- which(is.na(match(TP.genes, rownames(TP.weights))))
	if(length(not.found) > 0){
		TP.genes <- TP.genes[-not.found]
		}
		
	#dimension reduce the matrix of connections between the TP genes
	#(MP-related genes), and all genes in the genome if requested
	if(dimension.reduce){
		pheno.feat.mat = svd(TP.weights, nu = n.PC, nv = n.PC)
		pheno.feat = pheno.feat.mat$u
		rownames(pheno.feat) = rownames(TP.weights)
		}else{
		pheno.feat <- as.matrix(TP.weights)
		}
	
	#===============================================================		
	# Generate matrices that describe the gene-gene interactions
	# getween the TP genes and the TN genes
	#===============================================================		

	#Connections between genes sampled from outside the phenotype and 
	#all other genes
	
	#get the TP matrix from TP.weights
	TP.idx <- match(TP.genes, rownames(TP.weights))
	TP.mat <- as.matrix(TP.weights[TP.idx,])
	
	#For the TN matrix, select genes randomly from outside the MP term
	all.TN.genes <- setdiff(rownames(TP.weights), TP.genes)
	TN.genes <- sample(all.TN.genes, ncol(TP.weights))
	TN.idx <- match(TN.genes, rownames(TP.weights))
	#The TN matrix describes interactions between genes outside the MP term
	#and genes inside the MP term
	TN.mat <- as.matrix(TP.weights[TN.idx,])

	#Locate the TP and TN genes in the phenotype feature matrix
	pos.idx <- match(TP.genes, rownames(pheno.feat))
	# cbind(TP.genes, rownames(pheno.feat)[pos.idx])
	neg.idx <- match(TN.genes, rownames(pheno.feat))
	# cbind(TN.genes, rownames(pheno.feat)[neg.idx])
	
	#===============================================================		
	# Combine the gene-gene interaction matrices with the phenotype
	# information. The resulting matrix contains one row for each 
	# pair of genes in the TP.TP.mat or the TN.TN.mat. Each column
	# The columns describe how each gene pair is connected to each 
	# phenotype feature (either a TP gene or a PC of pheno.weights)
	# The final column is the network weight between the two genes.
	# These matrices thus describe how strong each gene-gene 
	# interaction is as well as how strongly each gene pair interacts 
	# with the phenotype
	#===============================================================		

	pos.feat = build.pairwise.features(A = TP.mat, X = as.matrix(pheno.feat[pos.idx, ]), interaction.effects, main.effects)
	neg.feat = build.pairwise.features(A = TN.mat, X = as.matrix(pheno.feat[neg.idx, ]), interaction.effects, main.effects)
		
	# Sub-sample pairs from each list if there are more pairs 
	# than the sample size
	if(max.pairs > nrow(pos.feat)){max.pairs <- nrow(pos.feat)}
	
	sub.pos.feat = as.matrix(pos.feat[sample(1:dim(pos.feat)[1], max.pairs), ])
	sub.neg.feat = as.matrix(neg.feat[sample(1:dim(neg.feat)[1], max.pairs), ])
	
	# Build training labels and examples
	train.lab = rbind(matrix(1, nrow = dim(sub.pos.feat)[1], 1), matrix(-1, nrow = dim(sub.neg.feat)[1], 1))
	train.examp = rbind(sub.pos.feat, sub.neg.feat)
	
	# Run SVM
	svm.out = cv.linear.svm(train.examp, train.lab, C.list = C.list)
	
	return(list("pheno.feat" = pheno.feat, "svm.out" = svm.out, "TN.genes" = TN.genes))
	
	}