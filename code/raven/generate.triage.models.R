#trait.genes must be the same ID type as the tissue
#network. For example, if the tissue network uses entrezgene
#ID, the trait genes need to be in entrezgene IDs
#fgn is the functional genomic network, can be downloaded
#from FNTM or GIANT or another source. Use download.tissue.net()
#to get a tissue-specific network from either FNTM or GIANT.
#cluster.modules indicates whether the gene list should be clustered
#into modules. If TRUE, the gene list is clustered using the fast greedy
#algorithm in igraph, and each module is run through the SVM pipeline
#separately. Functional enrichment of clusters can be acquired using
#net.module.enrichment(). If FALSE, the entire gene list will be run
#through the SVM in its entirety.
#max.cluster.size is used by iterative clustering to cluster the
#gene list into modules of a reasonable size for SVM. We recommend
#a maximum size of 500 genes. min.cluster.size denotes a cluster
#size below which the clusters will not be evaluated with the SVM.
#use.svd determines whether clusters are decomposed into eigentraits
#If TRUE, the SVD of the cluster is taken, and the optimal number
#of eigentraits to be selected is determined based which number of
#eigentraits gives the highest SVM accuracy.

generate.triage.models <- function(path = ".", project.name = "triage_project", 
trait.genes, fgn, n.trials = 100, cluster.modules = TRUE, cluster.threshold = 300, 
max.cluster.size = 500, min.cluster.size = 100, use.SVD = FALSE, verbose = TRUE, 
n.cores = 4){
	
	#create a new directory for the project
	project.dir <- file.path(path,  project.name)
	
	#if(!file.exists(project.dir)){system(paste("mkdir", project.dir))} does not work because of system
	if(!file.exists(project.dir)){dir.create(project.dir)}
	
	if(!cluster.modules){max.cluster.size = NULL}

	#=================================================================
	#generate a weighted adjacency matrix for the gene list
	#using the functional genomic network provided	
	#=================================================================
	if(verbose){cat("Generating FGN for specified genes.\n")}
	gene.fgn.file <- file.path(project.dir, "Gene.List.FGN.RData")
	if(!file.exists(gene.fgn.file)){
		gene.fgn <- tissue.adj.mat(tissue.net = fgn, gene.list = trait.genes, 
		inc.all.genes = TRUE, verbose = verbose)
		saveRDS(gene.fgn, gene.fgn.file)
		}else{
		if(verbose){cat("\tReading in existing file.\n")}
		gene.fgn <- readRDS(gene.fgn.file)
		}


	#=================================================================
	#If modules are specified, cluster the network into modules.
	#=================================================================
	
	module.file <- file.path(project.dir, "Module.Membership.RData")
	if(file.exists(module.file)){
		if(verbose){cat("\tReading in existing cluster file.\n")}
		mods <- readRDS(module.file)
		}else{		
		if(!cluster.modules || nrow(gene.fgn) < cluster.threshold){
			if(verbose){cat("Not performing clustering.\n")}
			mods <- rep(1, nrow(gene.fgn))
			u_mods <- 1
			}else{
			if(verbose){cat("Clustering functional network.\n")}
			net <- graph_from_adjacency_matrix(abs(gene.fgn[,1:nrow(gene.fgn)]), 
			weighted = TRUE, mode = "undirected")
			mods <- iter.cluster(net, max.mod.size = max.cluster.size)
			}
		saveRDS(mods, module.file)
		}
		

	u_mods = sort(unique(mods))		
		
	#=================================================================		
	#identify the modules with at least the minimum number of genes
	#=================================================================	
	mod.size <- unlist(lapply(u_mods, function(x) length(which(mods == x))))
	
	if(verbose){cat("Module sizes:\n", mod.size, "\n")}

	u_mods <- u_mods[which(mod.size > min.cluster.size)]
	if(verbose){cat("Analyzing", length(u_mods), "modules with more than", min.cluster.size, "genes.\n")}
		
	if(length(u_mods) == 0){
		stop("There are no modules greater than the minimum size specified by min.cluster.size.")
		}
		
	#======================================================================
	# go through the modules that are larger than the minimum size
	# and train 100 SVMs on each module.
	#======================================================================
	for(m in 1:length(u_mods)){
		module.dir <- file.path(project.dir, paste0("Module", m))
		if(!file.exists(module.dir)){dir.create(module.dir)}
		
		mod.locale <- which(mods == u_mods[m])
		sub.mat <- gene.fgn[mod.locale,]
		
		if(use.SVD){
			decomp.file <- file.path(module.dir, "Decomposition.RData")
			tuned.ev.file <- file.path(module.dir, "tuned.EV.RData")
			final.mat.file <- file.path(module.dir, "decomp.mat.RData")
			
			#decompose the module adjacency matrix
			if(!file.exists(decomp.file)){
				decomp <- svd(sub.mat)	
				saveRDS(decomp, decomp.file)
				}else{
				decomp <- readRDS(decomp.file)
				}	
	
			#Decide how many eigenvectors of the decomposed matrix to use
			if(!file.exists(tuned.ev.file)){
				orig.cost.list <- 10^seq(-5, 2, 1)
				#find the number of eigenvectors and 
				#cost list that maximize our classification 
				#accuracy
				tuned.ev <- tune.ev(full.mat = sub.mat, decomp.mat = decomp$v, ev.seq = c(seq(10, 40, 10), seq(50, min(300, nrow(sub.mat)), 50)))
				saveRDS(tuned.ev, tuned.ev.file)
				}else{
				tuned.ev <- readRDS(tuned.ev.file)
				}
	
			ev.to.take <- tuned.ev$best.ev.num
			tuned.cost.list <- tuned.ev$tuned.cost.list
			cost.fit.table <- tuned.ev$cost.fit.table

			pdf(file.path(module.dir, "Eigenvectors.v.Accuracy.pdf"), width = 7, height = 5)
			plot(as.numeric(rownames(cost.fit.table)), cost.fit.table[,2], type = "b", lwd = 3, xlab = "Number of Eiegenvectors", ylab = "Maximum Accuracy", axes = FALSE)
			axis(1, at = as.numeric(rownames(cost.fit.table)))
			axis(2)
			abline(v = ev.to.take)
			dev.off()
				
			if(!file.exists(final.mat.file)){
				final.mat <- t(decomp$v[,1:ev.to.take]) #tp.mat has genes in columns and eigenvectors in rows
				rownames(final.mat) <- paste("EV_", 1:ev.to.take)
				colnames(final.mat) <- colnames(sub.mat)
				# dim(final.mat)
				saveRDS(final.mat, final.mat.file)	
				}
			}else{ #end if use.SVD is TRUE
			final.mat.file <- file.path(module.dir, "non.decomp.mat.RData")
			final.mat <- sub.mat
			saveRDS(final.mat, final.mat.file)
			}			

		tuned.cost.file <- file.path(module.dir, "tuned.cost.list.RData")
		if(!file.exists(tuned.cost.file)){
			tuned.cost.list <- tune.cost.parameter(final.mat, n.trials = 1)
			saveRDS(tuned.cost.list, tuned.cost.file)
			}else{
			tuned.cost.list <- readRDS(tuned.cost.file)
			}
		

		#====================================================================
		# run the SVM models (first check to see which need to be run)
		#====================================================================
		
		svm.results <- list.files(path = module.dir, pattern = "SVM_Results")
		if(length(svm.results) > 0){
			svm.completed <- as.numeric(multi.strsplit(svm.results, c("SVM_Results", ".RData")))
			still.todo <- setdiff(1:n.trials, svm.completed)
			}else{
			still.todo <- 1:n.trials
			}

		if(length(still.todo) > 0){
			chunked.p <- chunkV(still.todo, n.cores)
			cl <- makeCluster(n.cores)
			registerDoParallel(cl)
			all.predictions <- foreach(p = 1:length(chunked.p), .export = c("tissue.svm", "cv.linear.svm")) %dopar% {
				tissue.svm(path = module.dir, tissue.mat = final.mat, num.tp = nrow(final.mat), trial.num = chunked.p[[p]], C.list = tuned.cost.list, verbose = TRUE)
				}				
			stopCluster(cl)
			}


		prediction.mat.file <- file.path(module.dir, "SVM.Prediction.Mat.RData")
		if(!file.exists(prediction.mat.file)){
			#prediction.mat holds the decision values for each gene in each model
			if(length(still.todo) == n.trials){
					prediction.mat <- Reduce("rbind", all.predictions)
					saveRDS(prediction.mat, prediction.mat.file)
					}else{
					#generate the prediction matrix from all models, 
					#and save as SVM.Prediction.Mat.RData
					genes.in.class <- predict.genes.in.class(path = module.dir, final.mat, trial.num = 1:n.trials)
					}
				}else{	
				prediction.mat <- readRDS(prediction.mat.file)
				}


					
	} #end looping through modules to build and evaluate SVM models


	
}