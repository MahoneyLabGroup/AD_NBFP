#This function does gene set expansion
#it takes a set of term genes and uses
#SVM to expand the list to other likely 
#candidates that have not been annotated 
#yet


expand.gene.set <- function(term.genes, edge.list, cluster.modules = FALSE, min.mod.size = 100, use.SVD = TRUE, max.mod.size = 300, num.trials = 100, vote.min = 0.5, path = ".", n.cores = 4){
	
	require(igraph)
	
	if(!file.exists(path)){system(paste("mkdir", path))}
	
	yeast.mat <- tissue.adj.mat(tissue.net = edge.list, gene.list = term.genes)
	tp.mat <- get.tp(yeast.mat)
	tp.net <- graph_from_adjacency_matrix(tp.mat, mode = "undirected", weighted = TRUE)


	if(cluster.modules){
		modules <- cluster_fast_greedy(tp.net)$membership
		# modules <- iter.cluster(tp.net, max.mod.size = 700)
		}else{
		modules <- rep(1, nrow(tp.mat))	
		}
		
	u_mods <- unique(modules)
	mod.size <- unlist(lapply(u_mods, function(x) length(which(modules == x))))
	
	all.mod.genes <- vector(mode = "list", length = length(u_mods))
	m.idx <- 1
	
	for(m in 1:length(u_mods)){
		if(mod.size[m] < min.mod.size){next()}
		
		module.dir <- paste0(path, '/Module', m)
		if(!file.exists(module.dir)){system(paste("mkdir", module.dir))}
		
		mod.locale <- which(modules == u_mods[m])
		mod.genes <- V(tp.net)$name[mod.locale]
				
		sub.net <- yeast.mat[mod.locale,]
	
		tuned.cost.file <- paste0(module.dir, "/tuned.cost.list.RData")
		prediction.mat.file <- paste0(module.dir, "/SVM.Prediction.Mat.RData")
		new.gene.file <- paste0(module.dir, "/SVM.Pheno.Genes.RData")
		
		if(use.SVD || nrow(sub.net) > max.mod.size){
			decomp.mat.file <- paste0(module.dir, "/decomp.mat.RData")
			decomp.file <- paste0(module.dir, "/Decomposition.RData")

			if(!file.exists(decomp.mat.file)){
				if(!file.exists(decomp.file)){
					decomp <- svd(sub.net)	
					saveRDS(decomp, decomp.file)
					}else{
					decomp <- readRDS(decomp.file)
					}	
		
				if(!file.exists(tuned.cost.file)){
					orig.cost.list <- 10^seq(-5, 2, 1)
					#find the number of eigenvectors and cost list that maximize our classification accuracy
					tuned.ev <- tune.ev(full.mat = sub.net, decomp.mat = decomp$v, ev.seq = c(seq(10, 40, 10), seq(50, min(300, nrow(sub.net)), 50)))
					saveRDS(tuned.ev, tuned.cost.file)
					}else{
					tuned.ev <- readRDS(tuned.cost.file)	
					}
	
				ev.to.take <- tuned.ev$best.ev.num
				tuned.cost.list <- tuned.ev$tuned.cost.list
				cost.fit.table <- tuned.ev$cost.fit.table
				pdf(paste0(module.dir, "/Eigenvectors.v.Accuracy.pdf"), width = 7, height = 5)
				plot(as.numeric(rownames(cost.fit.table)), cost.fit.table[,2], type = "b", lwd = 3, xlab = "Number of Eiegenvectors", ylab = "Maximum Accuracy", axes = FALSE)
				axis(1, at = as.numeric(rownames(cost.fit.table)))
				axis(2)
				abline(v = ev.to.take)
				dev.off()
					
				final.mat <- t(decomp$v[,1:ev.to.take]) #tp.mat has genes in columns and eigenvectors in rows
				rownames(final.mat) <- paste("EV_", 1:ev.to.take)
				colnames(final.mat) <- colnames(sub.net)
				# dim(final.mat)
				saveRDS(final.mat, decomp.mat.file)				
				}else{
				final.mat <- readRDS(decomp.mat.file)
				}	
			}else{
			#skip the svd if there are under 300 genes
			decomp.mat.file <- paste0(module.dir, "/non.decomp.mat.RData")
			final.mat <- sub.net
			saveRDS(final.mat, decomp.mat.file)
			if(!file.exists(tuned.cost.file)){
				tuned.cost.list <- tune.cost.parameter(final.mat, n.trials = 1)				
				saveRDS(tuned.cost.list, tuned.cost.file)
				}else{
				tuned.cost.list <- readRDS(tuned.cost.file)	
				}

			}
					
			svm.results <- list.files(path = module.dir, pattern = "SVM_Results", full.names = TRUE)
			if(length(svm.results) > 0){
				
				svm.completed <- as.numeric(matrix(multi.strsplit(svm.results, c("SVM_Results", ".RData")), ncol = 2, byrow = TRUE)[,2])
				still.todo <- setdiff(1:num.trials, svm.completed)
				}else{
				still.todo <- 1:num.trials
				}
	
		if(length(still.todo) > 0){	
	
				chunked.p <- chunkV(still.todo, n.cores)
				cl <- makeCluster(n.cores)
				registerDoParallel(cl)
				all.predictions <- foreach(p = 1:length(chunked.p), .export = c("tissue.svm", "cv.linear.svm")) %dopar% {
					tissue.svm(path = module.dir, tissue.mat = final.mat, num.tp = length(term.genes), trial.num = chunked.p[[p]], C.list = tuned.cost.list, verbose = TRUE)
					}				
				stopCluster(cl)
		
				#prediction.mat is the decision value for each gene in each model
				prediction.mat <- Reduce("rbind", all.predictions)
				saveRDS(prediction.mat, prediction.mat.file)
				}else{
				prediction.mat <- readRDS(prediction.mat.file)
				}
	
			all.votes <- apply(prediction.mat, 2, function(x) length(which(x > 0))/length(x))
			svm.pheno.genes <- colnames(prediction.mat)[which(all.votes > vote.min)]
	
			all.mod.genes[[m.idx]] <- svm.pheno.genes			
			names(all.mod.genes)[[m.idx]] <- paste0("Module", m) 
			svm.pheno.genes <- predict.genes.in.class(path = module.dir, tissue.mat = sub.net, prediction.mat = prediction.mat)
			saveRDS(svm.pheno.genes, new.gene.file)
			m.idx <- m.idx + 1
			
			}

		return(all.mod.genes)
	
	}