#This function gets the gene rankings for all submodules
#in a base directory
#Genes are in columns, and modules are in rows. In each 
#row, genes are listed from the highest ranking on the
#left to the lowest ranking on the right

get.genome.wide.scores <- function(base.dir){

    u_mods <- list.files(base.dir)
    mod.dir <- get.module.dir(base.dir, dir.table = TRUE)
    mod.names <- apply(mod.dir[[2]], 1, function(x) paste(x[1], x[2], sep = "_"))

    mod.scores <- vector(mode = "list", length = length(mod.dir[[1]]))
    names(mod.scores) <- mod.names

    for(i in 1:length(mod.dir[[1]])){
        prediction.mat <- readRDS(file.path(mod.dir[[1]][i], "SVM.Prediction.Mat.RData"))
        mean.svm <- colMeans(prediction.mat)
        mod.scores[[i]] <- mean.svm
    }

    u_genes <- unique(unlist(lapply(mod.scores,  names)))
    num.genes <- length(u_genes)
    score.mat <- matrix(NA, nrow = length(mod.scores), ncol = num.genes)
    rownames(score.mat) <- mod.names
    colnames(score.mat) <- u_genes

    for(i in 1:length(mod.scores)){
        score.mat[i,names(mod.scores[[i]])] <- mod.scores[[i]]
    }

    return(score.mat)
}