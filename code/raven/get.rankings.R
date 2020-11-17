#This function gets the gene rankings for all submodules
#in a base directory
#Genes are in columns, and modules are in rows. In each 
#row, genes are listed from the highest ranking on the
#left to the lowest ranking on the right

get.rankings <- function(base.dir){

    u_mods <- list.files(base.dir)
    mod.dir <- get.module.dir(base.dir, dir.table = TRUE)
    mod.names <- apply(mod.dir[[2]], 1, function(x) paste(x[1], x[2], sep = "_"))

    mod.rankings <- vector(mode = "list", length = length(mod.dir[[1]]))
    names(mod.rankings) <- mod.names

    for(i in 1:length(mod.dir[[1]])){
        prediction.mat <- readRDS(file.path(mod.dir[[1]][i], "SVM.Prediction.Mat.RData"))
        mean.svm <- colMeans(prediction.mat)
        mod.rankings[[i]] <- colnames(prediction.mat)[order(mean.svm, decreasing = TRUE)]
    }

    rank.mat <- Reduce("rbind", mod.rankings)
    rownames(rank.mat) <- mod.names
    
    return(rank.mat)
}