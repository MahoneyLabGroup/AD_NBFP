#This function is the same as compare.clustering() but does all
#pairwise comparisons, rather than just within module.
#This function compares submodule membership across two different
#clusterings. If the clusters were generated in two different 
#species, a conversion table must be provided. The conversion
#table is a two-column table that matches the IDs used in base.dir1
#in the first column to the IDs used in base.dir2 in the second.
#column
#For this function to work, the modules in the two base directories
#must be comparable and in the same order. For example, the DLPFCturquoise,
#DLPFCblue, and DLPFCbrown modules clustered in both human and mouse.
#or in different tissue networks, etc. The modules to be compared are
#subdirectories in base.dir1 and base.dir2. They will be compared in the
#order they are listed in the directory (alphabetical order)
#This function can be run after generate.triage.models()
#if condense.results is FALSE, the results are returned as a list of jaccard
#index matrices comparing submodules within modules.
#if condense.results is TRUE, these matrices are placed into a bigger matrix
#for easier comparison across all modules.

compare.clustering2 <- function(base.dir1, base.dir2, conversion.table = NULL, 
condense.results = FALSE){

    mod.dir1 <- get.module.dir(base.dir1, dir.table = TRUE)
    mod.dir2 <- get.module.dir(base.dir2, dir.table = TRUE)

    dir.names1 <- apply(mod.dir1[[2]], 1, function(x) paste(x[1], x[2], sep = "_"))
    dir.names2 <- apply(mod.dir2[[2]], 1, function(x) paste(x[1], x[2], sep = "_"))
    
    all.jaccard <- matrix(NA, nrow = length(mod.dir1[[1]]), ncol = length(mod.dir2[[1]]))
    rownames(all.jaccard) <- paste(basename(base.dir1), dir.names1, sep = "_")
    colnames(all.jaccard) <- paste(basename(base.dir2), dir.names2, sep = "_")

    #collect gene membership for each module    
    mod.entrez1 <- vector(mode = "list", length = length(mod.dir1[[1]]))
    for(i in 1:length(mod.dir1[[1]])){
        mod.mem <- read.csv(file.path(mod.dir1[[1]][i], "Module.Gene.Info.csv"), 
        stringsAsFactors = FALSE)
        mod.entrez1[[i]] <- mod.mem[,"entrezgene"]
    }

    mod.entrez2 <- vector(mode = "list", length = length(mod.dir2[[1]]))
    for(i in 1:length(mod.dir2[[1]])){
        mod.mem <- read.csv(file.path(mod.dir2[[1]][i], "Module.Gene.Info.csv"), 
        stringsAsFactors = FALSE)
        mod.entrez2[[i]] <- mod.mem[,"entrezgene"]
    }

    #if a conversion table is supplied, use it to translate module2 IDs to module1 IDs
    if(!is.null(conversion.table)){
        trans.list2 <- lapply(mod.entrez2, function(x) 
        conversion.table[which(conversion.table[,2] %in% x),1])
        mod.entrez2 <- trans.list2
    }

    
    #take out NAs from both lists
    mod.entrez1 <- lapply(mod.entrez1, function(x) x[which(!is.na(x))])
    mod.entrez2 <- lapply(mod.entrez2, function(x) x[which(!is.na(x))])

    #calculate jaccard indices between the gene IDs for all pairs of modules.
    for(i in 1:length(mod.entrez1)){
        for(j in 1:length(mod.entrez2)){
            all.jaccard[i,j] <- jaccard.ind(mod.entrez1[[i]], mod.entrez2[[j]])
        }
    }

    return(all.jaccard)

}