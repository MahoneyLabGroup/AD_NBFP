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

compare.clustering <- function(base.dir1, base.dir2, conversion.table = NULL, 
condense.results = FALSE){

    mod.dir1 <- list.files(base.dir1, full.names = TRUE)
    mod.dir2 <- list.files(base.dir2, full.names = TRUE)

    all.jaccard <- vector(mode = "list", length = length(mod.dir1))
    names(all.jaccard) <- basename(mod.dir1)

    for(i in 1:length(mod.dir1)){
        #get submodule membership for each module
        mod1.genes <- readRDS(file.path(mod.dir1[i], "Gene.List.FGN.RData"))
        mod1.mem <- readRDS(file.path(mod.dir1[i], "Module.Membership.RData"))
        u_mod1 <- sort(unique(mod1.mem))
        mem.list1 <- lapply(u_mod1, function(x) colnames(mod1.genes)[which(mod1.mem == x)])

        mod2.genes <- readRDS(file.path(mod.dir2[i], "Gene.List.FGN.RData"))
        mod2.mem <- readRDS(file.path(mod.dir2[i], "Module.Membership.RData"))
        u_mod2 <- sort(unique(mod2.mem))
        mem.list2 <- lapply(u_mod2, function(x) colnames(mod2.genes)[which(mod2.mem == x)])

        #if a conversion table has been provided, translate the genes in the 
        #second list to match the first
        if(!is.null(conversion.table)){
            trans.list2 <- lapply(mem.list2, function(x) 
            conversion.table[which(conversion.table[,2] %in% x),1])
            mem.list2 <- trans.list2
        }

        #calculate jaccard indices for each pair spanning the two lists
        pair.table <- cbind(
            rep(1:length(mem.list1), length(mem.list2)),
            rep(1:length(mem.list2), each = length(mem.list1)))

        jaccard.mat <- matrix(NA, nrow = length(mem.list1), ncol = length(mem.list2))
        rownames(jaccard.mat) <- paste(basename(base.dir1), basename(mod.dir1)[i], 1:length(mem.list1), sep = "_")
        colnames(jaccard.mat) <- paste(basename(base.dir2), basename(mod.dir2)[i], 1:length(mem.list2), sep = "_")
        for(j in 1:nrow(pair.table)){
            jaccard.mat[pair.table[j,1], pair.table[j,2]] <- 
            jaccard.ind(mem.list1[[pair.table[j,1]]], mem.list2[[pair.table[j,2]]])
            }
        all.jaccard[[i]] <- jaccard.mat
    }

    if(!condense.results){
        return(all.jaccard)
    }else{
        #put all of these matrices into one big matrix, so we can look
        #at them all together
        all.comp1 <- unlist(lapply(all.jaccard, rownames))
        all.comp2 <- unlist(lapply(all.jaccard, colnames))
        all.jaccard.mat <- matrix(NA, nrow = length(all.comp1), ncol = length(all.comp2))
        rownames(all.jaccard.mat) <- all.comp1
        colnames(all.jaccard.mat) <- all.comp2

        for(i in 1:length(all.jaccard)){
            all.jaccard.mat[rownames(all.jaccard[[i]]), colnames(all.jaccard[[i]])] <- all.jaccard[[i]]
        }
        return(all.jaccard.mat)
    }

}