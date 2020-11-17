#This function compares full genome gene rankings across two 
#different SVM runs. If the SVMs were trained in two different 
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

compare.ranking <- function(base.dir1, base.dir2, conversion.table = NULL, top.N = NULL){

    #get all rankings for each group
    rankings1 <- get.rankings(base.dir1)
    rankings2 <- get.rankings(base.dir2)

    #if a conversion table has been provided, translate the genes in the 
    #second list to match the first
    if(!is.null(conversion.table)){
        trans.rank2 <- t(apply(rankings2, 1, function(x) 
        conversion.table[match(x, conversion.table[,2]),1]))
        rankings2 <- trans.rank2
    }

    #compare pairs of rankings
    pair.cor.mat <- matrix(NA, nrow = nrow(rankings1), ncol = nrow(rankings2))
    rownames(pair.cor.mat) <- rownames(rankings1)
    colnames(pair.cor.mat) <- rownames(rankings2)

    pair.mat <- cbind(rep(1:nrow(rankings1), nrow(rankings2)),
                      rep(1:nrow(rankings2), each = nrow(rankings1)))

    for(i in 1:nrow(pair.mat)){
        ind1 <- pair.mat[i,1]
        ind2 <- pair.mat[i,2]
        rank.cor <- compare.ranked.lists(ranked.list1 = rankings1[ind1,], 
        ranked.list2 = rankings2[ind2,], top.N = top.N)
        pair.cor.mat[ind1, ind2] <- rank.cor$spearman.cor
    }

    return(pair.cor.mat)
}

