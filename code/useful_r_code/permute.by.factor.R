#This function returns a vector of permuted indices
#based on a vector of factors. For example, if your 
#experiment has multiple strains of mice, you should
#permute the strain labels such that the correlation
#structure within strains is preserved, but the correlation
#between strain and phenotype is removed.

permute.by.factor <- function(V.to.permute){

    u_levels <- unique(V.to.permute)
    new.levels <- sample(u_levels)
    
    permV <- rep(NA, length(V.to.permute))
    #level.key <- cbind(u_levels, new.levels)

    for(i in 1:length(u_levels)){
        orig.level.idx <- which(V.to.permute == u_levels[i])
        new.level.idx <- which(V.to.permute == new.levels[i])

        #When there are more than one or the other, bootstrap
        #and sample as appropriate
        if(length(orig.level.idx) > length(new.level.idx)){
            #sample the longer orig.level.idx to match the length of the new.level.idx
            sampled.orig <- sample(orig.level.idx, length(new.level.idx))
            sampled.new <- sample(new.level.idx, length(orig.level.idx), replace = TRUE)
        }

        if(length(new.level.idx) > length(orig.level.idx)){
            #sample the longer new.level.idx to match the length of the original.level.idx
            sampled.new <- sample(new.level.idx, length(orig.level.idx))
            sampled.orig <- sample(orig.level.idx, length(new.level.idx), replace = TRUE)
        }
    
        if(length(new.level.idx) == length(orig.level.idx)){
            sampled.new <- sample(new.level.idx, length(orig.level.idx))
            sampled.orig <- sample(orig.level.idx, length(new.level.idx))
        }
    
        #after sampling, swap the indices into the opposing idx's
        permV[orig.level.idx] <- sampled.new
        permV[new.level.idx] <- sampled.orig
    }

    return(permV)

}