#This function compares full genome gene SVM scores across two 
#different SVM runs. If the SVMs were trained in two different 
#species, a conversion table must be provided. The conversion
#table is a two-column table that matches the IDs used in base.dir1
#in the first column to the IDs used in base.dir2 in the second.
#column
#To compare only the tops of the lists, set top.N to a relatively
#small number, like 100. 
#This function identifies the pareto front of the two scores plotted
#against each other and returns the list of genes on the front.
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

compare.genome.wide.scores <- function(base.dir1, base.dir2, conversion.table = NULL, 
top.N = NULL, plot.results = FALSE){

    require(rPref)

    #get all rankings for each group
    cat("Retrieving genome-wide SVM scores...\n")
    scores1 <- get.genome.wide.scores(base.dir1)
    scores2 <- get.genome.wide.scores(base.dir2)

    #Get the training set for each module. These are used
    #to identify which genes were originally present
    #in each training set. This information is used
    #for plotting and is returned by the function.
    cat("Retrieving training sets...\n")
    tp1 <- get.training.set(base.dir1)
    tp2 <- get.training.set(base.dir2)


    convert.ids <- function(id.class2){
        id.locale <- match(id.class2, conversion.table[,2])
        new.id <- conversion.table[id.locale,1]
        return(new.id)
    }

    #if a conversion table has been provided, translate the genes in the 
    #second list to match the first
    if(!is.null(conversion.table)){   
        trans.score2 <- convert.ids(colnames(scores2))
        colnames(scores2) <- trans.score2

        #also convert training set names
        if(!is.null(tp2)){
            converted.tp2 <- lapply(tp2, convert.ids)
            tp2 <- converted.tp2
        }
    }

    #put both matrices in the same order
    common.ids <- intersect(colnames(scores1), colnames(scores2))
    common.locale1 <- match(common.ids, colnames(scores1))
    #head(cbind(common.ids, colnames(scores1)[common.locale1]))
    common.locale2 <- match(common.ids, colnames(scores2))
    #head(cbind(common.ids, colnames(scores2)[common.locale2]))
    
    matched.scores1 <- scores1[,common.locale1]
    matched.scores2 <- scores2[,common.locale2]

    if(is.null(top.N)){top.N = ncol(matched.scores1)}

    #compare pairs of scores
    pair.cor.mat <- matrix(NA, nrow = nrow(matched.scores1), ncol = nrow(matched.scores2))
    rownames(pair.cor.mat) <- rownames(matched.scores1)
    colnames(pair.cor.mat) <- rownames(matched.scores2)

    pair.mat <- cbind(rep(1:nrow(matched.scores1), nrow(matched.scores2)),
                      rep(1:nrow(matched.scores2), each = nrow(matched.scores1)))

    all.front.genes <- vector(mode = "list", length = nrow(pair.mat))
    comp.names <- rep(NA, length(nrow(pair.mat)))

    for(i in 1:nrow(pair.mat)){
        ind1 <- pair.mat[i,1]
        ind2 <- pair.mat[i,2]
        scores1.order <- order(matched.scores1[ind1,], decreasing = TRUE)
        score.cor <- cor.test(matched.scores1[ind1,scores1.order[1:top.N]], 
        matched.scores2[ind2,scores1.order[1:top.N]])
        
        cols <- c("red" = "#d7191c", "blue" = "#2c7bb6", "purple" = "#5e3c99")
        xlab = basename(base.dir1)
        ylab = basename(base.dir2)
        colv <- rep("darkgray", top.N)
        in1 <- which(colnames(matched.scores1)[scores1.order[1:top.N]] %in% tp1[[ind1]])
        in2 <- which(colnames(matched.scores2)[scores1.order[1:top.N]] %in% tp2[[ind2]])
        colv[in1] <- cols[1]
        colv[in2] <- cols[2]
        colv[intersect(in2, in2)] <- cols[3]
        scores.x <- matched.scores1[ind1,scores1.order[1:top.N]]
        scores.y <- matched.scores2[ind2,scores1.order[1:top.N]]
        
        comp1 <- paste(xlab, rownames(matched.scores1)[ind1], sep = "_")
        comp2 <- paste(ylab, rownames(matched.scores2)[ind2], sep = "_")
        comp.names[i] <- paste(comp1, comp2, sep = "_")

        if(plot.results){
            plot(scores.x, scores.y, col = colv, xlab = xlab, ylab = ylab, pch = 16,
            main = paste(gsub("Module", "", gsub("SVM_", "", comp1)), "vs\n", 
            gsub("Module", "", gsub("SVM_", "", comp2))))
            #model <- lm(scores.y~scores.x)
            #abline(model)
        }

        #plot the pareto front
        pref <- high(scores.x) * high(scores.y)
        data.df <- data.frame(cbind(scores.x, scores.y))
        rownames(data.df) <- colnames(matched.scores1)[scores1.order[1:top.N]]

        if(plot.results){
            plot_front(data.df, pref)
        }
        on.front <- psel(data.df, pref)

        #label the pareto front genes based on their membership in 
        #each training set
        in.mod <- rep("none", nrow(on.front))
        in1 <- which(rownames(on.front) %in% tp1[[ind1]])
        in2 <- which(rownames(on.front) %in% tp2[[ind2]])
        in.both <- intersect(in1, in2)
        in.mod[in1] <- xlab
        in.mod[in2] <- ylab
        in.mod[in.both] <- "both"
        on.front <- cbind(on.front, in.mod)
        all.front.genes[[i]] <- on.front
        
        pair.cor.mat[ind1, ind2] <- score.cor$estimate
    }

    #put in one legend not in the plotting area
    if(plot.results){
        plot.new()
        plot.window(xlim = c(0,1), ylim = c(0,1))
        legend(0, 0.5, pch = 16, col = c(cols, "darkgray"), 
        legend = c(paste(xlab, "TP"), paste(ylab, "TP"), "Both TP", "Neither"))
    }


    names(all.front.genes) <- comp.names
    final.results <- list(all.front.genes, pair.cor.mat)

    return(final.results)
}

