delete.candidate.results <- function(results.dir){

    really.do.it <- readline(prompt = "Are you sure you want to delete all results regarding 
    candidate genes (y/n)?\n")

    if(really.do.it == "n"){
        cat("Your results have not been touched.\n")
        return()
        }

    if(really.do.it == "y"){
        cat("Deleting all results regarding candidate genes.\n")
        cat("This does not affect trained SVM models.\n")
        
        top.file <- list.files(path = results.dir, pattern = "Candidate", full.names = TRUE)
        if(length(top.files) > 0){
            unlink(top.files)
        }

        locus.file <- list.files(path = results.dir, pattern = "Locus.Gene", full.names = TRUE)
        if(length(locus.file) > 0){
            unlink(locus.file)
            }

        dirs <- get.module.dir(results.dir)
        for(i in 1:length(dirs)){
            candidate.files <- list.files(path = dirs[i], pattern = "Candidate", full.names = TRUE)
            if(length(candidate.files) > 0){
                unlink(candidate.files)
            } 

        }
    }


}