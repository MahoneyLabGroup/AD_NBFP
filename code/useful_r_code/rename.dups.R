#This function appends text to duplicated values so you can
#use them as row.names
#the appended text can either be numeric or characters

rename.dups <- function(V, appended.text = c("character", "numeric")){
		appended.text <- appended.text[1]
		dups <- which(duplicated(V))
		
		if(length(dups) == 0){
			return(V)
			}
		
		for(i in 1:length(dups)){
			dup.locale <- which(V == V[dups[i]])
			if(appended.text == "character"){
				V[dup.locale] <- paste(V[dup.locale], letters[1:length(dup.locale)], sep = ".")
				}else{
				V[dup.locale] <- paste(V[dup.locale], 1:length(dup.locale), sep = ".")
				}
			}
		return(V)
		}
