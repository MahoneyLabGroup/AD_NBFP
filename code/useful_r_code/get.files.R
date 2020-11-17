#This script takes in a search pattern filenames and cleans
#it of unwanted patterns. For example, if you want
#eeg text files but not pdfs, use, .eeg as the find pattern
#and .pdf as the clean out pattern

get.files <- function(path = ".", want, dont.want = NULL, ignore.case = TRUE, full.names = FALSE){
		
	
	file.list <- list.files(path, full.names = full.names)
	
	all.want <- NULL
	for(i in 1:length(want)){
		if(i == 1){
			all.want <- c(all.want, grep(want[i], file.list, ignore.case = ignore.case))
			}else{
				all.want <- intersect(all.want, grep(want[i], file.list, ignore.case = ignore.case))
				}
			}
	
	#we only want files that are in all instances of all.want
	
	
	#If there are no files with what we want, just return NULL now
	if(length(all.want) == 0){
		return(NULL)
		}
	
	if(length(dont.want) > 0){
		all.dont.want <- NULL
		for(i in 1:length(dont.want)){
			all.dont.want <- c(all.dont.want, grep(dont.want[i], file.list, ignore.case = ignore.case))
			}
		
		#from the list of the potential files we want, remove those that we don't want
		final.file <- setdiff(all.want, all.dont.want)
		return(file.list[final.file])
		}else{
			return(file.list[all.want])
			}
	
	

	
}