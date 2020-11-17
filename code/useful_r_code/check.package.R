#This function checks that named packages are
#installed. If so it loads them. If not, it 
#installs them and loads them

check.package <- function(packages){

	all.packages <- rownames(installed.packages())	
	
	for(i in 1:length(packages)){
		package.locale <- which(all.packages == packages[i])
		if(length(package.locale) > 0){
			do.call("library", list(packages[i]))
			}else{
			install.packages(packages[i], repos = "http://cran.us.r-project.org")
			do.call("library", list(packages[i]))
			}
		}
	}
