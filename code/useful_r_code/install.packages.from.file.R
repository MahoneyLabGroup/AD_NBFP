#This function installs packages named in a file

install.packages.from.file <- function(filename = "R_packages_2015-07-20.txt"){
	
	package.table <- as.matrix(read.table(filename, sep = "\t", stringsAsFactors = FALSE, fill = TRUE))
	for(i in 1:dim(package.table)[1]){
		install.packages(package.table[i,1])
		}
	
}