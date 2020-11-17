#This function writes a list of all installed packages
#so if installing a new version of R deletes all your 
#packages, you can reinstall from a list.

write.installed.packages <- function(){
	
	all.packages <- installed.packages(.Library)[,c(1,3)]
	today <- Sys.Date()
	write.table(all.packages, file = paste("R_packages_", today, ".txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, na = "")
	
	
	
	
}