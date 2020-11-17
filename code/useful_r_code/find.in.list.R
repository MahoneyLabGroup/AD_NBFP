#This function finds an element in a list
#it returns the indices of the list that
#the element is in.


find.in.list <- function(val, listX){

    look.for <- lapply(listX, function(x) which(x == val))
    found <- which(sapply(look.for, length) > 0)
    if(length(found) == 0){
        return(NA)
        }else{
        return(found)
        }

}