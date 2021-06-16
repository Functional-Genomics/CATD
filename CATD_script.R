# This is the main script file to run the CATD benchmarking. 

source('framework.R')

#' @title get_name
#' @description
#' This function generates a result file name by joining all the input parameters using dot ".". 
#' @details
#' all the parameters are joined together. Only the ".rds" suffix of the rds data file is removed. 
#' @param 
#' params is a list of characters. 
#' @return
#' This function returns a string, which can be used as a result file name. 
#' @example
#' name = get_name(args)
get_name<-function(params){
    n = ''
    for( i in params){
        j = unlist(strsplit(i, "/"))
        j = gsub('.rds','',dplyr::last(j))
        n = paste0(n, j, sep = ".")
    }
    return(n)
}


# main function
args <- commandArgs(trailingOnly=TRUE)

if(args[1]=='s'){
	RESULTS = self_reference(args[2:length(args)])
}else if(args[1]=='c'){
	RESULTS = cross_reference(args[2:length(args)])
}else if(args[1]=='b'){
	RESULTS = bulk_2references(args[2:length(args)])
}else if(args[1] == 'p'){
	prepare_data(dataset, number_cells = 10000)
	
	exit(0)
}

name = get_name(args)
saveRDS(RESULTS, paste0("RDS/",name,"rds"))

print('Finished!')
