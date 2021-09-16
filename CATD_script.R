# This is the main script file to run the CATD benchmarking. 

source('framework.R')


# main function
args <- commandArgs(trailingOnly=TRUE)

if(args[1]=='s'){
	RESULTS = self_reference(args[2:length(args)])
}else if(args[1]=='c'){
	RESULTS = cross_reference(args[2:length(args)])
}else if(args[1]=='s1'){
	RESULTS = self_reference_pro(args[2:length(args)])
}else if(args[1]=='c1'){
	RESULTS = cross_reference_pro(args[2:length(args)])
}else if(args[1] == 't'){
	RESULTS = combine_reference_pro(args[2:length(args)])
}else if(args[1]=='b'){
	RESULTS = bulk_2references(args[2:length(args)])
}else if(args[1] == 'p'){
	RESULTS = prepare_self(args[2:length(args)])
}else if(args[1] == 'q'){
	RESULTS = prepare_cross(args[2:length(args)])
}

name = get_name(args)

if(args[1] %in% c('p','q')){
	saveRDS(RESULTS, paste0("RDS/",name,"rds"))
}else{
	saveRDS(RESULTS, paste0("RDS2/",name,"rds"))
}

print('Finished!')
