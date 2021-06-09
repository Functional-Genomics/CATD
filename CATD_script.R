args <- commandArgs(trailingOnly=TRUE)

if(args[1]=='s'){
	RESULTS = self_reference_pro(args[2:length(args)])
}else if(args[1]=='r'){
	RESULTS = self_reference(args[2:length(args)])
}else if(args[1]=='c'){
	RESULTS = cross_reference(args[2:length(args)])
}else if(args[1]=='b'){
	RESULTS = bulk_2references(args[2:length(args)])
}
name = get_name(args)
saveRDS(RESULTS, paste0("RDS/",name,"rds"))
