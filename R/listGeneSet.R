listGeneSet <- function(organism="hsapiens",hostName="http://www.webgestalt.org/"){

	json_data <- fromJSON(file=file.path(hostName,"data","genesetsummary.json"))
	idtype <- json_data[[organism]]
	name <- names(idtype)
	idList <- c()
	for(i in c(1:length(name))){
		x <- .getName(idtype[[i]])
		if(!is.null(x)){
			x <- paste(name[i],"_",x,sep="")
			idList <- c(idList,x)
		}
	}
	return(idList)
}

.getName <- function(idtype){
	return(unlist(lapply(idtype,function(e){return(e$name)})))
}