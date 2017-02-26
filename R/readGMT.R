readGMT <- function(gmtFile){
#####Change a gmt file to a three column matrix (gene set name, gene set description and genes)#######
	if(file_extension(gmtFile)!="gmt"){
		error <- "ERROR: Please upload a file with extension 'gmt'."
		cat(error)
		return(error)
	}else{
		data <- readLines(gmtFile)
		data <- strsplit(data,"\t")
		data <- lapply(data,.toList)
		
		if(is.null(data)){
			error <- "ERROR: GMT file does not contain any gene ids. Please check the format of the GMT file."
			cat(error)
			return(error)
		}else{
			data <- do.call("rbind",data)
			return(data)
		}
	}
}


.toList <- function(data){
	if(length(data)>2){
		data <- data[!is.na(data)]
		data1 <- cbind(rep(data[1],length(data)-2),rep(data[2],length(data)-2),data[c(-1,-2)])
		return(data1)
	}else{
		return(NULL)
	}
}