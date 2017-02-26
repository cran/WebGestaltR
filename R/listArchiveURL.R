listArchiveURL <- function(){
	archiveURL <- fread(input="http://www.webgestalt.org/archiveURL.txt",header=FALSE,sep="\t",stringsAsFactors=FALSE,colClasses="character",data.table=FALSE,showProgress=FALSE)
	
	for(i in c(1:nrow(archiveURL))){
		cat("Version:",archiveURL[i,1],"      URL:",archiveURL[i,2],sep="")
	}
	return(archiveURL)
}
