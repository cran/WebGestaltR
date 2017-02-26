listOrganism <- function(hostName="http://www.webgestalt.org/"){
	json_data <- fromJSON(file=paste(hostName,"/data/idtypesummary.json",sep=""))
	organisms <- names(json_data)
  return(organisms)
}