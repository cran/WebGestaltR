IDMapping <- function(organism="hsapiens",dataType="list",inputGeneFile=NULL,inputGene=NULL,sourceIdType,targetIdType, collapseMethod="mean",is_outputFile=FALSE, outputFileName="",methodType="R",hostName="http://www.webgestalt.org/"){
        
    standardId <- "entrezgene"
    goldIdType <- c("entrezgene","genesymbol","genename")
    
    largeIdList <- fread(input=paste(hostName,"/data/largeIdList.txt",sep=""),header=FALSE,sep="\t",stringsAsFactors=FALSE,colClasses="character",data.table=FALSE,showProgress=FALSE)
    largeIdList <- as.character(largeIdList[,1])
    
     #############Check organism########
    organisms <- listOrganism(hostName=hostName)
    
    if(length(which(organisms==organism))==0){
    	error <- paste("ERROR: ",organism," can not be supported.",sep="")
    	cat(error)
    	return(error)
    }
    
    ###########Check data type###########
    dataTypeA <- c("list","rnk","gmt")
    if(length(which(dataTypeA==dataType))==0){
    	error <- paste("ERROR: Data type ",dataType," can not be supported. Please select from 'list', 'rnk' and 'gmt'.",sep="")
    	cat(error)
    	return(error)
    }
    
     ############Check source id type#########
    idTypes <- listIDType(organism=organism,hostName=hostName)
    
    if(length(which(idTypes==sourceIdType))==0){
    	if(!is.null(inputGeneFile)){
    		error <- paste("ERROR: The ID type ",sourceIdType," of the uploaded file ",inputGeneFile," can not be supported for organism ",organism,".",sep="")
    	}else{
    		error <- paste("ERROR: The ID type ",sourceIdType," of the uploaded R object can not be supported for organism ",organism,".",sep="")
    	}
    	cat(error)
    	return(error)
    }
    
     ############Check target id type#########    
    if(length(which(idTypes==targetIdType))==0){
    	error <- paste("ERROR: The target ID type ",targetIdType," can not be supported for organism ",organism,".",sep="")
    	cat(error)
    	return(error)
    }
    
    
    ##########Check other parameters############
    
    collapseMethodList <- c("mean","median","min","max")
   	if(length(which(collapseMethodList==collapseMethod))==0){
   		error <- "ERROR: Please select the collapse method from mean, median, min and max."
   		cat(error)
   		return(error)
   	}    
   	
   	methodTypeA <- c("R","Python")
   	if(length(which(methodTypeA==methodType))==0){
   		error <- "ERROR: Please select 'R' or 'Python' to read the large ID mapping table (e.g. dbSNP)."
   		cat(error)
   		return(error)
   	}
    
    ###########Check input data type###############
    
    if(dataType=="gmt"){
    	if(!is.null(inputGeneFile)){
    		inputGene <- readGMT(inputGeneFile)
        if(is.character(inputGene) && length(inputGene)==1 && length(grep("ERROR:",inputGene))>0){
            return(inputGene)
        }
    	}else{
    		error <- paste("ERROR: For the functional database, please upload a 'gmt' file.",sep="")
        cat(error)
        return(error)
    	}
    }else{
    	inputGene <- formatCheck(dataType=dataType,inputGeneFile=inputGeneFile,inputGene=inputGene)
    	if(is.character(inputGene) && length(inputGene)==1 && length(grep("ERROR:",inputGene))>0){
    		return(inputGene)
    	}	
    }
    
    
	 ##########ID Mapping###############
	  
	  geneSymbol <- fread(input=paste(hostName,"/data/xref/",organism,"_genesymbol.table",sep=""),header=FALSE,sep="\t",stringsAsFactors=FALSE,colClasses="character",data.table=FALSE,showProgress=FALSE)
	 
	  geneName <- fread(input=paste(hostName,"/data/xref/",organism,"_genename.table",sep=""),header=FALSE,sep="\t",stringsAsFactors=FALSE,colClasses="character",data.table=FALSE,showProgress=FALSE)
	  
	  colnames(geneSymbol) <- c("entrezgeneS","genesymbol")
	  
	  colnames(geneName) <- c("entrezgeneS","genename")
	  
	  geneAnn <- merge(x=geneSymbol,y=geneName,by="entrezgeneS",all=TRUE)
	  
	  
	  re <- .processSourceIDMap(hostName=hostName,organism=organism,largeIdList=largeIdList,inputGeneFile=inputGeneFile,inputGene=inputGene,geneAnn=geneAnn,dataType=dataType,idType=sourceIdType,collapseMethod=collapseMethod,methodType=methodType)
		if(is.character(re) && length(re)==1 && length(grep("ERROR:",re))>0){
			return(re)
		}
		inputGene <- re$mapped
		unMapF <- re$unmapped
	  
	  ###########If the target ID type is not entrezgene, gene symbol or gene name#######
	  inputGene <- unique(inputGene)
	  
	  if(length(which(goldIdType==targetIdType))==0 && sourceIdType!=targetIdType){
	  
	  	x <- unique(inputGene[,"entrezgene"])
	  	targetF <- .IDMap(largeIdList=largeIdList,sourceIdType=targetIdType,hostName=hostName,organism=organism,inputGene=x,mapType="target")
	  	if(is.character(targetF) && length(targetF)==1 && length(grep("ERROR:",targetF))>0){
	  		return(targetF)
	  	}
	  	targetF <- targetF$mapF
	  	
	  	inputGene1 <- merge(x=inputGene,y=targetF,by="entrezgene")
	  	inputGene1 <- unique(inputGene1)
	  	unMapF2 <- setdiff(inputGene[,"userid"],inputGene1[,"userid"])
	  	unMapF <- unique(c(unMapF,unMapF2))
	  	inputGene <- inputGene1
	  	if(dataType=="list"){
	  		inputGene <- inputGene[,c(2,3,4,1,5)]
	  		colnames(inputGene) <- c("userid","genesymbol","genename","entrezgene","targetid")
	  	}
	  	if(dataType=="rnk"){
	  		inputGene <- inputGene[,c(2,3,4,1,6,5)]
	  		colnames(inputGene) <- c("userid","genesymbol","genename","entrezgene","targetid","score")
	  	}
	  	if(dataType=="gmt"){
	  		inputGene <- inputGene[,c(2:6,1,7)]
	  		colnames(inputGene) <- c("geneset","link","userid","genesymbol","genename","entrezgene","targetid")
	  	}
	  }
	  
	  #############Output#######################
	  
	  if(is_outputFile==TRUE){
	  	if(length(unMapF)>0){
	  		write.table(unMapF,file=paste(outputFileName,"_unmappedList.txt",sep=""),row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
	  	}
	  	if(dataType=="gmt"){
	  		if(length(which(goldIdType==targetIdType))==0){
	  			x <- tapply(inputGene[,ncol(inputGene)],inputGene[,1],paste,collapse="\t")
	  		}else{
	  			x <- tapply(inputGene[,targetIdType],inputGene[,1],paste,collapse="\t")
	  		}
	  		y <- unique(inputGene[,c(1,2)])
	  		x <- x[y[,1]]
	  		z <- data.frame(geneset=y[,1],link=y[,2],gene=x,stringsAsFactors=FALSE)
	  		write.table(z,file=paste(outputFileName,"_mappedList_from_",sourceIdType,"_to_",targetIdType,".gmt",sep=""),row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
	  	}else{
	  		if(dataType=="rnk"){
	  			write.table(inputGene,file=paste(outputFileName,"_mappedList_from_",sourceIdType,"_to_",targetIdType,".rnk",sep=""),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
	  		}else{
	  			write.table(inputGene,file=paste(outputFileName,"_mappedList_from_",sourceIdType,"_to_",targetIdType,".txt",sep=""),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
	  		}
	  	}
	 	}
	 	r <- list(mapped=inputGene,unmapped=unMapF)
		return(r)
}



.IDMap <- function(largeIdList,sourceIdType,hostName,organism,inputGene,mapType,methodType){
	if(length(which(largeIdList==sourceIdType))>0){
		if(methodType=="Python"){
	  	re <- .Python_IDMap(hostName,organism,sourceIdType,inputGene,mapType)
	  	if(is.character(re) && length(re)==1 && length(grep("ERROR:",re))>0){
	  		return(re)
	  	}
	  }else{
	  	re <- .R_IDMap(hostName,organism,sourceIdType,inputGene,mapType)
	  }
	}else{
		re <- .R_IDMap(hostName,organism,sourceIdType,inputGene,mapType)
	}
	return(re)
}

.R_IDMap <- function(hostName,organism,sourceIdType,inputGene,mapType){
#if mapType is source, inputGene is other ids. If mapType is target, inputGene is entrez gene ids

		  sourceFile <- fread(input=paste(hostName,"/data/xref/",organism,"_",sourceIdType,".table",sep=""),header=FALSE,sep="\t",stringsAsFactors=FALSE,colClasses="character",data.table=FALSE)
		  colnames(sourceFile) <- c("entrezgene",sourceIdType)
		  topF <- paste(sourceFile[1:5,2],collapse=",")
		  if(mapType=="source"){
				mapF <- sourceFile[sourceFile[,2] %in% inputGene,,drop=FALSE]
			}else{
				mapF <- sourceFile[sourceFile[,1] %in% inputGene,,drop=FALSE]
			}
			re <- list(topF=topF,mapF=mapF)
			return(re)
}

.Python_IDMap <- function(hostName,organism,sourceIdType,inputGene,mapType){

#for mac, install pip first sudo easy_install pip
#install pandas module, 
	re <- tryCatch(python.exec('import pandas as pd'),error=function(e){return(FALSE)})
	if(!is.null(re)){
		error <- "ERROR: Please install the python module 'pandas'."
		cat(error)
	  return(error)
	}
	python.assign('hostName',hostName)
	python.assign('organism',organism)
	python.assign('sourceIdType',sourceIdType)
	python.assign('inputGene',inputGene)
	
	python.exec('data = pd.read_csv(filepath_or_buffer=hostName+"/data/xref/"+organism+"_"+sourceIdType+".table",sep="\t",memory_map=True,header=None,names=["entrezgene",sourceIdType],dtype=str)')
	python.exec('top = data.head(n=5)')
	python.exec('top = top.to_dict()')
	topF <- python.get('top')
	topF <-do.call(cbind,topF)
	topF <- topF[,c("entrezgene",sourceIdType)]
	topF <- topF[,2]
	topF <- paste(topF,collapse=",")
	if(mapType=="source"){
		python.exec('map = data.loc[data[sourceIdType].isin(inputGene)]')
	}else{
		python.exec('map = data.loc[data["entrezgene"].isin(inputGene)]')
	}
	python.exec('map = map.to_dict()')
	mapF <- python.get('map')
	mapF <- do.call(cbind,mapF)
	mapF <- mapF[,c("entrezgene",sourceIdType)]
	re <- list(topF=topF,mapF=mapF)
	return(re)
}


.processSourceIDMap <- function(hostName,organism,largeIdList,inputGeneFile,inputGene,geneAnn,dataType,idType,collapseMethod,methodType){
	
	if(dataType=="list"){
		x <- unique(inputGene)
	}
	
	if(dataType=="rnk"){
		######Collapse the gene ids with multiple scores##########
	  x <- tapply(inputGene[,2],inputGene[,1],collapseMethod)
	  inputGene <- data.frame(id=names(x),score=as.numeric(x),stringsAsFactors=FALSE)
	  colnames(inputGene) <- c(idType,"score")
	  x <- inputGene[,1]
	}
	
	if(dataType=="gmt"){
		colnames(inputGene) <- c("geneset","link",idType)
		x <- unique(inputGene[,3])
	}
	
	mapR <- .IDMap(largeIdList=largeIdList,sourceIdType=idType,hostName=hostName,organism=organism,inputGene=x,mapType="source",methodType=methodType)
	if(is.character(mapR) && length(mapR)==1 && length(grep("ERROR:",mapR))>0){
		return(mapR)
	}
	topF <- mapR$topF
	mapF <- mapR$mapF
	
	if(nrow(mapF)==0){
		if(!is.null(inputGeneFile)){
	  		error <- paste("ERROR: The ID type of the uploaded file ",inputGeneFile," is not consistent with the input ID type ",idType,". Examples of the input ID type: ",topF,".",sep="")
	  }else{
	  		error <- paste("ERROR: The ID type of the uploaded R object is not consistent with the input ID type ",idType,". Examples of the input ID type: ",topF,".",sep="")
	  }
	  cat(error)
	  return(error)
	}
	
	unMapF <- setdiff(x,mapF[,2])
	
	colnames(mapF) <- c("entrezgeneS",idType)
	
	if(dataType=="list"){
		inputGene <- merge(x=mapF,y=geneAnn,by="entrezgeneS",all.x=TRUE)
	  inputGene <- inputGene[,c(2,3,4,1)]
	  colnames(inputGene) <- c("userid","genesymbol","genename","entrezgene")
	}
	
	if(dataType=="rnk"){
		inputGene <- merge(x=mapF,y=inputGene,by=idType,all.x=TRUE)
	  inputGene <- merge(x=inputGene,y=geneAnn,by="entrezgeneS",all.x=TRUE)
	  inputGene <- inputGene[,c(2,4,5,1,3)]
	  colnames(inputGene) <- c("userid","genesymbol","genename","entrezgene","score")
	}
	
	if(dataType=="gmt"){
	  inputGene <- merge(x=mapF,y=inputGene,by=idType,all.x=TRUE)
	  inputGene <- merge(x=inputGene,y=geneAnn,by="entrezgeneS",all.x=TRUE)
	  inputGene <- inputGene[,c(3,4,2,5,6,1)]
	  colnames(inputGene) <- c("geneset","link","userid","genesymbol","genename","entrezgene")
	}
	
	re <- list(mapped=inputGene,unmapped=unMapF)
	return(re)
}