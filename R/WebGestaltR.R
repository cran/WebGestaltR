WebGestaltR <- function(enrichMethod="ORA",organism="hsapiens",enrichDatabase="geneontology_Biological_Process",enrichDatabaseFile=NULL,enrichDatabaseType=NULL,enrichDatabaseDescriptionFile=NULL,interestGeneFile=NULL,interestGene=NULL,interestGeneType=NULL,collapseMethod="mean",referenceGeneFile=NULL,referenceGene=NULL,referenceGeneType=NULL,referenceSet=NULL,minNum=10,maxNum=500,fdrMethod="BH",sigMethod="fdr",fdrThr=0.05,topThr=10,dNum=20,perNum=1000,lNum=20,is.output=TRUE,outputDirectory=getwd(),projectName=NULL,keepGSEAFolder=FALSE,methodType="R",hostName="http://www.webgestalt.org/"){
		
		
		if(is.null(projectName)){
			timeStamp <- gsub("\\.","_",as.character(as.numeric(Sys.time())))
    }else{
    	timeStamp <- projectName
    }
   
    projectDir <- paste(outputDirectory,"/Project_",timeStamp,"/",sep="")
		    
    standardId <- "entrezgene"
    
    
		
		#########Web server will input "NULL" to the R package, thus, we need to change "NULL" to NULL########
		
		enrichDatabase <- .testNull(enrichDatabase)
		enrichDatabaseFile <- .testNull(enrichDatabaseFile)
		enrichDatabaseType <- .testNull(enrichDatabaseType)
		enrichDatabaseDescriptionFile <- .testNull(enrichDatabaseDescriptionFile)
		interestGeneFile <- .testNull(interestGeneFile)
		interestGene <- .testNull(interestGene)
		interestGeneType <- .testNull(interestGeneType)
		referenceGeneFile <- .testNull(referenceGeneFile)
		referenceGene <- .testNull(referenceGene)
		referenceGeneType <- .testNull(referenceGeneType)
		referenceSet <- .testNull(referenceSet)
		
		if(!is.logical(is.output)){
			error <- "ERROR: is.output should be an R logical object (TRUE or FALSE).\n"
			cat(error)
			return(error)
		}
		
		if(is.output==TRUE){
			if(!dir.exists(outputDirectory)){
				error <- paste("ERROR: The output directory ",outputDirectory," does not exist. please change another directory or create the directory.\n",sep="")
				cat(error)
				return(error)
			}
		}else{
			outputDirectory <- getwd()
		}
			
		if(!is.logical(keepGSEAFolder)){
			error <- "ERROR: keepGSEAFolder should be an R logical object (TRUE or FALSE)!\n"
			cat(error)
			return(error)
		}
		
		methodTypeA <- c("R","Python")
   	if(length(which(methodTypeA==methodType))==0){
   		error <- "ERROR: Please select 'R' or 'Python' to read the large ID mapping table (e.g. dbSNP)."
   		cat(error)
   		return(error)
   	}
		
    #########Check method############
    existingMethods <- c("ORA","GSEA")
		
    if(length(which(existingMethods==enrichMethod))==0){
			error <- paste("ERROR: The enrichment method ",enrichMethod, " is not supported.\n",sep="")
			cat(error)
			return(error)
		}
		
		organisms <- listOrganism(hostName=hostName)
		if(length(which(organisms==organism))==0){
				if(organism!="others"){
    			error <- paste("ERROR: The organism ",organism," is not supported.\n",sep="")
    			cat(error)
					return(error)
    		}
   	}
   	
   	
   	 ########Check enrichment parameters####################
   	 
   	collapseMethodList <- c("mean","median","min","max")
   	if(length(which(collapseMethodList==collapseMethod))==0){
   		error <- "ERROR: WebGesalt only supports the following collapse methods: mean, median, min and max."
    	cat(error)
			return(error)
   	}
       
    if(!.is.wholenumber(minNum) || minNum<1){
    	error <- "ERROR: The minimum number of genes annotated to the category should be an positive integer.\n"
    	cat(error)
			return(error)
    }
    
    if(!.is.wholenumber(maxNum) || maxNum<1){
    	error <- "ERROR: The maximum number of genes annotated to the category should be an positive integer.\n"
    	cat(error)
			return(error)
    }
    
    if(minNum>=maxNum){
    	error <- "ERROR: The minimum number of genes annotated to the category should be less than the maximum number.\n"
			cat(error)
			return(error)
    }
    
    fdrMethodList <- c("holm","hochberg","hommel","bonferroni","BH","BY")
    
    if(length(which(fdrMethodList==fdrMethod))==0){
    	error <- "ERROR: WebGesalt only supports the following FDR methods: holm, hochberg, hommel, bonferroni, BH and BY."
    	cat(error)
			return(error)
    }
    
    sigMethodList <- c("fdr","top")
    if(length(which(sigMethodList==sigMethod))==0){
    	error <- "ERROR: WebGesalt only supports two methods to identify the enriched categories: 'fdr' method identifies the categories with FDR less than 'fdrThr' and 'top' method ranks all categories based on the FDR and selects the 'topThr' categories."
    	cat(error)
			return(error)
    }
    
    if(!is.numeric(fdrThr) || fdrThr<0 || fdrThr>1){
    	error <- "ERROR: The FDR threshold should be a numberic from 0 to 1."
    	cat(error)
			return(error)
    }
    
    if(!.is.wholenumber(topThr) || topThr<1){
    	error <- "ERROR: The number of top categories should be a positive integer.\n"
    	cat(error)
			return(error)
    }
    
    if(!.is.wholenumber(dNum) || dNum<0 || dNum>100){
    	error <- "ERROR: The number of enriched categories shown in the final report should be a positive integer and less than 100.\n"
    	cat(error)
			return(error)
    }
    
    if(!.is.wholenumber(perNum) || perNum<0 || perNum>10000){
    	error <- "ERROR: The number of permutation for GSEA method should be a positive integer and less than 10000.\n"
    	cat(error)
			return(error)
    }
    
    
    if(!.is.wholenumber(lNum) || lNum<1){
    	error <- "ERROR: The number of categories with leading edge gene outputted from GSEA should be a positive integer.\n"
    	cat(error)
			return(error)
    }
    
    ####################################
   	
		
		#############Check enriched database#############
		cat("Uploading the functional categories...\n")
		geneSet <- NULL    ##gene sets
		geneSetDes <- NULL ##gene set description file
		geneSetDAG <- NULL ##gene set DAG file
		geneSetNet <- NULL ##gene set network file
		
		
		if(organism!="others"){
			geneSets <- listGeneSet(organism=organism,hostName=hostName)
			if(length(which(geneSets==enrichDatabase))==0){  ###if the upload gene set name is not in the database
				if(!is.null(enrichDatabase)){
					if(enrichDatabase=="others"){    ##if the user upload their own data
						#cat("Because 'enrichDatabase' is 'others', user can upload their own gene sets using GMT file and WebGestaltR will transform ids in the gene sets to entrez ids based on the parameter 'enrichDatabaseType'!\n")
						#######Read GMT File and transform id##########
						
						if(!is.null(enrichDatabaseFile)){
							geneSet <- IDMapping(organism=organism,dataType="gmt",inputGeneFile=enrichDatabaseFile,sourceIdType=enrichDatabaseType,targetIdType=standardId,collapseMethod=collapseMethod,is_outputFile=FALSE,methodType=methodType,hostName=hostName)
							if(is.character(geneSet) && length(geneSet)==1 && length(grep("ERROR:",geneSet))>0){
    							return(geneSet)
    					}
							geneSet <- geneSet$mapped
    					geneSet <- unique(geneSet[,c(1,2,6)])
    					if(!is.null(enrichDatabaseDescriptionFile)){     ##upload description file
    							geneSetDes <- .loadEnrichDatabaseDescriptionFile(geneSet,enrichDatabaseDescriptionFile)
    							if(is.character(geneSetDes) && length(geneSetDes)==1 && length(grep("ERROR:",geneSetDes))>0){
    									return(geneSetDes)
    							}
    					}
						}else{  #no enrichDatabaseFile for 'others' enrichDatabase
							error <- "ERROR: Please upload a 'gmt' file as the functional database."
    					cat(error)
							return(error)
						}
					}else{ ##input a wrong enrichDatabase name
    				error <- paste("ERROR: ",enrichDatabase," can not be supported for organism ",organism,".",sep="")
    				cat(error)
						return(error)
    			}
    		}else{  #enrichDatabase is NULL
    			error <- paste("ERROR: Please select the functional database or select 'others' to upload the functional database.",sep="")
    			cat(error)
					return(error)
    		}
   		}else{  #input a correct enrichDatabase
   		
   			#########Read GMT file from the existing database###########
   			geneSet <- readGMT(paste(hostName,"/data/geneset/",organism,"_",enrichDatabase,".gmt",sep=""))
   			if(is.character(geneSet) && length(geneSet)==1 && length(grep("ERROR:",geneSet))>0){
   				return(geneSet)
   			}
   			#########Read the description file#############
   			geneSetDesFile <- paste(hostName,"/data/geneset/",organism,"_",enrichDatabase,".des",sep="")
   			geneSetDes <- tryCatch(fread(input=geneSetDesFile,header=FALSE,sep="\t",stringsAsFactors=FALSE,colClasses="character",data.table=FALSE,showProgress=FALSE),warning=function(e){return(NULL)},error=function(e){return(NULL)})  #####read the des file. If no des file, return NULL. For the des file, First column is the category id and the second is the description
   			
   			###########Try to load the DAG file#################
   			geneSetDAGFile <- paste(hostName,"/data/geneset/",organism,"_",enrichDatabase,".dag",sep="")
   			geneSetDAG <- tryCatch(fread(input=geneSetDAGFile,header=FALSE,sep="\t",stringsAsFactors=FALSE,colClasses="character",data.table=FALSE,showProgress=FALSE),warning=function(e){return(NULL)},error=function(e){return(NULL)})  #####read the dag file. If no dag file, return NULL. For the dag file, First column is the parent term and the second is the child term
   			
   			
   			###########Try to load the network file if the gene sets are generated from the network##########
   			geneSetNetFile <- paste(hostName,"/data/geneset/",organism,"_",enrichDatabase,".net",sep="")
   			geneSetNet <- tryCatch(fread(input=geneSetNetFile,header=FALSE,sep="\t",stringsAsFactors=FALSE,colClasses="character",data.table=FALSE,showProgress=FALSE),warning=function(e){return(NULL)},error=function(e){return(NULL)})    ####Read the net file. If no net file, return NULL. For the net file, First column is the gene set name, the second is the network related to the gene set and the third column is the related functions.
   		}
		}else{
			#########Read GMT file for other orgnisms from user files###########

			if(!is.null(enrichDatabaseFile)){
					geneSet <- readGMT(enrichDatabaseFile)
					if(is.character(geneSet) && length(geneSet)==1 && length(grep("ERROR:",geneSet))>0){
   						return(geneSet)
   				}
					if(!is.null(enrichDatabaseDescriptionFile)){     ##upload description file
    					geneSetDes <- .loadEnrichDatabaseDescriptionFile(geneSet,enrichDatabaseDescriptionFile)
    					if(is.character(geneSetDes) & length(geneSetDes)==1 && length(grep("ERROR:",geneSetDes))>0){
    						return(geneSetDes)
    					}
    			}
			}else{  ##enrichDatabaseFile is NULL
				error <- "ERROR: Please upload a 'gmt' file as the functional database for 'others' organism."
    		cat(error)
				return(error)
			}
		}
    
   
    
    ###########Check input interesting gene list###############
    cat("Uploading the gene list...\n")
    interestGene_List <- NULL
    interestingGeneMap <- NULL
    
    
    if(is.null(interestGeneFile) && is.null(interestGene)){
    	if(enrichMethod==existingMethods[1]){
    		error <- paste("ERROR: Please upload a gene list.",sep="")
    	}else{
    		error <- paste("ERROR: Please upload a ranked gene list.",sep="")
    	}
    	cat(error)
			return(error)
    }else{
    	if(organism!="others"){
    		mapRe <- .uploadGene_existingOrganism(organism=organism,enrichMethod=enrichMethod,existingMethods=existingMethods,inputGeneFile=interestGeneFile,inputGene=interestGene,geneType=interestGeneType,standardId=standardId,collapseMethod=collapseMethod,geneSet=geneSet,methodType=methodType,hostName=hostName)
    		if(is.character(mapRe) && length(mapRe)==1 && length(grep("ERROR:",mapRe))>0){
    			return(mapRe)
    		}
    		interestGene_List <- mapRe$gene_List
    		interestingGeneMap <- mapRe$geneMap
    		
    	}else{
    		interestGene_List <- .uploadGene_Others(enrichMethod=enrichMethod,existingMethods=existingMethods,inputGeneFile=interestGeneFile,inputGene=interestGene,geneSet=geneSet)
    		if(is.character(interestGene_List) && length(interestGene_List)==1 && length(grep("ERROR:",interestGene_List))>0){
    			return(interestGene_List)
    		}
    	}
    }
    
    ###################load reference gene set for SEA method##############
    referenceGene_List <- NULL
    referenceGeneMap <- NULL
    
    if(enrichMethod==existingMethods[1]){
    	 cat("Uploading the reference gene list...\n")
    	 if(is.null(referenceGeneFile) && is.null(referenceGene) && is.null(referenceSet)){
    	 		error <- paste("ERROR: Please upload a reference gene list or select an existing reference set.",sep="")
    	 		cat(error)
					return(error)
    	 }else{
    	 		if(organism!="others"){
  					if(!is.null(referenceGeneFile) || !is.null(referenceGene)){
    					mapRe <- .uploadGene_existingOrganism(organism=organism,enrichMethod=enrichMethod,existingMethods=existingMethods,inputGeneFile=referenceGeneFile,inputGene=referenceGene,geneType=referenceGeneType,standardId=standardId,collapseMethod=collapseMethod,geneSet=geneSet,methodType=methodType,hostName=hostName)
    					if(is.character(mapRe) && length(mapRe)==1 && length(grep("ERROR:",mapRe))>0){
    						return(mapRe)
    					}
    					referenceGene_List <- mapRe$gene_List
    					referenceGeneMap <- mapRe$geneMap
    				}else{ ###referenceGeneFile and referenceGene are both NULL. But referenceSet is not NULL
    					refS <- listReferenceSet(organism=organism,hostName=hostName)
    					if(length(which(refS==referenceSet))==0){
    						   error <- paste("ERROR: Please select an existing reference set.",sep="")
									 cat(error)
									 return(error)
    					}
    					referenceGene_List <- fread(input=paste(hostName,"/data/reference/",organism,"_",referenceSet,".table",sep=""),header=FALSE,sep="\t",stringsAsFactors=FALSE,colClasses="character",data.table=FALSE)
      				referenceGene_List <- as.character(unique(referenceGene_List[,1]))
    				}
    			}else{ ##For other organisms
    				if(!is.null(referenceGeneFile) || !is.null(referenceGene)){
    						referenceGene_List <- .uploadGene_Others(enrichMethod=enrichMethod,existingMethods=existingMethods,inputGeneFile=referenceGeneFile,inputGene=referenceGene,geneSet=geneSet)
    						if(is.character(referenceGene_List) && length(referenceGene_List)==1 && length(grep("ERROR:",referenceGene_List))>0){
    							return(referenceGene_List)
    						}
    				}else{
    						error <- paste("ERROR: Please upload a reference gene list for 'others' organism.",sep="")
    	 					cat(error)
								return(error)
    				}
    			}
    	 }
    	 
    	 ##compare interest gene list and reference gene list
    	if(length(intersect(interestGene_List,referenceGene_List))==0){
    			error <- paste("ERROR: All genes in the interesting gene list are not included in the reference gene list.",sep="")
    			cat(error)
					return(error)
    	}
    	
    	if(length(intersect(interestGene_List,intersect(referenceGene_List,geneSet[,3])))==0){
    			error <- paste("ERROR: All genes in the interesting gene list are not included in the overlapping genes between genes in the reference list and genes annotated to the functional database.",sep="")
    			cat(error)
					return(error)
    	}
    	
    }
    
    
    	 
     ##########Create project folder##############
     if(is.output==TRUE){
     		dir.create(projectDir)
     
     		######Summary gene annotation based on the GOSlim###########
     
    		 if(organism!="others"){
     				cat("Summary the uploaded gene list by GO Slim data...\n")
     				if(enrichMethod==existingMethods[1]){
     					summaryGene <- interestGene_List
     				}else{
     					summaryGene <- interestGene_List[,1]
     				}
     		
     				goslim_output <- paste(projectDir,"/goslim_summary_",timeStamp,sep="")
     				re <- GOSlimSummary(organism=organism,genelist=summaryGene,outputFile=goslim_output,outputType="png",hostName=hostName)
     				if(is.character(re) && length(re)==1 && length(grep("ERROR:",re))>0){
     					return(re)
     				}
     		
						write.table(interestingGeneMap$mapped,file=paste(projectDir,"/interestingGene_Mappingtable_",timeStamp,".txt",sep=""),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
						write.table(interestingGeneMap$unmapped,file=paste(projectDir,"/interestingGene_unmappedList_",timeStamp,".txt",sep=""),row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
     		}else{
     				write.table(interestGene_List,file=paste(projectDir,"/interestList_",timeStamp,".txt",sep=""),row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
     		}
     }
     
     
     #############Run enrichment analysis###################
   
   
		cat("Perform the enrichment analysis...\n")
    if(enrichMethod==existingMethods[1]){
    	enrichedSig <- .ORAEnrichment(interestGene_List,referenceGene_List,geneSet,minNum=minNum,maxNum=maxNum,fdrMethod=fdrMethod,sigMethod=sigMethod,fdrThr=fdrThr,topThr=topThr)
    }else{
    	enrichedSig <- .GSEAEnrichment(hostName,outputDirectory,timeStamp,interestGene_List,geneSet,minNum=minNum,maxNum=maxNum,sigMethod=sigMethod,fdrThr=fdrThr,topThr=topThr,perNum=perNum,lNum=lNum,is.output=is.output,keepGSEAFolder=keepGSEAFolder)
    }
    
    if(is.character(enrichedSig) && length(enrichedSig)==1 && length(grep("ERROR:",enrichedSig))>0){
    	return(enrichedSig)
    }
    
    if(!is.null(enrichedSig)){
    	if(!is.null(geneSetDes)){ #######Add extra description information###########
    			colnames(geneSetDes) <- c("geneset","description")
    			enrichedSig <- merge(x=enrichedSig,y=geneSetDes,by="geneset",all.x=TRUE)
    			enrichedSig <- enrichedSig[,c(1,ncol(enrichedSig),2:(ncol(enrichedSig)-1))]
    			enrichedSig <- enrichedSig[order(enrichedSig[,"FDR"],enrichedSig[,"PValue"]),]
    	}
    	
    	if(!is.null(interestingGeneMap) && interestGeneType!=standardId){
    			enrichedSig <- .mapUserID(enrichedSig,enrichMethod,existingMethods,interestingGeneMap)
    	}
    	
    	if(is.output==TRUE){
    		write.table(enrichedSig,file=paste(projectDir,"/enrichment_results_",timeStamp,".txt",sep=""),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
    	}
    }
				
    if(is.output==TRUE){
    	
    	##############Create report##################
    
    	cat("Generate the final report...\n")
    	.createReport(hostName=hostName,outputDirectory=outputDirectory,organism=organism,timeStamp=timeStamp,enrichMethod=enrichMethod,existingMethods=existingMethods,geneSetDes=geneSetDes,geneSetDAG=geneSetDAG,geneSetNet=geneSetNet,standardId=standardId,interestGene_List=interestGene_List,interestingGeneMap=interestingGeneMap,referenceGene_List=referenceGene_List,referenceGeneMap=referenceGeneMap,enrichedSig=enrichedSig,enrichDatabase=enrichDatabase,enrichDatabaseFile=enrichDatabaseFile,enrichDatabaseType=enrichDatabaseType,enrichDatabaseDescriptionFile=enrichDatabaseDescriptionFile,interestGeneFile=interestGeneFile,interestGene=interestGene,interestGeneType=interestGeneType,collapseMethod=collapseMethod,referenceGeneFile=referenceGeneFile,referenceGene=referenceGene,referenceGeneType=referenceGeneType,referenceSet=referenceSet,minNum=minNum,maxNum=maxNum,fdrMethod=fdrMethod,sigMethod=sigMethod,fdrThr=fdrThr,topThr=topThr,dNum=dNum,perNum=perNum,lNum=lNum)
    	
    	comm <- paste("tar -C ",projectDir," -zcvf ",projectDir,"/Project_",timeStamp,".tar.gz .",sep="")
    	system(comm,ignore.stderr=TRUE,ignore.stdout=TRUE)
    	
    	cat("Results can be found in the ",projectDir,"!",sep="")
    }
   
    return(enrichedSig)
}

.is.wholenumber <- function(x, tol = .Machine$double.eps^0.5){
		re <- tryCatch(abs(x - round(x)) < tol,warning=function(e){return(FALSE)},error=function(e){return(FALSE)})
		return(re)
}



.loadEnrichDatabaseDescriptionFile <- function(geneSet,enrichDatabaseDescriptionFile){
	if(file_extension(enrichDatabaseDescriptionFile)!="des"){
		error <- "ERROR: The description file for the functional categetories should have a 'des' extension."
		cat(error)
		return(error)
	}else{
		geneSetDes <- fread(input=enrichDatabaseDescriptionFile,header=FALSE,sep="\t",stringsAsFactors=FALSE,colClasses="character",data.table=FALSE,showProgress=FALSE)
		if(ncol(geneSetDes)!=2){
			error <- "ERROR: The description file should contain two columns: the first column is the ID of the gene sets that should be the same with the uploaded enrichment gene sets, the second column is the description of the gene sets."
			cat(error)
			return(error)
		}else{
			if(length(intersect(unique(geneSet[,1]),geneSetDes[,1]))<0.6*length(unique(geneSet[,1]))){
					error <- "ERROR: The ID types of the uploaded functional database file and description file are different. Please check the uploaded files."
					cat(error)
					return(error)
			}else{
				colnames(geneSetDes) <- c("geneset","description")
				return(geneSetDes)
			}
		}
	}
}

.uploadGene_existingOrganism <- function(organism,enrichMethod,existingMethods,inputGeneFile,inputGene,geneType,standardId,collapseMethod,geneSet,methodType,hostName){
					
					if(enrichMethod==existingMethods[1]){
						geneMap <- IDMapping(organism=organism,dataType="list",inputGeneFile=inputGeneFile,inputGene=inputGene,sourceIdType=geneType,targetIdType=standardId,collapseMethod=collapseMethod,is_outputFile=FALSE,methodType=methodType,hostName=hostName)
					}else{
						geneMap <- IDMapping(organism=organism,dataType="rnk",inputGeneFile=inputGeneFile,inputGene=inputGene,sourceIdType=geneType,targetIdType=standardId,collapseMethod=collapseMethod,is_outputFile=FALSE,methodType=methodType,hostName=hostName)
					}
					if(is.character(geneMap) && length(geneMap)==1 && length(grep("ERROR:",geneMap))>0){
						return(geneMap)
					}
					
					geneMap_MappedList <- geneMap$mapped
					if(enrichMethod==existingMethods[1]){
    					gene_List <- as.character(unique(geneMap_MappedList[,standardId]))
    					ov <- intersect(gene_List,geneSet[,3])
    					if(length(ov)==0){
    						error <- "ERROR: All genes in the uploaded gene list are not annotated to any category of the functional database."
    						cat(error)
    						return(error)
    					}
    			}else{
    					gene_List <- geneMap_MappedList[,c(standardId,"score")]
    					a <- tapply(gene_List[,2],gene_List[,1],collapseMethod,na.rm=TRUE)
    					gene_List <- data.frame(geneid=names(a),score=as.numeric(a),stringsAsFactors=FALSE)
    					
    					ov <- intersect(gene_List[,1],geneSet[,3])
    					if(length(ov)==0){
    						error <- "ERROR: All genes in the uploaded ranked gene list are not annotated to any category of the functional database."
    						cat(error)
    						return(error)
    					}
    					
    					###Because if all genes are annotated to only one category, GSEA will return the error, we need to avoid this error by reporting the error in the R#
    					gL <- geneSet[geneSet[,3] %in% gene_List[,1],,drop=FALSE]
    					gL <- tapply(gL[,3],gL[,1],length)
    					if(length(gL)==1){
    						error <- "ERROR: All genes in the ranked gene list can only annotate to one category. GSEA can not support this type of data."
    						cat(error)
    						return(error)
    					}
    					     		
					}	
    			
    			re <- list(gene_List=gene_List,geneMap=geneMap)
    			return(re)
}


.uploadGene_Others <- function(enrichMethod,existingMethods,inputGeneFile,inputGene,geneSet){
		if(enrichMethod==existingMethods[2]){
			gene_List <- formatCheck(dataType="rnk",inputGeneFile=inputGeneFile,inputGene=inputGene)
		}else{
			gene_List <- formatCheck(dataType="list",inputGeneFile=inputGeneFile,inputGene=inputGene)
		}
		
		if(is.character(gene_List) && length(gene_List)==1 && length(grep("ERROR:",gene_List))>0){
			return(gene_List)
		}
		
		if(enrichMethod==existingMethods[1]){
    		
    		ov <- intersect(gene_List,geneSet[,3])
    		if(length(ov)==0){
    				error <- "ERROR: All genes in the uploaded gene list are not annotated to any category of the functional database."
    				cat(error)
    				return(error)
    		}
    }else{
    					
    		ov <- intersect(gene_List[,1],geneSet[,3])
    		if(length(ov)==0){
    				error <- "ERROR: All genes in the uploaded ranked gene list are not annotated to any category of the functional database."
    				cat(error)
    				return(error)
    		}
    					
    		###Because if all genes are annotated to only one category, GSEA will return the error, we need to avoid this error by reporting the error in the R#
    		gL <- geneSet[geneSet[,3] %in% gene_List[,1],,drop=FALSE]
    		gL <- tapply(gL[,3],gL[,1],length)
    		if(length(gL)==1){
    				error <- "ERROR: All genes in the ranked gene list can only annotate to one category. GSEA can not support this type of data."
    				cat(error)
    				return(error)
    		}
    					     		
		}	
    
    return(gene_List)
}	


.ORAEnrichment <- function(interestGene,referenceGene,geneSet,minNum=10,maxNum=500,fdrMethod="BH",sigMethod="fdr",fdrThr=0.05,topThr=10){
	
	#before running this code, the main code has checked the overlap among interestGene, referenceGene and geneSet.
	#And this three sets should have overlapping genes.
	
	interestGene <- as.character(interestGene)
	referenceGene <- as.character(referenceGene)
	
	geneSet[,3] <- as.character(geneSet[,3])
	
	referenceGene <- intersect(referenceGene,geneSet[,3])
	
 	
 	geneSet <- geneSet[geneSet[,3] %in% referenceGene,,drop=FALSE]
 	
	geneSetNum <- tapply(geneSet[,3],geneSet[,1],length)
	geneSetNum <- geneSetNum[geneSetNum>=minNum & geneSetNum<=maxNum]
	if(length(geneSetNum)==0){
		error <- paste("ERROR: The number of annotated genes for all functional categories are not from ",minNum," to ",maxNum," for the ORA enrichment method.",sep="")
		cat(error)
		return(error)
	}
	
	geneSet <- geneSet[geneSet[,1] %in% names(geneSetNum),,drop=FALSE]
 	
	
	interestGene <- intersect(interestGene,geneSet[,3])
	interestGene <- intersect(interestGene,referenceGene)
	
	
	
	###############Enrichment analysis###################
	ra <- length(interestGene)/length(referenceGene)
	
	
	
	g <- unique(geneSet[,c(1,2)])
	enrichedResult <- data.frame(geneset=g[,1],link=g[,2],stringsAsFactors=FALSE)
	
	refG <- data.frame(geneset=names(geneSetNum),C=geneSetNum,stringsAsFactors=FALSE)
	enrichedResult <- merge(enrichedResult,refG,by="geneset",all.x=TRUE)
	
	
	intG <- geneSet[geneSet[,3] %in% interestGene,,drop=FALSE]
	intGNum <- tapply(intG[,3],intG[,1],length)
	intGNum <- data.frame(geneset=names(intGNum),O=intGNum,stringsAsFactors=FALSE)
	enrichedResult <- merge(enrichedResult,intGNum,by="geneset",all.x=TRUE)
	enrichedResult[is.na(enrichedResult[,4]),4] <- 0
	
	enrichedResult[,5] <- enrichedResult[,3]*ra
	enrichedResult[,6] <- enrichedResult[,4]/enrichedResult[,5]
	enrichedResult[is.na(enrichedResult[,6]),6] <- NA
	enrichedResult[,7] <- 1-phyper(enrichedResult[,4]-1,length(interestGene),length(referenceGene)-length(interestGene),enrichedResult[,3],lower.tail = TRUE,log.p= FALSE)
	enrichedResult[,8] <- p.adjust(enrichedResult[,7],method="BH")
	colnames(enrichedResult)[5:8] <- c("E","R","PValue","FDR")
	
	
	intG <- tapply(intG[,3],intG[,1],paste,collapse=";")
	intG <- data.frame(geneset=names(intG),overlapGene=intG,stringsAsFactors=FALSE)
	enrichedResult <- merge(enrichedResult,intG,by="geneset",all.x=TRUE)
	enrichedResult[is.na(enrichedResult[,9]),9] <- NA
	

	enrichedResult <- enrichedResult[order(enrichedResult[,8],enrichedResult[,7]),]
	if(sigMethod=="fdr"){
		enrichedResult_sig <- enrichedResult[enrichedResult[,8]<fdrThr,]
		if(nrow(enrichedResult_sig)==0){
			cat("No significant gene set was identified based on FDR ",fdrThr,"!",sep="")
			return(NULL)
		}else{
			enrichedResult_sig <- enrichedResult_sig[order(enrichedResult_sig[,"FDR"],enrichedResult_sig[,"PValue"]),]
			return(enrichedResult_sig)
		}
	}else{
		#for the top method, we only select the terms with at least one annotated interesting gene
		x <- enrichedResult[!is.na(enrichedResult[,9]),]
		x <- x[order(x[,8],x[,7]),]
		if(nrow(x)>topThr){
			enrichedResult_sig <- x[1:topThr,]
		}else{
			enrichedResult_sig <- x
		}
		enrichedResult_sig <- enrichedResult_sig[order(enrichedResult_sig[,"FDR"],enrichedResult_sig[,"PValue"]),]
		return(enrichedResult_sig)	
	}
}

.GSEAEnrichment <- function(hostName,outputDirectory,projectName,geneRankList,geneSet,collapseMethod="mean",minNum=10,maxNum=500,sigMethod="fdr",fdrThr=0.05,topThr=10,perNum=1000,lNum=20,is.output=TRUE,keepGSEAFolder=FALSE){
	 
	  GSEAJarFile <- list.files(path=outputDirectory,pattern="gsea.*\\.jar$")
	  if(length(GSEAJarFile)==0){
				cat("No GSEA java jar file can be found in the current working directory. The function will download the GSEA java jar file to the ",outputDirectory,". The copyright of the GSEA java jar file belongs to the broad institute (http://software.broadinstitute.org/gsea/index.jsp).\n",sep="")
				GSEAJarFile <- paste(outputDirectory,"/gsea.jar",sep="")
				download.file(paste(hostName,"/gsea.jar",sep=''),GSEAJarFile,mode="wb")
		}else{
				GSEAJarFile <- paste(outputDirectory,"/",GSEAJarFile,sep="")
		}
		
		if(length(grep(" ",GSEAJarFile))>0){
			GSEAJarFile <- gsub(" ","\\ ",GSEAJarFile,fixed=TRUE)
		}
		
	 projectFolder <- paste(outputDirectory,"/Project_",projectName,"/",sep="")
	 if(!dir.exists(projectFolder)){
	 		dir.create(projectFolder)
	 }
	 
	 
	 geneRankList[,1] <- as.character(geneRankList[,1])
	 geneSet[,3] <- as.character(geneSet[,3])
	 
	 
	 x <- geneSet[geneSet[,3] %in% geneRankList[,1],,drop=FALSE]
 	
	 geneSetNum <- tapply(x[,3],x[,1],length)
	 geneSetNum <- geneSetNum[geneSetNum>=minNum & geneSetNum<=maxNum]
	 if(length(geneSetNum)==0){
			error <- paste("ERROR: The number of annotated genes for all functional categories are not from ",minNum," to ",maxNum," for the GSEA enrichment method.",sep="")
			cat(error)
			return(error)
		}
	 
	 a <- tapply(geneRankList[,2],geneRankList[,1],collapseMethod,na.rm=TRUE)
   geneRankList <- data.frame(geneid=names(a),score=a,stringsAsFactors=FALSE)
	 
		gsea_rnk <- paste(projectFolder,"/Project_",projectName,"_GSEA.rnk",sep="")
		write.table(geneRankList,file=gsea_rnk,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
		
		gsea_gmt <- paste(projectFolder,"/Project_",projectName,"_GSEA.gmt",sep="")
		x <- tapply(geneSet[,3],geneSet[,1],paste,collapse="\t")
		y <- unique(geneSet[,c(1,2)])
		x <- x[y[,1]]
		g <- data.frame(geneset=y[,1],link=y[,2],gene=x,stringsAsFactors=FALSE)
		write.table(g,file=gsea_gmt,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
		
		
		outputF <- paste(projectFolder,"/Project_",projectName,"_GSEA/",sep="")
		if(!dir.exists(outputF)){
			dir.create(outputF)
		}
						
		comm <- paste("java -Xmx4096m -cp ",GSEAJarFile," xtools.gsea.GseaPreranked -gmx ",gsea_gmt," -collapse false -mode Max_probe -norm meandiv -nperm ",perNum," -rnk ",gsea_rnk," -scoring_scheme weighted -rpt_label Project_",projectName," -include_only_symbols true -make_sets true -plot_top_x ",lNum," -rnd_seed timestamp -set_max ",maxNum," -set_min ",minNum," -zip_report false -out ",outputF," -gui false",sep="")
		
		system(comm)
		
		
		gseaR <- .readGSEA(outputF)
		gseaR_pos <- gseaR$positiveP
		gseaR_neg <- gseaR$negativeP
		if(is.null(gseaR_pos) && is.null(gseaR_neg)){
			error <- "No gene set is outputted from GSEA"
			cat(error)
			return(error)
		}else{
		
			if(sigMethod=="fdr"){
				if(!is.null(gseaR_pos)){
					gseaR_posSig <- gseaR_pos[gseaR_pos[,"FDR"]<fdrThr,]
				}else{
					gseaR_posSig <- NULL
				}
				
				if(!is.null(gseaR_neg)){
					gseaR_negSig <- gseaR_neg[gseaR_neg[,"FDR"]<fdrThr,]
				}else{
					gseaR_negSig <- NULL
				}
				
				sig <- rbind(gseaR_posSig,gseaR_negSig)
				if(nrow(sig)==0){
					cat("No significant gene set was identified based on FDR ",fdrThr,"!",sep="")
				  .removeFolder(projectFolder,is.output=is.output,keepGSEAFolder=keepGSEAFolder)
					return(NULL)
				}else{
							sig <- .mappingName(sig,geneSet)
							sig <- sig[order(sig[,"FDR"],sig[,"PValue"]),]
							.removeFolder(projectFolder,is.output=is.output,keepGSEAFolder=keepGSEAFolder)
							return(sig)
				}
			}else{
				if(!is.null(gseaR_pos)){
					gseaR_pos <- gseaR_pos[order(gseaR_pos[,"FDR"],gseaR_pos[,"PValue"]),]
					if(nrow(gseaR_pos)>topThr){
						gseaR_posSig <- gseaR_pos[1:topThr,]
					}else{
						gseaR_posSig <- gseaR_pos
					}
				}else{
					gseaR_posSig <- NULL
				}
				
				if(!is.null(gseaR_neg)){
					gseaR_neg <- gseaR_neg[order(gseaR_neg[,"FDR"],gseaR_neg[,"PValue"]),]
					if(nrow(gseaR_neg)>topThr){
						gseaR_negSig <- gseaR_neg[1:topThr,]
					}else{
						gseaR_negSig <- gseaR_neg
					}
				}else{
					gseaR_negSig <- NULL
				}
				
				sig <- rbind(gseaR_posSig,gseaR_negSig)
				sig <- .mappingName(sig,geneSet)
				sig <- sig[order(sig[,"FDR"],sig[,"PValue"]),]
	
				.removeFolder(projectFolder,is.output=is.output,keepGSEAFolder=keepGSEAFolder)
				return(sig)
			}
		}
}

.readGSEA <- function(gseaFolder){
	
	positiveP <- data.frame(geneset="",link="",Size=0,ES=0,NES=0,PValue=0,FDR=0,leadingEdgeNum=0,leadingEdgeGene="",stringsAsFactors=FALSE)
	rei <- 1
	
	negativeP <- data.frame(geneset="",link="",Size=0,ES=0,NES=0,PValue=0,FDR=0,leadingEdgeNum=0,leadingEdgeGene="",stringsAsFactors=FALSE)
	nei <- 1
	
	subFile <- list.files(gseaFolder,pattern="GseaPreranked")
	subF <- paste(gseaFolder,"/",subFile,sep="")
	
	
	positiveF <- list.files(subF,"gsea_report_for_na_pos_",full.names=TRUE)
	positiveF <- positiveF[grep("xls",positiveF,fixed=TRUE)]

	negativeF <- list.files(subF,"gsea_report_for_na_neg_",full.names=TRUE)
	negativeF <- negativeF[grep("xls",negativeF,fixed=TRUE)]	
	
	
			
	positiveD <- fread(input=positiveF,header=TRUE,sep="\t",stringsAsFactors=FALSE,data.table=FALSE,showProgress=FALSE)
	
	leadingF <- list.files(subF,pattern="xls")
	
	##if no gene is annotated to any category, GSEA will not output any category
	if(nrow(positiveD)>0){
		for(j in c(1:nrow(positiveD))){
			positiveP[j,1] <- positiveD[j,1]
			positiveP[j,2] <- positiveD[j,2]
			positiveP[j,3] <- positiveD[j,4]
			positiveP[j,4] <- positiveD[j,5]
			positiveP[j,5] <- positiveD[j,6]
			positiveP[j,6] <- positiveD[j,7]
			positiveP[j,7] <- positiveD[j,8]
			lF <- paste(positiveD[j,1],".xls",sep="")
			if(length(which(leadingF==lF))>0){
				x <- .readLeadingEdge(subF,positiveD[j,1])
				positiveP[j,8] <- x$geneListNum
				positiveP[j,9] <- x$genelist
			}else{
				positiveP[j,8] <- NA
				positiveP[j,9] <- NA
			}
		}
		positiveP <- positiveP[,c(1,2,3,8,4:7,9)]
	}else{
		positiveP <- NULL
	}
		
		
		
		
  negativeD <- fread(input=negativeF,header=TRUE,sep="\t",stringsAsFactors=FALSE,data.table=FALSE,showProgress=FALSE)
	if(nrow(negativeD)>0){
		for(j in c(1:nrow(negativeD))){
			negativeP[j,1] <- negativeD[j,1]
			negativeP[j,2] <- negativeD[j,2]
			negativeP[j,3] <- negativeD[j,4]
			negativeP[j,4] <- negativeD[j,5]
			negativeP[j,5] <- negativeD[j,6]
			negativeP[j,6] <- negativeD[j,7]
			negativeP[j,7] <- negativeD[j,8]
			lF <- paste(negativeD[j,1],".xls",sep="")
			if(length(which(leadingF==lF))>0){
				x <- .readLeadingEdge(subF,negativeD[j,1])
				negativeP[j,8] <- x$geneListNum
				negativeP[j,9] <- x$genelist
			}else{
				negativeP[j,8] <- NA
				negativeP[j,9] <- NA
			}
		}
		negativeP <- negativeP[,c(1,2,3,8,4:7,9)]
	}else{
		negativeP <- NULL
	}
	
	
	re <- list(positiveP=positiveP,negativeP=negativeP)
	
	return(re)
}

.readLeadingEdge <- function(dir,sigModule){
		
		moduleF <- fread(input=paste(dir,"/",sigModule,".xls",sep=""),header=TRUE,sep="\t",stringsAsFactors=FALSE,data.table=FALSE,showProgress=FALSE)
		
		sigGene <- moduleF[moduleF[,8]=="Yes",2]
		sigGeneNum <- length(sigGene)
		sigGeneL <- paste(sigGene,collapse=";")
		re <- list(geneListNum=sigGeneNum,genelist=sigGeneL)
		return(re)

}

.removeFolder <- function(dir,is.output=TRUE,keepGSEAFolder=FALSE){
	if(is.output==FALSE){
		comm <- paste("rm -rf ",dir,sep="")
		system(comm)
	}else{
		if(keepGSEAFolder==FALSE){
			all <- list.files(path=dir,full.names=TRUE)
			folder <- all[file.info(all)$isdir]
			for(i in c(1:length(folder))){
				comm <- paste("rm -rf ",folder[i],sep="")
				system(comm)
			}
		}
	}
}

.mappingName <- function(gseaResult,geneSet){
####GSEA will automatically change all ID to upper-case and change link to others. We need to change ID back and add the link
	geneSetName <- unique(geneSet[,c(1,2)])
	geneSetName <- cbind(geneSetName,toupper(geneSetName[,1]))
	colnames(geneSetName) <- c("geneset1","link1","geneset")
	mE <- merge(gseaResult,geneSetName,by="geneset")
	mE <- mE[,c("geneset1","link1",colnames(gseaResult)[3:ncol(gseaResult)])]
	colnames(mE)[1:2] <- c("geneset","link")
	mE[,1] <- as.character(mE[,1])
	mE[,2] <- as.character(mE[,2])
	return(mE)
}




.mapUserID <- function(enrichedSig,enrichMethod,existingMethod,interestingGeneMap){
	####map entrez gene back to the original user id and add one more column to the enrichedSig
	mapgene <- interestingGeneMap$mapped[,c("userid","entrezgene")]
	if(enrichMethod==existingMethod[1]){
		gene <- enrichedSig[,"overlapGene"]
	}else{
		gene <- enrichedSig[,"leadingEdgeGene"]
	}
	
	gene <- strsplit(gene,";")
	gene <- unlist(lapply(gene,.geneM,mapgene))
	if(enrichMethod==existingMethod[1]){
		enrichedSig <- data.frame(enrichedSig,OverlapGene_UserID=gene,stringsAsFactors=FALSE)
	}else{
		enrichedSig <- data.frame(enrichedSig,leadingEdgeGene_UserID=gene,stringsAsFactors=FALSE)
	}
	return(enrichedSig)
}

.geneM <- function(genelist,mappingtable){
	if(length(genelist)==1 && is.na(genelist)){
		###The categories outputted from GSEA may not have leading edge genes
		return(NA)
	}else{
		u <- mappingtable[mappingtable[,2] %in% genelist,1]
		u <- paste(u,collapse=";")
		return(u)
	}
}


.createReport <- function(hostName,outputDirectory,organism="hsapiens",timeStamp,enrichMethod="SEA",existingMethods,geneSetDes,geneSetDAG,geneSetNet,standardId,interestGene_List,interestingGeneMap,referenceGene_List,referenceGeneMap,enrichedSig,enrichDatabase="go_bp",enrichDatabaseFile=NULL,enrichDatabaseType=NULL,enrichDatabaseDescriptionFile=NULL,interestGeneFile=NULL,interestGene=NULL,interestGeneType=NULL,collapseMethod="mean",referenceGeneFile=NULL,referenceGene=NULL,referenceGeneType=NULL,referenceSet=NULL,minNum=10,maxNum=500,fdrMethod="BH",sigMethod="fdr",fdrThr=0.05,topThr=10,dNum=20,perNum=1000,lNum=20){

	 outputHtmlFile <- paste(outputDirectory,"/Project_",timeStamp,"/Report_",timeStamp,".html",sep="")
	 
	 cat("<html>\n",file=outputHtmlFile)
	 cat("<head>\n",file=outputHtmlFile,append=TRUE)
	 cat("<title>WebGestalt (WEB-based GEne SeT AnaLysis Toolkit)</title>\n",file=outputHtmlFile,append=TRUE)
	 cat('<script type="text/javascript" src="',paste(hostName,"/js/jquery.js",sep=""),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
	 cat('<script type="text/javascript" src="',paste(hostName,"/js/jquery-ui.js",sep=""),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
	 
	 if(organism!="others"){
			 cat('<script type="text/javascript" src="',paste(hostName,"/js/changeColor.js",sep=""),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
		   cat('<script type="text/javascript" src="',paste(hostName,"/js/tableExport.js",sep=""),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
		   cat('<script type="text/javascript" src="',paste(hostName,"/js/jquery.base64.js",sep=""),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
		   cat('<script type="text/javascript" src="',paste(hostName,"/js/html2canvas.js",sep=""),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
		   cat('<script type="text/javascript" src="',paste(hostName,"/js/jspdf/libs/sprintf.js",sep=""),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
		   cat('<script type="text/javascript" src="',paste(hostName,"/js/jspdf.js",sep=""),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
		   cat('<script type="text/javascript" src="',paste(hostName,"/js/jspdf/libs/base64.js",sep=""),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
		   cat('<script type="text/javascript" src="',paste(hostName,"/js/search.js",sep=""),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
		   
		   if(!is.null(geneSetNet)){
		  	 cat('<script type="text/javascript" src="',paste(hostName,"/js/cytoscapeweb/js/min/json2.min.js",sep=""),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
		   	 cat('<script type="text/javascript" src="',paste(hostName,"/js/cytoscapeweb/js/min/AC_OETags.min.js",sep=""),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
			   cat('<script type="text/javascript" src="',paste(hostName,"/js/cytoscapeweb/js/min/cytoscapeweb.min.js",sep=""),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
			   cat('<script type="text/javascript" src="',paste(hostName,"/js/cytoscapeWeb.js",sep=""),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
			 }
			 
			 if(!is.null(geneSetDAG)){
						 			
			 	cat('<script type="text/javascript" src="',paste(hostName,"/js/d3/d3.v3.min.js",sep=""),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
			 	cat('<script type="text/javascript" src="',paste(hostName,"/js/bootstrap.min.js",sep=""),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")	 
			 	cat('<script type="text/javascript" src="',paste(hostName,"/js/jquery.tipsy.js",sep=""),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
			 	cat('<script type="text/javascript" src="',paste(hostName,"/js/jquery.contextMenu.js",sep=""),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
			 	cat('<script type="text/javascript" src="',paste(hostName,"/js/jquery.floatThead.min.js",sep=""),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
		
				cat('<script type="text/javascript" src="',paste(hostName,"/js/dagre/dagre.js",sep=""),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
			 	cat('<script type="text/javascript" src="',paste(hostName,"/js/Minimap.js",sep=""),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
				 cat('<script type="text/javascript" src="',paste(hostName,"/js/MinimapZoom.js",sep=""),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
				 cat('<script type="text/javascript" src="',paste(hostName,"/js/DirectedAcyclicGraph.js",sep=""),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
				 cat('<script type="text/javascript" src="',paste(hostName,"/js/List.js",sep=""),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
				 cat('<script type="text/javascript" src="',paste(hostName,"/js/Selectable_svv.js",sep=""),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
				 cat('<script type="text/javascript" src="',paste(hostName,"/js/Graph.js",sep=""),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
				 cat('<script type="text/javascript" src="',paste(hostName,"/js/History.js",sep=""),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
				 cat('<script type="text/javascript" src="',paste(hostName,"/js/Tooltip.js",sep=""),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
			   cat('<script type="text/javascript" src="',paste(hostName,"/js/FileSaver.js",sep=""),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
			   cat('<script type="text/javascript" src="',paste(hostName,"/js/download.js",sep=""),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
			   cat('<script type="text/javascript" src="',paste(hostName,"/js/ContextMenu.js",sep=""),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
			   cat('<script type="text/javascript" src="',paste(hostName,"/js/searchView.js",sep=""),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
			   cat('<script type="text/javascript" src="',paste(hostName,"/js/canvg/rgbcolor.js",sep=""),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
			   cat('<script type="text/javascript" src="',paste(hostName,"/js/canvg/StackBlur.js",sep=""),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
			   cat('<script type="text/javascript" src="',paste(hostName,"/js/canvg/canvg.js",sep=""),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
			 
		   
			   #the following three js should be in the result folder#
			   cat('<script type="text/javascript" src="js/xtrace_utils.js"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
			   cat('<script type="text/javascript" src="js/xtrace_graph.js"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
			   cat('<script type="text/javascript" src="js/plotDAG.js"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
			 }
		}
	 
	 cat('<script type="text/javascript" src="',paste(hostName,"/js/tooltipster.bundle.min.js",sep=""),'"></script>\n',file=outputHtmlFile,append=TRUE,sep="")
	 cat("<script>$(document).ready(function(){$('.tooltip-container').tooltipster({functionInit: function (instance, helper) {var content = $(helper.origin).find('.tooltip_content').detach();instance.content(content);}});});</script>\n",file=outputHtmlFile,append=TRUE,sep="")

   cat('<link href="',paste(hostName,"/css/jquery-ui.css",sep=""),'" rel="stylesheet" type="text/css">\n',file=outputHtmlFile,append=TRUE,sep="")
   
   if(!is.null(geneSetDAG)){
	   cat('<link href="',paste(hostName,"/css/xtrace.css",sep=""),'" rel="stylesheet" type="text/css">\n',file=outputHtmlFile,append=TRUE,sep="")
	   cat('<link href="',paste(hostName,"/css/tipsy.css",sep=""),'" rel="stylesheet" type="text/css">\n',file=outputHtmlFile,append=TRUE,sep="")
	   cat('<link href="',paste(hostName,"/css/jquery.contextMenu.css",sep=""),'" rel="stylesheet" type="text/css">\n',file=outputHtmlFile,append=TRUE,sep="")
	   cat('<link href="http://maxcdn.bootstrapcdn.com/bootstrap/3.3.5/css/bootstrap-theme.min.css" rel="stylesheet" type="text/css">\n',file=outputHtmlFile,append=TRUE,sep="")
	   cat('<link href="http://maxcdn.bootstrapcdn.com/bootstrap/3.3.5/css/bootstrap.min.css" rel="stylesheet" type="text/css">\n',file=outputHtmlFile,append=TRUE,sep="")
	 	 cat('<link href="',paste(hostName,"/css/graph.css",sep=""),'" rel="stylesheet" type="text/css">\n',file=outputHtmlFile,append=TRUE,sep="")
   }
   cat('<link href="',paste(hostName,"/css/report.css",sep=""),'" rel="stylesheet" type="text/css">\n',file=outputHtmlFile,append=TRUE,sep="")
	 cat('<link href="',paste(hostName,"/css/zebra.css",sep=""),'" rel="stylesheet" type="text/css">\n',file=outputHtmlFile,append=TRUE,sep="")
	 cat('<link href="',paste(hostName,"/css/search.css",sep=""),'" rel="stylesheet" type="text/css">\n',file=outputHtmlFile,append=TRUE,sep="")
	 cat('<link href="',paste(hostName,"/css/tooltipster.bundle.min.css",sep=""),'" rel="stylesheet" type="text/css">\n',file=outputHtmlFile,append=TRUE,sep="")

	 cat("</head>\n",file=outputHtmlFile,append=TRUE)
	 
	 cat("<body>\n",file=outputHtmlFile,append=TRUE)
	 cat('<p align="center">\n',file=outputHtmlFile,append=TRUE)
	 cat('<iframe src="',paste(hostName,"/html/head.html",sep=""),'" style="border:0px;height:180px;width:100%"></iframe>\n',file=outputHtmlFile,append=TRUE,sep="")
    
	 
	 	if(organism!="others"){
				 ####################create the Tabs################################
				 cat('<div id ="tabs"> <ul>\n',file=outputHtmlFile,append=TRUE)
				 cat('<li><a href="#summary">Summary</a></li>\n',file=outputHtmlFile,append=TRUE)
				 cat('<li><a href="#interestGene">User ID Mapping Table</a></li>\n',file=outputHtmlFile,append=TRUE)
				 cat('<li><a href="#referenceGene">GOSlim Summary</a></li>\n',file=outputHtmlFile,append=TRUE)
				 if(!is.null(enrichedSig)){
				 		cat('<li><a href="#result">Enrichment Results</a></li>\n',file=outputHtmlFile,append=TRUE)
				 }
				 cat("</ul>\n",file=outputHtmlFile,append=TRUE)
				 
				 ##############Summary#################
				 cat('<div id="summary">\n',file=outputHtmlFile,append=TRUE)
				 #############summary for the analysis###################
	 			 cat('<h3>Summary&nbsp&nbsp&nbsp<a href="Project_',timeStamp,'.tar.gz" target="_blank"><font col="red">(Result Download)</font></a></h3>\n',file=outputHtmlFile,append=TRUE,sep="")
	 			 cat("<b>Enrich method:</b> ",enrichMethod,"<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
	 			 
	 			 cat("<b>Organism:</b>",organism,"<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
			
			###Introduce enrichdatabase
			
			if(enrichDatabase=="others"){
				if(is.null(enrichDatabaseDescriptionFile)){
					cat("<b>Enrichment Categories: </b>",enrichDatabaseFile," <b>ID Type:</b>",enrichDatabaseType,"<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
				}else{
					cat("<b>Enrichment Categories: </b>",enrichDatabaseFile," <b>ID Type:</b>",enrichDatabaseType," <b>Description File:</b>",enrichDatabaseDescriptionFile,"<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
				}
			}else{
				cat("<b>Enrichment Categories: </b>",enrichDatabase,"<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
			}
			
			###Introduce interesting gene list
			
			if(!is.null(interestGeneFile)){
    		cat("<b>Interesting gene list: </b>",basename(interestGeneFile),". <b>ID type: </b>",interestGeneType,"<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
   		}else{
    		cat("<b>Interesting gene list: </b> a R object. <b> ID type: </b>",interestGeneType,"<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
    	}
    	
   		cat("The interesting gene list contains <b>",nrow(interestingGeneMap$mapped)+length(interestingGeneMap$unmapped),"</b> user IDs in which <b>",nrow(interestingGeneMap$mapped),"</b> user IDs can unambiguously map to the unique Entrez Gene IDs and <b>",length(interestingGeneMap$unmapped),"</b> user IDs were mapped to multiple Entrez Gene IDs or could not be mapped to any Entrez Gene ID. The enrichment analysis and GO Slim summary were based upon the <b>",nrow(interestingGeneMap$mapped),"</b> unique Entrez Gene IDs.<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
			
			##Introduce reference gene list
			if(enrichMethod==existingMethods[1]){
				if(!is.null(referenceGeneFile)){
					cat("<b>Reference gene list: </b>",referenceGeneFile," <b>ID type: </b>",referenceGeneType,"<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
					cat("The reference gene list contains <b>",nrow(referenceGeneMap$mapped)+length(referenceGeneMap$unmapped),"</b> IDs in which <b>",nrow(referenceGeneMap$mapped),"</b> IDs can unambiguously map to the unique Entrez Gene IDs and <b>",length(referenceGeneMap$unmapped),"</b> IDs were mapped to multiple Entrez Gene IDs or could not be mapped to any Entrez Gene ID. <b>",nrow(referenceGeneMap$mapped),"</b> unique Entrez Gene IDs are used as the reference for the enrichment analysis.<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
				}else{
					if(!is.null(referenceGene)){
						cat("<b>Reference gene list: </b>a R object. <b>ID type: </b>",referenceGeneType,"<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
						cat("The reference gene list contains <b>",nrow(referenceGeneMap$mapped)+length(referenceGeneMap$unmapped),"</b> IDs in which <b>",nrow(referenceGeneMap$mapped),"</b> IDs can unambiguously map to the unique Entrez Gene IDs and <b>",length(referenceGeneMap$unmapped),"</b> IDs were mapped to multiple Entrez Gene IDs or could not be mapped to any Entrez Gene ID. <b>",nrow(referenceGeneMap$mapped),"</b> unique Entrez Gene IDs are used as the reference for the enrichment analysis.<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
					}else{
						cat("<b>Reference gene list: </b> all mapped Entrez Gene IDs from the selected platform ",referenceSet,"<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
						cat("The reference gene list contains <b>",length(referenceGene_List),"</b> IDs which are used as the reference for the enrichment analysis.<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
					}
				}
				
			}
			
			####Parameter summary sentences#####
			.parameterSummaryDescription(enrichedSig,enrichMethod,existingMethods,dNum,outputHtmlFile,minNum,maxNum,fdrMethod,sigMethod,fdrThr,topThr,perNum,lNum)
	 			 
				 cat('</div>\n',file=outputHtmlFile,append=TRUE)  ##Summary END
							 
				 ###############Plot the mapping table######################
				 
				 cat('<div id="interestGene">\n',file=outputHtmlFile,append=TRUE)
				 cat('<div id="mappedGene">\n',file=outputHtmlFile,append=TRUE)
				 cat('<h4>Mapped User IDs</h4>\n',file=outputHtmlFile,append=TRUE)
				 tableName <- colnames(interestingGeneMap$mapped)
				 cat('<table class="zebra"><thead><tr>\n',file=outputHtmlFile,append=TRUE)
				 for(i in c(1:length(tableName))){
				 		cat('<th>',tableName[i],'</th>\n',file=outputHtmlFile,append=TRUE,sep="")
				 }
				 cat('</tr></thead>\n',file=outputHtmlFile,append=TRUE,sep="")
				 cat('<tbody>\n',file=outputHtmlFile,append=TRUE,sep="")
			
				 for(i in c(1:nrow(interestingGeneMap$mapped))){
				 	cat('<tr>\n',file=outputHtmlFile,append=TRUE)
				 	a <- interestingGeneMap$mapped[i,]
				 	for(j in c(1:length(a))){
				 		if(names(a[j])==standardId){
				 			cat('<td><a target="entrez" href="http://www.ncbi.nlm.nih.gov/sites/entrez?db=gene&cmd=Retrieve&dopt=Graphics&list_uids=',as.vector(as.matrix(a[j])),'">',as.vector(as.matrix(a[j])),'</a></td>\n',file=outputHtmlFile,append=TRUE,sep="")
				 		}else{
				 			cat('<td>',as.vector(as.matrix(a[j])),'</td>\n',file=outputHtmlFile,append=TRUE,sep="")
				 		}
				 	}
				 	cat('</tr>\n',file=outputHtmlFile,append=TRUE)
				 }
				 
				 
				 cat('</tbody></table>\n',file=outputHtmlFile,append=TRUE)
				 cat("</div>\n",file=outputHtmlFile,append=TRUE)
				 
				 cat('<div id="empty">&nbsp&nbsp&nbsp&nbsp</div>\n',file=outputHtmlFile,append=TRUE)
				 
				 cat('<div id="unmappedGene">\n',file=outputHtmlFile,append=TRUE)
				 cat("<h4>User IDs mappted to multiple Entrtez IDs or not mapped</h4>\n",file=outputHtmlFile,append=TRUE)
				 
				 cat('<table class="zebra"><thead><tr>\n',file=outputHtmlFile,append=TRUE)
				 
				 cat('<th>userid</th>\n',file=outputHtmlFile,append=TRUE)
				 
				 cat('</tr></thead>\n',file=outputHtmlFile,append=TRUE,sep="")
				 cat('<tbody>\n',file=outputHtmlFile,append=TRUE,sep="")
			
				 for(i in c(1:length(interestingGeneMap$unmapped))){
				 	cat('<tr>\n',file=outputHtmlFile,append=TRUE)
				 	cat('<td>',interestingGeneMap$unmapped[i],'</td>\n',file=outputHtmlFile,append=TRUE,sep="")
				 	cat('</tr>\n',file=outputHtmlFile,append=TRUE)
				 }
				
				 cat('</tbody></table>\n',file=outputHtmlFile,append=TRUE)
				 cat("</div>\n",file=outputHtmlFile,append=TRUE)
				 cat("</div>\n",file=outputHtmlFile,append=TRUE)
				 
				 
				 ###########GOSlim summary#########################
				 
				 cat('<div id="referenceGene">\n',file=outputHtmlFile,append=TRUE)
				 cat('<h4>GOSlim summary for the user list genes</h4>\n',file=outputHtmlFile,append=TRUE)
				 cat('Each Biological Process, Cellular Component and Molecular Function category is represented by a red, blue and green bar, repectively.<br/>\n',file=outputHtmlFile,append=TRUE)
				 cat('The height of the bar represents the number of user list genes observed in the category.\n',file=outputHtmlFile,append=TRUE)
				 goslim <- paste("goslim_summary_",timeStamp,".png",sep="")
				 
				 cat('<img src="',goslim,'" width="100%" height="100%"/>',file=outputHtmlFile,append=TRUE,sep="")
				 #cat('<iframe src="',goslim,'" width="100%" height="100%">This brower does not support PDFs. Please download the PDF to view it: <a href="',goslim,'">Download PDF</a></iframe>\n',file=outputHtmlFile,append=TRUE,sep="")
				
				 cat("</div>\n",file=outputHtmlFile,append=TRUE)
				 
				 
				 ##########GO View####################
				 if(!is.null(enrichedSig)){
				 		
				 		enrichedSig_sub <- .extractSubSig(enrichMethod,existingMethods,enrichedSig,dNum)

				 		
				 		cat('<div id="result">\n',file=outputHtmlFile,append=TRUE)
				 		cat('<div id="dag">\n',file=outputHtmlFile,append=TRUE)
			
				 		if(!is.null(geneSetDAG)){
				 			
				 			#####Create DAG structure using GOView########
				 			####Create Json file#############
				 			jF <- .createJSONFile(enrichedSig_sub,enrichMethod,existingMethods,geneSetDes,geneSetDAG,reportDir=outputHtmlFile,outputDirectory=outputDirectory,timeStamp=timeStamp)
				 			
				 			###Download the necessary Javascript file to the folder
				 			jsD <- paste(outputDirectory,"/Project_",timeStamp,"/js/",sep="")
				 			if(!dir.exists(jsD)){
				 				dir.create(jsD)
				 			}
				 			
				 			download.file(paste(hostName,"/js/plotDAG.js",sep=''),paste(jsD,"/plotDAG.js",sep=""),mode="wb")
					 		download.file(paste(hostName,"/js/xtrace_utils.js",sep=''),paste(jsD,"/xtrace_utils.js",sep=""),mode="wb")
				 			download.file(paste(hostName,"/js/xtrace_graph.js",sep=''),paste(jsD,"/xtrace_graph.js",sep=""),mode="wb")

				 			#####This is the style for saving the DAG to the file. Because google chrome does not allow to 
				 			####read the cssRule from the local css by the javascript, we need to input this style to the javascript manually
				 			style <- read.table(paste(hostName,"/css/xtrace_style.txt",sep=""),header=FALSE,sep="\t",stringsAsFactors=FALSE)
				 			style <- as.vector(as.matrix(style))
				 			style <- paste(style,collapse="WJ")
				 			style <- paste("WJ",style,"WJ",sep="")
							
							if(enrichMethod==existingMethods[2]){
								cat("Categories with red color represent positive related categories while categories with blue color represent negative related categories.<br/>\n",file=outputHtmlFile,append=TRUE)
							}
				 			cat('<div id = "graphDAG" style="height:100%;width:100%"></div>\n',file=outputHtmlFile,append=TRUE)
 						}else{
				 			#########Create table to summary enriched results###########
				 			cat('<h4>Summary of the enriched categories</h4>\n',file=outputHtmlFile,append=TRUE)
				 			
				 			if(enrichMethod==existingMethods[1]){
				 				cat('This table lists the enriched categories, number of entrez genes in the user gene list and also in the categories and FDR.<br/><br/>\n',file=outputHtmlFile,append=TRUE)
				 			}else{
				 				cat('This table lists the enriched categories, number of entrez genes in the user gene list and also in the categories and FDR. Categories with red color represent the positive related categories while categories with blue color represent the positive related categories. <br/><br/>\n',file=outputHtmlFile,append=TRUE)
				 			}
					
				 			cat('<table class="zebra"><thead><tr>\n',file=outputHtmlFile,append=TRUE)
				 			if(!is.null(geneSetDes)){
				 				cat('<th>ID</th><th>Name</th><th>#Gene</th><th>FDR</th>\n',file=outputHtmlFile,append=TRUE)
				 				if(enrichMethod==existingMethods[1]){
				 					extractSig <- enrichedSig_sub[,c(1,2,5,9)]
				 				}else{
				 					extractSig <- enrichedSig_sub[,c(1,2,4,9)]
				 				}
				 		
				 			}else{
				 				cat('<th>ID</th><th>#Gene</th><th>FDR</th>\n',file=outputHtmlFile,append=TRUE)
				 				if(enrichMethod==existingMethods[1]){
				 					extractSig <- enrichedSig_sub[,c(1,4,8)]
				 				}else{
				 					extractSig <- enrichedSig_sub[,c(1,3,8)]
				 				}
				 			}
				 	
				 			cat('</tr></thead>\n',file=outputHtmlFile,append=TRUE,sep="")
				 			cat('<tbody>\n',file=outputHtmlFile,append=TRUE,sep="")
			
				 			for(i in c(1:nrow(enrichedSig_sub))){
				 				cat('<tr>\n',file=outputHtmlFile,append=TRUE)
				 				a <- extractSig[i,]
				 				
				 				###Set the row color for GSEA###
				 				if(enrichMethod==existingMethods[1]){
				 					color <- "black"
				 				}else{
				 					if(enrichedSig_sub[i,"NES"]>0){
				 						color <- "red"
				 					}else{
				 						color <- "blue"
				 					}
				 				}
				 				
				 				for(j in c(1:length(a))){
				 					if(j==1){
				 						cat('<td><a href="#',as.vector(as.matrix(a[j])),'" style="color:',color,'">',as.vector(as.matrix(a[j])),'</a></td>\n',file=outputHtmlFile,append=TRUE,sep="")
				 					}else{
				 						if(j==length(a)){
				 							cat('<td><font color="',color,'">',format(as.vector(as.matrix(a[j])),scientific=TRUE,digits=3),'</font></td>\n',file=outputHtmlFile,append=TRUE,sep="")
				 						}else{
					 						cat('<td><font color="',color,'">',as.vector(as.matrix(a[j])),'</font></td>\n',file=outputHtmlFile,append=TRUE,sep="")
					 					}
					 				}
				 				}
					 			cat('</tr>\n',file=outputHtmlFile,append=TRUE)
							}
				 
				 	 		cat('</tbody></table>\n',file=outputHtmlFile,append=TRUE)
				 		}
				 		
				 		cat("</div>\n",file=outputHtmlFile,append=TRUE)
				 
				 		cat('<div id="empty1">&nbsp&nbsp&nbsp&nbsp</div>\n',file=outputHtmlFile,append=TRUE)
				 
				 		#############################################
				 		cat('<div id="genelist_wrap">\n',file=outputHtmlFile,append=TRUE)
				 		cat('<div id="genelist_des">\n',file=outputHtmlFile,append=TRUE)
				 		
				 		.statisticDescription(enrichMethod,existingMethods,outputHtmlFile,fdrMethod)
				 		
				 		##########Search ID##############
				 		
				 		sigID <- unique(enrichedSig_sub[,1])
				 		cat('<div class="search">',file=outputHtmlFile,append=TRUE)
				 		cat('<input id="tags" type="text" name="q" placeholder="ID Search">',file=outputHtmlFile,append=TRUE)
				 		cat('<ul id="searchData">',file=outputHtmlFile,append=TRUE)
				 		for(i in c(1:nrow(enrichedSig_sub))){
				 			if(!is.null(geneSetDes)){
				 				cat('<li><a id="aid',i-1,'" href="#',enrichedSig_sub[i,1],'">',enrichedSig_sub[i,1],' ',enrichedSig_sub[i,2],'</a></li>',file=outputHtmlFile,append=TRUE,sep="")
				 			}else{
				 				cat('<li><a id="aid',i-1,'" href="#',enrichedSig_sub[i,1],'">',enrichedSig_sub[i,1],'</a></li>',file=outputHtmlFile,append=TRUE,sep="")
				 			}
				 		}
				 		cat("</ul></div>\n",file=outputHtmlFile,append=TRUE,sep="")
				 		cat('</div>\n',file=outputHtmlFile,append=TRUE)
            ###Show the detailed gene lists########
				 		cat('<div id="genelist">\n',file=outputHtmlFile,append=TRUE)

				 		for(i in c(1:nrow(enrichedSig_sub))){
				  
				  		id = paste("Table",i,sep="")
				  		
				  		if(!is.null(geneSetDes)){
								g <- enrichedSig_sub[i,10]
							}else{
								g <- enrichedSig_sub[i,9]
							}
							
							gNa <- 1
							if(!is.na(g)){  ##some enriched terms from GSEA may not have leading edge genes
								g <- unlist(strsplit(g,";"))
								x <- interestingGeneMap$mapped
								g <- x[x[,standardId] %in% g,,drop=FALSE]
								gNa <- 0
							}
							
				  		cat('<div><nav class="cl-effect-1">',file=outputHtmlFile,append=TRUE,sep="")
				  		
				  		#####Download table#######
				  		cat("<a style=\"cursor:pointer;color:red;\" onclick=\"$('#",id,"').tableExport({type:'csv',escape:'false'});\">Download Table</a>",file=outputHtmlFile,append=TRUE,sep="")
				  		
				  		if(!is.null(geneSetNet)){
				  			####Visualize Network######
				  		  if(enrichMethod==existingMethods[1]){
				  		  	cat("<a style=\"cursor:pointer;color:blue;\" onclick=\"cytoScapeWeb('",organism,"','",enrichedSig_sub[i,1],"','",paste(g[,"genesymbol"],collapse=";"),"','');\">Network View</a>",file=outputHtmlFile,append=TRUE,sep="")
				  		  }else{
				  		  	cat("<a style=\"cursor:pointer;color:blue;\" onclick=\"cytoScapeWeb('",organism,"','",enrichedSig_sub[i,1],"','",paste(g[,"genesymbol"],collapse=";"),"','",paste(g[,"score"],collapse=";"),"');\">Network View</a>",file=outputHtmlFile,append=TRUE,sep="")
				  		  }
				  		}
				  		
							cat('</nav></div>',file=outputHtmlFile,append=TRUE,sep="")
						

				 			cat('<table id="',id,'" class="zebra"><thead><tr>\n',file=outputHtmlFile,append=TRUE,sep="")
				 			cat('<th colspan="',ncol(interestingGeneMap$mapped),'">\n',file=outputHtmlFile,append=TRUE,sep="")
				 			
				 			if(enrichMethod==existingMethods[1]){
				 				elink <- .linkModification(enrichDatabase,enrichedSig_sub[i,"link"],enrichedSig_sub[i,"overlapGene"],interestingGeneMap)
				 			}else{
				 				elink <- .linkModification(enrichDatabase,enrichedSig_sub[i,"link"],enrichedSig_sub[i,"leadingEdgeGene"],interestingGeneMap)
				 			}
							
				 			if(!is.null(geneSetDes)){###Have description file
				 				cat('<a name="',enrichedSig_sub[i,1],'"><font color="black"><a href ="',elink,'" target="_blank"><b>ID:',enrichedSig_sub[i,1],'&nbsp&nbsp&nbsp&nbsp&nbsp&nbspName:',enrichedSig_sub[i,2],'</b></a></font></a>\n',file=outputHtmlFile,append=TRUE,sep="")
				 				cat('</th></tr>\n',file=outputHtmlFile,append=TRUE)
				 				cat('<tr><th colspan="',ncol(interestingGeneMap$mapped),'">\n',file=outputHtmlFile,append=TRUE,sep="")
				 				if(enrichMethod==existingMethods[1]){
				 					cat('C=',enrichedSig_sub[i,4],'; O=',enrichedSig_sub[i,5],'; E=',round(enrichedSig_sub[i,6],digits=2),'; R=',round(enrichedSig_sub[i,7],digits=2),'; PValue=',format(enrichedSig_sub[i,8],scientific=TRUE,digits=3),'; FDR=',format(enrichedSig_sub[i,9],scientific=TRUE,digits=3),'</th></tr>\n',file=outputHtmlFile,append=TRUE,sep="")
				 				}else{
				 					cat('Size=',enrichedSig_sub[i,4],'; L=',enrichedSig_sub[i,5],'; ES=',round(enrichedSig_sub[i,6],digits=2),'; NES=',round(enrichedSig_sub[i,7],digits=2),'; Pvalue=',format(enrichedSig_sub[i,8],scientific=TRUE,digits=3),'; FDR=',format(enrichedSig_sub[i,9],scientific=TRUE,digits=3),'</th></tr>\n',file=outputHtmlFile,append=TRUE,sep="")
				 				}
				 			}else{##No description file
				 				cat('<a name="',enrichedSig_sub[i,1],'"><font color="black"><a href ="',elink,'" target="_blank"><b>ID:',enrichedSig_sub[i,1],'</b></a></font></a>\n',file=outputHtmlFile,append=TRUE,sep="")
				 				cat('</th></tr>\n',file=outputHtmlFile,append=TRUE)
				 				cat('<tr><th colspan="',ncol(interestingGeneMap$mapped),'">\n',file=outputHtmlFile,append=TRUE,sep="")
				 				if(enrichMethod==existingMethods[1]){
				 					cat('C=',enrichedSig_sub[i,3],'; O=',enrichedSig_sub[i,4],'; E=',round(enrichedSig_sub[i,5],digits=2),'; R=',round(enrichedSig_sub[i,6],digits=2),'; PValue=',format(enrichedSig_sub[i,7],scientific=TRUE,digits=3),'; FDR=',format(enrichedSig_sub[i,8],scientific=TRUE,digits=3),'</th></tr>\n',file=outputHtmlFile,append=TRUE,sep="")
				 				}else{
				 					cat('Size=',enrichedSig_sub[i,3],'; L=',enrichedSig_sub[i,4],'; ES=',round(enrichedSig_sub[i,5],digits=2),'; NES=',round(enrichedSig_sub[i,6],digits=2),'; PValue=',format(enrichedSig_sub[i,7],scientific=TRUE,digits=3),'; FDR=',format(enrichedSig_sub[i,8],scientific=TRUE,digits=3),'</th></tr>\n',file=outputHtmlFile,append=TRUE,sep="")
				 				}
				 			}
				 	
				 			tableName <- colnames(interestingGeneMap$mapped)
				 			cat('<tr>\n',file=outputHtmlFile,append=TRUE)
				 			for(j in c(1:length(tableName))){
				 				cat('<th>',tableName[j],'</th>\n',file=outputHtmlFile,append=TRUE,sep="")
				 			}
				 			cat('</tr></thead>\n',file=outputHtmlFile,append=TRUE,sep="")
				 			cat('<tbody>\n',file=outputHtmlFile,append=TRUE,sep="")
							
					
							if(gNa==0){ ####GSEA method may have significant categories without leading edge genes. This may also happen when users use 'TOP' method to select the enriched categories. If all categories do not annotate to any interesting genes, all selcted 'TOP' categories should only have 'NA' genes
					
				 				for(j in c(1:nrow(g))){
				 					cat('<tr>\n',file=outputHtmlFile,append=TRUE)
				 					a <- g[j,]
				 					for(k in c(1:length(a))){
				 						if(names(a[k])==standardId){
				 							cat('<td><a target="entrez" href="http://www.ncbi.nlm.nih.gov/sites/entrez?db=gene&cmd=Retrieve&dopt=Graphics&list_uids=',as.vector(as.matrix(a[k])),'">',as.vector(as.matrix(a[k])),'</a></td>\n',file=outputHtmlFile,append=TRUE,sep="")
				 						}else{
				 							cat('<td>',as.vector(as.matrix(a[k])),'</td>\n',file=outputHtmlFile,append=TRUE,sep="")
				 						}
				 					}
				 					cat('</tr>\n',file=outputHtmlFile,append=TRUE)
				 				}
				 			}
				 			cat('</tbody></table><br/><br/>\n',file=outputHtmlFile,append=TRUE)
				 		}
				 
				 
				 		cat("</div>\n",file=outputHtmlFile,append=TRUE)				 
				 		cat("</div>\n",file=outputHtmlFile,append=TRUE)
				 		cat("</div>\n",file=outputHtmlFile,append=TRUE)
			
				 }
					
				#######################################################
					
				 cat("</div>\n",file=outputHtmlFile,append=TRUE)
				 if(!is.null(enrichedSig) && !is.null(geneSetDAG)){
				 		
				 		jID <- paste("Project_",timeStamp,"_.json",sep="")
				 		cat("<script>plotDAG('",jID,"','",jF,"','",style,"');</script>",file=outputHtmlFile,append=TRUE,sep="")
				 }
				 cat('<script>$( "#tabs" ).tabs();</script>',file=outputHtmlFile,append=TRUE)
			
				
		}else{
			###########Organism is others. No mapping information#############
			 #############summary for the analysis###################
	 		cat('<h3>Summary&nbsp&nbsp&nbsp<a href="Project_',timeStamp,'.tar.gz" target="_blank"><font col="red">(Result Download)</font></a></h3>\n',file=outputHtmlFile,append=TRUE,sep="")
	 		cat("<b>Enrich method:</b> ",enrichMethod,"<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
	 
			
			cat("<b>Organism</b>: Others<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
			if(is.null(enrichDatabaseDescriptionFile)){
				cat("<b>Enrichment Categories: </b>",enrichDatabaseFile,"<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
			}else{
				cat("<b>Enrichment Categories: </b>",enrichDatabaseFile," <b>Description File: </b>",enrichDatabaseDescriptionFile,"<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
			}
			
			if(!is.null(interestGeneFile)){
				cat("<b>Interesting List: </b>",basename(interestGeneFile)," <br/>\n",file=outputHtmlFile,append=TRUE,sep="")
			}else{
				cat("<b>Interesting List: </b> a R object <br/>\n",file=outputHtmlFile,append=TRUE,sep="")
			}
			
			if(enrichMethod==existingMethods[1]){
				cat("The file contains <b>",length(interestGene_List),"</b> user IDs (no ID mapping). All these IDs were used to perform the enrichment analysis.<br/>\n ",file=outputHtmlFile,append=TRUE,sep="")
				if(!is.null(referenceGeneFile)){
					cat("<b>Reference List: </b>",referenceGeneFile,"<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
				}else{
					cat("<b>Reference List: </b>a R object<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
				}
				cat("<b>Total number of reference IDs: </b>",length(referenceGene_List),"<br/>\n",file=outputHtmlFile,append=TRUE,sep="")

			}else{
				cat("The file contains <b>",nrow(interestGene_List),"</b> user IDs (no ID mapping). All these IDs were used to perform the enrichment analysis.<br/>\n ",file=outputHtmlFile,append=TRUE,sep="")
			}
			
			####Parameter summary sentences#####
			.parameterSummaryDescription(enrichedSig,enrichMethod,existingMethods,dNum,outputHtmlFile,minNum,maxNum,fdrMethod,sigMethod,fdrThr,topThr,perNum,lNum)
			
			##############Enrich Result################
			
			
			enrichedSig_sub <- .extractSubSig(enrichMethod,existingMethods,enrichedSig,dNum)

				
				.statisticDescription(enrichMethod,existingMethods,outputHtmlFile,fdrMethod)
				cat('<div width="100%">\n',file=outputHtmlFile,append=TRUE)
				cat('<table class="zebra"><thead><tr>\n',file=outputHtmlFile,append=TRUE)
				if(enrichMethod==existingMethods[1]){
				 		
				 		if(!is.null(geneSetDes)){
				 			cat('<th width="10%">ID</th><th width="20%">Name</th><th width="30%">Statistic</th><th width="40%">Genes</th></tr></thead>\n',file=outputHtmlFile,append=TRUE)
				 			extractSig <- data.frame(id=enrichedSig_sub[,1],name=enrichedSig_sub[,2],link=enrichedSig_sub[,3],statistic=paste("C=",enrichedSig_sub[,4],"; O=",enrichedSig_sub[,5],"; E=",round(enrichedSig_sub[,6],digits=2),"; R=",round(enrichedSig_sub[,7],digits=2),"; PValue=",format(enrichedSig_sub[,8],scientific=TRUE,digits=3),"; FDR=",format(enrichedSig_sub[,9],scientific=TRUE,digits=3),sep=""),genes=gsub(","," ",enrichedSig_sub[,10]),stringsAsFactors=FALSE)
				 		}else{
				 	  	cat('<th width="10%">ID</th width="30%"><th>Statistic</th><th width="60%">Genes</th></tr></thead>\n',file=outputHtmlFile,append=TRUE)
				 			extractSig <- data.frame(id=enrichedSig_sub[,1],link=enrichedSig_sub[,2],statistic=paste("C=",enrichedSig_sub[,3],"; O=",enrichedSig_sub[,4],"; E=",round(enrichedSig_sub[,5],digits=2),"; R=",round(enrichedSig_sub[,6],digits=2),"; PValue=",format(enrichedSig_sub[,7],scientific=TRUE,digits=3),"; FDR=",format(enrichedSig_sub[,8],scientific=TRUE,digits=3),sep=""),genes=gsub(","," ",enrichedSig_sub[,9]),stringsAsFactors=FALSE)
				 		}
				 	
				}else{
				 		cat("Categories with red color represent positive related categories while categories with blue color represent negative related categories.<br/>\n",file=outputHtmlFile,append=TRUE)
				 		if(!is.null(geneSetDes)){
				 			cat('<th width="10%">ID</th><th width="20%">Name</th><th width="30%">Statistic</th><th width="40%">Genes</th></tr></thead>\n',file=outputHtmlFile,append=TRUE)
				 			extractSig <- data.frame(id=enrichedSig_sub[,1],name=enrichedSig_sub[,2],link=enrichedSig_sub[,3],statistic=paste("Size=",enrichedSig_sub[,4],"; L=",enrichedSig_sub[,5],"; ES=",round(enrichedSig_sub[,6],digits=2),"; NES=",round(enrichedSig_sub[,7],digits=2),"; PValue=",format(enrichedSig_sub[,8],scientific=TRUE,digits=3),"; FDR=",format(enrichedSig_sub[,9],scientific=TRUE,digits=3),sep=""),genes=gsub(","," ",enrichedSig_sub[,10]),stringsAsFactors=FALSE)
				 		}else{
				 	  	cat('<th width="10%">ID</th width="30%"><th>Statistic</th><th width="60%">Genes</th></tr></thead>\n',file=outputHtmlFile,append=TRUE)
				 			extractSig <- data.frame(id=enrichedSig_sub[,1],link=enrichedSig_sub[,2],statistic=paste("Size=",enrichedSig_sub[,3],"; L=",enrichedSig_sub[,4],"; ES=",round(enrichedSig_sub[,5],digits=2),"; NES=",round(enrichedSig_sub[,6],digits=2),"; Pvalue=",format(enrichedSig_sub[,7],scientific=TRUE,digits=3),"; FDR=",format(enrichedSig_sub[,8],scientific=TRUE,digits=3),sep=""),genes=gsub(","," ",enrichedSig_sub[,9]),stringsAsFactors=FALSE)
				 		}
				 	
				}
				cat('<tbody>\n',file=outputHtmlFile,append=TRUE,sep="")

				for(i in c(1:nrow(enrichedSig_sub))){
						cat("<tr>\n",file=outputHtmlFile,append=TRUE)
						if(enrichMethod==existingMethods[1]){
							color <- "black"
						}else{
							if(enrichedSig_sub[i,"NES"]>0){
								color <- "red"
							}else{
								color <- "blue"
							}
						}
						id <- extractSig[i,1]
						if(!is.null(geneSetDes)){
							name <- extractSig[i,2]
							link <- extractSig[i,3]
							stat <- extractSig[i,4]
							genes <- gsub(";"," ",extractSig[i,5])
							cat('<td><a href="',link,'" target="_blank" style="color:',color,'">',id,'</a></font></td><td><font color="',color,'">',name,'</font></td><td><font color="',color,'">',stat,'</font></td><td><font color="',color,'">',genes,'</font></td>\n',file=outputHtmlFile,append=TRUE,sep="")
						}else{
							link <- extractSig[i,2]
							stat <- extractSig[i,3]
							genes <- gsub(";"," ",extractSig[i,4])
							cat('<td><a href="',link,'" target="_blank" style="color:',color,'">',id,'</a></font></td><td><font color="',color,'">',stat,'</font></td><td><font color="',color,'">',genes,'</font></td>\n',file=outputHtmlFile,append=TRUE,sep="")
						}
						cat("</tr>\n",file=outputHtmlFile,append=TRUE)
				}
				cat('</tbody></table></div><br/><br/>\n',file=outputHtmlFile,append=TRUE)
		}
		
		
	  cat('<iframe src="',paste(hostName,"/html/foot.html",sep=""),'" style="border:0px;height:130px;width:100%"></iframe>\n',file=outputHtmlFile,append=TRUE,sep="")

		cat("</body></html>",file=outputHtmlFile,append=TRUE)
}

.parameterSummaryDescription <- function(enrichedSig,enrichMethod,existingMethods,dNum,outputHtmlFile,minNum,maxNum,fdrMethod,sigMethod,fdrThr,topThr,perNum,lNum){

			cat("<b>Parameters for the enrichment analysis:</b>\n",file=outputHtmlFile,append=TRUE,sep="")
			cat("<ul><li><b>Minimum number of Entrez Gene IDs in the category:</b>",minNum,"</li>\n",file=outputHtmlFile,append=TRUE,sep="")
			cat("<li><b>Maximum number of Entrez Gene IDs in the category:</b>",maxNum,"</li>\n",file=outputHtmlFile,append=TRUE,sep="")
			if(enrichMethod==existingMethods[1]){
				cat("<li><b>FDR Method:</b>",fdrMethod,"</li>\n",file=outputHtmlFile,append=TRUE,sep="")
			}
			if(sigMethod=="fdr"){
				cat("<li><b>Significance Level:</b> FDR<",fdrThr,"</li>\n",file=outputHtmlFile,append=TRUE,sep="")
			}else{
				cat("<li><b>Significance Level:</b> Top",topThr,"</li>\n",file=outputHtmlFile,append=TRUE,sep="")
			}
			
			if(enrichMethod==existingMethods[2]){
				cat("<li><b>Number of permutation:</b>",perNum,"</li>\n",file=outputHtmlFile,append=TRUE,sep="")
				cat("<li><b>Number of categories with the outputted leading edge genes:</b>",lNum,"</li></ul>\n",file=outputHtmlFile,append=TRUE,sep="")
			}else{
				cat("</ul>\n",file=outputHtmlFile,append=TRUE,sep="")
			}

		if(!is.null(enrichedSig)){
				if(enrichMethod==existingMethods[1]){
					if(dNum>=nrow(enrichedSig)){
						cat("Based on the above parameters, <b>",nrow(enrichedSig),"</b> categories were identified as enriched categories and all were shown in this report.<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
					}else{
						cat("Based on the above parameters, <b>",nrow(enrichedSig),"</b> categories were identified as enriched categories, in which <b>",dNum,"</b> most significant categories were shown in this report. All significant categories can be downloaded from the 'Result Download' link.<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
					}
				}else{
					x <- enrichedSig[enrichedSig[,"NES"]>0,]
					y <- enrichedSig[enrichedSig[,"NES"]<0,]
					
					if(nrow(x)>0){
						if(dNum>=nrow(x)){
							cat("Based on the above parameters, <b>",nrow(x)," positive related </b>categories were identified as enriched categories and all were shown in this report.<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
						}else{
							cat("Based on the above parameters, <b>",nrow(x)," positive related </b>categories were identified as enriched categories, in which <b>",dNum,"</b> most significant categories were shown in this report. All positive related significant categories can be downloaded from the 'Result Download' link.<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
						}
					}else{
						cat("Based on the above parameters, <b>No positive related </b>category was identified as enriched category.<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
					}
					
					if(nrow(y)>0){
						if(dNum>=nrow(y)){
							cat("Based on the above parameters, <b>",nrow(y)," negative related </b>categories were identified as enriched categories and all were shown in this report.<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
						}else{
							cat("Based on the above parameters, <b>",nrow(y)," negative related </b>categories were identified as enriched categories, in which <b>",dNum,"</b> most significant categories were shown in this report. All positive related significant categories can be downloaded from the 'Result Download' link.<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
						}
					}else{
						cat("Based on the above parameters, <b>No negative related </b>category was identified as enriched category.<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
					}
					
				}
			}else{
				cat("Based on the above parameters, <b>No</b> category was identified as enriched category.<br/>\n",file=outputHtmlFile,append=TRUE,sep="")
			}
}

.extractSubSig <- function(enrichMethod,existingMethods,enrichedSig,dNum){
	###extract the dNum siganificant categories from the enrichedSig
		if(enrichMethod==existingMethods[1]){
				enrichedSig <- enrichedSig[order(enrichedSig[,"FDR"],enrichedSig[,"PValue"]),]
				if(nrow(enrichedSig)>dNum){
				 	enrichedSig_sub <- enrichedSig[1:dNum,]
				}else{
				 	enrichedSig_sub <- enrichedSig
				}
		}else{
				x <- enrichedSig[enrichedSig[,"NES"]>0,]
				y <- enrichedSig[enrichedSig[,"NES"]<0,]
				 			
				if(nrow(x)>0){
				 	x <- x[order(x[,"FDR"],x[,"PValue"]),]
				 	if(nrow(x)>dNum){
				 		x <- x[1:dNum,]
				 	}
				}
				 			
				if(nrow(y)>0){
				 	y <- y[order(y[,"FDR"],y[,"PValue"]),]
				 	if(nrow(y)>dNum){
				 		y <- y[1:dNum,]
				 	}
				}
				 			
				enrichedSig_sub <- rbind(x,y)
		}
		return(enrichedSig_sub)
}

.statisticDescription <- function(enrichMethod,existingMethods,outputHtmlFile,fdrMethod){
	####white the description of the statistic in the HTML page
	cat('<h4>Detailed information of the enriched categories</h4>\n',file=outputHtmlFile,append=TRUE)
	cat('<div class="tooltip_templates">The statistics\n',file=outputHtmlFile,append=TRUE)
	cat('<a class="tooltip-container"><img src="http://www.webgestalt.org/images/info.png" width="13" height="13"/>\n',file=outputHtmlFile,append=TRUE)
	cat('<span class="tooltip_content"><strong>',file=outputHtmlFile,append=TRUE)
	
	if(enrichMethod==existingMethods[1]){
		cat("<ul><li>C: the number of reference genes in the category</li>",file=outputHtmlFile,append=TRUE)
		cat("<li>O: the number of genes in the user gene list and also in the category</li>",file=outputHtmlFile,append=TRUE)
		cat("<li>E: The expected number in the category</li>",file=outputHtmlFile,append=TRUE)
		cat("<li>R: ratio of enrichment</li>",file=outputHtmlFile,append=TRUE)
		cat("<li>PValue: p value from hyergeometric test</li>",file=outputHtmlFile,append=TRUE)
		cat("<li>FDR: FDR from ",fdrMethod,"</li></ul>",file=outputHtmlFile,append=TRUE,sep="")
		cat("</strong></span></a>\n",file=outputHtmlFile,append=TRUE,sep="")
		cat(" for the enriched categories and the genes in the user gene list and also in the category are listed in the table. <br/>",file=outputHtmlFile,append=TRUE)
	}else{
		cat("<ul><li>Size: the number of genes in the user gene list and also in the category</li>",file=outputHtmlFile,append=TRUE)
		cat("<li>L: the number of leading edge genes</li>",file=outputHtmlFile,append=TRUE)
		cat("<li>ES: Enrichment Score</li>",file=outputHtmlFile,append=TRUE)
		cat("<li>NES: Normalized Enrichment Score</li>",file=outputHtmlFile,append=TRUE)
		cat("<li>PValue: p value from GSEA</li>",file=outputHtmlFile,append=TRUE)
		cat("<li>FDR: FDR from GSEA</li></ul>",file=outputHtmlFile,append=TRUE)
		cat("</strong></span></a>\n",file=outputHtmlFile,append=TRUE,sep="")
		cat(" for the enriched categories and the leading edge genes are listed in the table.<br/>",file=outputHtmlFile,append=TRUE)
	}
				 		
	cat("</div>\n",file=outputHtmlFile,append=TRUE)
}

.createJSONFile <- function(enrichedGO,enrichMethod,existingMethods,desFile,dagFile,reportDir,outputDirectory,timeStamp){
	
	sigGO <- unique(enrichedGO[,1])
	allGO <- unique(enrichedGO[,1])
	jsonF <- paste('[{"id":"Project_',timeStamp,'_.json","reports":[',sep="")
	
	while(length(allGO)>0){
		
		g <- allGO[1]
		gdes <- desFile[desFile[,1]==g,2]
		parents <- dagFile[dagFile[,2]==g,1]
		
		
		jsonF <- paste(jsonF,'{"gID":["',g,'"],"Edge":[',sep="")
		for(i in c(1:length(parents))){
			jsonF <- paste(jsonF,'"',parents[i],'",',sep="")
		}
		jsonF <- substring(jsonF,1,nchar(jsonF)-1)
		jsonF <- paste(jsonF,'],"Name":["',gdes,'"],"Sig":["',sep="")
		
		if(length(which(sigGO==g))>0){
		  l <- paste(reportDir,"#",g,sep="")
			if(enrichMethod==existingMethods[1]){
				gnum <- enrichedGO[enrichedGO[,1]==g,"O"]
				adjp <- format(enrichedGO[enrichedGO[,1]==g,"FDR"],scientific=TRUE,digits=3)
				sig <- "red"
			}else{
				gnum <- enrichedGO[enrichedGO[,1]==g,"leadingEdgeNum"]
				adjp <- format(enrichedGO[enrichedGO[,1]==g,"FDR"],scientific=TRUE,digits=3)
				if(enrichedGO[enrichedGO[,1]==g,"NES"]>0){
					sig <- "red"
				}else{
					sig <- "blue"
				}
			}
		}else{
			sig <- "white"
			l <- ""
			gnum <- ""
			adjp <- ""
		}
		
		jsonF <- paste(jsonF,sig,'"],"Link":["',l,'"],"gNum":["',gnum,'"],"FDR":["',adjp,'"]},',sep="")
		allGO <- allGO[-1]
		allGO <- unique(c(allGO,parents))
	}
	jsonF <- substring(jsonF,1,nchar(jsonF)-1)
	jsonF <- paste(jsonF,"]}]",sep="")
	outputFileName <- paste(outputDirectory,"/Project_",timeStamp,"/Project_",timeStamp,"_.json",sep="")
	write.table(jsonF,file=outputFileName,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
	return(jsonF)
}

.linkModification <- function(enrichDatabase,enrichpathwaylink,genelist,interestingGeneMap){
	#####Modify the link to highlight the genes in the pathways. Currently, we only have wikipathway and kegg pathways that need to modify the link########

	if(enrichDatabase=="pathway_KEGG"){
		link <- .keggLinkModification(enrichpathwaylink,genelist)
		return(link)
	}
	if(enrichDatabase=="pathway_Wikipathway"){
		link <- .wikiLinkModification(enrichpathwaylink,genelist,interestingGeneMap)
		return(link)
	}
	return(enrichpathwaylink)
	
}

.keggLinkModification <- function(enrichpathwaylink,genelist){
	genelist <- gsub(";","+",genelist)
	enrichpathwaylink <- paste(enrichpathwaylink,"+",genelist,sep="")
	return(enrichpathwaylink)
}

.wikiLinkModification <- function(enrichpathwaylink,genelist,interestingGeneMap){
	genemap <- interestingGeneMap$mapped
	genelist <- unlist(strsplit(genelist,";"))
	gene_symbol <- genemap[genemap[,"entrezgene"] %in% genelist,"genesymbol"]
	for(i in c(1:length(gene_symbol))){
		enrichpathwaylink <- paste(enrichpathwaylink,"&label[]=",gene_symbol[i],sep="")
	}
	enrichpathwaylink <- paste(enrichpathwaylink,"&colors=red",sep="")
	return(enrichpathwaylink)
}

.testNull <- function(parameter){
	if(!is.null(parameter)){
		if(length(parameter)==1 && length(which(parameter=="NULL"))==1){
			parameter <- NULL
		}
	}
	return(parameter)
}