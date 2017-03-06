#' Scans over a list of input files and annotates FAME standards
#'
#' @param inputFileList Vector of file paths to scan for FAME standards
#' @param FAME_frame A parsed FAME standard reference file (find at system.file("extdata", "FIND_FAME_FRAME.txt", package="CooperR2DGC"))
#' @param numCores The number of cores to use for parallel processing. Defaults to 1
#' @param RT1Penalty Penalty used for first retention time errors.  Defaults to 1.
#' @param RT2Penalty Penalty used for first retention time errors.  Defaults to 10.
#' @param similarityCutoffWarningThreshold Similarity score threshold at which to print a warning for user to double check FAME peak for a sample. Defaults to 85.

#' @return Will write out new files with properly annotated FAME standards with the extention _FAME_appended.txt added to the end of the file path
#' @import parallel
#' @export


Find_FAME_Standards<-function(inputFileList, FAME_Frame=system.file("extdata", "FIND_FAME_FRAME.txt", package="CooperR2DGC"), numCores=4, RT1Penalty=1, RT2Penalty=10, similarityCutoffWarningThreshold=80){
  FAMES<-read.table(FAME_Frame,sep="\t",header=T)
  FAMES[,2]<-as.character(FAMES[,2])
  RTSplit<-data.frame(strsplit(FAMES[,2], " , "), stringsAsFactors = F)
  RTSplit[1,]<-gsub("\"", "", RTSplit[1,])
  RTSplit[2,]<-gsub("\"", "", RTSplit[2,])
  FAMES[,"RT1"]<-as.numeric(t(RTSplit[1,]))
  FAMES[,"RT2"]<-as.numeric(t(RTSplit[2,]))
  FAMES[,3]<-as.character(FAMES[,3])

  FAMESFileSplit<-split(FAMES,1:nrow(FAMES))
  FAMESspectraSplit<-lapply(FAMESFileSplit, function(a) strsplit(a[[3]]," "))
  FAMESspectraSplit<-lapply(FAMESspectraSplit, function(b) lapply(b, function(c) strsplit(c,":")))
  FAMESspectraSplit<-lapply(FAMESspectraSplit, function(d) t(matrix(unlist(d),nrow=2)))
  FAMESspectraSplit<-lapply(FAMESspectraSplit, function(d) d[order(d[,1]),])
  FAMESspectraSplit<-lapply(FAMESspectraSplit, function(d) apply(d,2,as.numeric))

  for(File in inputFileList){
    print(File)
    currentRawFile<-read.table(File, sep="\t", fill=T, quote="",strip.white = T, stringsAsFactors = F,header=T)
    currentRawFile[,4]<-as.character(currentRawFile[,4])
    currentRawFile<-currentRawFile[which(!is.na(currentRawFile[,3])&nchar(currentRawFile[,4])!=0),]
    currentRawFile[,2]<-as.character(currentRawFile[,2])

    #Parse retention times
    RTSplit<-data.frame(strsplit(currentRawFile[,2], " , "), stringsAsFactors = F)
    RTSplit[1,]<-gsub("\"", "", RTSplit[1,])
    RTSplit[2,]<-gsub("\"", "", RTSplit[2,])
    currentRawFile[,"RT1"]<-as.numeric(t(RTSplit[1,]))
    currentRawFile[,"RT2"]<-as.numeric(t(RTSplit[2,]))

    #Remove identical metabolite rows
    uniqueIndex<-data.frame(paste(currentRawFile[,1], currentRawFile[,2], currentRawFile[,3]))
    currentRawFile<-currentRawFile[which(!duplicated(uniqueIndex)),]
    row.names(currentRawFile)<-c(1:nrow(currentRawFile))

    #Parse metabolite spectra into a list
    currentRawFileSplit<-split(currentRawFile,1:nrow(currentRawFile))
    spectraSplit<-lapply(currentRawFileSplit, function(a) strsplit(a[[4]]," "))
    spectraSplit<-lapply(spectraSplit, function(b) lapply(b, function(c) strsplit(c,":")))
    spectraSplit<-lapply(spectraSplit, function(d) t(matrix(unlist(d),nrow=2)))
    spectraSplit<-lapply(spectraSplit, function(d) d[order(d[,1]),])
    spectraSplit<-lapply(spectraSplit, function(d) apply(d,2,as.numeric))

    #CalcSims
    SimilarityScores<-mclapply(spectraSplit, function(e) lapply(FAMESspectraSplit, function(d) ((e[,2]%*%d[,2])/(sqrt(sum(e[,2]*e[,2]))*sqrt(sum(d[,2]*d[,2]))))*100), mc.cores=numCores)
    SimilarityMatrix<-matrix(unlist(SimilarityScores), nrow=nrow(FAMES))

    #Calculate pairwise RT penalties for each current file metabolite and seed file metabolites
    RT1Index<-matrix(unlist(lapply(currentRawFile[,6],function(x) abs(x-FAMES[,4])*RT1Penalty)),nrow=nrow(FAMES))
    RT2Index<-matrix(unlist(lapply(currentRawFile[,7],function(x) abs(x-FAMES[,5])*RT2Penalty)),nrow=nrow(FAMES))
    SimilarityMatrix<-SimilarityMatrix-RT1Index-RT2Index
    if(sum(apply(SimilarityMatrix,1,max)<similarityCutoffWarningThreshold)>1){
      message(paste0("Potential Problem Match: ", FAMES[which(apply(SimilarityMatrix,1,max)<similarityCutoffWarningThreshold),1], "  "))
    }
    currentRawFile[apply(SimilarityMatrix,1,which.max),1]<-as.character(FAMES[,1])
    write.table(currentRawFile[,c(1:5)], paste0(File,"_FAME_appended.txt"),sep="\t",row.names=F,quote=F)
  }
}
