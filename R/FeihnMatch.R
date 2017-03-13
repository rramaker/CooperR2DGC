#' Takes a metabolite info table from the ConsensusAlign output and annotates with the Feihn library

#' @param PeakInfoTable Metabolite info table (the second data frame in the output from ConsensusAlign)
#' @param FeihnLibrary Provide Feihn library generated from MakeReference function to ID metabolites with retention index. (Use data(FeihnStandards_030617) to find)
#' @param RT1_Standards Vector of standard names used to adjust first retention time. All names must be found in input files. Defaults to NULL.
#' @param RT1Penalty Penalty used for first retention time errors.  Defaults to 1.
#' @param numCores Number of cores used to parallelize alignment. See parallel package. Defaults to 4.
#' @param commonIons Provide a vector of ions to ignore from the FindProblemIons function. Defaults to empty vector.


#' @return A new metabolite info file with three columns appended to that are the top three Feihn library matches
#' @import parallel
#' @export
#'
#'

FeihnMatch<-function(PeakInfoTable, FeihnLibrary, numCores=4, commonIons=c(), RT1_Standards=NULL, RT1Penalty=1){

  #Parse seed file spectras
  peakSplit<-split(PeakInfoTable,1:nrow(PeakInfoTable))
  peakSpectraSplit<-lapply(peakSplit, function(a) strsplit(a[[8]]," "))
  peakSpectraSplit<-lapply(peakSpectraSplit, function(b) lapply(b, function(c) strsplit(c,":")))
  peakSpectraSplit<-lapply(peakSpectraSplit, function(d) t(matrix(unlist(d),nrow=2)))
  peakSpectraSplit<-lapply(peakSpectraSplit, function(d) d[order(d[,1]),])
  peakSpectraSplit<-lapply(peakSpectraSplit, function(d) d[which(!d[,1]%in%commonIons),])
  peakSpectraSplit<-lapply(peakSpectraSplit, function(d) apply(d,2,as.numeric))

  #Parse standard library spectras
  standardSplit<-split(FeihnLibrary,1:nrow(FeihnLibrary))
  standardSpectraSplit<-lapply(standardSplit, function(a) strsplit(a[[3]]," "))
  standardSpectraSplit<-lapply(standardSpectraSplit, function(b) lapply(b, function(c) strsplit(c,":")))
  standardSpectraSplit<-lapply(standardSpectraSplit, function(d) t(matrix(unlist(d),nrow=2)))
  standardSpectraSplit<-lapply(standardSpectraSplit, function(d) d[order(d[,1]),])
  standardSpectraSplit<-lapply(standardSpectraSplit, function(d) d[which(!d[,1]%in%commonIons),])
  standardSpectraSplit<-lapply(standardSpectraSplit, function(d) apply(d,2,as.numeric))

  #Calculate pairwise similarity scores
  SimilarityScores<-mclapply(peakSpectraSplit, function(e) lapply(standardSpectraSplit, function(d) ((e[,2]%*%d[,2])/(sqrt(sum(e[,2]*e[,2]))*sqrt(sum(d[,2]*d[,2]))))*100), mc.cores=numCores)
  SimilarityMatrix<-matrix(unlist(SimilarityScores), nrow=length(standardSpectraSplit))

  #Compute RT differences
  RT1Index<-matrix(unlist(lapply(PeakInfoTable[,6],function(x) abs(x-FeihnLibrary[,4])*RT1Penalty)),nrow=nrow(FeihnLibrary))

  #Use RT indexes to compute RT differences
  if(!is.null(RT1_Standards)){
    RT1_Length<-max(PeakInfoTable[which(PeakInfoTable[,4]%in%RT1_Standards),9])-min(PeakInfoTable[which(PeakInfoTable[,4]%in%RT1_Standards),9])
    RT1Index<-list()
    for(Standard in RT1_Standards){
      RT1Index[[Standard]]<-matrix(unlist(lapply(PeakInfoTable[,paste(Standard,"RT1",sep="_")],function(x) abs(x-FeihnLibrary[,paste(Standard,"RT1",sep="_")])*(RT1Penalty/length(RT1_Standards)))),nrow=nrow(FeihnLibrary))*RT1_Length
    }
    RT1Index<-Reduce("+",RT1Index)
  }

  #Subtract RT penalties from Similarity Scores
  SimilarityMatrix<-SimilarityMatrix-RT1Index
  row.names(SimilarityMatrix)<-FeihnLibrary[,1]

  #Append top three ID matches to each metabolite and scores to seedRaw for output
  PeakInfoTable<-cbind(t(apply(SimilarityMatrix,2,function(x) paste(names(x[order(-x)])[1:3],round(x[order(-x)][1:3],2),sep="_"))),PeakInfoTable)
  colnames(PeakInfoTable)[1:3]<-c("FeihnMatch_1","FeihnMatch_2","FeihnMatch_3")
  return(PeakInfoTable)
}
