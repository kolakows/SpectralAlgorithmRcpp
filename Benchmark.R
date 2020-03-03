
#install.packages("dendextend")
#install.packages("mclust")
library(dendextend)
library(mclust)




#assumming that all *.labels0.gz and *.data.gz files are in ./data directory

generateNames <- function(){
  DataNames <- list.files(path="./data",pattern="*.data.gz")
  lapply(DataNames,function(y){
    name <- gsub("^(.*?)\\..*", "\\1",y)
  })
}


DataNames <- generateNames()

#create folders for data

lapply(DataNames,function(name){
  dir.create("./allData",showWarnings = FALSE)
  dir.create(file.path("./allData",name),showWarnings = FALSE)
})


#put truth labels
lapply(DataNames,function(name){
  dirPath <- "./data/"
  labelsPath <- paste0(dirPath,name,".labels0.gz")
  labels <- data.frame(as.integer(read.table(labelsPath)[,1]))
  write.table(labels,paste0("./allData/",name,"/Truth.txt"),row.names = FALSE,col.names = FALSE)
})

#put data into folders, with scale
lapply(DataNames,function(name){
  dirPath <- "./data/"
  dataPath <- paste0("./data/",name,".data.gz")
  data <- scale(read.table(dataPath))
  write.table(data,paste0("./allData/",name,"/data.txt"),row.names = FALSE,col.names = FALSE)
})

#put data into folders, without scale
lapply(DataNames,function(name){
  dirPath <- "./data/"
  dataPath <- paste0("./data/",name,".data.gz")
  data <- read.table(dataPath)
  write.table(data,paste0("./allData/",name,"/data.txt"),row.names = FALSE,col.names = FALSE)
})


#generate number of clusters for each data file
createClusterK <- function(){
  DataNames <- generateNames()
  dirPath <- "./data/"
  saveDir <- "./ClusterK/"
  clusterK <- sapply(DataNames, function(name){
    labelsPath <- paste0(dirPath,name,".labels0.gz")
    K <- max(data.frame((read.table(labelsPath)[,1])))
    c(name,K)
  })
  colnames(clusterK) <- clusterK[1,]
  clusterK <- clusterK[2,]
  write.table(clusterK,paste0(saveDir,"K.txt"))
}
createClusterK()


#rootDir ./allData/ (with dir hierarchy)
hlustniem <- function(rootDir,methods){
  dataNames <- generateNames()
  refK <- read.table(".\\ClusterK\\K.txt")
  #pls dont tell anyone I used for loop
  for (meth in methods) {
    print(meth)
    lapply(dataNames,function(name){
      folder <- paste0(rootDir,name,"/")
      data <- read.table(paste0(folder,"data.txt"))
      h <- hclust(dist(data),meth)
      labels <- cutree(h,refK[name,])
      #createPlotFromData(data,as.factor(labels),folder,paste0(name," ",meth," clustering"),paste0(meth,"Plot"))
      write.table(labels,paste0(folder,meth,".txt"),row.names = FALSE,col.names = FALSE)
    })
  }
}
methods <- c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")
hlustniem("./allData/",methods)


#Generate Genie labels

GenerateGenie <- function(){
  library(genie)
  dirPath <- "./allData/"
  dataNames <- generateNames()
  refK <- read.table(".\\ClusterK\\K.txt")
  sapply(dataNames, function(name){
    dataPath <- paste0(dirPath,name,"/data.txt")
    data <- read.table(dataPath)
    h <- hclust2(objects = as.matrix(data))
    labels <- cutree(h,refK[name,])
    write.table(labels,paste0("./allData/",name,"/Genie.txt"),row.names = FALSE,col.names = FALSE);
  })
}
GenerateGenie()


methods <- c("Genie", "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")

#rootDir ./allData/ (with dir hierarchy) 
CalcIndexes <- function(rootDir,methodNames){
  dataNames <- generateNames()
  for (name in dataNames) {
    truth <- read.table(paste0(rootDir,name,"/Truth.txt"))
    benchm <- lapply(methodNames, function(meth){
      labels <- read.table(paste0(rootDir,name,"/",meth,".txt"))
      Rand <- adjustedRandIndex(truth$V1,labels$V1)
      FM <- FM_index(truth$V1,labels$V1)
      c(Rand=Rand,FM=FM)
    })
    benchm <- data.frame(benchm)
    colnames(benchm) <- methodNames
    write.csv(benchm,paste0(rootDir,name,"/benchm.csv"))
  }
}

CalcIndexes("./allData/",methods)


ReadBenchm<- function(rootDir){
  dataNames <- generateNames()
  lbenchm <- lapply(dataNames,function(name){
    read.csv(paste0(rootDir,name,"/benchm.csv"))
  })
  names(lbenchm) <- dataNames
  lbenchm
}
b <- ReadBenchm("./allData/")


CalcMeanIndexesAndWriteCSV <- function(b){
  Randbenchm <- lapply(b,function(x){
    x[1,]
  })
  
  RandStats <- do.call(rbind,Randbenchm)
  RandStats$X <- NULL
  RandIndexForData <- colMeans(RandStats)
  RandIndexForData <- RandIndexForData[order(-RandIndexForData)]
  RandIndexForData$Index <- "Rand"
  RandIndexForData <- data.frame(RandIndexForData)
  
  FMbenchm <- lapply(b,function(x){
    x[2,]
  })
  
  FMStats <- do.call(rbind,FMbenchm)
  FMStats$X <- NULL
  FMIndexForData <- colMeans(FMStats)
  FMIndexForData <- FMIndexForData[order(-FMIndexForData)]
  FMIndexForData$Index <- "FM"
  FMIndexForData <- data.frame(FMIndexForData)
  
  library(dplyr)
  
  FMIndexForData <- FMIndexForData %>% select(Index,everything())
  RandIndexForData <- RandIndexForData %>% select(Index,everything())
  write.csv(RandIndexForData,"./RandIndex.csv")
  write.csv(FMIndexForData,"./FMIndex.csv")
}

CalcMeanIndexesAndWriteCSV(b)




#drawing plots 

GenerateTruthPlots2d <- function(){
  DataNames <- list.files(path="./data",pattern="*.data.gz")
  LabelNames <- list.files(path="./data",pattern="*.labels0.gz")
  stats <- sapply(DataNames,function(x){
    dim(read.table(paste0("./data/",x)))
  })
  statsdf <- as.data.frame(stats)
  dim2 <- colnames(statsdf[,statsdf[2,]==2])
  dirPath <- "./data/"
  saveDir <- "./2d/"
  sapply(dim2, function(dataFile){
    name <- gsub("^(.*?)\\..*", "\\1",dataFile)
    labelsPath <- paste0(dirPath,name,".labels0.gz")
    dataPath <- paste0(dirPath,dataFile)
    createPlot(dataPath,labelsPath,paste0(saveDir,name,"/"),paste0(name," Ground Truth"),"GroundTruth")
  })
  dim2
}
GenerateTruthPlots2d()


#Genie plots with simultaneous writing labels
GenerateGeniePlot2d <- function(){
  library(genie)
  dirPath <- "./data/"
  saveDir <- "./allData/"
  dim2 <- generateNames()
  refK <- read.table(".\\ClusterK\\K.txt")
  sapply(dim2, function(dataFile){
    name <- gsub("^(.*?)\\..*", "\\1",dataFile)
    dataPath <- paste0(dirPath,dataFile)
    data <- read.table(dataPath)
    h <- hclust2(objects = as.matrix(data))
    labels <- cutree(h,refK[name,])
    createPlotFromData(data,as.factor(labels),paste0(saveDir,name,"/"),paste0(name," Genie clustering"),"GeniePlot")
    write.table(labels,paste0("./allData/",name,"/","Genie.txt"),row.names = FALSE,col.names = FALSE);
  })
}

# general function for drawing plots
createPlotFromData <- function(data,labels,saveDir,plotTitle,saveName){
  Alldata <- cbind(data,labels)
  colnames(Alldata) <- c("V1","V2","labels")
  ggplot(Alldata,aes(x=V1,y=V2)) + 
    geom_point(aes(col=labels),size=1.5) + 
    labs(title=plotTitle) + 
    theme(axis.title = element_blank()) + theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      legend.position = "none")
  ggsave(paste0(saveName,".png"),path=saveDir)
}

