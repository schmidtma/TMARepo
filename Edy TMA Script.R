# Script to analyze Cell Profiler results for Edy's TMA and hypothetically
# eventually maybe Nikki/Ally's CD166 work as well but it relies on (x,y)
# TMA coordinates to cobble together image sets for spots.
#
# Certain columns of data are more pertinent than others. In total there may 
# be more than 200 columns worth of data but only about 20-30 will be used in
# this script. The rest are useful features for use in machine learning tools
# such as Cell Profiler Analyst.
# 
# Important columns (names) from Cell-specific analysis spreadsheets:
#   ObjectNumber, Metadata_ELEMENTID, Metadata_ImageN, 
#   Metdata_X_INDEX, Metadata_Y_INDEX, Classify_CKPos_CD166Neg_CD98Neg
#   Classify_CKPos_CD166Pos_CD98Neg, Classify_CKPos_CD166Neg_CD98Pos
#   Classify_CKPos_CD166Pos_CD98Pos
# 
# Important columns (names) from Image spreadsheets:
#   Classify_CKPos_CD166Neg_CD98Neg_NumObjectsPerBin,
#   Classify_CKPos_CD166Neg_CD98Neg_PctObjectsPerBin,	
#   Classify_CKPos_CD166Neg_CD98Pos_NumObjectsPerBin,	
#   Classify_CKPos_CD166Neg_CD98Pos_PctObjectsPerBin,
#   Classify_CKPos_CD166Pos_CD98Neg_NumObjectsPerBin,
#   Classify_CKPos_CD166Pos_CD98Neg_PctObjectsPerBin,
#   Classify_CKPos_CD166Pos_CD98Pos_NumObjectsPerBin,
#   Classify_CKPos_CD166Pos_CD98Pos_PctObjectsPerBin,
#   Count_CD166Pos,  Count_CD166PosCells,	Count_CD166WithNuclei,
#   Count_CD98Pos,	Count_CD98PosCells,	Count_CD98WithNuclei,
#   Count_CKPos,	Count_CKPosCells,	Count_CKWithNuclei,	Count_Nuclei
#   Metadata_X_INDEX,  Metadata_Y_INDEX
# 
# The Metadata_X_INDEX and _Y_INDEX can be used via a match/merge function
# to relate them to TMA spot position, this has already been done so I'm
# pulling that code from another script I've written to be used below.
# 
# In totality this will just run and produce graphs from a given Cell Profiler
# run and output them to a new directory for pdf's of graphs.

#setwd("E:/Cell Profiler/20150807 Edy TMA 2-0-1")

setwd("E:/Cell Profiler/20150825 Images for Full Spot Analysis/Outputs from 20150825 Images/Compiled data")

ImageData <- read.csv("MyExpt_Image.csv", header=TRUE)
# ImageData <- read.csv("MyImage_master.csv", header=TRUE)
Data <- ImageData[,c(
  "Classify_CKPos_CD166Neg_CD98Neg_NumObjectsPerBin",
  "Classify_CKPos_CD166Neg_CD98Neg_PctObjectsPerBin",  
  "Classify_CKPos_CD166Neg_CD98Pos_NumObjectsPerBin",	
  "Classify_CKPos_CD166Neg_CD98Pos_PctObjectsPerBin",
  "Classify_CKPos_CD166Pos_CD98Neg_NumObjectsPerBin",
  "Classify_CKPos_CD166Pos_CD98Neg_PctObjectsPerBin",
  "Classify_CKPos_CD166Pos_CD98Pos_NumObjectsPerBin",
  "Classify_CKPos_CD166Pos_CD98Pos_PctObjectsPerBin",
  "Count_CD166Pos",  "Count_CD166PosCells",	"Count_CD166WithNuclei",
  "Count_CD98Pos",	"Count_CD98PosCells",	"Count_CD98WithNuclei",
  "Count_CKPos",	"Count_CKPosCells",	"Count_CKWithNuclei",	"Count_Nuclei",
  "Metadata_X_INDEX",  "Metadata_Y_INDEX")]

rm(ImageData)

attach(Data)

# ## Make TMA map by SpotID:
# t(matrix(48:1, 6,8))
# [,1] [,2] [,3] [,4] [,5] [,6]
# [1,]   48   47   46   45   44   43
# [2,]   42   41   40   39   38   37
# [3,]   36   35   34   33   32   31
# [4,]   30   29   28   27   26   25
# [5,]   24   23   22   21   20   19
# [6,]   18   17   16   15   14   13
# [7,]   12   11   10    9    8    7
# [8,]    6    5    4    3    2    1

SpotID <- c(t(matrix(48:1, 6,8)))

### For simplicity CD98Pos and CD166Pos by default mean also CK Positive, this is a condition of the CP Pipeline.
DblNegPct <- c(); CD98PosPct <- c(); CD166PosPct <- c(); DblPosPct <- c()
DblNegCnt <- c(); CD98PosCnt <- c(); CD166PosCnt <- c(); DblPosCnt <- c()
NucleiTotal <- c()

# note the values 0:5 and 0:7 may need to be modified for other TMAs,
# these values and the SpotID map above are dependent on the shape and size
# of the given TMA. These values can be inferred from the Metadata_x_index
# value but only if given a full data set ie: 0:max(x) and 0:max(y) are 
# functionalized inputs but only work if the maximum number of rows/columns
# from the TMA are taken as inputs to Cell Profiler.

for( x in c(0:5) )  {
  
  for( y in c(0:7) ) {
    
    DblNegPct <- append(DblNegPct, mean(Classify_CKPos_CD166Neg_CD98Neg_PctObjectsPerBin[which(Metadata_X_INDEX == x & Metadata_Y_INDEX == y)]), length(DblNegPct))
    names(DblNegPct)[length(DblNegPct)] <- paste("(",x, ", ", y,")", sep="")
    
    CD98PosPct <-append(CD98PosPct, mean(Classify_CKPos_CD166Neg_CD98Pos_PctObjectsPerBin[which(Metadata_X_INDEX == x & Metadata_Y_INDEX == y)]), length(CD98PosCnt))
    names(CD98PosPct)[length(CD98PosPct)] <- paste("(",x, ", ", y,")", sep='')
    
    CD166PosPct <-append(CD166PosPct, mean(Classify_CKPos_CD166Pos_CD98Neg_PctObjectsPerBin[which(Metadata_X_INDEX == x & Metadata_Y_INDEX == y)]), length(CD166PosPct))
    names(CD166PosPct)[length(CD166PosPct)] <- paste("(",x, ", ", y,")", sep='')
    
    DblPosPct <-append(DblPosPct, mean(Classify_CKPos_CD166Pos_CD98Pos_PctObjectsPerBin[which(Metadata_X_INDEX == x & Metadata_Y_INDEX == y)]), length(DblPosPct))
    names(DblPosPct)[length(DblPosPct)] <- paste("(",x, ", ", y,")", sep='')
    
    ### Percent values were misleading due to some controls having low 
    ### expression of CK but high colocalization with CK.
    
    DblNegCnt <- append(DblNegCnt, sum(Classify_CKPos_CD166Neg_CD98Neg_NumObjectsPerBin[which(Metadata_X_INDEX == x & Metadata_Y_INDEX == y)]), length(DblNegCnt))
    names(DblNegCnt)[length(DblNegCnt)] <- paste("(",x, ", ", y,")", sep="")
    
    CD98PosCnt <-append(CD98PosCnt, sum(Classify_CKPos_CD166Neg_CD98Pos_NumObjectsPerBin[which(Metadata_X_INDEX == x & Metadata_Y_INDEX == y)]), length(CD98PosCnt))
    names(CD98PosCnt)[length(CD98PosCnt)] <- paste("(",x, ", ", y,")", sep='')
    
    CD166PosCnt <-append(CD166PosCnt, sum(Classify_CKPos_CD166Pos_CD98Neg_NumObjectsPerBin[which(Metadata_X_INDEX == x & Metadata_Y_INDEX == y)]), length(CD166PosCnt))
    names(CD166PosCnt)[length(CD166PosCnt)] <- paste("(",x, ", ", y,")", sep='')
    
    DblPosCnt <-append(DblPosCnt, sum(Classify_CKPos_CD166Pos_CD98Pos_NumObjectsPerBin[which(Metadata_X_INDEX == x & Metadata_Y_INDEX == y)]), length(DblPosCnt))
    names(DblPosCnt)[length(DblPosCnt)] <- paste("(",x, ", ", y,")", sep='')
    
    NucleiTotal <- append(NucleiTotal, sum(Count_Nuclei[which(Metadata_X_INDEX == x & Metadata_Y_INDEX == y)]), length(NucleiTotal))
    names(NucleiTotal)[length(NucleiTotal)] <- paste("(",x, ", ", y, ")", sep='')
    
  }
}

SpotInfo <- read.csv("SpotInfo.csv", header=TRUE)
SpotStage <- read.csv("SpotStage.csv", row.names=c(1), header=TRUE)

AllInfo <- merge(SpotInfo, 
                 data.frame(DblNegPct, CD166PosPct, CD98PosPct, DblPosPct,
                                      DblNegCnt, CD166PosCnt, CD98PosCnt, DblPosCnt,
                                      NucleiTotal, SpotID, names(DblNegPct)))

### Object cleanup
rm(list=colnames(AllInfo))


colnames(AllInfo)[16] <- 'xy'

AllInfo$pT <- as.factor(substr(AllInfo$pTNM, 1, 2))

### Box and whisker plots for different ways of comparing the data

### By percentage of total CK+ cells
pdf("Plots of Expression by Percent.pdf")

plot(AllInfo$DblPosPct ~ AllInfo$Stage, main=c('% CK+/CD166+/CD98+ Cells'),
     xlab=c('Stage'), ylab=c('% CK+/CD166+/CD98+ Cells'))

plot(AllInfo$CD166PosPct ~ AllInfo$Stage, main=c('% CK+/CD166+/CD98- Cells'),
     xlab=c('Stage'), ylab=c('% CK+/CD166+/CD98- Cells'))

plot(AllInfo$CD98PosPct ~ AllInfo$Stage, main=c('% CK+/CD166-/CD98+ Cells'),
     xlab=c('Stage'), ylab=c('% CK+/CD166-/CD98+ Cells'))

plot(AllInfo$DblNegPct ~ AllInfo$Stage, main=c('% CK+/CD166-/CD98- Cells'),
     xlab=c('Stage'), ylab=c('% CK+/CD166+/CD98+ & CD166+/Cd98- Cells'))

plot(AllInfo$DblPosPct + AllInfo$CD166PosPct ~ AllInfo$Stage, main=c('% CK+/CD166+ cells'),
     xlab=c('Stage'), ylab=c('% CK+/CD166-/CD98- Cells'))
dev.off()

write.csv(AllInfo, "Data used for Graphs.csv", row.names=FALSE)

# sessionInfo()
# R version 3.2.0 (2015-04-16)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 7 x64 (build 7601) Service Pack 1
# 
# locale:
#   [1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252   
# [3] LC_MONETARY=English_United States.1252 LC_NUMERIC=C                          
# [5] LC_TIME=English_United States.1252    
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# loaded via a namespace (and not attached):
#   [1] tools_3.2.0

#########################################################################

# The following is the script that acted as a template for this.

#########################################################################
########### BEGIN TEMPLATE SCRIPT #######################################
#########################################################################

# setwd("C:/CP Outputs/Excel Files/")
# 
# ### Quick aside: By taking 5 2048^2px images from each TMA spot I'm effectively surveying 25% of the total area of the spot
# ### This can be addressed simply by taking more representative snapshots from within each TMA spot.
# ### (2048^2)*5/(((2048*5)^2)*(pi*5^2)/100) = ~25%. Each TMA spot is approximately 5x5 2048 images, however it's a square
# ### circumscribing a circle so only ~78% of the area circumscribed is covered in tissue.
# 
# ImageData <- read.csv("MyExpt_Image.csv", header=T)
# 
# ###THIS WILL VARY AND IS NOT FUNCTIONALLY PROGRAMMED, PULL COLUMNS AS NECESSARY FOR YOUR DATA
# ###First set of 8 + 18 correspond to the classify categories and total nuclei count per image.
# 
# FrameData <- ImageData[,c(1:8, 18,
#                           which(colnames(ImageData) == "Metadata_ImageN"),
#                           which(colnames(ImageData) == "Metadata_X_INDEX"),
#                           which(colnames(ImageData) == "Metadata_Y_INDEX"))]
# 
# attach(FrameData)
# SpotID <- 48:1
# 
# ### For simplicity CD98Pos and CD166Pos by default mean also CK Positive, this is a condition of the CP Pipeline.
# DblNegPct <- c(); CD98PosPct <- c(); CD166PosPct <- c(); DblPosPct <- c()
# DblNegCnt <- c(); CD98PosCnt <- c(); CD166PosCnt <- c(); DblPosCnt <- c()
# NucleiTotal <- c()
# 
# ### Should include some graphs to help look at QC on a per image basis too. What values are good for this?
# ### Maybe a ratio of primary objects in RGPink versus DAPI nuclei? In any image it should always be DAPI>>other colors
# ### Maybe pull in one of the cell object sheets and make a table to look at the most common child number. If for 
# ### any CP pipeline run CK+ cells have more than 1 child in more than 3-5% of cases then I'd call the run garbage.
# ### 
# 
# ### Using two for loops, make per image means of percentage positive cells for each TMA spot by coordinate.
# 
# for(x in min(unique(Metadata_X_INDEX)):max(unique(Metadata_X_INDEX))){
#   
#   for(y in min(unique(Metadata_Y_INDEX)):max(unique(Metadata_Y_INDEX))){
#     
#     DblNegPct <- append(DblNegPct, mean(Classify_CKPos_CD166Neg_CD98Neg_PctObjectsPerBin[which(Metadata_X_INDEX == x & Metadata_Y_INDEX == y)]), length(DblNegPct))
#     names(DblNegPct)[length(DblNegPct)] <- paste("(",x, ", ", y,")", sep="")
#     
#     CD98PosPct <-append(CD98PosPct, mean(Classify_CKPos_CD166Neg_CD98Pos_PctObjectsPerBin[which(Metadata_X_INDEX == x & Metadata_Y_INDEX == y)]), length(CD98PosCnt))
#     names(CD98PosPct)[length(CD98PosPct)] <- paste("(",x, ", ", y,")", sep='')
#     
#     CD166PosPct <-append(CD166PosPct, mean(Classify_CKPos_CD166Pos_CD98Neg_PctObjectsPerBin[which(Metadata_X_INDEX == x & Metadata_Y_INDEX == y)]), length(CD166PosPct))
#     names(CD166PosPct)[length(CD166PosPct)] <- paste("(",x, ", ", y,")", sep='')
#     
#     DblPosPct <-append(DblPosPct, mean(Classify_CKPos_CD166Pos_CD98Pos_PctObjectsPerBin[which(Metadata_X_INDEX == x & Metadata_Y_INDEX == y)]), length(DblPosPct))
#     names(DblPosPct)[length(DblPosPct)] <- paste("(",x, ", ", y,")", sep='')
#     
#     ### Percent values were misleading due to some controls having low 
#     ### expression of CK but high colocalization with CK.
#     
#     DblNegCnt <- append(DblNegCnt, sum(Classify_CKPos_CD166Neg_CD98Neg_NumObjectsPerBin[which(Metadata_X_INDEX == x & Metadata_Y_INDEX == y)]), length(DblNegCnt))
#     names(DblNegCnt)[length(DblNegCnt)] <- paste("(",x, ", ", y,")", sep="")
#     
#     CD98PosCnt <-append(CD98PosCnt, sum(Classify_CKPos_CD166Neg_CD98Pos_NumObjectsPerBin[which(Metadata_X_INDEX == x & Metadata_Y_INDEX == y)]), length(CD98PosCnt))
#     names(CD98PosCnt)[length(CD98PosCnt)] <- paste("(",x, ", ", y,")", sep='')
#     
#     CD166PosCnt <-append(CD166PosCnt, sum(Classify_CKPos_CD166Pos_CD98Neg_NumObjectsPerBin[which(Metadata_X_INDEX == x & Metadata_Y_INDEX == y)]), length(CD166PosCnt))
#     names(CD166PosCnt)[length(CD166PosCnt)] <- paste("(",x, ", ", y,")", sep='')
#     
#     DblPosCnt <-append(DblPosCnt, sum(Classify_CKPos_CD166Pos_CD98Pos_NumObjectsPerBin[which(Metadata_X_INDEX == x & Metadata_Y_INDEX == y)]), length(DblPosCnt))
#     names(DblPosCnt)[length(DblPosCnt)] <- paste("(",x, ", ", y,")", sep='')
#     
#     NucleiTotal <- append(NucleiTotal, sum(Count_Nuclei[which(Metadata_X_INDEX == x & Metadata_Y_INDEX == y)]), length(NucleiTotal))
#     names(NucleiTotal)[length(NucleiTotal)] <- paste("(",x, ", ", y, ")", sep='')
#     
#   }
# }
# 
# ## Per image means Mapped as heatmap? Then place text of cancer stage over it?
# ## Easy enough to do but Edy is less interested in that.
# 
# SpotInfo <- read.csv("SpotInfo.csv", header=T)
# SpotStage <- read.csv("SpotStage.csv", row.names=c(1), header=T)
# 
# AllInfo <- merge(SpotInfo, data.frame(DblNegPct, CD166PosPct, CD98PosPct, DblPosPct,
#                                       DblNegCnt, CD166PosCnt, CD98PosCnt, DblPosCnt,
#                                       NucleiTotal, SpotID, names(DblNegPct)))
# 
# ### Object cleanup
# rm(list=colnames(AllInfo))
# 
# colnames(AllInfo)[16] <- 'xy'
# 
# AllInfo$pT <- as.factor(substr(AllInfo$pTNM, 1, 2))
# 
# ### Box and whisker plots for different ways of comparing the data
# 
# ### By percentage of total CK+ cells
# pdf("Plots of Expression by Percent.pdf")
# 
# plot(AllInfo$DblPosPct ~ AllInfo$Stage, main=c('% CK+/CD166+/CD98+ Cells'),
#      xlab=c('Stage'), ylab=c('% CK+/CD166+/CD98+ Cells'))
# 
# plot(AllInfo$CD166PosPct ~ AllInfo$Stage, main=c('% CK+/CD166+/CD98- Cells'),
#      xlab=c('Stage'), ylab=c('% CK+/CD166+/CD98- Cells'))
# 
# plot(AllInfo$CD98PosPct ~ AllInfo$Stage, main=c('% CK+/CD166-/CD98+ Cells'),
#      xlab=c('Stage'), ylab=c('% CK+/CD166-/CD98+ Cells'))
# 
# plot(AllInfo$DblNegPct ~ AllInfo$Stage, main=c('% CK+/CD166-/CD98- Cells'),
#      xlab=c('Stage'), ylab=c('% CK+/CD166+/CD98+ & CD166+/Cd98- Cells'))
# 
# plot(AllInfo$DblPosPct + AllInfo$CD166PosPct ~ AllInfo$Stage, main=c('% CK+/CD166+ cells'),
#      xlab=c('Stage'), ylab=c('% CK+/CD166-/CD98- Cells'))
# dev.off()
# 
# ### By absolute count of cells positive for maker (must also be CK+)
# pdf('Plots of Expression by Count.pdf')
# 
# plot(AllInfo$DblPosCnt ~ AllInfo$Stage, main=c('Total CK+/CD166+/CD98+ Cells'), 
#      xlab=c('Stage'), ylab=c('Total CK+/CD166+/CD98+ Cells'))
# 
# plot(AllInfo$CD166PosCnt ~ AllInfo$Stage, main=c('Total CK+/CD166+/CD98- Cells'),
#      xlab=c('Stage'), ylab=c('Total CK+/CD166+/CD98- Cells'))
# 
# plot(AllInfo$CD98PosCnt ~ AllInfo$Stage, main=c('Total CK+/CD166-/CD98+ Cells'),
#      xlab=c('Stage'), ylab=c('Total CK+/CD166-/CD98+ Cells'))
# 
# plot(AllInfo$DblNegCnt ~ AllInfo$Stage, main=c('Total CK+/CD166-/CD98- Cells'),
#      xlab=c('Stage'), ylab=c('Total CK+/CD166+ Cells'))
# 
# plot(AllInfo$DblPosCnt + AllInfo$CD166PosCnt ~ AllInfo$Stage, main=c('Total CK+/CD166+ cells'),
#      xlab=c('Stage'), ylab=c('Total CK+/CD166+/CD98+ & CK+/CD166+/CD98- Cells'))
# dev.off()
# 
# ### By pTNM T stage (T1 versus T2 and control)
# ### Percent
# pdf("Plots of Expression by Percent Against pTNM.pdf")
# 
# plot(AllInfo$DblPosPct ~ as.factor(AllInfo$pT), main=c('% CK+/CD166+/CD98+ Cells'),
#      xlab=c('T Parameter'), ylab=c('% CK+/CD166+/CD98+ Cells'))
# 
# plot(AllInfo$CD166PosPct ~ as.factor(AllInfo$pT), main=c('% CK+/CD166+/CD98- Cells'),
#      xlab=c('T Parameter'), ylab=c('% CK+/CD166+/CD98- Cells'))
# 
# plot(AllInfo$CD98PosPct ~ as.factor(AllInfo$pT), main=c('% CK+/CD166-/CD98+ Cells'),
#      xlab=c('T Parameter'), ylab=c('% CK+/CD166-/CD98+ Cells'))
# 
# plot(AllInfo$DblNegPct ~ as.factor(AllInfo$pT), main=c('% CK+/CD166-/CD98- Cells'),
#      xlab=c('T Parameter'), ylab=c('% CK+/CD166+/CD98+ & CD166+/Cd98- Cells'))
# 
# plot(AllInfo$DblPosPct + AllInfo$CD166PosPct ~ as.factor(AllInfo$pT), main=c('% CK+/CD166+ cells'),
#      xlab=c('T Parameter'), ylab=c('% CK+/CD166-/CD98- Cells'))
# dev.off()
# 
# pdf("Plots of Expression by Count Against pTNM.pdf")
# 
# plot(AllInfo$DblPosCnt ~ as.factor(AllInfo$pT), main=c('Total CK+/CD166+/CD98+ Cells'),
#      xlab=c('T Parameter'), ylab=c('% CK+/CD166+/CD98+ Cells'))
# 
# plot(AllInfo$CD166PosCnt ~ as.factor(AllInfo$pT), main=c('Total CK+/CD166+/CD98- Cells'),
#      xlab=c('T Parameter'), ylab=c('% CK+/CD166+/CD98- Cells'))
# 
# plot(AllInfo$CD98PosCnt ~ as.factor(AllInfo$pT), main=c('Total CK+/CD166-/CD98+ Cells'),
#      xlab=c('T Parameter'), ylab=c('% CK+/CD166-/CD98+ Cells'))
# 
# plot(AllInfo$DblNegCnt ~ as.factor(AllInfo$pT), main=c('Total CK+/CD166-/CD98- Cells'),
#      xlab=c('T Parameter'), ylab=c('% CK+/CD166+/CD98+ & CD166+/Cd98- Cells'))
# 
# plot(AllInfo$DblPosCnt + AllInfo$CD166PosCnt ~ as.factor(AllInfo$pT), main=c('Total CK+/CD166+ cells'),
#      xlab=c('T Parameter'), ylab=c('% CK+/CD166-/CD98- Cells'))
# dev.off()
# 
# pdf("Plots of Expression by Count Over Total Cells Against pTNM.pdf")
# 
# plot(AllInfo$DblPosCnt/AllInfo$NucleiTotal * 100 ~ as.factor(AllInfo$pT), 
#      main=c('% CK+/CD166+/CD98+ Over Total Cells'),
#      xlab=c('T Parameter'), ylab=c('% CK+/CD166+/CD98+ Cells'))
# 
# plot(AllInfo$CD166PosCnt/AllInfo$NucleiTotal * 100 ~ as.factor(AllInfo$pT), 
#      main=c('% CK+/CD166+/CD98- Over Total Cells'),
#      xlab=c('T Parameter'), ylab=c('% CK+/CD166+/CD98- Cells'))
# 
# plot(AllInfo$CD98PosCnt/AllInfo$NucleiTotal * 100 ~ as.factor(AllInfo$pT), 
#      main=c('% CK+/CD166-/CD98+ Over Total Cells'),
#      xlab=c('T Parameter'), ylab=c('% CK+/CD166-/CD98+ Cells'))
# 
# plot(AllInfo$DblNegCnt/AllInfo$NucleiTotal * 100 ~ as.factor(AllInfo$pT), 
#      main=c('% CK+/CD166-/CD98- Over Total Cells'),
#      xlab=c('T Parameter'), ylab=c('% CK+/CD166+/CD98+ & CD166+/Cd98- Cells'))
# 
# plot((AllInfo$DblPosCnt + AllInfo$CD166PosCnt)/AllInfo$NucleiTotal * 100 ~ as.factor(AllInfo$pT), 
#      main=c('% CK+/CD166+ Over Total Cells'),
#      xlab=c('T Parameter'), ylab=c('% CK+/CD166-/CD98- Cells'))
# dev.off()
# 
# pdf("Plots of Expression by Count Over Total Cells Against Stage.pdf")
# 
# plot(AllInfo$DblPosCnt/AllInfo$NucleiTotal * 100 ~ AllInfo$Stage, main=c('% CK+/CD166+/CD98+ Cells'),
#      xlab=c('Stage'), ylab=c('% CK+/CD166+/CD98+ Cells'))
# 
# plot(AllInfo$CD166PosCnt/AllInfo$NucleiTotal * 100 ~ AllInfo$Stage, main=c('% CK+/CD166+/CD98- Cells'),
#      xlab=c('Stage'), ylab=c('% CK+/CD166+/CD98- Cells'))
# 
# plot(AllInfo$CD98PosCnt/AllInfo$NucleiTotal * 100 ~ AllInfo$Stage, main=c('% CK+/CD166-/CD98+ Cells'),
#      xlab=c('Stage'), ylab=c('% CK+/CD166-/CD98+ Cells'))
# 
# plot(AllInfo$DblNegCnt/AllInfo$NucleiTotal * 100 ~ AllInfo$Stage, main=c('% CK+/CD166-/CD98- Cells'),
#      xlab=c('Stage'), ylab=c('% CK+/CD166+/CD98+ & CD166+/Cd98- Cells'))
# 
# plot((AllInfo$DblPosCnt + AllInfo$CD166PosCnt)/AllInfo$NucleiTotal * 100 ~ AllInfo$Stage,
#      main=c('% CK+/CD166+ cells'), xlab=c('Stage'), ylab=c('% CK+/CD166-/CD98- Cells'))
# dev.off()
# 
# ### PDF plots of % positive and total positive Versus HPV status.
# 
# pdf("Plots of Expression by Count Over Total Cells Against HPV Status.pdf")
# 
# plot(AllInfo$DblPosCnt/AllInfo$NucleiTotal * 100 ~ AllInfo$HPV, main=c('% CK+/CD166+/CD98+ Cells'),
#      xlab=c('HPV'), ylab=c('% CK+/CD166+/CD98+ Cells'))
# 
# plot(AllInfo$CD166PosCnt/AllInfo$NucleiTotal * 100 ~ AllInfo$HPV, main=c('% CK+/CD166+/CD98- Cells'),
#      xlab=c('HPV'), ylab=c('% CK+/CD166+/CD98- Cells'))
# 
# plot(AllInfo$CD98PosCnt/AllInfo$NucleiTotal * 100 ~ AllInfo$HPV, main=c('% CK+/CD166-/CD98+ Cells'),
#      xlab=c('HPV'), ylab=c('% CK+/CD166-/CD98+ Cells'))
# 
# plot(AllInfo$DblNegCnt/AllInfo$NucleiTotal * 100 ~ AllInfo$HPV, main=c('% CK+/CD166-/CD98- Cells'),
#      xlab=c('HPV'), ylab=c('% CK+/CD166+/CD98+ & CD166+/Cd98- Cells'))
# 
# plot((AllInfo$DblPosCnt + AllInfo$CD166PosCnt)/AllInfo$NucleiTotal * 100 ~ AllInfo$HPV,
#      main=c('% CK+/CD166+ cells'), xlab=c('HPV'), ylab=c('% CK+/CD166+/CD98+ & CK+/CD166+/CD98-  Cells'))
# dev.off()
# 
# pdf("Plots of Expression by Total Count Against HPV Status.pdf")
# 
# plot(AllInfo$DblPosCnt ~ AllInfo$HPV, main=c('Total CK+/CD166+/CD98+ Cells'),
#      xlab=c('HPV'), ylab=c('Total CK+/CD166+/CD98+ Cells'))
# 
# plot(AllInfo$CD166PosCnt ~ AllInfo$HPV, main=c('Total CK+/CD166+/CD98- Cells'),
#      xlab=c('HPV'), ylab=c('Total CK+/CD166+/CD98- Cells'))
# 
# plot(AllInfo$CD98PosCnt ~ AllInfo$HPV, main=c('Total CK+/CD166-/CD98+ Cells'),
#      xlab=c('HPV'), ylab=c('Total CK+/CD166-/CD98+ Cells'))
# 
# plot(AllInfo$DblNegCnt ~ AllInfo$HPV, main=c('Total CK+/CD166-/CD98- Cells'),
#      xlab=c('HPV'), ylab=c('Total CK+/CD166+/CD98+ & CD166+/Cd98- Cells'))
# 
# plot((AllInfo$DblPosCnt + AllInfo$CD166PosCnt) ~ AllInfo$HPV,
#      main=c('Total CK+/CD166+ cells'), xlab=c('HPV'), ylab=c('Total CK+/CD166+/CD98+ & CK+/CD166+/CD98- Cells'))
# dev.off()
# 
# pdf("Plots of Expression by Percentage Positive Against HPV Status.pdf")
# 
# plot(AllInfo$DblPosPct ~ AllInfo$HPV, main=c('% CK+/CD166+/CD98+ Cells'),
#      xlab=c('HPV'), ylab=c('% CK+/CD166+/CD98+ Cells'))
# 
# plot(AllInfo$CD166PosPct ~ AllInfo$HPV, main=c('% CK+/CD166+/CD98- Cells'),
#      xlab=c('HPV'), ylab=c('% CK+/CD166+/CD98- Cells'))
# 
# plot(AllInfo$CD98PosPct ~ AllInfo$HPV, main=c('% CK+/CD166-/CD98+ Cells'),
#      xlab=c('HPV'), ylab=c('% CK+/CD166-/CD98+ Cells'))
# 
# plot(AllInfo$DblNegPct ~ AllInfo$HPV, main=c('% CK+/CD166-/CD98- Cells'),
#      xlab=c('HPV'), ylab=c('% CK+/CD166+/CD98+ & CD166+/Cd98- Cells'))
# 
# plot((AllInfo$DblPosPct + AllInfo$CD166PosPct) ~ AllInfo$HPV,
#      main=c('% CK+/CD166+ cells'), xlab=c('HPV'), ylab=c('% CK+/CD166+/CD98+ & CK+/CD166+/CD98- Cells'))
# dev.off()
# 
# ### Data used for graphs
# write.csv(AllInfo, "Data used for Graphs.csv", row.names=F)


