# Script for analysis of significance in EdyTMA data
# 
# Step 1: import data
# Step 2: append smoker status, re-write data (change file name)
# Step 3: Table pTNM status. pTNM status may be used as category
# (not just pT)
# Step 4: use Wilcoxon Ranked test of Unity to compare HPV +/- by
# marker presence (percentile). Pairwise across all markers.
# Consider ANOVA for the multivariate analysis across pTNM? 

setwd("E:/Cell Profiler/20150825 Images for Full Spot Analysis/3rd Run outputs/")
# TMAdat <- read.csv("Data used for Graphs.csv", header=T)
# map <- read.csv("OPSCC TMA Map.csv", header=T)
# TMAdat$Smoker <- map$Smoker
# write.csv(TMAdat, "20151110 Data used for analysis.csv", row.names=F)

# setwd("~/disk1/EdyTMA")
TMAdat <- read.csv("2016 Run (hard cluster size) - Data used for Graphs - bad spots dropped.csv", stringsAsFactors=F)

# Caveat to tabling: no smoker status for controls. Additionally one pt is
# ambiguous. 
table(TMAdat$pTNM, TMAdat$Smoker)
table(TMAdat$pTNM, TMAdat$HPV)
table(TMAdat$pTNM, TMAdat$HPV, TMAdat$Smoker)

# wilcox.test will be used for nonparametric testing with modifications as
# necessary.
# colnames(TMAdat)
# [1] "SpotID"         "Donor.Block.ID" "Age"            "HPV"           
# [5] "pTNM"           "Stage"          "DblNegPct"      "CD166PosPct"   
# [9] "CD98PosPct"     "DblPosPct"      "DblNegCnt"      "CD166PosCnt"   
# [13] "CD98PosCnt"     "DblPosCnt"      "NucleiTotal"    "xy"            
# [17] "pT"             "Smoker"  

# find numeric columns
numcol = logical(ncol(TMAdat))
# cutoff for max proportion of NAs in a numeric column  
NAcut = .10
for(i in 1:ncol(TMAdat) ){
  NA0 = sum(is.na(TMAdat[,i]))
  NA1 = sum(is.na(as.numeric(TMAdat[,i])))
  if( (NA1-NA0)/(nrow(TMAdat)-NA0)<NAcut ){ numcol[i] = TRUE }
}

# masks for HPV+ and HPV-, omitting controls
# indirectly, testing for unknown smoking status removes controls.  
#  it also removes one incompletely characterized patient
HPVpos = TMAdat$Smoker != "?" & TMAdat$HPV == "+"
HPVneg = TMAdat$Smoker != "?" & TMAdat$HPV == "-"

# loop over numeric columns to run Wilcoxon tests
for(i in which(numcol) ){
  med.pos = median(TMAdat[HPVpos,i],na.rm=T)
  med.neg = median(TMAdat[HPVneg,i],na.rm=T)
  w.obj = wilcox.test(TMAdat[HPVpos,i],TMAdat[HPVneg,i],conf.int=TRUE)
  cat(names(TMAdat)[i], signif(c(w.obj$p.value,w.obj$estimate,
                                 med.pos,med.neg),3), w.obj$statistic,
                                  "\n", sep="\t")
}

# Pulled ugly and incorrect code (ran tests that were inclusive of
# controls in HPV- subgroup, incorrectly defined CD166+). See below
# for corrections on the Control group comparison.

# No real success by pairwise testing - CD166PosPct by Smoker or 
# by HPV yields p close to 0.05 (0.14 and 0.06 respectively)...
# Also I just realized in my tests I was including Controls in 
# the HPV- group.

# This list approach for collecting wilcox results is ugly as sin but 
# it works...for now.

markers <- c("DblPosPct", "CD166PosPct", "CD98PosPct", "DblNegPct"); 
byHPV <- c(); bySmoker <- c();

for(i in markers){
  byHPV <- list(byHPV, wilcox.test(get(i) ~ HPV, 
                 data=TMAdat[c(which(TMAdat$pTNM != "Ctrl")),]))
}

for(i in markers){
  bySmoker <- list(bySmoker, wilcox.test(get(i) ~ Smoker, 
                    data=TMAdat[c(which(TMAdat$Smoker != "?")),]))
}

## Using Julja's code to iterate over Smoker status by column and look for significance.

Smopos = TMAdat$Smoker == "+"
Smoneg = TMAdat$Smoker == "-"

for(i in which(numcol) ){
  med.pos = median(TMAdat[Smopos,i],na.rm=T)
  med.neg = median(TMAdat[Smoneg,i],na.rm=T)
  w.obj = wilcox.test(TMAdat[Smopos,i],TMAdat[Smoneg,i],conf.int=TRUE)
  cat(names(TMAdat)[i], signif(c(w.obj$p.value,w.obj$estimate,
                                 med.pos,med.neg),3), w.obj$statistic,
      "\n", sep="\t")
}

### Sophie asked me to compare Smoker status versus HPV and to the best of my understanding I believe
### this means to essentially create a fourway comparison: HPVneg&Smokerneg Vs. HPVNeg&SmokerPos and 
### HPVPos&SmokerNeg vs. HPVPos&SmokerPos. So that's what I'll do below by using boolean completion in
### Julja's code to winnow for the above categories.


for(i in which(numcol) ){
  med.pos = median(TMAdat[Smopos & HPVneg,i],na.rm=T)
  med.neg = median(TMAdat[Smoneg & HPVneg,i],na.rm=T)
  w.obj = wilcox.test(TMAdat[Smopos & HPVneg,i],TMAdat[Smoneg & HPVneg,i],conf.int=TRUE)
  cat(names(TMAdat)[i], signif(c(w.obj$p.value,w.obj$estimate,
                                 med.pos,med.neg),3), w.obj$statistic,
      "\n", sep="\t")
}

for(i in which(numcol) ){
  med.pos = median(TMAdat[Smopos & HPVpos,i],na.rm=T)
  med.neg = median(TMAdat[Smoneg & HPVpos,i],na.rm=T)
  w.obj = wilcox.test(TMAdat[Smopos & HPVpos,i],TMAdat[Smoneg & HPVpos,i],conf.int=TRUE)
  cat(names(TMAdat)[i], signif(c(w.obj$p.value,w.obj$estimate,
                                 med.pos,med.neg),3), w.obj$statistic,
      "\n", sep="\t")
}

### Comparison against controls
for(i in which(numcol) ){
  med.pos = median(TMAdat[Smopos & HPVpos,i],na.rm=T)
  med.neg = median(TMAdat[Smoneg & HPVpos,i],na.rm=T)
  w.obj = wilcox.test(TMAdat[which(TMAdat$Stage != "Ctrl"),i],TMAdat[which(TMAdat$Stage == "Ctrl"),i],conf.int=TRUE)
  cat(names(TMAdat)[i], signif(c(w.obj$p.value,w.obj$estimate,
                                 med.pos,med.neg),3), w.obj$statistic,
      "\n", sep="\t")
}



pdf("CD166 Positivity in HNSCC versus Normal.pdf")
boxplot(data=TMAdat, log(CD166PosCnt) ~ Cancer,
        ylab="log10 Count CD166+", xlab="Cancer Status", 
        main="CD166 Positivity in HNSCC versus Normal")
stripchart(data=TMAdat, log(CD166PosCnt) ~ Cancer, 
           vertical=TRUE,method="jitter", pch=21,
           col= "maroon", bg = "bisque", add=TRUE)
dev.off()
pdf("CD98 Positivity in HNSCC versus Normal.pdf")
boxplot(data=TMAdat, log(CD98PosCnt) ~ Cancer,
        ylab="log10 Count CD98+", xlab="Cancer Status", 
        main="CD98 Positivity in HNSCC versus Normal")
stripchart(data=TMAdat, log(CD98PosCnt) ~ Cancer, 
           vertical=TRUE, method="jitter", pch=21,
           col= "maroon", bg = "bisque", add=TRUE)
dev.off()

boxplot(data=TMAdat, log(CD166PosCnt) ~ HPVpos,
        ylab="log10 Count CD166+", xlab="HPV Status", 
        main="CD166 Log Counts by HPV Status")
stripchart(data=TMAdat, log(CD166PosCnt) ~ HPVneg, 
           vertical=TRUE,method="jitter", pch=21,
           col= "maroon", bg = "bisque", add=TRUE)

C166 <- (TMAdat$CD166PosPct+TMAdat$DblPosPct)
names(C166) <- TMAdat$Stage


boxplot(data=TMAdat, C166[1:38] ~ HPV[1:38],
        ylab="log10 Count CD166+", xlab="HPV Status", 
        main="CD166 Log Counts by HPV Status")
stripchart(data=TMAdat, (CD166PosPct+DblPosPct)[1:38] ~ HPVneg[1:38], 
           vertical=TRUE,method="jitter", pch=21,
           col= "maroon", bg = "bisque", add=TRUE)

hist(data=TMAdat, C166)
summary(C166)

# FDR correction for controls
FDRtest <- c(.249, .392, .91, .987, .017, .0115, .00156, .153, .367)
p.adjust(FDRtest, method="BH")

# FDR Correction for smoker positive versus negative
FDRtestSmokers <- c(.412, .0319, .549,1,.921,.046,.73,.289,.443)
p.adjust(FDRtestSmokers, method="BH")
