# Alvarez, Robertson, et al. 2019
# Code for working with DWH epiGBS data
# Aug 23, 2019

# Set working directory here
setwd("/Volumes/ANALYSIS3/2016SpAlt_epiGBS/RedoAnalysis")
library(tidyverse)

# This function formats the SNP table from the VCF file to create allele frequencies. It also implements a depth cutoff of 10x
PrepGenetic<-function(Depth=10){
  #bcftools query -f '%CHROM %POS  %REF  %ALT [ %SAMPLE=%RO] [ %SAMPLE=%DP]\n' snp.vcf.gz > AD.table
  
  ## Preliminary steps to get data organized and ready to convert
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(stringr)
  
  allele.depth2 <-fread("AD.table", header=FALSE,data.table = FALSE)
  
  temp<-data.frame(t(allele.depth2[1,-c(1:4)]))
  temp<-data.frame(str_split_fixed(temp$X1,"=",2))
  colnames(allele.depth2)<-c("Chr","Pos","Ref","Alt",as.character(temp$X1))
  
  allele.count <- cbind(allele.depth2[1:4],allele.depth2[,c(5:52)])
  allele.depth <- cbind(allele.depth2[1:4],allele.depth2[,c(53:ncol(allele.depth2))])
  fwrite(allele.depth,"depthMatrix.txt",sep="\t",quote=FALSE,row.names = FALSE)
  
  for(i in 5:ncol(allele.count)){
    print(i)
    allele.count[,i]<-gsub(paste(colnames(allele.count)[i],"=",sep=""),"",allele.count[,i])
    allele.depth[,i]<-gsub(paste(colnames(allele.depth)[i],"=",sep=""),"",allele.depth[,i])
  }
  
  #Depth=10
  
  tempmat<-data.matrix(allele.depth[,-c(1:4)])
  tempmat[tempmat < Depth] <- NA
  
  AlleleFreq<-data.frame(allele.count[,c(1:4)],(data.matrix(allele.count[,-c(1:4)])/tempmat))
  
  fwrite(AlleleFreq,paste("Allele.frequency.",Depth,".txt",sep=""),quote=FALSE,row.names = FALSE,sep="\t")
  
}
# This function filters loci and samples with too much missing data, then imputes the rest.
FilterGenetic<-function(){
  library(data.table)
  library(tidyr)
  ######
  
  AlleleFreq<-fread("Allele.frequency.10.txt",data.table = FALSE)
  
  AlleleFreq$LocusMissing<-rowSums(is.na(AlleleFreq[,-c(1:4)]))/ncol(AlleleFreq[,-c(1:4)])
  
  AlleleFreq<-dplyr::filter(AlleleFreq,LocusMissing<1)
  AlleleFreq$LocusMissing<-NULL
  
  InvariantFilter<-function(){
    # Filter out fully methylated invariant sites
    print("Filtering 100% methylated....")
    AlleleFreq2<-AlleleFreq[,-c(1:4)]
    AlleleFreq2[is.na(AlleleFreq2)] <- 1
    AlleleFreq$Invariant<-rowSums(AlleleFreq2==1)/ncol(AlleleFreq2)
    rm(AlleleFreq2)
    AlleleFreq<-dplyr::filter(AlleleFreq,!Invariant==1)
    AlleleFreq$Invariant<-NULL
    
    # Filter out not methylated invariant sites
    print("Filtering 0% methylated...")
    AlleleFreq2<-AlleleFreq[,-c(1:4)]
    AlleleFreq2[is.na(AlleleFreq2)] <- 0
    AlleleFreq$Invariant<-rowSums(AlleleFreq2==0)/ncol(AlleleFreq2)
    rm(AlleleFreq2)
    AlleleFreq<-dplyr::filter(AlleleFreq,!Invariant==1)
    AlleleFreq$Invariant<-NULL
    return(AlleleFreq)
  }
  AlleleFreq<-InvariantFilter()
  
  SampleFilter<-function(MissingMax=0.80){
    SampMissing<-c()
    for(i in 5:ncol(AlleleFreq)){
      SampMissing<-c(SampMissing,sum(is.na(AlleleFreq[,i]))/nrow(is.na(AlleleFreq)))
    }
    SampMissing<-c(0,0,0,0,SampMissing)
    SampMissing
    sum(SampMissing>.80)
    AlleleFreq<-AlleleFreq[,!(SampMissing>MissingMax)]
    return(AlleleFreq)
  }
  AlleleFreq<-SampleFilter(0.8)
  
  LocusFilter<-function(MissingMax=0.5){
    AlleleFreq$LocusMissing<-rowSums(is.na(AlleleFreq[,-c(1:4)]))/ncol(AlleleFreq[,-c(1:4)])
    summary(AlleleFreq$LocusMissing)
    AlleleFreq<-dplyr::filter(AlleleFreq,LocusMissing<=MissingMax)
    AlleleFreq$LocusMissing<-NULL
    return(AlleleFreq)
  }
  AlleleFreq<-LocusFilter(0.5)
  
  SampleFilter<-function(MissingMax=0.80){
    SampMissing<-c()
    for(i in 5:ncol(AlleleFreq)){
      SampMissing<-c(SampMissing,sum(is.na(AlleleFreq[,i]))/nrow(is.na(AlleleFreq)))
    }
    SampMissing<-c(0,0,0,0,SampMissing)
    SampMissing
    sum(SampMissing>.80)
    AlleleFreq<-AlleleFreq[,!(SampMissing>MissingMax)]
    return(AlleleFreq)
  }
  AlleleFreq<-SampleFilter(0.65)
  
  LocusFilter<-function(MissingMax=0.5){
    AlleleFreq$LocusMissing<-rowSums(is.na(AlleleFreq[,-c(1:4)]))/ncol(AlleleFreq[,-c(1:4)])
    summary(AlleleFreq$LocusMissing)
    AlleleFreq<-dplyr::filter(AlleleFreq,LocusMissing<=MissingMax)
    AlleleFreq$LocusMissing<-NULL
    return(AlleleFreq)
  }
  #AlleleFreq<-LocusFilter(0.2)
  
  AlleleFreq<-InvariantFilter()
  
  PlotMissing<-function(){
    PlotDat<-AlleleFreq
    PlotDat<-data.frame(ChrPos=paste(PlotDat$Chr,PlotDat$Pos,sep=":"),PlotDat[,-c(1:4)])
    library(reshape2)
    PlotDat2<-reshape2::melt(PlotDat,id.vars=colnames(PlotDat)[1])
    PlotDat2$variable<-gsub("SPALT_","",PlotDat2$variable)
    PlotDat2$value<-is.na(PlotDat2$value)
    PlotDat2$value<-gsub(FALSE,"Present",PlotDat2$value)
    PlotDat2$value<-gsub(TRUE,"Imputed",PlotDat2$value)
    PlotDat3<-data.frame(table(PlotDat2$variable,PlotDat2$value))
    PlotDat4 <- PlotDat3  %>%
      group_by(Var1, Var2) %>%
      summarise(n = sum(Freq)) %>%
      mutate(percentage = n / sum(n))
    
    library(ggplot2)
    library(ggthemes)
    library(ggsci)
    library(tidyverse)
    
    #pdf("CrazyPlot.pdf",height=10,width=8)
    #ggplot(PlotDat2,aes(y=as.numeric(as.factor(ChrPos)),x=variable,color=value)) + 
    #  geom_point(shape=15,alpha=0.5) + #coord_polar() +
    #  theme_bw() + scale_color_npg(name="Value handling") +
    #  theme(axis.title.y=element_blank(),
    #        axis.text.y=element_blank(),
    #        axis.ticks.y=element_blank())
    #dev.off()
    
    # Samples
    #pdf("SampleCoverage.pdf",height=8,width=8)
    GenSamp<-ggplot(PlotDat4,aes(x=Var1,y=percentage,fill=Var2)) + 
      geom_bar(stat="identity") + labs(x="Sample",y="Percentage of loci") +
      geom_hline(yintercept = 0.2,linetype = "dashed",color="white")+
      theme_minimal() + scale_fill_hc(name="Data handling") + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = "bottom") +
      theme(axis.title.x=element_text(vjust=-2)) +
      theme(axis.title.y=element_text(angle=90, vjust=-0.5)) +
      theme(plot.title=element_text(size=15, vjust=3)) +
      theme(plot.margin = unit(c(1,1,1,1), "cm"))
    #dev.off()
    
    PlotDat3<-data.frame(table(PlotDat2$ChrPos,PlotDat2$value))
    PlotDat4 <- PlotDat3  %>%
      group_by(Var1, Var2) %>%
      summarise(n = sum(Freq)) %>%
      mutate(percentage = n / sum(n))
    
    #pdf("LocusCoverage.pdf",height=8,width=8)
    GenLoc<-ggplot(PlotDat4[1:10000,],aes(x=as.numeric(as.factor(Var1)),y=percentage,fill=Var2)) + 
      geom_area() + labs(x="Locus",y="Percentage of data at a locus") +
      geom_hline(yintercept = 0.5,linetype = "dashed",color="white") +
      theme_minimal() + scale_fill_hc(guide=FALSE) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      theme(axis.title.x=element_text(vjust=-2)) +
      theme(axis.title.y=element_text(angle=90, vjust=-0.5)) +
      theme(plot.title=element_text(size=15, vjust=3)) +
      theme(plot.margin = unit(c(1,1,1,1), "cm"))
    #dev.off()
    
    library(cowplot)
    #pdf("FigS2.pdf",height=8,width=8)
    #print(plot_grid(Samp,Loc,rel_heights = c(1.2, 1),
    #                #rel_widths = c(1, 1.2),
    #                nrow = 2,labels="AUTO"))
    #dev.off()
    save(GenLoc,GenSamp,file="GenPlots.RData")
    
  }
  PlotMissing()
  
  
  library(impute)
  tempmat<-data.matrix(AlleleFreq[,-c(1:4)])
  imputed<-impute.knn(tempmat,colmax=0.9)
  imputed<-data.frame(AlleleFreq[,c(1:4)],imputed$data)
  
  fwrite(imputed,"AlleleFreq.txt",quote=FALSE,sep="\t",row.names = FALSE)
  
}
# This function uses STAMPP to get the genetic relatedness matrix (based on Nei's d).
PopGenPrep_stampp<-function(){
  ###
  
  library(StAMPP)
  library(data.table)
  freq<-fread("AlleleFreq.txt")
  #freq$SPALT_MS02<-NULL
  #freq$SPALT_MS02
  
  stamped<-data.frame(Sample=colnames(freq)[-c(1:4)],Pop=NA,Ploidy=rep(6,ncol(freq)-4),Format="freq",t(freq[,-c(1:4)]))
  #stamped$Pop<-c(1,2,5,4,6,1,4,6,1,2,5,3,4,6,1,2,5,6,1,3,4,6,1,2,5,4,6,1,2,5,3,4,6,1,2,5,4,6)
  stamped$Pop<-c(1,2,5,4,6,1,4,6,1,2,5,3,4,6,1,2,5,6,1,3,4,6,1,2,5,4,6,1,2,5,3,4,6,1,2,5,4,6)
  stamped$Pop<-as.character(stamped$Pop)
  # 1= GIN1, 2=GIN2, 5=MSN, 4=GIO2, 6=MSO, 3=GIO1
  
  
  allele.freq <- stamppConvert(stamped, type='r')
  allele.D.ind <- stamppNeisD(allele.freq, pop=FALSE)
  allele.D.pop <- stamppNeisD(allele.freq, pop=TRUE)
  allele.relationship<-stamppGmatrix(allele.freq)
  fwrite(allele.relationship,"Gmatrix.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names = FALSE)
  
  allele.D.ind<-data.frame(allele.D.ind)
  colnames(allele.D.ind) <- c(row.names(allele.D.ind))
  #names(allele.D.pop)[1:45] <- c(row.names(allele.D.pop))
  fwrite(allele.D.ind, "genetic.distance.ind.txt", quote=FALSE, sep="\t", row.names = TRUE, col.names = FALSE)
  
  allele.D.pop<-data.frame(allele.D.pop)
  colnames(allele.D.pop)<-c("GIN1","GIN2","GIO1","GIO2","MSN","MSO")
  fwrite(allele.D.pop, "genetic.distance.pop.txt", quote=FALSE, sep="\t", row.names = FALSE, col.names = TRUE)
  
  #Amova <- stamppAmova(allele.D.ind, allele.freq, 10000)
  
  allele.fst <- stamppFst(allele.freq, nboots=1, percent=95, nclusters=6)
  save(list = c("allele.D.ind","allele.D.pop","allele.fst"), file = "GeneticDistanceAF.RData",compress = TRUE)
  load("GeneticDistanceAF.RData")
}
# This function formats the methylation table from the BED file to create methylation frequencies. It also implements a depth cutoff of 10x
PrepMeth<-function(Depth=10){
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(stringr)
  allele.depth2<-fread("methylation.bed",data.table = FALSE)
  
  allele.count <- cbind(allele.depth2[1:4],dplyr::select(allele.depth2,contains("methylated")))
  allele.depth <- cbind(allele.depth2[1:4],dplyr::select(allele.depth2,contains("total")))
  
  
  
  tempmat<-data.matrix(allele.depth[,-c(1:4)])
  tempmat[tempmat < Depth] <- NA
  
  AlleleFreq<-data.frame(allele.count[,c(1:4)],(data.matrix(allele.count[,-c(1:4)])/tempmat))
  
  fwrite(AlleleFreq,paste("Meth.frequency.",Depth,".txt",sep=""),quote=FALSE,row.names = FALSE,sep="\t")
  
  
  
}
# This function filters loci and samples with too much missing data, then imputes the rest.
FilterMeth<-function(){
  ######
  
  library(data.table)
  library(dplyr)
  library(stringr)
  
  AlleleFreq<-fread("Meth.frequency.10.txt",data.table = FALSE)
  Genetic<-fread("AlleleFreq.txt",data.table = FALSE)
  
  
  AlleleFreq$LocusMissing<-rowSums(is.na(AlleleFreq[,-c(1:4)]))/ncol(AlleleFreq[,-c(1:4)])
  
  AlleleFreq<-dplyr::filter(AlleleFreq,LocusMissing<1)
  AlleleFreq$LocusMissing<-NULL
  
  InvariantFilter<-function(){
    # Filter out fully methylated invariant sites
    print("Filtering 100% methylated....")
    AlleleFreq2<-AlleleFreq[,-c(1:4)]
    AlleleFreq2[is.na(AlleleFreq2)] <- 1
    AlleleFreq$Invariant<-rowSums(AlleleFreq2==1)/ncol(AlleleFreq2)
    rm(AlleleFreq2)
    AlleleFreq<-dplyr::filter(AlleleFreq,!Invariant==1)
    AlleleFreq$Invariant<-NULL
    
    # Filter out not methylated invariant sites
    print("Filtering 0% methylated...")
    AlleleFreq2<-AlleleFreq[,-c(1:4)]
    AlleleFreq2[is.na(AlleleFreq2)] <- 0
    AlleleFreq$Invariant<-rowSums(AlleleFreq2==0)/ncol(AlleleFreq2)
    rm(AlleleFreq2)
    AlleleFreq<-dplyr::filter(AlleleFreq,!Invariant==1)
    AlleleFreq$Invariant<-NULL
    return(AlleleFreq)
  }
  AlleleFreq<-InvariantFilter()
  
  colnames(AlleleFreq)<-gsub("_methylated","",colnames(AlleleFreq))
  
  AlleleFreq<-data.frame(AlleleFreq[,c(1:4)],dplyr::select(AlleleFreq, one_of(colnames(Genetic)[-c(1:4)])))
  
  SampMissing<-c()
  for(i in 5:ncol(AlleleFreq)){
    SampMissing<-c(SampMissing,sum(is.na(AlleleFreq[,i]))/nrow(is.na(AlleleFreq)))
  }
  SampMissing<-c(0,0,0,0,SampMissing)
  SampMissing
  sum(SampMissing>.80)
  #AlleleFreq<-AlleleFreq[,!(SampMissing>.80)]
  
  AlleleFreq$LocusMissing<-rowSums(is.na(AlleleFreq[,-c(1:4)]))/ncol(AlleleFreq[,-c(1:4)])
  summary(AlleleFreq$LocusMissing)
  AlleleFreq<-dplyr::filter(AlleleFreq,LocusMissing<=.5)
  AlleleFreq$LocusMissing<-NULL
  
  AlleleFreq<-InvariantFilter()
  
  AlleleFreq2<-data.frame(ID=paste(AlleleFreq$chr,AlleleFreq$pos,sep=":"),AlleleFreq)
  
  PlotMissing<-function(){
    PlotDat<-AlleleFreq2[,-c(2:5)]
    #PlotDat<-data.frame(ChrPos=paste(PlotDat$Chr,PlotDat$Pos,sep=":"),PlotDat[,-c(1:4)])
    library(reshape2)
    PlotDat2<-reshape2::melt(PlotDat,id.vars=colnames(PlotDat)[1])
    PlotDat2$variable<-gsub("SPALT_","",PlotDat2$variable)
    PlotDat2$value<-is.na(PlotDat2$value)
    PlotDat2$value<-gsub(FALSE,"Present",PlotDat2$value)
    PlotDat2$value<-gsub(TRUE,"Imputed",PlotDat2$value)
    PlotDat3<-data.frame(table(PlotDat2$variable,PlotDat2$value))
    PlotDat4 <- PlotDat3  %>%
      group_by(Var1, Var2) %>%
      summarise(n = sum(Freq)) %>%
      mutate(percentage = n / sum(n))
    
    library(ggplot2)
    library(ggthemes)
    library(ggsci)
    
    #pdf("CrazyPlot.pdf",height=10,width=8)
    #ggplot(PlotDat2,aes(y=as.numeric(as.factor(ChrPos)),x=variable,color=value)) + 
    #  geom_point(shape=15,alpha=0.5) + #coord_polar() +
    #  theme_bw() + scale_color_npg(name="Value handling") +
    #  theme(axis.title.y=element_blank(),
    #        axis.text.y=element_blank(),
    #        axis.ticks.y=element_blank())
    #dev.off()
    
    # Samples
    #pdf("SampleCoverage.pdf",height=8,width=8)
    MethSamp<-ggplot(PlotDat4[1:10000,],aes(x=Var1,y=percentage,fill=Var2)) + 
      geom_bar(stat="identity") + labs(x="Sample",y="Percentage of loci") +
      geom_hline(yintercept = 0.2,linetype = "dashed",color="white")+
      theme_minimal() + scale_fill_hc(name="Data handling") + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = "bottom") +
      theme(axis.title.x=element_text(vjust=-2)) +
      theme(axis.title.y=element_text(angle=90, vjust=-0.5)) +
      theme(plot.title=element_text(size=15, vjust=3)) +
      theme(plot.margin = unit(c(1,1,1,1), "cm"))
    #dev.off()
    
    PlotDat3<-data.frame(table(PlotDat2$ID,PlotDat2$value))
    PlotDat4 <- PlotDat3  %>%
      group_by(Var1, Var2) %>%
      summarise(n = sum(Freq)) %>%
      mutate(percentage = n / sum(n))
    PlotDat4<-data.frame(PlotDat4)
    
    #pdf("LocusCoverage.pdf",height=8,width=8)
    MethLoc<-ggplot(PlotDat4[1:10000,],aes(x=as.numeric(as.factor(Var1)),y=percentage,fill=Var2)) + 
      geom_area() + labs(x="Locus",y="Percentage of data at a locus") +
      geom_hline(yintercept = 0.5,linetype = "dashed",color="white") +
      theme_minimal() + scale_fill_hc(guide=FALSE) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      theme(axis.title.x=element_text(vjust=-2)) +
      theme(axis.title.y=element_text(angle=90, vjust=-0.5)) +
      theme(plot.title=element_text(size=15, vjust=3)) +
      theme(plot.margin = unit(c(1,1,1,1), "cm"))
    #dev.off()
    
    library(cowplot)
    #pdf("FigS2.pdf",height=8,width=8)
    #print(plot_grid(Samp,Loc,rel_heights = c(1.2, 1),
    #                #rel_widths = c(1, 1.2),
    #                nrow = 2,labels="AUTO"))
    #dev.off()
    
    save(MethLoc,MethSamp,file="MethPlots.RData")
    
  }
  PlotMissing()
  
  out<-data.frame()
  library(impute)
  
  while (nrow(AlleleFreq2)>0){
    print("looping")
    
    if(nrow(AlleleFreq2)>round(nrow(AlleleFreq)/10)){
      tempdf<-dplyr::sample_n(AlleleFreq2,round(nrow(AlleleFreq)/10))
    } else {
      tempdf<-AlleleFreq2
    }
    
    tempmat<-data.matrix(tempdf[,-c(1:5)])
    
    imputed<-impute.knn(tempmat,colmax=0.8,maxp=4000)
    imputed<-data.frame(tempdf[,c(1:5)],imputed$data)
    out<-rbind(out,imputed)
    AlleleFreq2<-dplyr::filter(AlleleFreq2,!ID %in% out$ID)
  }
  
  imputed<-out
  
  fwrite(imputed,"MethFreq.txt",quote=FALSE,sep="\t",row.names = FALSE)
  
  # macau
  
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(stringr)
  allele.depth2<-fread("methylation.bed",data.table = FALSE)
  
  #freq<-fread("MethFreq.txt",data.table = FALSE)[,-c(1:2)]
  freq<-fread("MethFreq.txt",data.table = FALSE)
  allele.depth <- cbind(allele.depth2[1:4],dplyr::select(allele.depth2,contains("total")))
  #rm(allele.depth2)
  colnames(allele.depth)<-gsub("_total","",colnames(allele.depth))
  colnames(allele.depth)<-gsub("-",".",colnames(allele.depth),fixed = TRUE)
  #colnames(allele.depth)
  allele.depth<-dplyr::select(allele.depth, one_of(colnames(freq)))
  allele.depth$chrpos<-paste(allele.depth$chr,allele.depth$pos,allele.depth$context,sep=":")
  #allele.depth<-data.frame(chrpos=paste(allele.depth$chr,allele.depth$pos,allele.depth$context,sep=":"),allele.depth[,-c(1:4)])
  freq$chrpos<-paste(freq$chr,freq$pos,freq$context,sep=":")
  #freq<-data.frame(chrpos=paste(freq$chr,freq$pos,freq$context,sep=":"),freq[,-c(1:4)])
  allele.depth<-dplyr::filter(allele.depth, chrpos %in% freq$chrpos)
  #allele.depth$chrpos<-NULL
  
  tmp<-data.matrix(allele.depth[,-c(1:4,ncol(allele.depth))])
  allele.depth[,-c(1:4,ncol(allele.depth))]<-tmp
  rm(tmp)
  #for(i in 2:ncol(allele.depth)){
  #  print(i)
  #  allele.depth[,i]<-as.numeric(as.character(allele.depth[,i]))
  #}
  
  #allele.depth<-allele.depth %>% mutate_all(~replace(., is.na(.), 10))
  allele.depth[is.na(allele.depth)]<-10
  allele.depth<-allele.depth[order(match(allele.depth$chrpos,freq$chrpos)),]
  
  allele.depth<-allele.depth[,-c(1:4)]
  allele.depth<-allele.depth[,c(ncol(allele.depth),1:(ncol(allele.depth)-1))]
  freq<-freq[,-c(1:5)]
  freq<-freq[,c(ncol(freq),1:(ncol(freq)-1))]
  
  fwrite(allele.depth,"counts.macau.txt",quote=FALSE,row.names = FALSE,sep = "\t")
  
  tempmat<-data.matrix(freq[,-1])
  tempdep<-data.matrix(allele.depth[,-1])
  tempmeth<-data.matrix(tempmat*tempdep)
  tempmeth<-round(tempmeth)
  #sum(tempmeth>tempdep)
  freq[,-1]<-tempmeth
  fwrite(freq,"mcounts.macau.txt",quote=FALSE,row.names = FALSE,sep = "\t")
  
}
# This function makes Figure S2
MakeFilterPlots<-function(){
  load("GenPlots.RData")
  load("MethPlots.RData")
  
  library(cowplot)
  pdf("FigS2.pdf",height=8,width=12)
  print(plot_grid(GenSamp,MethSamp,GenLoc,MethLoc,nrow = 2,ncol=2,labels="AUTO",rel_heights = c(1.2,1.2,1,1)))
  dev.off()
}


# Load and format analysis files - SNP matrix, methylation frequencies, phenotypes

library(data.table)
allele<-fread("AlleleFreq.txt",data.table = FALSE)
pheno<-fread("env_spartina2.csv",data.table=FALSE)

pheno$treatment<-gsub("nooil","no oil",pheno$treatment)


colnames(allele)<-gsub(".","-",colnames(allele),fixed = TRUE)
pheno<-dplyr::filter(pheno,sampleID %in% colnames(allele))
pheno<-pheno[match(colnames(allele)[-c(1:4)], pheno$sampleID),]
meth<-fread("MethFreq.txt",data.table = FALSE)[,-1]

GetLFs<-function(df=allele,K=3){
  library(lfmm)
  Y <- df[,-c(1:4)]
  pc <- prcomp(t(Y))
  plot((pc$sdev[1:20])^2, xlab = 'PC', ylab = "Variance explained")
  points(4,(pc$sdev[4])^2, type = "h", lwd = 3, col = "blue")
  # two latent factors
  
  Y<-data.frame(t(df[,-c(1:4)]))
  mod.lfmm <- lfmm_ridge(Y = Y, 
                         X = as.numeric(as.factor(pheno$treatment)), 
                         K = K)
  return(mod.lfmm$U)
}
LF.allele<-GetLFs(df=allele,K=5) 
LF.meth<-GetLFs(df=meth,K=4)

# Depth PCA
DepthPCA<-function(){
  depth<-fread("depthMatrix.txt",data.table = FALSE)
  depth$chrpos<-paste(depth$Chr,depth$Pos,sep=":")
  allele$chrpos<-paste(allele$Chr,allele$Pos,sep=":")
  depth<-dplyr::filter(depth,chrpos %in% allele$chrpos)
  depth<-dplyr::select(depth, one_of(colnames(allele)))
  depth$chrpos<-NULL
  depth[,-c(1:4)]<-data.matrix(depth[,-c(1:4)])
  depth[is.na(depth)]<-10
  
  topc<-data.matrix(t(depth[,-c(1:4)]))
  pcdata<-prcomp(topc)
  sum<-summary(pcdata)
  sum<-data.frame((sum$importance))
  sum<-data.frame(var=rownames(sum),sum)[2,]
  sum<-reshape2::melt(sum,id.vars="var")
  pcdata<-data.frame(pheno,pcdata$x)
  pcdata$site<-gsub("MSO9","MSO",pcdata$site)
  
  library(RColorBrewer)
  colourCount = length(unique(sum$variable))
  getPalette = colorRampPalette(brewer.pal(8, "Set1"))
  
  library(ggplot2)
  library(ggthemes)
  PC<-ggplot(pcdata,aes(x=PC1,y=PC2,color=site,shape=treatment)) + 
    geom_point() + theme_bw() + scale_color_few() + guides(color=guide_legend(title="Site")) +
    guides(shape=guide_legend(title="Exposure"))
  PVE<-ggplot(sum,aes(x=variable,y=value,fill=variable)) + 
    geom_bar(stat="identity",aes(factor(variable)),fill=getPalette(colourCount)) +
    theme_bw() + labs(x="Principal component",y="Proportion of variance explained")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  library(cowplot)
  
  pdf("DepthPCA.pdf",height=8,width=8)
  print(plot_grid(PC,PVE,nrow=2,labels="AUTO"))
  dev.off()
}
#DepthPCA()

# This function formats and writes out methylation data for running MACAU on the command line
MACAU<-function(){
  genomacau<-fread("genetic.distance.ind.txt",data.table = FALSE,header=FALSE)
  genomacau$V1<-gsub(".","-",genomacau$V1,fixed = TRUE)
  colnames(genomacau)[1]<-"sampleID"
  genomacau<-genomacau[order(match(genomacau$sampleID,pheno$sampleID)),]
  fwrite(genomacau[,-1],"geno.macau.txt",quote=FALSE,row.names = FALSE,col.names = FALSE,sep="\t")
  cov<-data.frame(LF.meth)
  fwrite(cov,"cov.macau.txt",quote=FALSE,row.names = FALSE,col.names = FALSE,sep="\t")
  fwrite(data.frame(as.numeric(as.factor(pheno$treatment))),"pheno.macau.txt",quote=FALSE,row.names = FALSE,col.names = FALSE,sep="\t")
  
}

# This function contains the analyses performed for the genetic information: RDA, outlier tests, and principal components analysis
GeneticAnalysis<-function(){
  
  # PCA for visualization
  
  library(flashpcaR)
  
  Y <- allele[,-c(1:4)]
  pc <- prcomp(t(Y))
  
  pheno$sampleID<-gsub("SPALT_","",pheno$sampleID)
  
  
  library(vegan)
  PCdata<-data.frame(pheno,pc$x)
  library(ggplot2)
  library(ggthemes)
  library(ggrepel)
  
  PCA<-ggplot(PCdata,aes(x=PC1,y=PC2,color=treatment,label=sampleID, shape=site)) + geom_point() +
    #geom_text(size=2,position = position_nudge(y = -0.5)) +
    geom_text_repel(size=2.5) + labs(title="A") +
    scale_color_fivethirtyeight() + theme_bw() +
    guides(color=guide_legend(title="Exposure"),shape=guide_legend(title="Site"))
  
  
  # Generate neutral distribution for tests via pruning
  
  #allele.in<-allele[,-c(2:4)]
  
  GetPruned<-function(data=allele,avgR2level=0.1){
    chrtable<-data.frame(table(data$Chr))
    chrtable<-dplyr::filter(chrtable,Freq>=3)
    chrs<-unique(as.numeric(as.character(chrtable[,1])))
    GetCors<-function(chrnum){
      data2<-dplyr::filter(data,Chr == chrs[chrnum])
      GetAvgCor<-function(dataIn=data2){
        tmp<-data.matrix(t(dataIn[,-c(1:4)]))
        colnames(tmp)<-paste(dataIn$Chr,data2$Pos,sep=":")
        tmp<-data.frame(as.table(cor(tmp)^2))
        Cors<-tmp %>% group_by(Var1) %>% dplyr::summarize(Mean = mean(Freq, na.rm=TRUE))
        #Corout<-data.frame(dataIn[i,c(1:4)],CorMean=mean(Cors,na.rm=TRUE))
        #return(dplyr::filter(Cors,Mean==max(Mean)))
        return(dplyr::filter(Cors,Mean<avgR2level))
      }
      #print(paste(chrnum,"out of",length(chrs)))
      return(GetAvgCor())
      #
      #return(out)
    }
    library(parallel)
    out<-mclapply(1:length(chrs),GetCors,mc.cores=7)
    out<-rbindlist(out)
    colnames(out)<-c("ChrPos","Mean")
    allele$ChrPos<-paste(allele$Chr,allele$Pos,sep=":")
    LDprun<-dplyr::filter(allele, ChrPos %in% out$ChrPos)
    LDprun$ChrPos<-NULL
    return(LDprun)
  }
  #LDprun<-GetPruned()
  
  # RDA as both a test for differentiation and a genome scan
  
  Y <- data.frame(t(allele[,-c(1:4)]))
  fit<-rda(Y~treatment+Condition(LF.allele),data=pheno)
  sc<-scores(fit,choices=1:5,display="sites")
  scdata<-data.frame(pheno,RDA1=sc[,1],pc$x)
  
  # Permutation test - takes a few min
  anova(fit,parallel=3,permutations = 999,by="terms")
  
  # differences in dispersion
  tempdist<-dist((Y))
  #colnames(genomacau)<-colnames(allele[,-c(1:4)])
  beta.test<-betadisper(d=tempdist,group=as.factor(pheno$treatment),bias.adjust = FALSE)
  anova(beta.test)
  
  #pdf("RDA.genetic.pdf",height=6,width=6)
  RDA<-ggplot(scdata,aes(y=RDA1,x=PC1,color=treatment,label=sampleID, shape=site)) + geom_point() +
    #geom_text(size=2,position = position_nudge(y = -0.5)) +
    geom_text_repel(size=2.5) + labs(title="B") +
    scale_color_fivethirtyeight() + theme_bw() + guides(shape=FALSE,color=FALSE)
  #dev.off()
  
  # call significance
  sc<-scores(fit,choices=1:5,display="species")
  scdata<-data.frame(Locus=1:nrow(allele),allele[,c(1:4)],sc)
  std<-sd(scdata$RDA1)
  scdata$sig<-abs(scdata$RDA1)>(std*3)
  sum(scdata$sig==TRUE)
  
  
  #fitPrun<-rda(data.frame(t(LDprun[,-c(1:4)]))~treatment+Condition(LF.allele),data=pheno)
  #scPrun<-scores(fitPrun,choices=1:5,display="species")
  #perc<-quantile(abs(scPrun[,1]),probs = seq(0, 1, 0.01))
  #Upper<-perc["99%"]
  #scdata$sig<-abs(scdata$RDA1)>=Upper
  #sum(scdata$sig==TRUE)
  
  #SigLoci<-dplyr::filter(scdata, sig==TRUE)
  SigLoci<-scdata
  fwrite(dplyr::filter(scdata, sig==TRUE),"Genetic.sig.txt",quote=FALSE,row.names = FALSE,sep="\t")
  
  # get change through simple regression
  GetChange<-function(){
    #atemp<-allele[,-c(1:4)]
    #aNo<-dplyr::select(atemp,contains("N"))
    #aO<-dplyr::select(atemp,-contains("N"))
    #atemp$No<-rowMeans(data.matrix(aNo))
    #atemp$O<-rowMeans(data.matrix(aO))
    #atemp$change<-atemp$No-atemp$O
    #atemp<-data.frame(allele,Change=atemp$change)
    #atemp<-atemp[,c(1:2,ncol(atemp))]
    
    atemp<-data.frame(t(allele[,-c(1:4)]))
    atemp<-data.frame(sampleID=rownames(atemp),atemp)
    atemp$sampleID<-gsub("SPALT_","",atemp$sampleID)
    pheno$sampleID<-gsub("SPALT_","",pheno$sampleID)
    atemp<-merge(pheno,atemp,by="sampleID")
    atemp<-atemp[,-c(1,3,4)]
    X=model.matrix(~atemp[,1])
    Loop<-function(loc){
      #lmfit<-lm(atemp[,i]~treatment,data=atemp)
      lmfit<-.lm.fit(x=X,y=atemp[,loc])
      #lmfit$coefficients
      #lmfit2$coefficients
      #cor(atemp[,i],as.numeric(as.factor(atemp$treatment)))
      return(lmfit$coefficients[2])
    }
    library(parallel)
    Avg<-sapply(2:ncol(atemp),Loop)
    atemp<-data.frame(allele,Change=Avg)
    atemp<-atemp[,c(1:2,ncol(atemp))]
    return(atemp)
    
  }
  atemp<-GetChange()
  
  
  SigLoci$ID<-paste(SigLoci$Chr,SigLoci$Pos,sep=":")
  atemp$ID<-paste(atemp$Chr,atemp$Pos,sep=":")
  
  SigLoci<-merge(SigLoci,atemp,by="ID",all.x=TRUE)
  
  # Plot
  GenPlot<-ggplot(SigLoci,aes(x=(RDA1),y=Change,color=sig)) + geom_point(alpha=0.25) + 
    scale_color_hc(name="Significant?") +# coord_polar() +
    geom_hline(yintercept = 0.2) + geom_hline(yintercept = -0.2) +
    theme_bw() +
    #theme(axis.text.x=element_blank(),
    #      axis.ticks.x=element_blank()) +
    theme(legend.position = "bottom") +
    labs(y="Change due to oil exposure (no oil - oil)",
         x="Locus score", title="C")
  
  # just counting stuff
  sum(abs(SigLoci$Change)>0.05 & SigLoci$sig == TRUE)
  sum(abs(SigLoci$Change)>0.2 & SigLoci$sig == TRUE)
  sum(SigLoci$Change>0 & SigLoci$sig == TRUE)
  sum(SigLoci$Change<0 & SigLoci$sig == TRUE)
  sum(SigLoci$Change>0 & SigLoci$sig == TRUE & abs(SigLoci$Change)>0.05)
  sum(SigLoci$Change<0 & SigLoci$sig == TRUE & abs(SigLoci$Change)>0.05)
  sum(SigLoci$Change>0 & SigLoci$sig == TRUE & abs(SigLoci$Change)>0.2)
  sum(SigLoci$Change<0 & SigLoci$sig == TRUE & abs(SigLoci$Change)>0.2)
  
  library(cowplot)
  library(gridExtra)
  pdf("Genetic.plots.pdf",height=7,width=10)
  #plot_grid(PCA,RDA,rel_widths = c(1.2,1))
  grid.arrange(arrangeGrob(PCA, RDA, ncol = 2,widths=c(1.2,1)), # Second row with 2 plots in 2 different columns
               GenPlot,                             # First row with one plot spaning over 2 columns
               nrow = 2)                       # Number of rows
  dev.off()
  
  popdist<-fread("genetic.distance.pop.txt",data.table = FALSE)
  LatLong<-fread("coordinates.txt",data.table=FALSE)
  colnames(LatLong)[1]<-"site"
  LatLong<-LatLong[c(4,5,1,2,6,3),]
  
  PPheno<-pheno[,-1]
  PPheno<-PPheno[!duplicated(PPheno),]
  PPheno$site<-gsub("9","",PPheno$site)
  PPheno<-merge(LatLong,PPheno,by="site")
  
  fit0<-capscale(popdist~Lat*Lon,data=PPheno)
  anova(fit0,by="term")
  plot(fit0)
  
}

# This function contains the analyses performed for the methylation information: RDA (both conditional on genetic structure and unconditional) and principal components analysis
MethylAnalysis<-function(){
  
  library(vegan)
  library(ggplot2)
  library(ggthemes)
  library(ggrepel)
  
  Y <- meth[,-c(1:4)]
  pc <- prcomp(t(Y))
  PCdata<-data.frame(pheno,pc$x)
  
  #pdf("PCA.methylation.pdf",height=6,width=6)
  PCA<-ggplot(PCdata,aes(x=PC1,y=PC2,color=treatment,label=sampleID, shape=site)) + geom_point() +
    #geom_text(size=2,position = position_nudge(y = -0.5)) +
    geom_text_repel(size=2.5) +
    scale_color_fivethirtyeight() + theme_bw() + 
    guides(color=guide_legend(title="Exposure"),shape=guide_legend(title="Site"))
  #dev.off()
  
  genpc<-prcomp(allele[,-c(1:4)])
  genpc<-genpc$rotation
  
  Y <- data.frame(t(meth[,-c(1:4)]))
  fit.uncontrolled<-rda(Y~treatment+Condition(LF.meth),data=pheno)
  fit<-rda(Y~treatment+Condition(LF.meth+genpc[,1:5]),data=pheno)
  sc<-scores(fit,choices=1:5,display="sites")
  scdata<-data.frame(pheno,sc)
  
  #anova(fit.uncontrolled,permutations = 999,by="terms",parallel=2) # varExp = 16.746, F=1.9647, P<0.001
  #anova(fit,permutations = 999,by="terms",parallel=2) # varExp = 9.91, F=1.199, P=0.1
  
  #pdf("RDA.meth.pdf",height=6,width=6)
  RDA<-ggplot(scdata,aes(x=RDA1,y=PC1,color=treatment,label=sampleID, shape=site)) + geom_point() +
    #geom_text(size=2,position = position_nudge(y = -0.5)) +
    geom_text_repel(size=2.5) +
    scale_color_fivethirtyeight() + theme_bw()+ guides(shape=FALSE,color=FALSE)
  #dev.off()
  
  library(cowplot)
  pdf("Methyl.plots.pdf",height=5,width=10)
  plot_grid(PCA,RDA,rel_widths = c(1.2,1))
  dev.off()
  
  #sc<-scores(fit,choices=1:5,display="species")
  #scdata<-data.frame(Locus=1:nrow(meth),meth[,c(1:4)],sc)
  #std<-sd(scdata$RDA1)
  #scdata$sig<-abs(scdata$RDA1)>(std*3)
  
  #SigLoci<-dplyr::filter(scdata, sig==TRUE)
  fwrite(scdata,"Methylation.RDA.csv",quote=FALSE,row.names = FALSE)
  
  # look for context overrepresentation
  
  #M <- data.frame(table(meth$context))
  #M<-cbind(M,Obs=data.frame(table(SigLoci$))[,2])
  #M$Freq<-as.numeric(as.character(M$Freq))
  #M$Obs<-as.numeric(as.character(M$Obs))
  #M2<-data.frame(t(M))
  #colnames(M2)<-c("Unknown","CG","CHG","CHH")
  #M2<-M2[-1,]
  #for(i in 1:ncol(M2)){
  #  M2[,i]<-as.numeric(as.character(M2[,i]))
  #}
  #M <- as.table(as.numeric(M2[1,]),as.numeric(M2[2,]))
  #rownames(M2)<-c("Expected","Observed")
  #chisq.test((M))
}

# Intermediate step: run MACAU

# This function contains the assessment/qvalue correction for the methylation DMP analysis (performed separately using MACAU)
# This function also contains the chi-squared test used to look for overrepresentation of contexts among DMPs
Methyl_Locus<-function(){
  library(stringr)
  Mac<-fread("spalt_epigbs.assoc.txt",data.table = FALSE)
  tempdf<-data.frame(str_split_fixed(Mac$id,":",3))
  colnames(tempdf)<-c("chr","pos","context")
  Mac<-data.frame(tempdf,Mac)
  hist(Mac$pvalue,breaks=50)
  
  library(qvalue)
  qs<-qvalue(p=Mac$pvalue)
  Mac$qval<-qs$qvalues
  nrow(dplyr::filter(Mac, qval<=0.05))
  sigs<-dplyr::filter(Mac, qval<=0.05)
  sigs$chr<-as.numeric(as.character(sigs$chr))
  sigs$pos<-as.numeric(as.character(sigs$pos))
  rm(tempdf)
  fwrite(sigs,"Methylation.sig.csv",quote=FALSE,row.names = FALSE)
  
  # test for overrepresentation
  
  M <- data.frame(table(Mac$context))
  M<-cbind(M,Obs=data.frame(table(sigs$context))[,2])
  M$Freq<-as.numeric(as.character(M$Freq))
  M$Obs<-as.numeric(as.character(M$Obs))
  M2<-data.frame(t(M))
  colnames(M2)<-c("Unknown","CG","CHG","CHH")
  M2<-M2[-1,]
  for(i in 1:ncol(M2)){
    M2[,i]<-as.numeric(as.character(M2[,i]))
  }
  M <- as.table(as.numeric(M2[1,]),as.numeric(M2[2,]))
  #rownames(M2)<-c("Expected","Observed")
  chisq.test((M)) # this suggests that the distributions are significantly different. the plot indicates that CG is overrepresented
  
  #CG overrrepresentation
  binom.test(x=125,n=240,p=M2[1,2]/nrow(Mac))
  
  #CHG overrrepresentation
  binom.test(x=57,n=240,p=M2[1,3]/nrow(Mac))
  
  #CHH overrrepresentation
  binom.test(x=57,n=240,p=M2[1,4]/nrow(Mac))
  
  # percent changes
  
  
  GetChange<-function(allele=meth){
    #atemp<-allele[,-c(1:4)]
    #aNo<-dplyr::select(atemp,contains("N"))
    #aO<-dplyr::select(atemp,-contains("N"))
    #atemp$No<-rowMeans(data.matrix(aNo))
    #atemp$O<-rowMeans(data.matrix(aO))
    #atemp$change<-atemp$No-atemp$O
    #atemp<-data.frame(allele,Change=atemp$change)
    #atemp<-atemp[,c(1:2,ncol(atemp))]
    
    atemp<-data.frame(t(allele[,-c(1:4)]))
    atemp<-data.frame(sampleID=rownames(atemp),atemp)
    atemp$sampleID<-gsub("SPALT_","",atemp$sampleID)
    pheno$sampleID<-gsub("SPALT_","",pheno$sampleID)
    atemp<-merge(pheno,atemp,by="sampleID")
    atemp<-atemp[,-c(1,3,4)]
    X=model.matrix(~atemp[,1])
    Loop<-function(loc){
      #lmfit<-lm(atemp[,i]~treatment,data=atemp)
      lmfit<-.lm.fit(x=X,y=atemp[,loc])
      #lmfit$coefficients
      #lmfit2$coefficients
      #cor(atemp[,i],as.numeric(as.factor(atemp$treatment)))
      return(lmfit$coefficients[2])
    }
    library(parallel)
    Avg<-sapply(2:ncol(atemp),Loop)
    atemp<-data.frame(allele,Change=Avg)
    atemp<-atemp[,c(1:3,ncol(atemp))]
    return(atemp)
    
  }
  mtemp<-GetChange()
  
  
  Meth2<-mtemp
  Meth2$ChrPos<-paste("C",Meth2$chr,Meth2$pos,Meth2$context,sep=":")
  sigs$ChrPos<-paste("C",sigs$chr,sigs$pos,sigs$context,sep=":")
  Meth2<-Meth2[,-c(1:3)]
  
  Mac2<-Mac
  Mac2$ChrPos<-paste("C",Mac2$chr,Mac2$pos,Mac2$context,sep=":")
  sigs2<-merge(Mac2,Meth2,by="ChrPos")
  fwrite(sigs2,"Methylation.all.csv",quote=FALSE,row.names = FALSE)
  
  
  sigs<-merge(sigs,Meth2,by="ChrPos")
  #sigs$Change<-out$change
  sigs$ChrPos<-NULL
  
  fwrite(sigs,"Methylation.sig.csv",quote=FALSE,row.names = FALSE)
  
  # counting
  mean(sigs$Change)
  sum(sigs$Change>0.05)
  sum(sigs$Change>0.2)
  
  
  # visualization
  
  # first, get BLAST
  
  GetBlast<-function(){
    library(data.table)
    library(stringr)
    setwd("/Volumes/Analysis3/2016SpAlt_epiGBS/blastn_results_DG2/blastn_vs_Transcriptome_Salt/")
    
    # first, blasted the consensus reads against the published transcriptome
    # now looking for overlaps
    
    # Get key
    epiVStrans<-fread("/Volumes/Analysis3/2016SpAlt_epiGBS/blastn_results_DG2/blastn_vs_Transcriptome_Salt/res_blastn_ConsEpiGBS_vs_TransSalt_BH.txt",data.table = FALSE)
    colnames(epiVStrans)<-c("queryID","seqID","pctIdent","alignmentLength",
                            "mismatch","gapopen","queryStart","queryEnd","seqStart","seqEnd","evalue","bitscore")
    
    epiVStrans<-data.frame(data.frame(str_split_fixed(epiVStrans$queryID,";",2)),epiVStrans)
    epiVStrans$X1<-paste("C",as.character(epiVStrans$X1),sep=":")
    #epiVStrans2<-data.frame(table(epiVStrans$X1,epiVStrans$seqID))
    #epiVStrans2<-dplyr::filter(epiVStrans2,!Freq==0)
    
    epiVStrans2<-epiVStrans[,c(1,4,9:12)]
    setwd("/Volumes/Analysis3/2016SpAlt_epiGBS/RedoAnalysis/")
    return(epiVStrans2)
  }
  blast<-GetBlast()
  colnames(blast)[1]<-"chr"
  
  library(ggplot2)
  library(ggthemes)
  library(ggsci)
  library(ggrepel)
  
  sigs$chr<-paste("C",sigs$chr,sep=":")
  sigs$big<-abs(sigs$Change)>=0.2
  sigs<-merge(sigs,blast,by="chr",all.x=TRUE)
  
  Plot<- ggplot(sigs,aes(x=pvalue,y=abs(Change),color=big,label=seqID)) 
  Plot<-Plot + geom_point() + 
    geom_text_repel(data=dplyr::filter(sigs,big==TRUE)) + 
    labs(x="P-value",y="Absolute change due to oil treatment") +
    geom_hline(yintercept = 0.2) + scale_color_npg(guide=FALSE) + theme_bw() +
    theme(legend.position = "bottom")
  
  pdf("Meth.sigChange.pdf",height=6,width=6)
  print(Plot)
  dev.off()
}

# This function examines the relationship between the outlier SNPs/DMPs and the previously-generated transcriptome
Overlaps<-function(){
  library(data.table)
  library(stringr)
  setwd("/Volumes/Analysis3/2016SpAlt_epiGBS/blastn_results_DG2/blastn_vs_Transcriptome_Salt/")
  
  
  # first, blasted the consensus reads against the published transcriptome
  # now looking for overlaps
  
  # Get key
  epiVStrans<-fread("/Volumes/Analysis3/2016SpAlt_epiGBS/blastn_results_DG2/blastn_vs_Transcriptome_Salt/res_blastn_ConsEpiGBS_vs_TransSalt_BH.txt",data.table = FALSE)
  colnames(epiVStrans)<-c("queryID","seqID","pctIdent","alignmentLength",
                          "mismatch","gapopen","queryStart","queryEnd","seqStart","seqEnd","evalue","bitscore")
  
  epiVStrans<-data.frame(data.frame(str_split_fixed(epiVStrans$queryID,";",2)),epiVStrans)
  epiVStrans$X1<-paste("C",as.character(epiVStrans$X1),sep=":")
  #epiVStrans2<-data.frame(table(epiVStrans$X1,epiVStrans$seqID))
  #epiVStrans2<-dplyr::filter(epiVStrans2,!Freq==0)
  
  epiVStrans2<-epiVStrans[,c(1,4,9:12)]
  setwd("/Volumes/Analysis3/2016SpAlt_epiGBS/RedoAnalysis/")
  fwrite(epiVStrans2,"TranscriptomeOverlap.csv",quote=FALSE,row.names = FALSE)
  
  # 5' overlap
  sum(epiVStrans2$seqStart<250)
  
  
  # VS microarray
  
  GetBLAST<-function(x="epiGBSblast.txt"){
    Data<-fread(x, header=FALSE, sep="\t", data.table = TRUE)
    colnames(Data)<-c("queryId", "subjectId", "percIdentity", "alnLength", "mismatchCount",
                      "gapOpenCount", "queryStart", "queryEnd", "subjectStart", 
                      "subjectEnd", "eVal", "bitScore")
    
    print("File loaded. Filtering out duplicate entries, keeping value with highest bitScore...")
    setkey(Data,"queryId")
    Data<-Data[, .SD[which.max(abs(bitScore))], by=queryId]
    Data<-filter(Data, percIdentity > 90)
    return(Data)
  }
  mic<-GetBLAST(x="/Volumes/Analysis3/2016SpAlt_epiGBS/Spart_array_probes_vs_Salt_2.txt")
  #mic<-fread("/Volumes/Analysis3/2016SpAlt_epiGBS/Spart_array_probes_vs_Salt_2.txt",data.table = FALSE)
  nrow(dplyr::filter(epiVStrans2, seqID %in% mic[,2]))
  length(unique(epiVStrans2$seqID) %in% unique(mic[,2]))
  
  micAll<-mic
  
  # This creates a table of only significantly expressed loci 
  micsig<-fread("pctlnorm_amr_sig.txt",data.table = FALSE)
  colnames(micsig)[28]<-"OilSig"
  
  # Filter to only the loci that were significantly differentially expressed
  micsig<-dplyr::filter(micsig, OilSig==1)
  mic<-dplyr::filter(mic,queryId %in% micsig$contig)
  rm(micsig)
  
  MicTest<-mic
  colnames(MicTest)[2]<-"seqID"
  
  # This creates a harmonized list of labels that appear in the epiGBS, the microarray, and the transcriptome
  MicTest<-merge(MicTest[,1:2],epiVStrans2,by="seqID",all=FALSE)
  
  
  # Now look for overlaps with methylation and genetics
  BLAST<-fread("TranscriptomeOverlap.csv",data.table = FALSE)
  
  MethSig<-fread("Methylation.sig.csv",data.table = FALSE)
  length(unique(MethSig$chr))
  
  GenSig<-fread("Genetic.sig.txt",data.table = FALSE)
  
  MethSig$chr<-paste("C",as.character(MethSig$chr),sep=":")
  GenSig$Chr<-paste("C",as.character(GenSig$Chr),sep=":")
  
  # just counting here
  library(stringr)
  
  allele2<-allele
  allele2$Chr<-paste("C",allele2$Chr,sep=":")
  
  # Overlapping transcriptome
  nrow(dplyr::filter(allele2, Chr %in% unique(BLAST$X1)))
  #length(unique(allele2$Chr) %in% unique(BLAST$X1)
  #length(intersect(GenSig$Chr,BLAST$X1))
  
  # Number of significant genetic loci in coding regions
  nrow(dplyr::filter(GenSig, Chr %in% unique(BLAST$X1)))
  # Number of significant methylation loci in coding regions
  nrow(dplyr::filter(MethSig, chr %in% unique(BLAST$X1)))
  # Number of significant genetic loci among microarray hits
  SigNum<-dplyr::filter(GenSig, sig==TRUE)
  length(unique(dplyr::filter(SigNum, Chr %in% unique(MicTest$X1))$Chr))
  fwrite(dplyr::filter(SigNum, Chr %in% unique(MicTest$X1)),"GenMicOverlap.csv",quote=FALSE,row.names = FALSE)
  Micsigss<-(nrow(dplyr::filter(SigNum, Chr %in% unique(MicTest$X1))))
  
  
  # bootstrap test to see if genetic loci are more significant than you'd expect by chance
  # make sure it's just the values and a column of loci
  allele.in<-allele[,-c(2:4)]
  
  BootTest<-function(Reps=9999,input=allele.in){
    
    #incor<-cor(data.matrix(t(input[,-c(1)])))
    
    Test<-function(i,data=input){
      # first, randomly draw a number of row equal to the number of significant loci
      sub<-sample_n(data,nrow(SigNum),replace = TRUE)
      sub[,1]<-paste("C",sub[,1],sep=":")
      colnames(sub)[1]<-"Chr"
      # then, see how many overlap the significant microarray data
      (OL<-nrow(dplyr::filter(sub, Chr %in% unique(MicTest$X1))))
      return(OL)
    }
    
    
    library(parallel)
    
    ## start cluster
    parallelCluster <- parallel::makeCluster(detectCores()-1, type = "FORK")
    print(parallelCluster)
    
    result<-parSapply(parallelCluster,1:Reps,Test)
    
    # Shutdown cluster neatly
    if(!is.null(parallelCluster)) {
      parallel::stopCluster(parallelCluster)
      parallelCluster <- c()
    }
    
    return((sum(result>=Micsigss)+1)/Reps)
  }
  BootTest(input=allele) # P>0.5
  
  # Number of significant methylation loci among microarray hits
  SigNum<-MethSig
  fwrite(dplyr::filter(SigNum, chr %in% unique(MicTest$X1)),"MethMicOverlap.csv",quote=FALSE,row.names = FALSE)
  Micsigss<-(nrow(dplyr::filter(MethSig, chr %in% unique(MicTest$X1))))
  #colnames(meth)[1]<-"Chr"
  BootTest(input=meth) # P<0.05
  
  # Before doing anything else, format significant loci for publication
  FormatsSig<-function(){
    Gens<-fread("Genetic.sig.txt",data.table = FALSE)
    Gens<-Gens[,2:3]
    Meths<-fread("Methylation.sig.csv",data.table=FALSE)
    Meths<-Meths[,c(1:3)]
    Gens$Type<-"SNP"
    colnames(Meths)<-colnames(Gens)
    
    Comb<-rbind(Gens,Meths)
    colnames(epiVStrans2)[1]<-"Chr"
    Comb$Chr<-paste("C",Comb$Chr,sep=":")
    
    Comb<-merge(Comb,epiVStrans2[,1:2],by="Chr",all.x = TRUE)
    
    colnames(MicTest)[3]<-"Chr"
    
    Comb<-merge(Comb,MicTest[,c(2:3)],by="Chr",all.x = TRUE)
    
    
    OS<-fread("OrSatBLAST.csv",data.table = FALSE,na.strings = "")
    colnames(OS)[1]<-"seqID"
    Comb<-merge(Comb,OS,by="seqID",all.x = TRUE)
    Comb<-Comb[,c(2:5,1,8)]
    colnames(Comb)<-c("Contig","Position","Type","Microarray ID","Transcriptome ID","Oryza sativa ID")
    fwrite(Comb,"TableS1.csv",quote=FALSE,row.names = FALSE) 
  }
  FormatsSig()
  
  # Oryza overlap
  
  Oryza<-function(){
    # now looking for overlaps with oryza annotations
    OS<-fread("OrSatBLAST.csv",data.table = FALSE,na.strings = "")
    
    colnames(OS)[1]<-"seqID"
    #colnames(epiVStrans2)[1]<-"contig"
    annot<-merge(epiVStrans2,OS,by="seqID",all.x = TRUE)
    
    # writing out results: gen
    annot.gensig<-dplyr::filter(annot,X1 %in% GenSig$Chr)
    nrow(dplyr::filter(annot.gensig, !is.na(V4)))
    fwrite(annot.gensig,"Annot.gensig.csv",quote=FALSE,row.names = FALSE)
    
    nrow(dplyr::filter(annot.gensig, seqID %in% mic$V3))
    
    # meth
    annot.methsig<-dplyr::filter(annot,X1 %in% MethSig$chr)
    nrow(dplyr::filter(annot.methsig, !is.na(V4)))
    fwrite(annot.gensig,"Annot.methsig.csv",quote=FALSE,row.names = FALSE)
    
    colnames(annot.methsig)[2]<-"chr"
    annot.methsig<-merge(annot.methsig,MethSig,by="chr",all=TRUE)
    annot.methsig$abschange<-abs(annot.methsig$Change)
    
    nrow(dplyr::filter(annot.methsig, seqID %in% mic$V3))
  }
  
  # Distance comparison
  
  Comparison<-function(){
    Evalues<-fread("/Volumes/Analysis3/2016SpAlt_epiGBS/medsummclean.txt",data.table = FALSE)
    #Key<-fread("MicroarrayKey.csv",data.table = FALSE)
    #Key<-dplyr::filter(Key, TargetID %in% micAll$queryId)
    
    
    # Load meth data, separate contig
    library(stringr)
    Meth<-fread("MethFreq.txt",data.table = FALSE)
    #Meth.contig<-data.frame(str_split_fixed(Meth$ChrPos,":",2))
    #Meth<-data.frame(Meth.contig$X1,Meth)
    Meth[,2]<-Meth$chr
    colnames(Meth)[2]<-"epiGBS"
    Meth<-Meth[,-c(3:6)]
    colnames(Meth)[1]<-"ChrPos"
    
    # recreate pools using averaging
    GetPool<-function(samples,name){
      Pool<-data.frame(Meth)
      Pool<-Pool[,-c(1:ncol(Pool))]
      for(i in samples){
        temp<-dplyr::select(Meth,contains(i))
        Pool<-cbind(Pool,temp)
        rm(temp)
      }
      Pool$mean<-rowMeans(Pool)
      Pool<-data.frame(Pool$mean)
      colnames(Pool)<-name
      return(Pool)
    }
    
    ReturnMethylation<-function(){
      colnames(Meth)
      
      samples<-c("GIO1.6","GIO1.7")
      name<-"P1"
      
      # oil
      P1<-GetPool(samples=c("GIO1.6","GIO1.7","GIO1.5"),name="P1")
      P2<-GetPool(samples=c("GIO1.9","GIO1.4","GIO1.10"),name="P2")
      P3<-GetPool(samples=c("GIO1.1","GIO1.2","GIO1.3"),name="P3")
      P4<-GetPool(samples=c("GIO2.8","GIO2.6","GIO2.5"),name="P4")
      P5<-GetPool(samples=c("GIO2.7","GIO2.1","GIO2.4"),name="P5")
      P6<-GetPool(samples=c("GIO2.9","GIO2.3","GIO2.2"),name="P6")
      P7<-GetPool(samples=c("MSO6","MSO10","MSO4"),name="P7")
      P8<-GetPool(samples=c("MSO9","MSO1","MSO8"),name="P8")
      P9<-GetPool(samples=c("MSO2","MSO3","MSO5"),name="P9")
      
      oil<-cbind(Meth$ChrPos,P1,P2,P3,P4,P5,P6,P7,P8,P9)
      colnames(oil)[1]<-"ChrPos"
      rm(P1,P2,P3,P4,P5,P6,P7,P8,P9)
      
      # no oil
      P1<-GetPool(samples=c("GIN1.8","GIN1.4","GIN1.3"),name="P10")
      P2<-GetPool(samples=c("GIN1.9","GIN1.1","GIN1.2"),name="P11")
      P3<-GetPool(samples=c("GIN1.6","GIN1.7","GIN1.5"),name="P12")
      
      P4<-GetPool(samples=c("GIN2.9","GIN2.7","GIN2.11"),name="P13")
      P5<-GetPool(samples=c("GIN2.3","GIN2.4","GIN2.2"),name="P14")
      P6<-GetPool(samples=c("GIN2.1","GIN2.6","GIN2.12"),name="P15")
      
      P7<-GetPool(samples=c("MSN5","MSN8","MSN7"),name="P16")
      P8<-GetPool(samples=c("MSN4","MSN2","MSN3"),name="P17")
      P9<-GetPool(samples=c("MSN9","MSN1","MSN6"),name="P18")
      
      NOoil<-cbind(Meth$ChrPos,P1,P2,P3,P4,P5,P6,P7,P8,P9)
      colnames(NOoil)[1]<-"ChrPos"
      rm(P1,P2,P3,P4,P5,P6,P7,P8,P9)
      
      oil$ChrPos<-as.character(oil$ChrPos)
      NOoil$ChrPos<-as.character(NOoil$ChrPos)
      
      Meth.pools<-merge(oil,NOoil,by="ChrPos")
      Meth.pools$ChrPos<-as.character(Meth.pools$ChrPos)
      rm(oil,NOoil)
      return(Meth.pools)
    }
    Meth.pools<-ReturnMethylation()
    
    ## Genetic data
    ReturnGenetic<-function(){
      #Gen<-read.table("SPALT_combined.output.raw.snps.AD.table")
      #Gen2<-fread("allele_freq.txt",data.table = FALSE)
      Gen3<-fread("AlleleFreq.txt",data.table = FALSE,fill=TRUE)
      #Gen3<-data.frame(ChrPos=paste(Gen2$Chr,Gen2$Pos,sep=":"),Gen2)
      
      GetPool.gen<-function(samples,name){
        Pool<-data.frame(Gen3)
        Pool<-Pool[,-c(1:ncol(Pool))]
        for(i in samples){
          temp<-dplyr::select(Gen3,contains(i))
          Pool<-cbind(Pool,temp)
          rm(temp)
        }
        Pool$mean<-rowMeans(Pool)
        Pool<-data.frame(Pool$mean)
        colnames(Pool)<-name
        return(Pool)
      }
      
      # oil
      P1<-GetPool.gen(samples=c("GIO1.6","GIO1.7","GIO1.5"),name="P1")
      P2<-GetPool.gen(samples=c("GIO1.9","GIO1.4","GIO1.10"),name="P2")
      P3<-GetPool.gen(samples=c("GIO1.1","GIO1.2","GIO1.3"),name="P3")
      P4<-GetPool.gen(samples=c("GIO2.8","GIO2.6","GIO2.5"),name="P4")
      P5<-GetPool.gen(samples=c("GIO2.7","GIO2.1","GIO2.4"),name="P5")
      P6<-GetPool.gen(samples=c("GIO2.9","GIO2.3","GIO2.2"),name="P6")
      P7<-GetPool.gen(samples=c("MSO6","MSO10","MSO4"),name="P7")
      P8<-GetPool.gen(samples=c("MSO9","MSO1","MSO8"),name="P8")
      P9<-GetPool.gen(samples=c("MSO2","MSO3","MSO5"),name="P9")
      
      oil<-cbind(Gen3[,1:2],P1,P2,P3,P4,P5,P6,P7,P8,P9)
      colnames(oil)[1:2]<-c("Chr","Pos")
      rm(P1,P2,P3,P4,P5,P6,P7,P8,P9)
      
      # no oil
      P1<-GetPool.gen(samples=c("GIN1.8","GIN1.4","GIN1.3"),name="P10")
      P2<-GetPool.gen(samples=c("GIN1.9","GIN1.1","GIN1.2"),name="P11")
      P3<-GetPool.gen(samples=c("GIN1.6","GIN1.7","GIN1.5"),name="P12")
      
      P4<-GetPool.gen(samples=c("GIN2.9","GIN2.7","GIN2.11"),name="P13")
      P5<-GetPool.gen(samples=c("GIN2.3","GIN2.4","GIN2.2"),name="P14")
      P6<-GetPool.gen(samples=c("GIN2.1","GIN2.6","GIN2.12"),name="P15")
      
      P7<-GetPool.gen(samples=c("MSN5","MSN8","MSN7"),name="P16")
      P8<-GetPool.gen(samples=c("MSN4","MSN2","MSN3"),name="P17")
      P9<-GetPool.gen(samples=c("MSN9","MSN1","MSN6"),name="P18")
      
      NOoil<-cbind(Gen3[,1:2],P1,P2,P3,P4,P5,P6,P7,P8,P9)
      colnames(NOoil)[1:2]<-c("Chr","Pos")
      rm(P1,P2,P3,P4,P5,P6,P7,P8,P9)
      
      ChrPos<-paste(oil$Chr,oil$Pos,sep = ":")
      oil<-data.frame(ChrPos,oil[,-c(1:2)])
      NOoil<-data.frame(ChrPos,NOoil[,-c(1:2)])
      
      oil$ChrPos<-as.character(oil$ChrPos)
      NOoil$ChrPos<-as.character(NOoil$ChrPos)
      
      Gen.pools<-merge(oil,NOoil,by="ChrPos")
      Gen.pools$ChrPos<-as.character(Gen.pools$ChrPos)
      rm(oil,NOoil)
      return(Gen.pools)
    }
    Gen.pools<-ReturnGenetic()
    
    # Unify ChrPos formats
    
    library(tidyr)
    Meth.ChrPos<-data.frame(str_split_fixed(Meth.pools$ChrPos,":",2))
    Meth.pools<-data.frame(Meth.ChrPos,Meth.pools)
    colnames(Meth.pools)[1:2]<-c("Chr","Pos")
    for(i in 1:2){
      Meth.pools[,i]<-as.character(Meth.pools[,i])
    }
    
    Gen.ChrPos<-data.frame(str_split_fixed(Gen.pools$ChrPos,":",2))
    Gen.pools<-data.frame(Gen.ChrPos,Gen.pools)
    colnames(Gen.pools)[1:2]<-c("Chr","Pos")
    for(i in 1:2){
      Gen.pools[,i]<-as.character(Gen.pools[,i])
    }
    
    temp<-data.frame(str_split_fixed(colnames(Evalues),"_",5))[-1,]
    colnames(Evalues)[-1]<-paste("P",as.character(temp[,3]),sep="")
    E.dist<-dist(t(Evalues[,-1]))
    M.dist<-dist(t(Meth.pools[,-c(1:3)]))
    G.dist<-dist(t(Gen.pools[,-c(1:3)]))
    
    library(vegan)
    mantel.partial(E.dist,M.dist,G.dist)
    mantel(E.dist,G.dist)
    mantel(E.dist,M.dist)
  }
  
  
  
  
  # methylation and genetics overlap - no changes to trinucleotide
  
  overlap<-dplyr::filter(MethSig, chr %in% GenSig$Chr)
  overlap.gen<-dplyr::filter(GenSig, Chr %in% MethSig$chr)
  
  out<-data.frame()
  
  contigs<-unique(overlap$chr)
  
  for( i in 1:length(contigs)){
    tempdf1<-dplyr::filter(overlap,chr == contigs[i])
    tempdf2<-dplyr::filter(overlap.gen,Chr == contigs[i])
    
    contextchk<-nrow(dplyr::filter(tempdf2, Pos %in% tempdf1[,2] || Pos %in% (tempdf1[,2]+1) || Pos %in% (tempdf1[,2] +2)))
    
    toout<-data.frame(Chr=contigs[i],Num=contextchk)
    out<-rbind(out,toout)
  }
  
}

# This function contains the methylation figures
Figures<-function(){
  library(ggthemes)
  library(ggsci)
  ### Methylation
  
  MethSig<-fread("methylation.all.csv",data.table = FALSE)
  MethRDA<-fread("methylation.rda.csv",data.table = FALSE)
  
  library(flashpcaR)
  PCmeth<-flashpca(t(meth[,-c(1:4)]),stand="none")
  PCmeth$vectors
  PCdata<-data.frame(pheno,PCmeth$vectors)
  
  FixPoly<-function(){
    temp<-dplyr::filter(MethSig,!context %in% c("CG","CHG","CHH"))
    temp[,4]<-"Polymorphic"
    temp2<-dplyr::filter(MethSig,context %in% c("CG","CHG","CHH"))
    temp<-rbind(temp,temp2)
    return(temp)
  }
  MethSig<-FixPoly()
  
  scdata<-data.frame(pheno,RDA1=MethRDA$RDA1,PCmeth$vectors)
  
  library(ggrepel)
  PCA<-ggplot(PCdata,aes(x=X1,y=X2,color=treatment,label=sampleID, shape=site)) + geom_point() +
    #geom_text(size=2,position = position_nudge(y = -0.5)) +
    geom_text_repel(size=2.5) + labs(x="PC1",y="PC2") +
    scale_color_fivethirtyeight() + theme_bw() + 
    guides(color=guide_legend(title="Exposure"),shape=guide_legend(title="Site"))
  
  RDA<-ggplot(scdata,aes(x=X1,y=RDA1,color=treatment,label=sampleID, shape=site)) + geom_point() +
    #geom_text(size=2,position = position_nudge(y = -0.5)) +
    geom_text_repel(size=2.5) + labs(x="PC1",y="RDA1") +
    scale_color_fivethirtyeight() + theme_bw()+ guides(shape=FALSE,color=FALSE)
  
  library(ggplot2)
  
  contexts_frags<-data.frame(table(MethSig$chr,MethSig$context))
  contexts_frags$Var2<-as.character(contexts_frags$Var2)
  #contexts_frags[1:6279,2]<-"Polymorphic"
  #contexts_frags$Var2<-gsub(".","Polymorphic",contexts_frags,fixed = TRUE)
  #contexts_frags<-dplyr::filter(contexts_frags, Var2 %in% c("CG","CHG","CHH"))
  
  # context distribution
  Methy_cxt<-ggplot(contexts_frags,aes(y=as.numeric(as.factor(Var1)),x=1,fill=Freq)) +
    geom_raster() + scale_fill_distiller(type="seq",palette = 7,direction=1,name="Frequency") +
    facet_grid(~as.factor(Var2)) +
    theme_bw() + 
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank()
    ) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "bottom") +
    labs(y="Consensus fragment number",x="Trinucleotide context",title="A")
  
  # Methylation outliers
  MethSig$sig<-MethSig$qval<=0.05
  MethyOut<-ggplot(MethSig,aes(y=-log10(pvalue),x=Change,color=sig)) + 
    geom_point(alpha=0.5) + theme_bw() +
    facet_wrap(~context) + 
    geom_vline(xintercept = 0.2,linetype="dashed",alpha=0.5) + 
    geom_vline(xintercept = -0.2,linetype="dashed",alpha=0.5) +
    #geom_vline(xintercept = 0) +
    labs(x="Change in methylation frequency",y="Negative log10 P-value")+
    theme(legend.position = "bottom") +
    scale_color_hc(name="Significant?")
  
  library(gridExtra)
  library(cowplot)
  
  pdf("Methylation.fig.pdf",height=9,width=9)
  #plot_grid(Methy_cxt,MethyOut,align="v",nrow=2)
  plot1<-plot_grid(PCA,RDA,ncol=2,labels = "AUTO", rel_widths = c(1.2,1))
  plot2<-plot_grid(plot1,MethyOut,nrow=2,labels=c("","C"))
  print(plot2)
  dev.off()
  
}