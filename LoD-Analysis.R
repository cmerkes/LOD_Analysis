## Run LoD-Calculator script first to get all the data in.
## Load additional R package required for making the figures:
library(cowplot)

## Adjust assay names to simplify figures:
TarMatch <- data.frame(Target=as.character(Targets),
                       Number=c(36,35,1,6,5,32,7,10,9,11,31,30,29,25,26,27,28,33,34,
                                13,12,2,3,4,15,14,17,16,18,19,20,21,22,23,24,8))
TarMatch <- TarMatch[c(24,23,1:22,25:36),]
for(i in 1:nrow(TarMatch)) {
  DAT$Target <- gsub(TarMatch[i,1],TarMatch[i,2],DAT$Target)
  DAT2$Target <- gsub(TarMatch[i,1],TarMatch[i,2],DAT2$Target)
  DAT3$Assay <- gsub(TarMatch[i,1],TarMatch[i,2],DAT3$Assay)
}
Targets <- c(1:36)
LOD.list2 <- LOD.list2[c(1,4,23,24,25,6,5,8,37,10,9,11,22,21,27,26,29,28,30,31,32,33,
                         34,35,36,15,16,17,18,14,13,12,7,19,20,3,2)]
LOD.list3 <- LOD.list3[c(1,4,23,24,25,6,5,8,37,10,9,11,22,21,27,26,29,28,30,31,32,33,
                         34,35,36,15,16,17,18,14,13,12,7,19,20,3,2)]
LOQ.list <- LOQ.list[c(1,4,23,24,25,6,5,8,37,10,9,11,22,21,27,26,29,28,30,31,32,33,
                       34,35,36,15,16,17,18,14,13,12,7,19,20,3,2)]
curve.list <- curve.list[c(1,4,23,24,25,6,5,8,37,10,9,11,22,21,27,26,29,28,30,31,32,33,
                       34,35,36,15,16,17,18,14,13,12,7,19,20,3,2)]

## Generate LOD figures:
png(filename="LOD-figure.png",width=6.5,height=8,units="in",res=300)
par(mfrow=c(4,2),mar=c(3,3,2.5,0.5),mgp=c(1.5,0.75,0))
for(i in 1:8) {
  if(!is.na(LOD.list2[i+1])) {
    DAT4 <- rbind(ED(get(LOD.list2[i+1]),0.95,interval="delta",type="absolute"),
                  ED(get(LOD.list2[i+1]),1-sqrt(0.05),interval="delta",type="absolute"),
                  ED(get(LOD.list2[i+1]),1-0.05^(1/3),interval="delta",type="absolute"),
                  ED(get(LOD.list2[i+1]),1-0.05^0.25,interval="delta",type="absolute"),
                  ED(get(LOD.list2[i+1]),1-0.05^0.2,interval="delta",type="absolute"),
                  ED(get(LOD.list2[i+1]),1-0.05^0.125,interval="delta",type="absolute"))
    if(substr(LOD.list3[i+1],1,3)=="LL2") {
      DAT4 <- exp(DAT4)
    }
    DAT4 <- data.frame(DAT4,LoD=c("1rep.LOD","2rep.LOD","3rep.LOD","4rep.LOD",
                                  "5rep.LOD","8rep.LOD"),
                       Assay=rep(Targets[i],nrow(DAT4)))
    DAT4$Assay <- as.character(DAT4$Assay)
    if(i==1) {
      LOD.CI <- DAT4
    }
    if(i>1) {
      LOD.CI <- rbind(LOD.CI,DAT4)
    }
    if(sum(!is.na(DAT4[,3])&DAT4[,3]<=0)>0) { #Unable to plot negative lower limits, converting any lower limit values to 0.0001
      DAT4[!is.na(DAT4[,3])&DAT4[,3]<=0,3] <- 0.0001
    }
    plot(get(LOD.list2[i+1]),main=paste0("LoD Plot for assay ",Targets[i]),
         ylab="Detection Probability",xlab="Standard concentrations (Copies / Reaction)",
         xlim=c(min(DAT4[,1:4],na.rm=T),max(DAT$SQ,na.rm=T)),cex=0.6,cex.axis=0.6,
         cex.lab=0.6)
    LODS <- sum(!is.na(DAT4[,1]))
    COLS <- c(rgb(0.8,0.47,0.65),rgb(0,0.45,0.7),rgb(0.94,0.89,0.26),
              rgb(0.84,0.37,0),rgb(0,0.62,0.45),rgb(0.90,0.62,0))
    PNTS <- c(15,16,17,18,25,3)
    YS <- c(0.95,1-sqrt(0.05),1-0.05^(1/3),1-0.05^0.25,1-0.05^0.2,1-0.05^0.125)
    LODS2 <- c("Limit of Detection","2 Replicates LoD","3 Replicates LoD",
               "4 Replicates LoD","5 Replicates LoD","8 Replicates LoD")
    if(LODS<6) {
      LODS2[(LODS+1):6] <- gsub("licates LoD","s: Insufficient Data",LODS2[(LODS+1):6])
    }
    points(x=DAT4[1:LODS,1],y=YS[1:LODS],pch=PNTS,col=COLS,cex=0.9)
    for(j in 1:LODS) {
      lines(x=DAT4[j,3:4],y=rep(YS[j],2),col=COLS[j],lwd=1)
      lines(x=rep(DAT4[j,3],2),y=c(YS[j]-0.02,YS[j]+0.02),lwd=1,col=COLS[j])
      lines(x=rep(DAT4[j,4],2),y=c(YS[j]-0.02,YS[j]+0.02),lwd=1,col=COLS[j])
    }
    legend("bottomright",legend=LODS2,pch=PNTS,col=COLS,text.col=COLS,cex=0.7)
    Pval <- modelFit(get(LOD.list2[i+1]))[[5]][2]
    mtext(paste0("FCT used: ",LOD.list3[i+1],"    Lack of fit test: p = ",Pval),side=3,
          cex=0.5)
  }
  if(is.na(LOD.list2[i+1])) {
    plot(DAT2$Rate[DAT2$Target==Targets[i]]~log10(DAT2$Standards[DAT2$Target==Targets[i]]),
         ylim=c(0,1),ylab="Detection Probability",
         xlab=expression("Log of standard concentrations (Log"[10]*"Copies / Reaction)"),
         main=paste0("LoD for assay ",Targets[i]," unsolvable"),cex=0.6,cex.axis=0.6,
         cex.lab=0.6)
  }
}
dev.off()

png(filename="LOD-figure2.png",width=6.5,height=8,units="in",res=300)
par(mfrow=c(4,2),mar=c(3,3,2.5,0.5),mgp=c(1.5,0.75,0))
for(i in 9:16) {
  if(!is.na(LOD.list2[i+1])) {
    DAT4 <- rbind(ED(get(LOD.list2[i+1]),0.95,interval="delta",type="absolute"),
                  ED(get(LOD.list2[i+1]),1-sqrt(0.05),interval="delta",type="absolute"),
                  ED(get(LOD.list2[i+1]),1-0.05^(1/3),interval="delta",type="absolute"),
                  ED(get(LOD.list2[i+1]),1-0.05^0.25,interval="delta",type="absolute"),
                  ED(get(LOD.list2[i+1]),1-0.05^0.2,interval="delta",type="absolute"),
                  ED(get(LOD.list2[i+1]),1-0.05^0.125,interval="delta",type="absolute"))
    if(substr(LOD.list3[i+1],1,3)=="LL2") {
      DAT4 <- exp(DAT4)
    }
    DAT4 <- data.frame(DAT4,LoD=c("1rep.LOD","2rep.LOD","3rep.LOD","4rep.LOD",
                                  "5rep.LOD","8rep.LOD"),
                       Assay=rep(Targets[i],nrow(DAT4)))
    DAT4$Assay <- as.character(DAT4$Assay)
    if(i==1) {
      LOD.CI <- DAT4
    }
    if(i>1) {
      LOD.CI <- rbind(LOD.CI,DAT4)
    }
    if(sum(!is.na(DAT4[,3])&DAT4[,3]<=0)>0) { #Unable to plot negative lower limits, converting any lower limit values to 0.0001
      DAT4[!is.na(DAT4[,3])&DAT4[,3]<=0,3] <- 0.0001
    }
    plot(get(LOD.list2[i+1]),main=paste0("LoD Plot for assay ",Targets[i]),
         ylab="Detection Probability",xlab="Standard concentrations (Copies / Reaction)",
         xlim=c(min(DAT4[,1:4],na.rm=T),max(DAT$SQ,na.rm=T)),cex=0.6,cex.axis=0.6,
         cex.lab=0.6)
    LODS <- sum(!is.na(DAT4[,1]))
    COLS <- c(rgb(0.8,0.47,0.65),rgb(0,0.45,0.7),rgb(0.94,0.89,0.26),
              rgb(0.84,0.37,0),rgb(0,0.62,0.45),rgb(0.90,0.62,0))
    PNTS <- c(15,16,17,18,25,3)
    YS <- c(0.95,1-sqrt(0.05),1-0.05^(1/3),1-0.05^0.25,1-0.05^0.2,1-0.05^0.125)
    LODS2 <- c("Limit of Detection","2 Replicates LoD","3 Replicates LoD",
               "4 Replicates LoD","5 Replicates LoD","8 Replicates LoD")
    if(LODS<6) {
      LODS2[(LODS+1):6] <- gsub("licates LoD","s: Insufficient Data",LODS2[(LODS+1):6])
    }
    points(x=DAT4[1:LODS,1],y=YS[1:LODS],pch=PNTS,col=COLS,cex=0.9)
    for(j in 1:LODS) {
      lines(x=DAT4[j,3:4],y=rep(YS[j],2),col=COLS[j],lwd=1)
      lines(x=rep(DAT4[j,3],2),y=c(YS[j]-0.02,YS[j]+0.02),lwd=1,col=COLS[j])
      lines(x=rep(DAT4[j,4],2),y=c(YS[j]-0.02,YS[j]+0.02),lwd=1,col=COLS[j])
    }
    legend("bottomright",legend=LODS2,pch=PNTS,col=COLS,text.col=COLS,cex=0.7)
    Pval <- modelFit(get(LOD.list2[i+1]))[[5]][2]
    mtext(paste0("FCT used: ",LOD.list3[i+1],"    Lack of fit test: p = ",Pval),side=3,
          cex=0.5)
  }
  if(is.na(LOD.list2[i+1])) {
    plot(DAT2$Rate[DAT2$Target==Targets[i]]~log10(DAT2$Standards[DAT2$Target==Targets[i]]),
         ylim=c(0,1),ylab="Detection Probability",
         xlab=expression("Log of standard concentrations (Log"[10]*"Copies / Reaction)"),
         main=paste0("LoD for assay ",Targets[i]," unsolvable"),cex=0.6,cex.axis=0.6,
         cex.lab=0.6)
  }
}
dev.off()

png(filename="LOD-figure3.png",width=6.5,height=8,units="in",res=300)
par(mfrow=c(4,2),mar=c(3,3,2.5,0.5),mgp=c(1.5,0.75,0))
for(i in 17:24) {
  if(!is.na(LOD.list2[i+1])) {
    DAT4 <- rbind(ED(get(LOD.list2[i+1]),0.95,interval="delta",type="absolute"),
                  ED(get(LOD.list2[i+1]),1-sqrt(0.05),interval="delta",type="absolute"),
                  ED(get(LOD.list2[i+1]),1-0.05^(1/3),interval="delta",type="absolute"),
                  ED(get(LOD.list2[i+1]),1-0.05^0.25,interval="delta",type="absolute"),
                  ED(get(LOD.list2[i+1]),1-0.05^0.2,interval="delta",type="absolute"),
                  ED(get(LOD.list2[i+1]),1-0.05^0.125,interval="delta",type="absolute"))
    if(substr(LOD.list3[i+1],1,3)=="LL2") {
      DAT4 <- exp(DAT4)
    }
    DAT4 <- data.frame(DAT4,LoD=c("1rep.LOD","2rep.LOD","3rep.LOD","4rep.LOD",
                                  "5rep.LOD","8rep.LOD"),
                       Assay=rep(Targets[i],nrow(DAT4)))
    DAT4$Assay <- as.character(DAT4$Assay)
    if(i==1) {
      LOD.CI <- DAT4
    }
    if(i>1) {
      LOD.CI <- rbind(LOD.CI,DAT4)
    }
    if(sum(!is.na(DAT4[,3])&DAT4[,3]<=0)>0) { #Unable to plot negative lower limits, converting any lower limit values to 0.0001
      DAT4[!is.na(DAT4[,3])&DAT4[,3]<=0,3] <- 0.0001
    }
    plot(get(LOD.list2[i+1]),main=paste0("LoD Plot for assay ",Targets[i]),
         ylab="Detection Probability",xlab="Standard concentrations (Copies / Reaction)",
         xlim=c(min(DAT4[,1:4],na.rm=T),max(DAT$SQ,na.rm=T)),cex=0.6,cex.axis=0.6,
         cex.lab=0.6)
    LODS <- sum(!is.na(DAT4[,1]))
    COLS <- c(rgb(0.8,0.47,0.65),rgb(0,0.45,0.7),rgb(0.94,0.89,0.26),
              rgb(0.84,0.37,0),rgb(0,0.62,0.45),rgb(0.90,0.62,0))
    PNTS <- c(15,16,17,18,25,3)
    YS <- c(0.95,1-sqrt(0.05),1-0.05^(1/3),1-0.05^0.25,1-0.05^0.2,1-0.05^0.125)
    LODS2 <- c("Limit of Detection","2 Replicates LoD","3 Replicates LoD",
               "4 Replicates LoD","5 Replicates LoD","8 Replicates LoD")
    if(LODS<6) {
      LODS2[(LODS+1):6] <- gsub("licates LoD","s: Insufficient Data",LODS2[(LODS+1):6])
    }
    points(x=DAT4[1:LODS,1],y=YS[1:LODS],pch=PNTS,col=COLS,cex=0.9)
    for(j in 1:LODS) {
      lines(x=DAT4[j,3:4],y=rep(YS[j],2),col=COLS[j],lwd=1)
      lines(x=rep(DAT4[j,3],2),y=c(YS[j]-0.02,YS[j]+0.02),lwd=1,col=COLS[j])
      lines(x=rep(DAT4[j,4],2),y=c(YS[j]-0.02,YS[j]+0.02),lwd=1,col=COLS[j])
    }
    legend("bottomright",legend=LODS2,pch=PNTS,col=COLS,text.col=COLS,cex=0.7)
    Pval <- modelFit(get(LOD.list2[i+1]))[[5]][2]
    mtext(paste0("FCT used: ",LOD.list3[i+1],"    Lack of fit test: p = ",Pval),side=3,
          cex=0.5)
  }
  if(is.na(LOD.list2[i+1])) {
    plot(DAT2$Rate[DAT2$Target==Targets[i]]~log10(DAT2$Standards[DAT2$Target==Targets[i]]),
         ylim=c(0,1),ylab="Detection Probability",
         xlab=expression("Log of standard concentrations (Log"[10]*"Copies / Reaction)"),
         main=paste0("LoD for assay ",Targets[i]," unsolvable"),cex=0.6,cex.axis=0.6,
         cex.lab=0.6)
  }
}
dev.off()

png(filename="LOD-figure4.png",width=6.5,height=8,units="in",res=300)
par(mfrow=c(4,2),mar=c(3,3,2.5,0.5),mgp=c(1.5,0.75,0))
for(i in 25:32) {
  if(!is.na(LOD.list2[i+1])) {
    DAT4 <- rbind(ED(get(LOD.list2[i+1]),0.95,interval="delta",type="absolute"),
                  ED(get(LOD.list2[i+1]),1-sqrt(0.05),interval="delta",type="absolute"),
                  ED(get(LOD.list2[i+1]),1-0.05^(1/3),interval="delta",type="absolute"),
                  ED(get(LOD.list2[i+1]),1-0.05^0.25,interval="delta",type="absolute"),
                  ED(get(LOD.list2[i+1]),1-0.05^0.2,interval="delta",type="absolute"),
                  ED(get(LOD.list2[i+1]),1-0.05^0.125,interval="delta",type="absolute"))
    if(substr(LOD.list3[i+1],1,3)=="LL2") {
      DAT4 <- exp(DAT4)
    }
    DAT4 <- data.frame(DAT4,LoD=c("1rep.LOD","2rep.LOD","3rep.LOD","4rep.LOD",
                                  "5rep.LOD","8rep.LOD"),
                       Assay=rep(Targets[i],nrow(DAT4)))
    DAT4$Assay <- as.character(DAT4$Assay)
    if(i==1) {
      LOD.CI <- DAT4
    }
    if(i>1) {
      LOD.CI <- rbind(LOD.CI,DAT4)
    }
    if(sum(!is.na(DAT4[,3])&DAT4[,3]<=0)>0) { #Unable to plot negative lower limits, converting any lower limit values to 0.0001
      DAT4[!is.na(DAT4[,3])&DAT4[,3]<=0,3] <- 0.0001
    }
    plot(get(LOD.list2[i+1]),main=paste0("LoD Plot for assay ",Targets[i]),
         ylab="Detection Probability",xlab="Standard concentrations (Copies / Reaction)",
         xlim=c(min(DAT4[,1:4],na.rm=T),max(DAT$SQ,na.rm=T)),cex=0.6,cex.axis=0.6,
         cex.lab=0.6)
    LODS <- sum(!is.na(DAT4[,1]))
    COLS <- c(rgb(0.8,0.47,0.65),rgb(0,0.45,0.7),rgb(0.94,0.89,0.26),
              rgb(0.84,0.37,0),rgb(0,0.62,0.45),rgb(0.90,0.62,0))
    PNTS <- c(15,16,17,18,25,3)
    YS <- c(0.95,1-sqrt(0.05),1-0.05^(1/3),1-0.05^0.25,1-0.05^0.2,1-0.05^0.125)
    LODS2 <- c("Limit of Detection","2 Replicates LoD","3 Replicates LoD",
               "4 Replicates LoD","5 Replicates LoD","8 Replicates LoD")
    if(LODS<6) {
      LODS2[(LODS+1):6] <- gsub("licates LoD","s: Insufficient Data",LODS2[(LODS+1):6])
    }
    points(x=DAT4[1:LODS,1],y=YS[1:LODS],pch=PNTS,col=COLS,cex=0.9)
    for(j in 1:LODS) {
      lines(x=DAT4[j,3:4],y=rep(YS[j],2),col=COLS[j],lwd=1)
      lines(x=rep(DAT4[j,3],2),y=c(YS[j]-0.02,YS[j]+0.02),lwd=1,col=COLS[j])
      lines(x=rep(DAT4[j,4],2),y=c(YS[j]-0.02,YS[j]+0.02),lwd=1,col=COLS[j])
    }
    legend("bottomright",legend=LODS2,pch=PNTS,col=COLS,text.col=COLS,cex=0.7)
    Pval <- modelFit(get(LOD.list2[i+1]))[[5]][2]
    mtext(paste0("FCT used: ",LOD.list3[i+1],"    Lack of fit test: p = ",Pval),side=3,
          cex=0.5)
  }
  if(is.na(LOD.list2[i+1])) {
    plot(DAT2$Rate[DAT2$Target==Targets[i]]~log10(DAT2$Standards[DAT2$Target==Targets[i]]),
         ylim=c(0,1),ylab="Detection Probability",
         xlab=expression("Log of standard concentrations (Log"[10]*"Copies / Reaction)"),
         main=paste0("LoD for assay ",Targets[i]," unsolvable"),cex=0.6,cex.axis=0.6,
         cex.lab=0.6)
  }
}
dev.off()

png(filename="LOD-figure5.png",width=6.5,height=4,units="in",res=300)
par(mfrow=c(2,2),mar=c(3,3,2.5,0.5),mgp=c(1.5,0.75,0))
for(i in 33:36) {
  if(!is.na(LOD.list2[i+1])) {
    DAT4 <- rbind(ED(get(LOD.list2[i+1]),0.95,interval="delta",type="absolute"),
                  ED(get(LOD.list2[i+1]),1-sqrt(0.05),interval="delta",type="absolute"),
                  ED(get(LOD.list2[i+1]),1-0.05^(1/3),interval="delta",type="absolute"),
                  ED(get(LOD.list2[i+1]),1-0.05^0.25,interval="delta",type="absolute"),
                  ED(get(LOD.list2[i+1]),1-0.05^0.2,interval="delta",type="absolute"),
                  ED(get(LOD.list2[i+1]),1-0.05^0.125,interval="delta",type="absolute"))
    if(substr(LOD.list3[i+1],1,3)=="LL2") {
      DAT4 <- exp(DAT4)
    }
    DAT4 <- data.frame(DAT4,LoD=c("1rep.LOD","2rep.LOD","3rep.LOD","4rep.LOD",
                                  "5rep.LOD","8rep.LOD"),
                       Assay=rep(Targets[i],nrow(DAT4)))
    DAT4$Assay <- as.character(DAT4$Assay)
    if(i==1) {
      LOD.CI <- DAT4
    }
    if(i>1) {
      LOD.CI <- rbind(LOD.CI,DAT4)
    }
    if(sum(!is.na(DAT4[,3])&DAT4[,3]<=0)>0) { #Unable to plot negative lower limits, converting any lower limit values to 0.0001
      DAT4[!is.na(DAT4[,3])&DAT4[,3]<=0,3] <- 0.0001
    }
    plot(get(LOD.list2[i+1]),main=paste0("LoD Plot for assay ",Targets[i]),
         ylab="Detection Probability",xlab="Standard concentrations (Copies / Reaction)",
         xlim=c(min(DAT4[,1:4],na.rm=T),max(DAT$SQ,na.rm=T)),cex=0.6,cex.axis=0.6,
         cex.lab=0.6)
    LODS <- sum(!is.na(DAT4[,1]))
    COLS <- c(rgb(0.8,0.47,0.65),rgb(0,0.45,0.7),rgb(0.94,0.89,0.26),
              rgb(0.84,0.37,0),rgb(0,0.62,0.45),rgb(0.90,0.62,0))
    PNTS <- c(15,16,17,18,25,3)
    YS <- c(0.95,1-sqrt(0.05),1-0.05^(1/3),1-0.05^0.25,1-0.05^0.2,1-0.05^0.125)
    LODS2 <- c("Limit of Detection","2 Replicates LoD","3 Replicates LoD",
               "4 Replicates LoD","5 Replicates LoD","8 Replicates LoD")
    if(LODS<6) {
      LODS2[(LODS+1):6] <- gsub("licates LoD","s: Insufficient Data",LODS2[(LODS+1):6])
    }
    points(x=DAT4[1:LODS,1],y=YS[1:LODS],pch=PNTS,col=COLS,cex=0.9)
    for(j in 1:LODS) {
      lines(x=DAT4[j,3:4],y=rep(YS[j],2),col=COLS[j],lwd=1)
      lines(x=rep(DAT4[j,3],2),y=c(YS[j]-0.02,YS[j]+0.02),lwd=1,col=COLS[j])
      lines(x=rep(DAT4[j,4],2),y=c(YS[j]-0.02,YS[j]+0.02),lwd=1,col=COLS[j])
    }
    legend("bottomright",legend=LODS2,pch=PNTS,col=COLS,text.col=COLS,cex=0.7)
    Pval <- modelFit(get(LOD.list2[i+1]))[[5]][2]
    mtext(paste0("FCT used: ",LOD.list3[i+1],"    Lack of fit test: p = ",Pval),side=3,
          cex=0.5)
  }
  if(is.na(LOD.list2[i+1])) {
    plot(DAT2$Rate[DAT2$Target==Targets[i]]~log10(DAT2$Standards[DAT2$Target==Targets[i]]),
         ylim=c(0,1),ylab="Detection Probability",
         xlab=expression("Log of standard concentrations (Log"[10]*"Copies / Reaction)"),
         main=paste0("LoD for assay ",Targets[i]," unsolvable"),cex=0.6,cex.axis=0.6,
         cex.lab=0.6)
  }
}
dev.off()

## Generate LOQ plots:
POLYS <- data.frame(Target=Targets,Poly=rep(0,length(Targets)))
for(i in 1:length(Targets)) {
  if(is.na(LOQ.list[i+1])==F) {
    ## Re-generate prediction data for the model:
    newData <- data.frame(Standards = seq(1, 10000))
    newData$Cq.CV <- predict(get(LOQ.list[i+1]), newData)
    ## Define LOQ polygon coordinates:
    PDAT <- data.frame(x=c(min(DAT2$Standards[DAT2$Target==Targets[i]&!is.na(DAT2$Cq.CV)]),
                           min(DAT2$Standards[DAT2$Target==Targets[i]&!is.na(DAT2$Cq.CV)]),
                           DAT3$LOQ[DAT3$Assay==Targets[i]],
                           DAT3$LOQ[DAT3$Assay==Targets[i]]),
                       y=c(min(c(DAT2$Cq.CV[DAT2$Target==Targets[i]],newData$Cq.CV[newData$Standards<=max(DAT2$Standards[DAT2$Target==Targets[i]],na.rm=T)&newData$Standards>=min(DAT2$Standards[DAT2$Target==Targets[i]],na.rm=T)]),na.rm=T)*0.9,
                           newData$Cq.CV[newData$Standards==DAT3$LOQ[DAT3$Assay==Targets[i]]],
                           newData$Cq.CV[newData$Standards==DAT3$LOQ[DAT3$Assay==Targets[i]]],
                           min(c(DAT2$Cq.CV[DAT2$Target==Targets[i]],newData$Cq.CV[newData$Standards<=max(DAT2$Standards[DAT2$Target==Targets[i]],na.rm=T)&newData$Standards>=min(DAT2$Standards[DAT2$Target==Targets[i]],na.rm=T)]),na.rm=T)*0.9))
    if(DAT3$LOQ[DAT3$Assay==Targets[i]]!=floor(DAT3$LOQ[DAT3$Assay==Targets[i]])) {
      PDAT$y[2:3] <- LOQ.Threshold
    }
  }
  
  Decay.Plot <- ggplot(DAT2[DAT2$Target==Targets[i]&!is.na(DAT2$Cq.CV),], aes(x= Standards, y = Cq.CV)) +
    geom_point(size=2) +
    scale_x_continuous(trans = 'log10') +
    ylab("Coefficient of variation for Cq-Values") +
    xlab("Standard concentrations (Copies / Reaction)") +
    geom_vline(xintercept=DAT3$LOD[DAT3$Assay==Targets[i]],color="red") +
    theme(legend.position="none") +
    theme(plot.title=element_text(hjust=0.5,size=9))
  
  if(is.na(LOQ.list[i+1])==F) {
    if(DAT3$LOQ[DAT3$Assay==Targets[i]]<=min(DAT2$Standards[DAT2$Target==Targets[i]])) {
      PDAT$x[3:4] <- NA
      Decay.Plot <- Decay.Plot + 
        annotate("text",y=max(DAT2$Cq.CV[DAT2$Target==Targets[i]],na.rm=T)*0.99,
                 x=median(DAT2$Standards[DAT2$Target==Targets[i]])/2,
                 label="LOQ may be outside tested range.",hjust=0,size=2)
    }
    assign(paste0("PDAT",i),PDAT)
    #Decay.Plot <- Decay.Plot + geom_polygon(data=get(paste0("PDAT",i)),
                                            #aes(x=x,y=y,alpha=0.5))
    
    if(as.character(get(LOQ.list[i+1])$call)[1]=="nls") {
      Decay.Plot <- Decay.Plot + 
        stat_smooth(method = "nls", formula = y ~ SSasymp(x, Asym, R0, lrc), se = FALSE) +
        ggtitle(paste0("Exponential Decay LOQ model for: ",Targets[i]))
    }
    
    if(as.character(get(LOQ.list[i+1])$call)[1]=="lm") {
      if(grepl("poly",as.character(get(LOQ.list[i+1])$call)[2])==T) {
        POLYS[i,2] <- length(get(LOQ.list[i+1])$coefficients)-1
        Decay.Plot <- Decay.Plot +
          ggtitle(paste0(POLYS[i,2],"-order polynomial LOQ model for: ",Targets[i]))
      }
      else {
        Decay.Plot <- Decay.Plot +
          stat_smooth(method = "lm", formula = y ~ x,se=F) +
          ggtitle(paste0("Linear LOQ model for: ",Targets[i]))
      }
    } 
  }
  
  if(is.na(LOQ.list[i+1])==T) {
    Decay.Plot <- Decay.Plot + 
      ggtitle(paste0("LOQ model for: ",Targets[i]," not solvable."))
  }
  Decay.Plot <- Decay.Plot +
    geom_polygon(data=get(paste0("PDAT",i)),aes(x=x,y=y,alpha=0.5)) +
    theme(text=element_text(size=7)) +
    annotate("text",y=max(DAT2$Cq.CV[DAT2$Target==Targets[i]],na.rm=T)*0.98,
             x=DAT3$LOD[DAT3$Assay==Targets[i]]*0.8,angle=90,label="LOD",
             color="red",size=2) +
    theme(axis.text=element_text(size=7))
  assign(paste0("LOQ.Plot",i),Decay.Plot)
}
LOQ.Plot1 <- LOQ.Plot1 + stat_smooth(method = "lm",formula = y ~ poly(x,4),se=F)
LOQ.Plot2 <- LOQ.Plot2 + stat_smooth(method = "lm",formula = y ~ poly(x,2),se=F)
LOQ.Plot3 <- LOQ.Plot3 + stat_smooth(method = "lm",formula = y ~ poly(x,4),se=F)
LOQ.Plot4 <- LOQ.Plot4 + stat_smooth(method = "lm",formula = y ~ poly(x,5),se=F)
LOQ.Plot6 <- LOQ.Plot6 + stat_smooth(method = "lm",formula = y ~ poly(x,4),se=F)
LOQ.Plot7 <- LOQ.Plot7 + stat_smooth(method = "lm",formula = y ~ poly(x,6),se=F)
LOQ.Plot8 <- LOQ.Plot8 + stat_smooth(method = "lm",formula = y ~ poly(x,3),se=F)
LOQ.Plot14 <- LOQ.Plot14 + stat_smooth(method = "lm",formula = y ~ poly(x,5),se=F)
LOQ.Plot15 <- LOQ.Plot15 + stat_smooth(method = "lm",formula = y ~ poly(x,3),se=F)
LOQ.Plot16 <- LOQ.Plot16 + stat_smooth(method = "lm",formula = y ~ poly(x,6),se=F)
LOQ.Plot17 <- LOQ.Plot17 + stat_smooth(method = "lm",formula = y ~ poly(x,3),se=F)
LOQ.Plot18 <- LOQ.Plot18 + stat_smooth(method = "lm",formula = y ~ poly(x,6),se=F)
LOQ.Plot20 <- LOQ.Plot20 + stat_smooth(method = "lm",formula = y ~ poly(x,3),se=F)
LOQ.Plot21 <- LOQ.Plot21 + stat_smooth(method = "lm",formula = y ~ poly(x,6),se=F)
LOQ.Plot22 <- LOQ.Plot22 + stat_smooth(method = "lm",formula = y ~ poly(x,6),se=F)
LOQ.Plot23 <- LOQ.Plot23 + stat_smooth(method = "lm",formula = y ~ poly(x,6),se=F)
LOQ.Plot24 <- LOQ.Plot24 + stat_smooth(method = "lm",formula = y ~ poly(x,4),se=F)
LOQ.Plot25 <- LOQ.Plot25 + stat_smooth(method = "lm",formula = y ~ poly(x,6),se=F)
LOQ.Plot26 <- LOQ.Plot26 + stat_smooth(method = "lm",formula = y ~ poly(x,6),se=F)
LOQ.Plot27 <- LOQ.Plot27 + stat_smooth(method = "lm",formula = y ~ poly(x,6),se=F)
LOQ.Plot28 <- LOQ.Plot28 + stat_smooth(method = "lm",formula = y ~ poly(x,4),se=F)
LOQ.Plot29 <- LOQ.Plot29 + stat_smooth(method = "lm",formula = y ~ poly(x,6),se=F)
LOQ.Plot30 <- LOQ.Plot30 + stat_smooth(method = "lm",formula = y ~ poly(x,6),se=F)
LOQ.Plot31 <- LOQ.Plot31 + stat_smooth(method = "lm",formula = y ~ poly(x,6),se=F)
LOQ.Plot33 <- LOQ.Plot33 + stat_smooth(method = "lm",formula = y ~ poly(x,6),se=F)
LOQ.Plot34 <- LOQ.Plot34 + stat_smooth(method = "lm",formula = y ~ poly(x,6),se=F)

## Apply a higher polynomial model for assay 23 LOQ:
LOQ.mod37 <- lm(Cq.CV~poly(log10(Standards),7),data=DAT2[DAT2$Target==23,])

## Check to make sure 7th order polynomial is a better fit:
summary(LOQ.mod37)$sigma
summary(get(LOQ.list[24]))$sigma

## Generate model data and get LOQ:
newData <- data.frame(Standards = seq(1, 10000))
newData$Cq.CV <- predict(LOQ.mod37, newData)
A <- max(DAT2$Standards[DAT2$Target==23])
B <- LOQ.Threshold
C <- max(newData$Standards[newData$Cq.CV<=B&newData$Standards<=A])
D <- max(newData$Standards[newData$Cq.CV>B&newData$Standards<C])
DAT3$LOQ[DAT3$Assay==23] <- D+1

## Generate new LOQ plot:
PDAT37 <- data.frame(x=c(min(DAT2$Standards[DAT2$Target==23&!is.na(DAT2$Cq.CV)]),
                       min(DAT2$Standards[DAT2$Target==23&!is.na(DAT2$Cq.CV)]),
                       32,32),
                   y=c(min(c(DAT2$Cq.CV[DAT2$Target==23],newData$Cq.CV[newData$Standards<=max(DAT2$Standards[DAT2$Target==23],na.rm=T)&newData$Standards>=min(DAT2$Standards[DAT2$Target==23],na.rm=T)]),na.rm=T)*0.9,
                       newData$Cq.CV[newData$Standards==32],
                       newData$Cq.CV[newData$Standards==32],
                       min(c(DAT2$Cq.CV[DAT2$Target==23],newData$Cq.CV[newData$Standards<=max(DAT2$Standards[DAT2$Target==23],na.rm=T)&newData$Standards>=min(DAT2$Standards[DAT2$Target==23],na.rm=T)]),na.rm=T)*0.9))
LOQ.Plot37 <- ggplot(DAT2[DAT2$Target==23&!is.na(DAT2$Cq.CV),], aes(x= Standards, y = Cq.CV)) +
  geom_point(size=2) +
  scale_x_continuous(trans = 'log10') +
  ylab("Coefficient of variation for Cq-Values") +
  xlab("Standard concentrations (Copies / Reaction)") +
  geom_vline(xintercept=DAT3$LOD[DAT3$Assay==23],color="red") +
  theme(legend.position="none") +
  theme(plot.title=element_text(hjust=0.5,size=9)) +
  geom_polygon(data=PDAT37,aes(x=x,y=y,alpha=0.75)) +
  theme(text=element_text(size=7)) +
  annotate("text",y=max(DAT2$Cq.CV[DAT2$Target==23],na.rm=T)*0.98,
           x=DAT3$LOD[DAT3$Assay==23]*0.8,angle=90,label="LOD",
           color="red",size=2) +
  theme(axis.text=element_text(size=7)) +
  ggtitle("7-order polynomial LOQ model for: 23") +
  stat_smooth(method = "lm",formula = y ~ poly(x,7),se=F)

## Adjust precision threshold and recalculate LOQ for assay 31:
newData <- data.frame(Standards = seq(1, 10000))
newData$Cq.CV <- predict(get(LOQ.list[32]), newData)
B <- min(DAT2$Cq.CV[DAT2$Target==31],na.rm=T)*1.5
C <- max(newData$Standards[newData$Cq.CV<=B&newData$Standards<=A])
D <- max(newData$Standards[newData$Cq.CV>B&newData$Standards<C])
DAT3$LOQ[DAT3$Assay==31] <- D+1

## Generate new LOQ plot:
PDAT38 <- data.frame(x=c(min(DAT2$Standards[DAT2$Target==31&!is.na(DAT2$Cq.CV)]),
                         min(DAT2$Standards[DAT2$Target==31&!is.na(DAT2$Cq.CV)]),
                         239,239),
                     y=c(min(c(DAT2$Cq.CV[DAT2$Target==31],newData$Cq.CV[newData$Standards<=max(DAT2$Standards[DAT2$Target==31],na.rm=T)&newData$Standards>=min(DAT2$Standards[DAT2$Target==31],na.rm=T)]),na.rm=T)*0.9,
                         newData$Cq.CV[newData$Standards==239],
                         newData$Cq.CV[newData$Standards==239],
                         min(c(DAT2$Cq.CV[DAT2$Target==31],newData$Cq.CV[newData$Standards<=max(DAT2$Standards[DAT2$Target==31],na.rm=T)&newData$Standards>=min(DAT2$Standards[DAT2$Target==31],na.rm=T)]),na.rm=T)*0.9))
LOQ.Plot38 <- ggplot(DAT2[DAT2$Target==31&!is.na(DAT2$Cq.CV),], aes(x= Standards, y = Cq.CV)) +
  geom_point(size=2) +
  scale_x_continuous(trans = 'log10') +
  ylab("Coefficient of variation for Cq-Values") +
  xlab("Standard concentrations (Copies / Reaction)") +
  geom_vline(xintercept=DAT3$LOD[DAT3$Assay==31],color="red") +
  theme(legend.position="none") +
  theme(plot.title=element_text(hjust=0.5,size=9)) +
  geom_polygon(data=PDAT38,aes(x=x,y=y,alpha=0.75)) +
  theme(text=element_text(size=7)) +
  annotate("text",y=max(DAT2$Cq.CV[DAT2$Target==31],na.rm=T)*0.98,
           x=DAT3$LOD[DAT3$Assay==31]*0.8,angle=90,label="LOD",
           color="red",size=2) +
  theme(axis.text=element_text(size=7)) +
  ggtitle("6-order polynomial LOQ model for: 31") +
  stat_smooth(method = "lm",formula = y ~ poly(x,6),se=F)

## Generate LOQ figures:
png(filename="LOQ-figure.png",width=6.5,height=8,units="in",res=300)
plot_grid(LOQ.Plot1,LOQ.Plot2,LOQ.Plot3,LOQ.Plot4,LOQ.Plot5,LOQ.Plot6,LOQ.Plot7,LOQ.Plot8,
          ncol=2,nrow=4)
dev.off()
png(filename="LOQ-figure2.png",width=6.5,height=8,units="in",res=300)
plot_grid(LOQ.Plot9,LOQ.Plot10,LOQ.Plot11,LOQ.Plot12,LOQ.Plot13,LOQ.Plot14,LOQ.Plot15,
          LOQ.Plot16,ncol=2,nrow=4)
dev.off()
png(filename="LOQ-figure3.png",width=6.5,height=8,units="in",res=300)
plot_grid(LOQ.Plot17,LOQ.Plot18,LOQ.Plot19,LOQ.Plot20,LOQ.Plot21,LOQ.Plot22,LOQ.Plot23,
          LOQ.Plot37,ncol=2,nrow=4)
dev.off()
png(filename="LOQ-figure4.png",width=6.5,height=8,units="in",res=300)
plot_grid(LOQ.Plot24,LOQ.Plot25,LOQ.Plot26,LOQ.Plot27,LOQ.Plot28,LOQ.Plot29,LOQ.Plot30,
          LOQ.Plot32,ncol=2,nrow=4)
dev.off()
png(filename="LOQ-figure5.png",width=6.5,height=6,units="in",res=300)
plot_grid(LOQ.Plot31,LOQ.Plot38,LOQ.Plot33,LOQ.Plot34,LOQ.Plot35,LOQ.Plot36,
          ncol=2,nrow=3)
dev.off()



######
######
##
## Compare LoD models
##
######
######

## Run LoD-Calculator script with LOD.FCT specified to each sigmoidal function and
##   run this code to pull out data between each run:
## (Do it in order: LL.2, LL.3, LL.3u, LL.4, LL.5, W1.2, W1.3, W1.4, W2.2, W2.3,
##                  W2.4, AR.2, AR.3, MM.2, MM.3)

## After running LL.2():
COMP.DAT <- DAT3

## After running the rest (IN ORDER!!! run repeatedly for each function):
COMP.DAT <- rbind(COMP.DAT,DAT3)

## Add in function details:
COMP.DAT$Function <- as.factor(c(rep("LL.2",36),rep("LL.3",36),rep("LL.3u",36),
                                 rep("LL.4",36),rep("LL.5",36),rep("W1.2",36),
                                 rep("W1.3",36),rep("W1.4",36),rep("W2.2",36),
                                 rep("W2.3",36),rep("W2.4",36),rep("AR.2",36),
                                 rep("AR.3",36),rep("MM.2",36),rep("MM.3",36)))
COMP.DAT$Type <- as.factor(c(rep("Log-Logistic",180),rep("Weibull Type I",108),
                             rep("Weibull Type II",108),rep("Asymptotic Regression",72),
                             rep("Michaelis-Menten",72)))
COMP.DAT$NoParm <- c(rep(2,36),rep(3,72),rep(4,36),rep(5,36),rep(2,36),rep(3,36),
                     rep(4,36),rep(2,36),rep(3,36),rep(4,36),rep(2,36),rep(3,36),
                     rep(2,36),rep(3,36))
COMP.DAT$Solved <- as.numeric(!is.na(COMP.DAT$LOD))

## Export Comp Data in case I need to close R and come back to it:
write.csv(COMP.DAT,file="Comparison-Data.csv",row.names=F)
COMP.DAT <- read.csv(file="Comparison-Data.csv",header=T)

## Rename assays by number:
for(i in 1:nrow(TarMatch)) {
  COMP.DAT$Assay <- gsub(TarMatch[i,1],TarMatch[i,2],COMP.DAT$Assay)
}
COMP.DAT$Assay <- factor(COMP.DAT$Assay,levels=sort(unique(as.numeric(COMP.DAT$Assay))))

## Exploring - Do function type or number of parameters have any consistent effects?:
COMP.DAT$Type <- factor(COMP.DAT$Type,levels=c("Asymptotic Regression",
                        "Log-Logistic","Michaelis-Menten","Weibull Type I",
                        "Weibull Type II"))
REP1.mod <- lm(LOD~Assay+Type+NoParm,data=COMP.DAT)
summary(REP1.mod)
COMP.DAT$Type <- factor(COMP.DAT$Type,levels=c("Weibull Type II",
                        "Log-Logistic","Michaelis-Menten","Weibull Type I",
                        "Asymptotic Regression"))
REP1.mod2 <- lm(LOD~Assay+Type+NoParm,data=COMP.DAT)
summary(REP1.mod2)
COMP.DAT$Type <- factor(COMP.DAT$Type,levels=c("Log-Logistic",
                        "Michaelis-Menten","Weibull Type I","Weibull Type II",
                        "Asymptotic Regression"))
REP1.mod3 <- lm(LOD~Assay+Type+NoParm,data=COMP.DAT)
summary(REP1.mod3)
COMP.DAT$Type <- factor(COMP.DAT$Type,levels=c("Weibull Type I",
                        "Michaelis-Menten","Weibull Type II","Log-Logistic",
                        "Asymptotic Regression"))
REP1.mod4 <- lm(LOD~Assay+Type+NoParm,data=COMP.DAT)
summary(REP1.mod4)

## Plot LOD vs Function by assay:
COMP.DAT$Function <- factor(COMP.DAT$Function,levels=c("W2.2","W2.3","W2.4",
                       "AR.2","AR.3","LL.2","LL.3","LL.3u","LL.4","LL.5",
                       "MM.2","MM.3","W1.2","W1.3","W1.4"))
LODvFCT <- ggplot(data=COMP.DAT,aes(x=Function,y=LOD,col=Assay,group=Assay)) +
  geom_point(stat="summary",fun.y=sum) +
  stat_summary(fun.y=sum,geom="line") +
  scale_y_log10() +
  theme(text=element_text(size=8),
        axis.text=element_text(size=8),
        axis.title=element_text(size=10),
        legend.key.size=unit(0.7,"lines"))
ggsave(filename="LODvFCT.png",plot=LODvFCT,width=6.5,height=4,units="in",dpi=300)

## Plot LOD vs Function Type:
COMP.DAT$Type <- factor(COMP.DAT$Type,levels=c("Weibull Type II",
                   "Asymptotic Regression","Log-Logistic","Michaelis-Menten",
                   "Weibull Type I"))
LODvTYP <- ggplot(data=COMP.DAT,aes(x=Type,y=LOD)) +
  geom_bar(stat="summary",fun.y="mean") +
  ylab("Mean LOD") + xlab("Function Type") +
  coord_cartesian(ylim=c(0,60)) +
  theme(axis.text.x=element_text(angle=45,hjust=1))

## Sum and plot solved data:
COMP.DAT2 <- data.frame(Function=unique(COMP.DAT$Function),
                        Solved=rep(0,length(unique(COMP.DAT$Function))))
for(i in 1:nrow(COMP.DAT2)) {
  COMP.DAT2$Solved[i] <- sum(COMP.DAT$Solved[COMP.DAT$Function==COMP.DAT2$Function[i]])
}
ggplot(data=COMP.DAT2,aes(x=Function,y=Solved)) +
  geom_bar(stat="identity") +
  scale_y_continuous(limits=c(0,35),breaks=c(0,10,20,30)) +
  ylab("Number of Assays Solvable (out of 36)")
print(COMP.DAT2)

## Determine the range of LOD's and LOQ's:
DAT3[DAT3$LOD==min(DAT3$LOD,na.rm=T),]
DAT3[DAT3$LOD==max(DAT3$LOD,na.rm=T),]
mean(DAT3$LOD,na.rm=T)
sd(DAT3$LOD,na.rm=T)
DAT3[DAT3$LOQ==min(DAT3$LOQ),]
DAT3[DAT3$LOQ==max(DAT3$LOQ),]
mean(DAT3$LOQ)
sd(DAT3$LOQ)

## Rearrange data and plot Effective LOD vs number of reps:
DAT5 <- data.frame(Assay=rep(DAT3$Assay,6),LOD=c(DAT3$LOD,DAT3$rep2.LOD,DAT3$rep3.LOD,
                    DAT3$rep4.LOD,DAT3$rep5.LOD,DAT3$rep8.LOD),Reps=c(rep(1,36),
                    rep(2,36),rep(3,36),rep(4,36),rep(5,36),rep(8,36)))
DAT5$Assay <- factor(DAT5$Assay,levels=sort(unique(as.numeric(as.character(DAT5$Assay)))))
LODvREP <- ggplot(data=DAT5,aes(x=Reps,y=LOD,color=Assay,group=Assay)) +
  geom_point(stat="summary",fun.y=sum) +
  stat_summary(fun.y=sum,geom="line") +
  scale_y_log10() +
  ylab("Effective LOD") +
  xlab("Number of Replicates Analyzed") +
  theme(text=element_text(size=8),
        axis.text=element_text(size=8),
        axis.title=element_text(size=10),
        legend.key.size=unit(0.7,"lines"))
ggsave(filename="LODvREP.png",plot=LODvREP,width=4,height=4,units="in",dpi=300)

## Get the relative drops:
REL.DROP <- data.frame(Assay=as.character(c(1:36)),Drop=rep(NA,36))
for(i in 1:36) {
  REL.DROP[i,2] <- DAT5$LOD[DAT5$Assay==as.character(i)&DAT5$Reps==1]/DAT5$LOD[DAT5$Assay==as.character(i)&DAT5$Reps==8]
}
mean(REL.DROP$Drop,na.rm=T)
sd(REL.DROP$Drop,na.rm=T)

## Get example code plots for BH1:
i=27
LOD.smoother <- plot(get(LOD.list2[i+1]))
LOD.smoother2 <- data.frame(LOD.smoother)
colnames(LOD.smoother2) <- c("Standards","Rate")
DAT4 <- rbind(ED(get(LOD.list2[i+1]),0.95,interval="delta",type="absolute"),
              ED(get(LOD.list2[i+1]),1-sqrt(0.05),interval="delta",type="absolute"),
              ED(get(LOD.list2[i+1]),1-0.05^(1/3),interval="delta",type="absolute"),
              ED(get(LOD.list2[i+1]),1-0.05^0.25,interval="delta",type="absolute"),
              ED(get(LOD.list2[i+1]),1-0.05^0.2,interval="delta",type="absolute"),
              ED(get(LOD.list2[i+1]),1-0.05^0.125,interval="delta",type="absolute"))
DAT4 <- data.frame(DAT4,LoD=c("Limit of Detection","2 Replicates LoD",
                  "3 Replicates LoD","4 Replicates LoD","5 Replicates LoD",
                  "8 Replicates LoD"),Rate=c(0.95,1-sqrt(0.05),1-0.05^(1/3),
                  1-0.05^0.25,1-0.05^0.2,1-0.05^0.125))
colnames(DAT4) <- c("Point","SE","Standards","xend","LoD","Rate")
DAT4$LoD <- factor(DAT4$LoD,levels=c("Limit of Detection","2 Replicates LoD",
                   "3 Replicates LoD","4 Replicates LoD","5 Replicates LoD",
                   "8 Replicates LoD"))
Pval <- modelFit(get(LOD.list2[i+1]))[[5]][2]
BH1.LOD <- ggplot(data=DAT2[DAT2$Target==27,],aes(x=Standards,y=Rate)) +
  geom_line(data=LOD.smoother2) +
  scale_x_log10() +
  geom_point(shape=1,size=1) +
  ylab("Detection Probability") +
  xlab("Standard concentrations (Copies / Reaction)") +
  ggtitle(label=paste0("LoD Plot for: ",Targets[i])) +
  labs(subtitle=paste0("FCT used: ",LOD.list3[i+1],"    Lack of fit test: p = ",Pval)) +
  geom_segment(data=DAT4,aes(x=Standards,xend=xend,y=Rate,yend=Rate,col=LoD)) +
  geom_point(data=DAT4,aes(x=Point,y=Rate,col=LoD,shape=LoD),size=1.5) +
  geom_segment(data=DAT4,aes(x=Standards,xend=Standards,y=Rate-0.01,
                             yend=Rate+0.01,col=LoD)) +
  geom_segment(data=DAT4,aes(x=xend,xend=xend,y=Rate-0.01,yend=Rate+0.01,col=LoD)) +
  theme(plot.title=element_text(hjust=0.5,size=9),
        plot.subtitle=element_text(hjust=0.5,size=7),
        text=element_text(size=7),axis.text=element_text(size=7),
        legend.position=c(0.99,0),legend.justification=c(1,0),
        legend.title=element_blank(),legend.text=element_text(size=4),
        legend.key.size=unit(0.5,"lines"))

ggOut <- ggplot(data=DAT[DAT$Target==Targets[i]&is.na(DAT$SQ)==F,],
                aes(x=SQ,y=Cq,color=factor(Mod),shape=factor(Mod),size=factor(Mod))) + 
  geom_jitter(width=0.1,alpha=0.75) + 
  scale_shape_manual("",values=c(3,20),guide=F) +
  scale_size_manual("",values=c(0.3,1)) +
  scale_x_log10() +
  scale_color_manual("",values=c("blue", "black")) +
  xlab("Standard Concentrations (Copies / Reaction)") +
  ylab("Cq-value") +
  geom_abline(intercept=coef(get(curve.list[i+1]))[1],
              slope=coef(get(curve.list[i+1]))[2]) +
  geom_vline(xintercept=DAT3$LOD[i],linetype=2) +
  geom_vline(xintercept=DAT3$LOQ[i],colour="red") +
  annotate("text",y=max(DAT$Cq[DAT$Target==Targets[i]],na.rm=T)*0.99,
           x=DAT3$LOD[i]*0.8,angle=90,label="LOD",size=2) +
  annotate("text",y=max(DAT$Cq[DAT$Target==Targets[i]],na.rm=T)*0.94,color="red",
           x=DAT3$LOQ[i]*0.8,angle=90,label="LOQ",size=2) +
  theme_bw() + theme(legend.justification=c(1,1),legend.position=c(1,0.99)) +
  ggtitle(paste0("Standard curve for: ",Targets[i])) +
  theme(plot.title=element_text(hjust=0.5,size=9),
        axis.title=element_text(size=7),text=element_text(size=7)) +
  theme(legend.title=element_blank(),legend.key.size=unit(0.5,"lines"),
        legend.text=element_text(size=6)) +
  annotate("text",y=min(DAT$Cq[DAT$Target==Targets[i]&is.na(DAT$SQ)==F],na.rm=T)*1.05,
           x=min(DAT$SQ[DAT$Target==Targets[i]&is.na(DAT$SQ)==F],na.rm=T)*1.01,hjust=0,
           label=(paste0("R-squared: ",DAT3$R.squared[i],"\ny = ",DAT3$Slope[i],"x + ",DAT3$Intercept[i])),
           size=2)

png(filename="BH1-example.png",width=6.5,height=4,units="in",res=300)
plot_grid(ggOut,BH1.LOD,LOQ.Plot27,ncol=2,nrow=2)
dev.off()

## Rearrange and export Table 1:
TarMatch <- TarMatch[order(TarMatch[,2]),]
TarMatch$Modeled.LOD <- DAT3$LOD[match(TarMatch[,2],DAT3$Assay)]
TarMatch$Threshold.LOD <- DAT3$Low.95[match(TarMatch[,2],DAT3$Assay)]
TarMatch$LOQ <- DAT3$LOQ[match(TarMatch[,2],DAT3$Assay)]
TarMatch$Lab <- c(rep("CERC",4),rep("ERDC",2),"MNRF","NWRC",rep("UMESC",5),
                  rep("UVIC",11),rep("WGL",10),rep("WSU",2))
TarMatch$Modeled.LOD <- signif(TarMatch$Modeled.LOD,digits=3)
TarMatch$Threshold.LOD <- signif(TarMatch$Threshold.LOD,digits=3)
TarMatch$LOQ <- signif(TarMatch$LOQ,digits=3)
TarMatch <- TarMatch[,c(2,1,3:6)]
write.csv(TarMatch,file="Table1.csv",row.names=F)
