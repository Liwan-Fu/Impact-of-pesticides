library(metafor)
library(xlsx)
library(Hmisc)
library(export)
setwd('...')
combined=read.csv('Meta combined data.csv',header=T)

########## Figures 2a-d ###########

Pesticide_category=c('Pesticide',
                     unique(combined$pesticide_by_target_organisms))
dat_list1=lapply(unique(combined$Taxonomic_group_response),
                   function(x){subset(combined,Taxonomic_group_response==x)})
dat_fig2=matrix(NA,length(Pesticide_category)*length(dat_list1),3)
for (i in 1:length(Pesticide_category)) {
  for (j in 1:length(dat_list1)){
    mydat=subset(dat_list1[[j]],
                   pesticide_by_target_organisms==Pesticide_category[i])
    res=rma.mv(yi, vi, 
                 random = list(~ 1 | factor(Insecticide.name),~1|Code),
                 data=mydat)
    dat_fig2[(i-1)*length(dat_list1)+j,]=c(res$b,res$ci.lb,res$ci.ub)
  }
}
colnames(dat_fig2)<-c("Effect size","lb","ub")

### Figure 2
fig2a=dat_fig2
errbar(1:nrow(fig2a),fig2a[,1], fig2a[,2], fig2a[,3],pch=19,
       ann = FALSE,lwd=4,cex=2,
       col=c(c('blue','green','black','red'),
             rep(c('blue','green','red'),2)),
       errbar.col=c(c('blue','green','black','red'),
                    rep(c('blue','green','red'),2)))
mtext("Effect size",side=4,line=3,cex=1.5,family="sans")
abline(h=0,lty=2)
graph2jpg(file='Figure 2a-d',height=10,width=8.5,dpi=300)


########## Extended Data Figs. 2a-d ###########
dat_list2=c(lapply(unique(combined$Primary.classification.of.animal......一级分类)[c(2,1,4,3,5:6)],
                   function(x){
                     subset(combined,Primary.classification.of.animal......一级分类==x)
                   }))
dat_list3=c(lapply(unique(dat_list2[[1]]$Taxonomic_group_response),
                   function(x) {
                     subset(dat_list2[[1]],Taxonomic_group_response==x)
                   }),
            lapply(unique(dat_list2[[2]]$Taxonomic_group_response),
                   function(x) {
                     subset(dat_list2[[2]],Taxonomic_group_response==x)
                   }),
            lapply(unique(dat_list2[[3]]$Taxonomic_group_response),
                   function(x) {
                     subset(dat_list2[[3]],Taxonomic_group_response==x)
                   }),
            lapply(unique(dat_list2[[4]]$Taxonomic_group_response),
                   function(x) {
                     subset(dat_list2[[4]],Taxonomic_group_response==x)
                   }),
            lapply(unique(dat_list2[[5]]$Taxonomic_group_response),
                   function(x) {
                     subset(dat_list2[[5]],Taxonomic_group_response==x)
                   }),
            lapply(unique(dat_list2[[6]]$Taxonomic_group_response),
                   function(x) {
                     subset(dat_list2[[6]],Taxonomic_group_response==x)
                   }))
Pesticide_category=c('Pesticide',
                     unique(combined$pesticide_by_target_organisms))
dat_EDF1=matrix(NA,length(Pesticide_category)*length(dat_list3),3)
for (i in 1:length(Pesticide_category)) {
  for (j in 1:length(dat_list3)){
      mydat=subset(dat_list3[[j]],
                   pesticide_by_target_organisms==Pesticide_category[i])
    if (nrow(mydat)<=1) {
      dat_EDF1[(i-1)*length(dat_list3)+j,]=rep(NA,3)
    } else {
      res=rma.mv(yi, vi, 
                 random = list(~ 1 | factor(Insecticide.name),~1|Code),
                 data=mydat)
      dat_EDF1[(i-1)*length(dat_list3)+j,]=c(res$b,res$ci.lb,res$ci.ub)
    }
  }
}
colnames(dat_EDF1)<-c("Effect size","lb","ub")

### Extended Data Fig 2a-d
dat_EDF1a=dat_EDF1
errbar(1:nrow(dat_EDF1a),dat_EDF1a[,1], dat_EDF1a[,2], dat_EDF1a[,3],
       pch=19,
       ann = FALSE,lwd=4,cex=2,
       col=c(rep(c('blue','green','black','red'),2),
             rep(c('blue','green','red'),4)),
       errbar.col=c(rep(c('blue','green','black','red'),2),
                    rep(c('blue','green','red'),4)))
mtext("Effect size",side=4,line=3,cex=1.5,family="sans") #Arial font
abline(h=0,lty=2)
graph2jpg(file='Extended Data Fig 2a-d',height=10,width=12,dpi=300)


########## Extended Data Figs. 3a-k ###########
Pesticide_category=c('Pesticide',
                     unique(combined$pesticide_by_target_organisms))
dat_list1=lapply(unique(combined$Taxonomic_group_response),
                 function(x){subset(combined,Taxonomic_group_response==x)})
S6_category=c(unique(combined$pesticide_by_source),
              unique(combined$Primary.classification.of.insecticide........一级分类))[c(1,3,2,4,6,5,7:11)]
dat_EDF2=matrix(NA,length(S6_category)*length(dat_list1),3)
for (i in 1:length(S6_category)) {
  for (j in 1:length(dat_list1)){
    if (i<=3) {
      mydat=subset(dat_list1[[j]],
                   pesticide_by_source==S6_category[i])
      
    } else {
      mydat=subset(dat_list1[[j]],
                   Primary.classification.of.insecticide........一级分类==S6_category[i])
      if (nrow(mydat)<=1) {
        dat_EDF2[(i-1)*length(dat_list1)+j,]=rep(NA,3)
      } else {
        res=rma.mv(yi, vi, 
                   random = list(~ 1 | factor(Insecticide.name),~1|Code),
                   data=mydat)
        dat_EDF2[(i-1)*length(dat_list1)+j,]=c(res$b,res$ci.lb,res$ci.ub)
      }      
    }    
  }
}
dat_EDF2<-as.data.frame(dat_EDF2,row.names<-Rnames)
colnames(dat_EDF2)<-c("Effect size","lb","ub")

### Extended Data Fig 3a-k
dat_EDF2a=dat_EDF2
errbar(1:nrow(dat_EDF2a),dat_EDF2a[,1], dat_EDF2a[,2], dat_EDF2a[,3],
       pch=19,
       ann = FALSE,lwd=4,cex=2,
       col=c(c('blue','green','black','red'),
             rep(c('blue','green','red'),2)),
       errbar.col=c(c('blue','green','black','red'),
                    rep(c('blue','green','red'),2)))
mtext("Effect size",side=4,line=3,cex=1.5,family="sans") #Arial font
abline(h=0,lty=2)
graph2jpg(file='Extended Data Fig 3a-k',height=10,width=8.5,dpi=300)


########## Extended Data Figs. 4a-h ###########
Pesticide_category=c('Pesticide',
                     unique(combined$pesticide_by_target_organisms))
dat_list1=lapply(unique(combined$Taxonomic_group_response),
                 function(x){subset(combined,Taxonomic_group_response==x)})
Experiment=unique(combined$Experiment.type)
dat_EDF3=matrix(NA,length(Experiment)*length(Pesticide_category)*length(dat_list1),3)
for (k in 1:length(Experiment)) {
  for (i in 1:length(Pesticide_category)) {
    for (j in 1:length(dat_list1)){
        mydat=subset(dat_list1[[j]],
                     Experiment.type==Experiment[k]&pesticide_by_target_organisms==Pesticide_category[i])
      if (nrow(mydat)<=1) {
        dat_EDF3[(k-1)*length(Pesticide_category)*length(dat_list1)+(i-1)*length(dat_list1)+j,]=rep(NA,3)
      } else {
        res=rma.mv(yi, vi, 
                   random = list(~ 1 | factor(Insecticide.name),~1|Code),
                   data=mydat)
        dat_EDF3[(k-1)*length(Pesticide_category)*length(dat_list1)+(i-1)*length(dat_list1)+j,]=c(res$b,res$ci.lb,res$ci.ub)
      }
    }
  }
}

dat_EDF3<-as.data.frame(dat_EDF3,row.names<-Rnames)
colnames(dat_EDF3)<-c("#Effect","lb","ub")

### Extended Data Fig 4a-h
dat_EDF3a=dat_EDF3
errbar(1:nrow(dat_EDF3a),dat_EDF3a[,1], dat_EDF3a[,2], dat_EDF3a[,3],
       pch=19,
       ann = FALSE,lwd=4,cex=2,
       col=c(c('blue','green','black','red'),
             rep(c('blue','green','red'),2)),
       errbar.col=c(c('blue','green','black','red'),
                    rep(c('blue','green','red'),2)))
mtext("Effect size",side=4,line=3,cex=1.5,family="sans")
abline(h=0,lty=2)
graph2jpg(file='Extended Data Fig 4a-h',height=10,width=8.5,dpi=300)


########## Extended Data Figs. 5a-h ###########
Pesticide_category=c('Pesticide',
                     unique(combined$pesticide_by_target_organisms))
dat_list1=lapply(unique(combined$Taxonomic_group_response),
                 function(x){subset(combined,Taxonomic_group_response==x)})
Experiment=unique(combined$Major.climatic.zones)
dat_EDF4=matrix(NA,length(Experiment)*length(Pesticide_category)*length(dat_list1),3)
for (k in 1:length(Experiment)) {
  for (i in 1:length(Pesticide_category)) {
    for (j in 1:length(dat_list1)){
        mydat=subset(dat_list1[[j]],
                     Major.climatic.zones==Experiment[k]&pesticide_by_target_organisms==Pesticide_category[i]&Experiment.type=='Field')
      if (nrow(mydat)<=1) {
        dat_EDF4[(k-1)*length(Pesticide_category)*length(dat_list1)+(i-1)*length(dat_list1)+j,]=rep(NA,3)
      } else {
        res=rma.mv(yi, vi, 
                   random = list(~ 1 | factor(Insecticide.name),~1|Code),
                   data=mydat)
        dat_EDF4[(k-1)*length(Pesticide_category)*length(dat_list1)+(i-1)*length(dat_list1)+j,]=c(res$b,res$ci.lb,res$ci.ub)
      }
    }
  }
}
colnames(dat_EDF4)<-c("#Effect","lb","ub")

### Extended Data Fig 5a-h
dat_EDF4a=dat_EDF4
errbar(1:nrow(dat_EDF4a),dat_EDF4a[,1], dat_EDF4a[,2], dat_EDF4a[,3],
       pch=19,
       ann = FALSE,lwd=4,cex=2,
       col=c(c('blue','green','black','red'),
             rep(c('blue','green','red'),2)),
       errbar.col=c(c('blue','green','black','red'),
                    rep(c('blue','green','red'),2)))
mtext("Effect size",side=4,line=3,cex=1.5,family="sans") #Arial font
abline(h=0,lty=2)
graph2jpg(file='Extended Data Fig 5a-h',height=10,width=8.5,dpi=300)


########## Extended Data Figs. 6a-h ###########
Pesticide_category=c('Pesticide',
                     unique(combined$pesticide_by_target_organisms))
dat_list1=lapply(unique(combined$Taxonomic_group_response),
                 function(x){subset(combined,Taxonomic_group_response==x)})
Experiment=unique(combined$Types.of.organism.exposure.to.pesticides)[c(2,1)]
dat_EDF5=matrix(NA,length(Experiment)*length(Pesticide_category)*length(dat_list1),3)
for (k in 1:length(Experiment)) {
  for (i in 1:length(Pesticide_category)) {
    for (j in 1:length(dat_list1)){
        mydat=subset(dat_list1[[j]],
                     Types.of.organism.exposure.to.pesticides==Experiment[k]&pesticide_by_target_organisms==Pesticide_category[i])
      if (nrow(mydat)<=1) {
        dat_EDF5[(k-1)*length(Pesticide_category)*length(dat_list1)+(i-1)*length(dat_list1)+j,]=rep(NA,3)
      } else {
        res=rma.mv(yi, vi, 
                   random = list(~ 1 | factor(Insecticide.name),~1|Code),
                   data=mydat)
        dat_EDF5[(k-1)*length(Pesticide_category)*length(dat_list1)+(i-1)*length(dat_list1)+j,]=c(res$b,res$ci.lb,res$ci.ub)
      }
    }
  }
}
colnames(dat_EDF5)<-c("#Effect","lb","ub")

### Extended Data Fig 6a-h
dat_EDF5a=dat_EDF5
errbar(1:nrow(dat_EDF5a),dat_EDF5a[,1], dat_EDF5a[,2], dat_EDF5a[,3],
       pch=19,
       ann = FALSE,lwd=4,cex=2,
       col=c(c('blue','green','black','red'),
             rep(c('blue','green','red'),2)),
       errbar.col=c(c('blue','green','black','red'),
                    rep(c('blue','green','red'),2)))
mtext("Effect size",side=4,line=3,cex=1.5,family="sans") #Arial font
abline(h=0,lty=2)
graph2jpg(file='Extended Data Fig 6a-h',height=10,width=8.5,dpi=300)


########## Relationship plot ###########
require(nlme)
require(effects)

linear_plot=function(dat,columname,roname,number){
  
  if (columname %in% c('Animal','Invertebrate animal','Vertebrate animal')) {
    mylabel=c('growth','reproduction','behavior','biomarker')
  } else {
    mylabel=c('growth','reproduction','biomarker')
  }
  for(i in 1:length(mylabel)) {
    subdat=subdat[!is.na(subdat$xi),]
    if (nrow(subdat)>2) {
      x <- subdat[,'xi']
      y <- subdat[,'yi']
      mydata=data.frame(x=x,y=y,Code=subdat[,'Code'],
                        pesticide_by_target_organisms=subdat[,'pesticide_by_target_organisms'],
                        Insecticide.name=subdat[,'Insecticide.name'],vi=subdat[,'vi'])
      dat_lis3[[i]]=subset(mydata)
      # for gls and Effect, mydata needs to be in global environment (renamed as "mydata2"):    
      newx<-log2(seq(min(2^x),max(2^x),length=500))
      assign("mydata2", mydata, envir=.GlobalEnv) 
      fit<-try(lme(y~x,
                   random =list(~ 1|Insecticide.name,~ 1|Code),
                   weights = varFixed(~ vi),data=mydata2,control = lmeControl(sigma = 1)))
      P[i]=signif(summary(fit)$tTable[2,5],4)
      e1=Effect("x",fit,xlevels=list(x=newx))
      remove("mydata2", envir=.GlobalEnv) # remove from global environment
      pred.df=as.data.frame(e1)
      pred.c=pred.df[,c(2,4,5)]
      dat_lis1[[i]]=fit
      #colnames(dat_lis1[[i]])=c('a','b')
      #dat_lis1[[i]]=subset(dat_lis1[[i]],b>-22&b<14)
      dat_lis2[[i]]=as.data.frame(cbind(rev(newx),newx,rev(pred.c[ ,3]),pred.c[ ,2]))
      #colnames(dat_lis2[[i]])=c('a','b','c','d')
      #dat_lis2[[i]]=subset(dat_lis2[[i]],c>-20&d>-22&c<20&d<22)
    }
  }
  
  #poly_index=do.call('rbind',dat_lis2)
  #x_lim=range(c(range(poly_index[,1]),range(poly_index[,2])))
  #y_lim=range(c(range(poly_index[,3]),range(poly_index[,4])))
  
  y_lim=c(-20,20)
  
  plot(x,y,type = "n",
       xlab=paste0('Log2 (',roname,' application rate added over the control)'),
       las=1,cex.axis=1.2,cex.lab=1.2,
       ylab=paste0(columname,' effect size'),family="sans")
  abline(h=0,lty=2,lwd=1.5)
  box(which = "plot", lwd = 2)
  
  
  col1=c(rgb(180,82,205,maxColorValue = 255),
         rgb(0,0,0,maxColorValue = 255))
  col2=c(rgb(180,82,205,50,maxColorValue = 255),
         rgb(0,0,0,50,maxColorValue = 255))
  if (length(mylabel)==3) {
    col1=col1[c(1:2,4)]
    col2=col2[c(1:2,4)]
  }
  
  for (i in 1:length(mylabel)) {
    points(dat_lis3[[i]]$x,dat_lis3[[i]]$y)
    if (!is.na(P[i])&P[i]<0.05) {
      polygon(c(dat_lis2[[i]][,1], dat_lis2[[i]][,2]), c(dat_lis2[[i]][,3], dat_lis2[[i]][,4]),
              col = col2[i], border = NA)
      #lines(dat_lis1[[i]][,1],dat_lis1[[i]][,2],lwd=2.5,col=col1[i])
      #abline(summary(dat_lis1[[i]])$tTable[1,1],summary(dat_lis1[[i]])$tTable[2,1],lwd=2.5,col=col1[i])
      a=summary(dat_lis1[[i]])$tTable[1,1]
      b=summary(dat_lis1[[i]])$tTable[2,1]
      segments(x0=0,y0=a,x1=max(dat_lis2[[i]][,2]),
               lwd=2.5,col=col1[i])
    }
  }
  
  if (number=='a') {
    legend('topright',inset = c(0.02,0),bty="n",
           pch = 21,col = col1,
           ncol = 1, y.intersp =0.9)
  }
}
  
  if (columname %in% c('Animal','Invertebrate animal','Vertebrate animal')) {
    mylabel=c('growth','reproduction','behavior','biomarker')
  } else {
    mylabel=c('growth','reproduction','biomarker')
  }
  for(i in 1:length(mylabel)) {
    subdat=subset(dat,Taxonomic_group_response==grep(mylabel[i],unique(dat$Taxonomic_group_response),value = T))
    subdat=subdat[!is.na(subdat$xi),]
    if (nrow(subdat)>2) {
      x <- subdat[,'xi']
      y <- subdat[,'yi']
      mydata=data.frame(x=x,y=y,Code=subdat[,'Code'],
                        pesticide_by_target_organisms=subdat[,'pesticide_by_target_organisms'],
                        Insecticide.name=subdat[,'Insecticide.name'],vi=subdat[,'vi'])
      dat_lis3[[i]]=subset(mydata)
      # for gls and Effect, mydata needs to be in global environment (renamed as "mydata2"):    
      newx<-log2(seq(min(2^x),max(2^x),length=500))
      assign("mydata2", mydata, envir=.GlobalEnv) 
      fit<-try(lme(y~x,
                   random =list(~ 1|Insecticide.name,~ 1|Code),
                   weights = varFixed(~ vi),data=mydata2,control = lmeControl(sigma = 1)))
      P[i]=signif(summary(fit)$tTable[2,5],4)
      e1=Effect("x",fit,xlevels=list(x=newx))
      remove("mydata2", envir=.GlobalEnv) # remove from global environment
      pred.df=as.data.frame(e1)
      pred.c=pred.df[,c(2,4,5)]
      dat_lis1[[i]]=fit
      #colnames(dat_lis1[[i]])=c('a','b')
      #dat_lis1[[i]]=subset(dat_lis1[[i]],b>-22&b<14)
      dat_lis2[[i]]=as.data.frame(cbind(rev(newx),newx,rev(pred.c[ ,3]),pred.c[ ,2]))
      #colnames(dat_lis2[[i]])=c('a','b','c','d')
      #dat_lis2[[i]]=subset(dat_lis2[[i]],c>-20&d>-22&c<20&d<22)
    }
  }
  
  
  #poly_index=do.call('rbind',dat_lis2)
  #x_lim=range(c(range(poly_index[,1]),range(poly_index[,2])))
  #y_lim=range(c(range(poly_index[,3]),range(poly_index[,4])))
  
  maxindex=do.call('rbind',dat_lis2)
  
  plot(x,y,type = "n",
       xlab=paste0('Log2 (',roname,' application rate added over the control)'),
       las=1,cex.axis=1.2,cex.lab=1.2,mgp=c(2.5,1,0),
       ylab=paste0(columname,' effect size'),family="sans")
  abline(h=0,lty=2,lwd=1.5)
  box(which = "plot", lwd = 2)
  
  
  col1=c(rgb(180,82,205,maxColorValue = 255),rgb(255,0,0,maxColorValue = 255))
  col2=c(rgb(180,82,205,50,maxColorValue = 255),rgb(255,0,0,50,maxColorValue = 255))
  if (length(mylabel)==3) {
    col1=col1[c(1:2,4)]
    col2=col2[c(1:2,4)]
  }
  
  for (i in 1:length(mylabel)) {
    points(dat_lis3[[i]]$x,dat_lis3[[i]]$y,col = col1[i])
    if (!is.na(P[i])&P[i]<0.05) {
      polygon(c(dat_lis2[[i]][,1], dat_lis2[[i]][,2]), c(dat_lis2[[i]][,3], dat_lis2[[i]][,4]),
              col = col2[i], border = NA)
      #lines(dat_lis1[[i]][,1],dat_lis1[[i]][,2],lwd=2.5,col=col1[i])
      #abline(summary(dat_lis1[[i]])$tTable[1,1],summary(dat_lis1[[i]])$tTable[2,1],lwd=2.5,col=col1[i])
      a=summary(dat_lis1[[i]])$tTable[1,1]
      b=summary(dat_lis1[[i]])$tTable[2,1]
      segments(x0=0,y0=a,x1=max(dat_lis2[[i]][,2]),y1=a+max(dat_lis2[[i]][,2])*b,
               lwd=2.5,col=col1[i])
    }
  }
  
  if (number=='a') {
      legend('topright',inset = c(0.02,0),bty="n",
             stringr::str_to_title(mylabel),
             pch = 21,col = col1,
             ncol = 1, y.intersp =0.9)
  }
}


###### Figure 3 ######

### Figure 3a
dat=combined[combined$Taxonomic_group=='animals',]
columname='Animal'
roname='pesticide'
number='a'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Figure 3a',height=6,width=8,dpi=300)

### Figure 3b
dat=combined[combined$Taxonomic_group=='plants',]
columname='Plant'
roname='pesticide'
number='b'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Figure 3b',height=6,width=8,dpi=300)

### Figure 3c
dat=combined[combined$Taxonomic_group=='microorganisms',]
columname='Microorganism'
roname='pesticide'
number='c'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Figure 3c',height=6,width=8,dpi=300)

### Figure 3d
dat=subset(combined,Taxonomic_group=='animals'&pesticide_by_target_organisms=='insecticides')
columname='Animal'
roname='insecticide'
number='d'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Figure 3d',height=6,width=8,dpi=300)

### Figure 3e
dat=subset(combined,Taxonomic_group=='plants'&pesticide_by_target_organisms=='insecticides')
columname='Plant'
roname='insecticide'
number='e'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Figure 3e',height=6,width=8,dpi=300)

### Figure 3f
dat=subset(combined,Taxonomic_group=='microorganisms'&pesticide_by_target_organisms=='insecticides')
columname='Microorganism'
roname='insecticide'
number='f'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Figure 3f',height=6,width=8,dpi=300)

### Figure 3g
dat=subset(combined,Taxonomic_group=='animals'&pesticide_by_target_organisms=='fungicides')
columname='Animal'
roname='fungicide'
number='g'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Figure 3g',height=6,width=8,dpi=300)

### Figure 3h
dat=subset(combined,Taxonomic_group=='plants'&pesticide_by_target_organisms=='fungicides')
columname='Plant'
roname='fungicide'
number='h'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Figure 3h',height=6,width=8,dpi=300)

### Figure 3i
dat=subset(combined,Taxonomic_group=='microorganisms'&pesticide_by_target_organisms=='fungicides')
columname='Microorganism'
roname='fungicide'
number='i'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Figure 3i',height=6,width=8,dpi=300)

### Figure 3j
dat=subset(combined,Taxonomic_group=='animals'&pesticide_by_target_organisms=='herbicides')
columname='Animal'
roname='herbicides'
number='j'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Figure 3j',height=6,width=8,dpi=300)

### Figure 3k
dat=subset(combined,Taxonomic_group=='plants'&pesticide_by_target_organisms=='herbicides')
columname='Plant'
roname='herbicides'
number='k'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Figure 3k',height=6,width=8,dpi=300)

### Figure 3l
dat=subset(combined,Taxonomic_group=='microorganisms'&pesticide_by_target_organisms=='herbicides')
columname='Microorganism'
roname='herbicides'
number='l'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Figure 3l',height=6,width=8,dpi=300)

####### Supplementary Fig 1 ######
dat_F6=subset(combined,Experiment.type=='Laboratory')

### Supplementary Fig 1a
dat=dat_F6[dat_F6$Taxonomic_group=='animals',]
columname='Animal'
roname='pesticide'
number='a'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Supplementary Fig 1a',height=6,width=8,dpi=300)

### Supplementary Fig 1b
dat=dat_F6[dat_F6$Taxonomic_group=='plants',]
columname='Plant'
roname='pesticide'
number='b'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Supplementary Fig 1b',height=6,width=8,dpi=300)

### Supplementary Fig 1c
dat=dat_F6[dat_F6$Taxonomic_group=='microorganisms',]
columname='Microorganism'
roname='pesticide'
number='c'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Supplementary Fig 1c',height=6,width=8,dpi=300)

### Supplementary Fig 1d
dat=subset(dat_F6,Taxonomic_group=='animals'&pesticide_by_target_organisms=='insecticides')
columname='Animal'
roname='insecticide'
number='d'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Supplementary Fig 1d',height=6,width=8,dpi=300)

### Supplementary Fig 1e
dat=subset(dat_F6,Taxonomic_group=='plants'&pesticide_by_target_organisms=='insecticides')
columname='Plant'
roname='insecticide'
number='e'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Supplementary Fig 1e',height=6,width=8,dpi=300)

### Supplementary Fig 1f
dat=subset(dat_F6,Taxonomic_group=='microorganisms'&pesticide_by_target_organisms=='insecticides')
columname='Microorganism'
roname='insecticide'
number='f'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Supplementary Fig 1f',height=6,width=8,dpi=300)

### Supplementary Fig 1g
dat=subset(dat_F6,Taxonomic_group=='animals'&pesticide_by_target_organisms=='fungicides')
columname='Animal'
roname='fungicide'
number='g'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Supplementary Fig 1g',height=6,width=8,dpi=300)

### Supplementary Fig 1h
dat=subset(dat_F6,Taxonomic_group=='plants'&pesticide_by_target_organisms=='fungicides')
columname='Plant'
roname='fungicide'
number='h'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Supplementary Fig 1h',height=6,width=8,dpi=300)

### Supplementary Fig 1i
dat=subset(dat_F6,Taxonomic_group=='microorganisms'&pesticide_by_target_organisms=='fungicides')
columname='Microorganism'
roname='fungicide'
number='i'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Supplementary Fig 1i',height=6,width=8,dpi=300)

### Supplementary Fig 1j
dat=subset(dat_F6,Taxonomic_group=='animals'&pesticide_by_target_organisms=='herbicides'&Taxonomic_group_response!='animal reproduction')
columname='Animal'
roname='herbicides'
number='j'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Supplementary Fig 1j',height=6,width=8,dpi=300)

### Supplementary Fig 1k
dat=subset(dat_F6,Taxonomic_group=='plants'&pesticide_by_target_organisms=='herbicides')
columname='Plant'
roname='herbicides'
number='k'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Supplementary Fig 1k',height=6,width=8,dpi=300)

### Supplementary Fig 1l
dat=subset(dat_F6,Taxonomic_group=='microorganisms'&pesticide_by_target_organisms=='herbicides')
columname='Microorganism'
roname='herbicides'
number='l'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Supplementary Fig 1l',height=6,width=8,dpi=300)


####### Supplementary Fig 2 ######
dat_F7=subset(combined,Experiment.type=='Field')

### Supplementary Fig 2a
dat=dat_F7[dat_F7$Taxonomic_group=='animals',]
columname='Animal'
roname='pesticide'
number='a'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Supplementary Fig 2a',height=6,width=8,dpi=300)

### Supplementary Fig 2b
dat=dat_F7[dat_F7$Taxonomic_group=='plants',]
columname='Plant'
roname='pesticide'
number='b'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Supplementary Fig 2b',height=6,width=8,dpi=300)

### Supplementary Fig 2c
dat=dat_F7[dat_F7$Taxonomic_group=='microorganisms',]
columname='Microorganism'
roname='pesticide'
number='c'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Supplementary Fig 2c',height=6,width=8,dpi=300)

### Supplementary Fig 2d
dat=subset(dat_F7,Taxonomic_group=='animals'&pesticide_by_target_organisms=='insecticides')
columname='Animal'
roname='insecticide'
number='d'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Supplementary Fig 2d',height=6,width=8,dpi=300)

### Supplementary Fig 2e
dat=subset(dat_F7,Taxonomic_group=='plants'&pesticide_by_target_organisms=='insecticides')
columname='Plant'
roname='insecticide'
number='e'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Supplementary Fig 2e',height=6,width=8,dpi=300)

### Supplementary Fig 2f
dat=subset(dat_F7,Taxonomic_group=='microorganisms'&pesticide_by_target_organisms=='insecticides')
columname='Microorganism'
roname='insecticide'
number='f'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Supplementary Fig 2f',height=6,width=8,dpi=300)

### Supplementary Fig 2g
dat=subset(dat_F7,Taxonomic_group=='animals'&pesticide_by_target_organisms=='fungicides'&Taxonomic_group_response!='animal biomarker')
columname='Animal'
roname='fungicide'
number='g'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Supplementary Fig 2g',height=6,width=8,dpi=300)

### Supplementary Fig 2h
dat=subset(dat_F7,Taxonomic_group=='plants'&pesticide_by_target_organisms=='fungicides')
columname='Plant'
roname='fungicide'
number='h'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Supplementary Fig 2h',height=6,width=8,dpi=300)

### Supplementary Fig 2i
dat=subset(dat_F7,Taxonomic_group=='microorganisms'&pesticide_by_target_organisms=='fungicides'&Taxonomic_group_response!='microorganism growth')
columname='Microorganism'
roname='fungicide'
number='i'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Supplementary Fig 2i',height=6,width=8,dpi=300)

### Supplementary Fig 2j
dat=subset(dat_F7,Taxonomic_group=='animals'&pesticide_by_target_organisms=='herbicides'&Taxonomic_group_response!='animal behavior')
columname='Animal'
roname='herbicides'
number='j'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Supplementary Fig 2j',height=6,width=8,dpi=300)

### Supplementary Fig 2k
dat=subset(dat_F7,Taxonomic_group=='plants'&pesticide_by_target_organisms=='herbicides')
columname='Plant'
roname='herbicides'
number='k'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Supplementary Fig 2k',height=6,width=8,dpi=300)

### Supplementary Fig 2l
dat=subset(dat_F7,Taxonomic_group=='microorganisms'&pesticide_by_target_organisms=='herbicides')
columname='Microorganism'
roname='herbicides'
number='l'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Supplementary Fig 2l',height=6,width=8,dpi=300)


####### Supplementary Fig 3 ######
dat_F8=subset(combined,Experiment.type=='Field'&Major.climatic.zones=='Temperate')

### Supplementary Fig 3a-l
dat=dat_F8[dat_F8$Taxonomic_group=='animals',]
columname='Animal'
roname='pesticide'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Supplementary Fig 3a-l',height=6,width=8,dpi=300)


####### Supplementary Fig 4 ######
dat_F9=subset(combined,Experiment.type=='Field'&Major.climatic.zones=='Tropical')

### Supplementary Fig 4a-l
dat=dat_F9[dat_F9$Taxonomic_group=='animals',]
columname='Animal'
roname='pesticide'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Supplementary Fig 4a-l',height=6,width=8,dpi=300)


####### Supplementary Fig 5 ######
dat_F10=subset(combined,Types.of.organism.exposure.to.pesticides=='Aquatic')

### Supplementary Fig 5a-l
dat=dat_F10[dat_F10$Taxonomic_group=='animals',]
columname='Animal'
roname='pesticide'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Supplementary Fig 5a-l',height=6,width=8,dpi=300)


####### Supplementary Fig 6 ######
dat_F11=subset(combined,Types.of.organism.exposure.to.pesticides=='Terrestrial')

### Supplementary Fig 6a-l
dat=dat_F11[dat_F11$Taxonomic_group=='animals',]
columname='Animal'
roname='pesticide'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Supplementary Fig 6a-l',height=6,width=8,dpi=300)


####### Supplementary Fig 7 ######
dat_F12=subset(combined,pesticide_by_source=='chemical')

### Supplementary Fig 7a-l
dat=dat_F12[dat_F12$Taxonomic_group=='animals',]
columname='Animal'
roname='chemical synthetic pesticide'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Supplementary Fig 7a-l',height=6,width=8,dpi=300)


####### Supplementary Fig 8 ######
dat_F13=subset(combined,pesticide_by_source=='mineral-based')

### Supplementary Fig 8a-l
dat=dat_F13[dat_F13$Taxonomic_group=='animals',]
columname='Animal'
roname='mineral-based pesticide'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Supplementary Fig 8a-l',height=6,width=8,dpi=300)


####### Supplementary Fig 9 ######
dat_F14=subset(combined,pesticide_by_source=='biogenic')

### Supplementary Fig 9a-l
dat=dat_F14[dat_F14$Taxonomic_group=='animals',]
columname='Animal'
roname='biogenic pesticide'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Supplementary Fig 9a-l',height=6,width=8,dpi=300)


####### Supplementary Fig 10 ######
dat_F15_I=subset(combined,Primary.classification.of.animal......一级分类=='Invertebrate')
dat_F15_V=subset(combined,Primary.classification.of.animal......一级分类=='Vertebrate')

### Supplementary Fig 10a-l
dat=dat_F15_I
columname='Invertebrate animal'
roname='Pesticide'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Supplementary Fig 10a-l',height=6,width=8,dpi=300)


####### Supplementary Fig 11 ######
dat_F16_I=subset(combined,Primary.classification.of.animal......一级分类=='Invertebrate'&pesticide_by_target_organisms=='insecticides')
dat_F16_V=subset(combined,Primary.classification.of.animal......一级分类=='Vertebrate'&pesticide_by_target_organisms=='insecticides')

### Supplementary Fig 11a-l
dat=dat_F16_I
columname='Invertebrate animal'
roname='insecticide'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Supplementary Fig 11a-l',height=6,width=8,dpi=300)


####### Supplementary Fig 12 ######
dat_F17_I=subset(combined,Primary.classification.of.animal......一级分类=='Invertebrate'&pesticide_by_target_organisms=='fungicides')
dat_F17_V=subset(combined,Primary.classification.of.animal......一级分类=='Vertebrate'&pesticide_by_target_organisms=='fungicides')

### Supplementary Fig 12a-h
dat=dat_F17_I
columname='Invertebrate animal'
roname='fungicide'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Supplementary Fig 12a-l',height=6,width=8,dpi=300)


####### Supplementary Fig 13 ######
dat_F18_I=subset(combined,Primary.classification.of.animal......一级分类=='Invertebrate'&pesticide_by_target_organisms=='herbicides')
dat_F18_V=subset(combined,Primary.classification.of.animal......一级分类=='Vertebrate'&pesticide_by_target_organisms=='herbicides')

### Supplementary Fig 13a-h
dat=dat_F18_I
columname='Invertebrate animal'
roname='herbicide'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Supplementary Fig 13a-h',height=6,width=8,dpi=300)


####### Supplementary Fig 14 ######
dat_F19_1=subset(combined,Primary.classification.of.animal......一级分类=='Seed plant')
dat_F19_2=subset(combined,Primary.classification.of.animal......一级分类=='Spore-producing  plant')

### Supplementary Fig 14a-h
dat=dat_F19_1
columname='Seed plant'
roname='Pesticide'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Supplementary Fig 14a-h',height=6,width=8,dpi=300)


####### Supplementary Fig 15 ######
dat_F20_1=subset(combined,Primary.classification.of.animal......一级分类=='Seed plant' & pesticide_by_target_organisms=='insecticides')
dat_F20_2=subset(combined,Primary.classification.of.animal......一级分类=='Spore-producing  plant' & pesticide_by_target_organisms=='insecticides')

### Supplementary Fig 15a-h
dat=dat_F20_1
columname='Seed plant'
roname='insecticide'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Supplementary Fig 15a-h',height=6,width=8,dpi=300)


####### Supplementary Fig 16 ######
dat_F21_1=subset(combined,Primary.classification.of.animal......一级分类=='Seed plant' & pesticide_by_target_organisms=='fungicides')
dat_F21_2=subset(combined,Primary.classification.of.animal......一级分类=='Spore-producing  plant' & pesticide_by_target_organisms=='fungicides')

### Supplementary Fig 16a-h
dat=dat_F21_1
columname='Seed plant'
roname='fungicide'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Supplementary Fig 16a-h',height=6,width=8,dpi=300)


####### Supplementary Fig 17 ######
dat_F22_1=subset(combined,Primary.classification.of.animal......一级分类=='Seed plant' & pesticide_by_target_organisms=='herbicides')
dat_F22_2=subset(combined,Primary.classification.of.animal......一级分类=='Spore-producing  plant' & pesticide_by_target_organisms=='herbicides')

### Supplementary Fig 17a-h
dat=dat_F22_1
columname='Seed plant'
roname='herbicide'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Supplementary Fig 17a-h',height=6,width=8,dpi=300)


####### Supplementary Fig 18 ######
dat_F23_1=subset(combined,Primary.classification.of.animal......一级分类=='Bacteria')
dat_F23_2=subset(combined,Primary.classification.of.animal......一级分类=='Fungus')

### Supplementary Fig 18a-h
dat=dat_F23_1
columname='Bacteria'
roname='Pesticide'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Supplementary Fig 18a-h',height=6,width=8,dpi=300)


####### Supplementary Fig 19 ######
dat_F24_1=subset(combined,Primary.classification.of.animal......一级分类=='Bacteria' & pesticide_by_target_organisms=='insecticides')
dat_F24_2=subset(combined,Primary.classification.of.animal......一级分类=='Fungus' & pesticide_by_target_organisms=='insecticides')

### Supplementary Fig 19a-h
dat=dat_F24_1
columname='Bacteria'
roname='insecticide'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Supplementary Fig 19a-h',height=6,width=8,dpi=300)


####### Supplementary Fig 20 ######
dat_F25_1=subset(combined,Primary.classification.of.animal......一级分类=='Bacteria' & pesticide_by_target_organisms=='fungicides')
dat_F25_2=subset(combined,Primary.classification.of.animal......一级分类=='Fungus' & pesticide_by_target_organisms=='fungicides')

### Supplementary Fig 20a-h
dat=dat_F25_1
columname='Bacteria'
roname='fungicide'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Supplementary Fig 20a-h',height=6,width=8,dpi=300)


####### Supplementary Fig 21 ######
dat_F26_1=subset(combined,Primary.classification.of.animal......一级分类=='Bacteria' & pesticide_by_target_organisms=='herbicides')
dat_F26_2=subset(combined,Primary.classification.of.animal......一级分类=='Fungus' & pesticide_by_target_organisms=='herbicides')

### Supplementary Fig 21a-h
dat=dat_F26_1
columname='Bacteria'
roname='herbicide'
linear_plot(dat,columname,roname,number)
graph2jpg(file='Supplementary Fig 21a-h',height=6,width=8,dpi=300)
