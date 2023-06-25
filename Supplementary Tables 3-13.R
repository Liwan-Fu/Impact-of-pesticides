library(metafor)
library(xlsx)
library(optimParallel)
setwd('...')
combined=read.csv('Meta combined data.csv',header=T)

# combined=na.omit(combined)
#  combined=combined[sample(nrow(combined))[1:1000],]

########## Supplementary Table 3 ###########
Pesticide_category=c('Pesticide',
                     unique(combined$pesticide_by_target_organisms))
dat_list1=lapply(unique(combined$Taxonomic_group_response),
                 function(x){subset(combined,Taxonomic_group_response==x)})
S4=matrix(NA,length(Pesticide_category)*length(dat_list1),8)
for (i in 1:length(Pesticide_category)) {
  for (j in 1:length(dat_list1)){
      mydat=subset(dat_list1[[j]],
                   pesticide_by_target_organisms==Pesticide_category[i])
      res=rma.mv(yi, vi, 
                 random = list(~ 1 | factor(Insecticide.name),~1|Code),
                 data=mydat)
    S4[(i-1)*length(dat_list1)+j,]=c(nrow(mydat),length(unique(mydat$Code)),
                                     res$b,res$zval,nrow(mydat)-1,
                                     res$pval,res$ci.lb,res$ci.ub)
  }
}
colnames(S4)<-c("#Obs","#Studies","Effect size","T-value","Df",
                "P-value","CL.lb","CL.ub")
write.xlsx(S4,"Supplementary Table 3.xls")

########## Supplementary Table 4 ###########
dat_list2=c(lapply(unique(combined$Primary.classification.of.animal......一级分类),
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
S5=matrix(NA,length(Pesticide_category)*length(dat_list3),8)
for (i in 1:length(Pesticide_category)) {
  for (j in 1:length(dat_list3)){
      mydat=subset(dat_list3[[j]],
                   pesticide_by_target_organisms==Pesticide_category[i])
    if (nrow(mydat)<=1) {
      S5[(i-1)*length(dat_list3)+j,]=c(nrow(mydat),length(unique(mydat$Code)),
                                       rep(NA,6))
    } else {
      res=rma.mv(yi, vi, 
                 random = list(~ 1 | factor(Insecticide.name),~1|Code),
                 data=mydat)
      S5[(i-1)*length(dat_list3)+j,]=c(nrow(mydat),length(unique(mydat$Code)),
                                       res$b,res$zval,nrow(mydat)-1,
                                       res$pval,res$ci.lb,res$ci.ub)
    }
  }
}
colnames(S5)<-c("#Obs","#Studies","Effect size","T-value","Df",
                "P-value","CL.lb","CL.ub")
write.xlsx(S5,"Supplementary Table 4.xls")

########## Supplementary Table 5 ###########
Pesticide_category=c('Pesticide',
                     unique(combined$pesticide_by_target_organisms))
dat_list1=lapply(unique(combined$Taxonomic_group_response),
                   function(x){subset(combined,Taxonomic_group_response==x)})
S6_category=c(unique(combined$pesticide_by_source),
              unique(combined$Primary.classification.of.insecticide........一级分类))[c(1,3,2,4,7,10,6,8,5,9,11)]
S6=matrix(NA,length(S6_category)*length(dat_list1),8)
for (i in 1:length(S6_category)) {
  for (j in 1:length(dat_list1)){
      mydat=subset(dat_list1[[j]],
                   Primary.classification.of.insecticide........一级分类==S6_category[i])
      if (nrow(mydat)<=1) {
        S6[(i-1)*length(dat_list1)+j,]=c(nrow(mydat),length(unique(mydat$Code)),
                                         rep(NA,6))
      } else {
        res=rma.mv(yi, vi, 
                   random = list(~ 1 | factor(Insecticide.name),~1|Code),
                   data=mydat)
        S6[(i-1)*length(dat_list1)+j,]=c(nrow(mydat),length(unique(mydat$Code)),
                                         res$b,res$zval,nrow(mydat)-1,
                                         res$pval,res$ci.lb,res$ci.ub)
      }

    }

  }
colnames(S6)<-c("#Obs","#Studies","Effect size","T-value","Df",
                "P-value","CL.lb","CL.ub")
write.xlsx(S6,"Supplementary Table 5.xls")

########## Supplementary Table 6 ########### 
Pesticide_category=c('Pesticide',
                     unique(combined$pesticide_by_target_organisms))
dat_list1=lapply(unique(combined$Taxonomic_group_response),
                   function(x){subset(combined,Taxonomic_group_response==x)})
Experiment=unique(combined$Experiment.type)
S7=matrix(NA,length(Experiment)*length(Pesticide_category)*length(dat_list1),8)
for (k in 1:length(Experiment)) {
  for (i in 1:length(Pesticide_category)) {
    for (j in 1:length(dat_list1)){
        mydat=subset(dat_list1[[j]],
                     Experiment.type==Experiment[k]&pesticide_by_target_organisms==Pesticide_category[i])
      if (nrow(mydat)<=1) {
        S7[(k-1)*40+(i-1)*10+j,]=c(nrow(mydat),length(unique(mydat$Code)),
                                   rep(NA,6))
      } else {
        res=rma.mv(yi, vi, 
                   random = list(~ 1 | factor(Insecticide.name),~1|Code),
                   data=mydat)
        S7[(k-1)*40+(i-1)*10+j,]=c(nrow(mydat),length(unique(mydat$Code)),
                                   res$b,res$zval,nrow(mydat)-1,res$pval,res$ci.lb,res$ci.ub)
      }
    }
  }
}

colnames(S7)<-c("#Obs","#Studies","Effect size","T-value","Df",
                "P-value","CL.lb","CL.ub")
write.xlsx(S7,"Supplementary Table 6.xls")

########## Supplementary Table 7 ###########
Pesticide_category=c('Pesticide',
                     unique(combined$pesticide_by_target_organisms))
dat_list1=lapply(unique(combined$Taxonomic_group_response),
                   function(x){subset(combined,Taxonomic_group_response==x)})
Experiment=unique(combined$Major.climatic.zones)
S8=matrix(NA,length(Experiment)*length(Pesticide_category)*length(dat_list1),8)
for (k in 1:length(Experiment)) {
  for (i in 1:length(Pesticide_category)) {
    for (j in 1:length(dat_list1)){
        mydat=subset(dat_list1[[j]],
                     Major.climatic.zones==Experiment[k]&pesticide_by_target_organisms==Pesticide_category[i]&Experiment.type=='Field')
      if (nrow(mydat)<=1) {
        S8[(k-1)*40+(i-1)*10+j,]=c(nrow(mydat),length(unique(mydat$Code)),
                                   rep(NA,6))
      } else {
        res=rma.mv(yi, vi, 
                   random = list(~ 1 | factor(Insecticide.name),~1|Code),
                   data=mydat)
        S8[(k-1)*40+(i-1)*10+j,]=c(nrow(mydat),length(unique(mydat$Code)),
                                   res$b,res$zval,nrow(mydat)-1,res$pval,res$ci.lb,res$ci.ub)
      }
    }
  }
}

colnames(S8)<-c("#Obs","#Studies","Effect size","T-value","Df",
                "P-value","CL.lb","CL.ub")
write.xlsx(S8,"Supplementary Table 7.xls")

########## Supplementary Table 8 ###########
Pesticide_category=c('Pesticide',
                     unique(combined$pesticide_by_target_organisms))
dat_list1=lapply(unique(combined$Taxonomic_group_response),
                   function(x){subset(combined,Taxonomic_group_response==x)})
Experiment=unique(combined$Types.of.organism.exposure.to.pesticides)[c(2,1)]
S9=matrix(NA,length(Experiment)*length(Pesticide_category)*length(dat_list1),8)
for (k in 1:length(Experiment)) {
  for (i in 1:length(Pesticide_category)) {
    for (j in 1:length(dat_list1)){
        mydat=subset(dat_list1[[j]],
                     Types.of.organism.exposure.to.pesticides==Experiment[k]&pesticide_by_target_organisms==Pesticide_category[i])
      if (nrow(mydat)<=1) {
        S9[(k-1)*40+(i-1)*10+j,]=c(nrow(mydat),length(unique(mydat$Code)),
                                   rep(NA,6))
      } else {
        res=rma.mv(yi, vi, 
                   random = list(~ 1 | factor(Insecticide.name),~1|Code),
                   data=mydat)
        S9[(k-1)*40+(i-1)*10+j,]=c(nrow(mydat),length(unique(mydat$Code)),
                                   res$b,res$zval,nrow(mydat)-1,res$pval,res$ci.lb,res$ci.ub)
      }
    }
  }
}

colnames(S9)<-c("#Obs","#Studies","Effect size","T-value","Df",
                "P-value","CL.lb","CL.ub")
write.xlsx(S9,"Supplementary Table 8.xls")

########## Supplementary Table 9 ###########
animal_sed=sort(read.xlsx('inter.xls',sheetName = 'animal-1',header = T)[,1])
plant_sed=sort(read.xlsx('inter.xls',sheetName = 'plant-1',header = T)[,1])
microo_sed=sort(read.xlsx('inter.xls',sheetName = 'microo-1',header = T)[,1])
dat_list10=c(lapply(unique(combined$Taxonomic_group),
                   function(x){subset(combined,Taxonomic_group==x)}),
             lapply(unique(combined$Taxonomic_group_response),
                   function(x){subset(combined,Taxonomic_group_response==x)}))
mycat=unique(combined$Taxonomic_group_response)
for(i in 1:length(mycat)) {
      dat_list10=c(dat_list10,lapply(plant_sed,
                                     function(x){subset(combined,Taxonomic_group_response==mycat[i]&Secondary.classification.of.animal.........二级分类==x)}))
      dat_list10=c(dat_list10,lapply(microo_sed,
                                     function(x){subset(combined,Taxonomic_group_response==mycat[i]&Secondary.classification.of.animal.........二级分类==x)}))
}

for (i in 1:length(mycat)) {
  if (i<=4) {
    Nontarget_category=c(Nontarget_category,paste(plant_sed,mycat[i],sep = ' '))
    Nontarget_category=c(Nontarget_category,paste(microo_sed,mycat[i],sep = ' '))
}
S10=list()
for (k in 1:length(insecti)) {
  insecti_cate=sort(unique(subset(combined,
                                  pesticide_by_target_organisms=='insecticides'&pesticide_by_source==insecti[k])$Secondary.classification.of.insecticide........二级分类))
  S10[[k]]=matrix(NA,length(insecti_cate)*length(dat_list10),11)
  for (i in 1:length(insecti_cate)) {
    for (j in 1:length(dat_list10)){
      mydat=subset(dat_list10[[j]],
                   pesticide_by_target_organisms=='insecticides'&pesticide_by_source==insecti[k]&Secondary.classification.of.insecticide........二级分类==insecti_cate[i])
      if (nrow(mydat)<=1) {
        S10[[k]][(i-1)*length(dat_list10)+j,]=c(paste(insecti[k],'insecticides',sep = '_'),
                                                insecti_cate[i],Nontarget_category[j],
                                                nrow(mydat),length(unique(mydat$Code)),rep(NA,6))
      } else {
        res=rma.mv(yi, vi, 
                   random = list(~ 1 | factor(Insecticide.name),~1|Code),
                   data=mydat)
        S10[[k]][(i-1)*length(dat_list10)+j,]=c(paste(insecti[k],'insecticides',sep = '_'),
                                                insecti_cate[i],Nontarget_category[j],
                                                nrow(mydat),length(unique(mydat$Code)),
                                                res$b,res$zval,nrow(mydat)-1,
                                                res$pval,res$ci.lb,res$ci.ub)
      }
    }
  }
}
S10=do.call('rbind',S10)
colnames(S10)<-c('Insecticide','Insecticide category','Non-target organism category',
                "#Obs","#Studies","Effect size","T-value","Df",
                "P-value","CL.lb","CL.ub")
write.xlsx(S10,"Supplementary Table 9.xls",row.names = F)

########## Supplementary Table 10 ###########
animal_sed=sort(read.xlsx('inter.xls',sheetName = 'animal-1',header = T)[,1])
plant_sed=sort(read.xlsx('inter.xls',sheetName = 'plant-1',header = T)[,1])
microo_sed=sort(read.xlsx('inter.xls',sheetName = 'microo-1',header = T)[,1])
dat_list10=c(lapply(unique(combined$Taxonomic_group),
                    function(x){subset(combined,Taxonomic_group==x)}),
             lapply(unique(combined$Taxonomic_group_response),
                    function(x){subset(combined,Taxonomic_group_response==x)}))
mycat=unique(combined$Taxonomic_group_response)
for(i in 1:length(mycat)) {
      dat_list10=c(dat_list10,lapply(plant_sed,
                                     function(x){subset(combined,Taxonomic_group_response==mycat[i]&Secondary.classification.of.animal.........二级分类==x)}))
      dat_list10=c(dat_list10,lapply(microo_sed,
                                     function(x){subset(combined,Taxonomic_group_response==mycat[i]&Secondary.classification.of.animal.........二级分类==x)}))
}

for (i in 1:length(mycat)) {
    Nontarget_category=c(Nontarget_category,paste(plant_sed,mycat[i],sep = ' '))
    Nontarget_category=c(Nontarget_category,paste(microo_sed,mycat[i],sep = ' '))
}


S11=list()
for (k in 1:length(fungici)) {
  fungici_cate=sort(unique(subset(combined,
                                  pesticide_by_target_organisms=='fungicides'&pesticide_by_source==fungici[k])$Secondary.classification.of.insecticide........二级分类))
  S11[[k]]=matrix(NA,length(fungici_cate)*length(dat_list10),11)
  for (i in 1:length(fungici_cate)) {
    for (j in 1:length(dat_list10)){
      mydat=subset(dat_list10[[j]],
                   pesticide_by_target_organisms=='fungicides'&pesticide_by_source==fungici[k]&Secondary.classification.of.insecticide........二级分类==fungici_cate[i])
      if (nrow(mydat)<=1) {
        S11[[k]][(i-1)*length(dat_list10)+j,]=c(paste(fungici[k],'fungicides',sep = '_'),
                                                fungici_cate[i],Nontarget_category[j],
                                                nrow(mydat),length(unique(mydat$Code)),rep(NA,6))
      } else {
        res=rma.mv(yi, vi, 
                   random = list(~ 1 | factor(Insecticide.name),~1|Code),
                   data=mydat)
        S11[[k]][(i-1)*length(dat_list10)+j,]=c(paste(fungici[k],'fungicides',sep = '_'),
                                                fungici_cate[i],Nontarget_category[j],
                                                nrow(mydat),length(unique(mydat$Code)),
                                                res$b,res$zval,nrow(mydat)-1,
                                                res$pval,res$ci.lb,res$ci.ub)
      }
    }
  }
}
S11=do.call('rbind',S11)
colnames(S11)<-c('fungicides','fungicides category','Non-target organism category',
                 "#Obs","#Studies","Effect size","T-value","Df",
                 "P-value","CL.lb","CL.ub")
write.xlsx(S11,"Supplementary Table 10.xls",row.names = F)

########## Supplementary Table 11 ###########
animal_sed=sort(read.xlsx('inter.xls',sheetName = 'animal-1',header = T)[,1])
plant_sed=sort(read.xlsx('inter.xls',sheetName = 'plant-1',header = T)[,1])
microo_sed=sort(read.xlsx('inter.xls',sheetName = 'microo-1',header = T)[,1])
dat_list10=c(lapply(unique(combined$Taxonomic_group),
                    function(x){subset(combined,Taxonomic_group==x)}),
             lapply(unique(combined$Taxonomic_group_response),
                    function(x){subset(combined,Taxonomic_group_response==x)}))
mycat=unique(combined$Taxonomic_group_response)
for(i in 1:length(mycat)) {
      dat_list10=c(dat_list10,lapply(plant_sed,
                                     function(x){subset(combined,Taxonomic_group_response==mycat[i]&Secondary.classification.of.animal.........二级分类==x)}))
      dat_list10=c(dat_list10,lapply(microo_sed,
                                     function(x){subset(combined,Taxonomic_group_response==mycat[i]&Secondary.classification.of.animal.........二级分类==x)}))
}
for (i in 1:length(mycat)) {
  Nontarget_category=c(Nontarget_category,paste(plant_sed,mycat[i],sep = ' '))
    Nontarget_category=c(Nontarget_category,paste(microo_sed,mycat[i],sep = ' '))
}



S12=list()
for (k in 1:length(herbici)) {
  herbici_cate=sort(unique(subset(combined,
                                  pesticide_by_target_organisms=='herbicides'&pesticide_by_source==herbici[k])$Secondary.classification.of.insecticide........二级分类))
  S12[[k]]=matrix(NA,length(herbici_cate)*length(dat_list10),11)
  for (i in 1:length(herbici_cate)) {
    for (j in 1:length(dat_list10)){
      mydat=subset(dat_list10[[j]],
                   pesticide_by_target_organisms=='herbicides'&pesticide_by_source==herbici[k]&Secondary.classification.of.insecticide........二级分类==herbici_cate[i])
      if (nrow(mydat)<=1) {
        S12[[k]][(i-1)*length(dat_list10)+j,]=c(paste(herbici[k],'herbicides',sep = '_'),
                                                herbici_cate[i],Nontarget_category[j],
                                                nrow(mydat),length(unique(mydat$Code)),rep(NA,6))
      } else {
        res=rma.mv(yi, vi, 
                   random = list(~ 1 | factor(Insecticide.name),~1|Code),
                   data=mydat)
        S12[[k]][(i-1)*length(dat_list10)+j,]=c(paste(herbici[k],'herbicides',sep = '_'),
                                                herbici_cate[i],Nontarget_category[j],
                                                nrow(mydat),length(unique(mydat$Code)),
                                                res$b,res$zval,nrow(mydat)-1,
                                                res$pval,res$ci.lb,res$ci.ub)
      }
    }
  }
}
S12=do.call('rbind',S12)
colnames(S12)<-c('herbicides','herbicides category','Non-target organism category',
                 "#Obs","#Studies","Effect size","T-value","Df",
                 "P-value","CL.lb","CL.ub")
write.xlsx(S12,"Supplementary Table 11.xls",row.names = F)

############### Supplementary Table 12 ##############
datlis_S13=list(combined[combined$Taxonomic_group=='animals',],
                combined[combined$Primary.classification.of.animal......一级分类=='Invertebrate',],
                combined[combined$Primary.classification.of.animal......一级分类=='Vertebrate',],
                combined[combined$Taxonomic_group=='plants',],
                combined[combined$Primary.classification.of.animal......一级分类=='Seed plant',],
                combined[combined$Primary.classification.of.animal......一级分类=='Spore-producing  plant',],
                combined[combined$Taxonomic_group=='microorganisms',],
                combined[combined$Primary.classification.of.animal......一级分类=='Bacteria',],
                combined[combined$Primary.classification.of.animal......一级分类=='Fungus',])
S13=list()
for (i in 1:length(datlis_S13)) {
  S13[[i]]=matrix(NA,length(Indicator),9)
  for (j in 1:length(Indicator)) {
    mydat=subset(datlis_S13[[i]])
    res1=rma(yi, vi, data = mydat,)
    res2=trimfill(res1,control=list(stepadj=0.5,maxiter=200),verbose=T)
    S13_group=unique(datlis_S13[[i]]$Primary.classification.of.animal......一级分类)
    if (length(S13_group)>1) {
      S13_group=unique(datlis_S13[[i]]$Taxonomic_group)
    }
    S13[[i]][j,]=c(S13_group,Indicator[j],res2$k0,
                   res2$b,res2$zval,nrow(mydat)-1,
                   res2$pval,res2$ci.lb,res2$ci.ub)
  }
}
S13=do.call('rbind',S13)
colnames(S13)<-c('Phylogenetic group','Indicator',"Missing studies",
                 "Effect size","T-value","Df",
                 "P-value","CL.lb","CL.ub")
write.xlsx(S13,"Supplementary Table 12.xls",row.names = F)

############### Supplementary Table for relationship of figures ##############
require(nlme)
require(effects)
library(MuMIn)

linear_table=function(dat,Fig_Source,columname,roname){
  
  if (columname %in% c('Animal','Invertebrate animal','Vertebrate animal')) {
    mylabel=c('growth','reproduction','behavior','biomarker')
  } else {
    mylabel=c('growth','reproduction','biomarker')
  }
  
  res14=matrix(NA,length(mylabel),14)
  for(i in 1:length(mylabel)) {
    subdat=subset(dat,Taxonomic_group_response==grep(mylabel[i],unique(dat$Taxonomic_group_response),value = T))
    subdat=subdat[!is.na(subdat$xi),]
    x <- subdat[,'xi']
    y <- subdat[,'yi']
    mydata=data.frame(x=x,y=y,Code=subdat[,'Code'],
                      pesticide_by_target_organisms=subdat[,'pesticide_by_target_organisms'],
                      Insecticide.name=subdat[,'Insecticide.name'],vi=subdat[,'vi'])
    assign("mydata2", mydata, envir=.GlobalEnv) 
    fit<-try(lme(y~x,
                 random =list(~ 1|Insecticide.name,~ 1|Code),
                 weights = varFixed(~ vi),data=mydata2,control = lmeControl(sigma = 1)))
      parameter=signif(summary(fit)$tTable[2,-1],4)
      con95=round(intervals(fit, 'fixed', level=0.95)[[1]][2,c(1,3)],4)
      if (Intercept_x[1]<0) {
        equation=paste0('Y = ',Intercept_x[2],'X - ',abs(Intercept_x[1]))
      } else {
        equation=paste0('Y = ',Intercept_x[2],'X + ',abs(Intercept_x[1]))
      }
      Threshhold=round(-Intercept_x[1]/Intercept_x[2],4)
      R2=round(r.squaredGLMM(fit)[2],4)
      res14[i,]=c(Fig_Source,paste(columname,mylabel[i],'response',sep = ' '),
                paste(roname,'application rate added over the control',sep = ' '),
                nrow(subdat),length(unique(subdat$Code)),
                equation,con95,parameter,Threshhold,R2)
  }
  return(as.data.frame(res14))
}

###### Supplementary Table 13 ######
res=list()

###### Figure 3 ######
### Figure 3a
dat=combined[combined$Taxonomic_group=='animals',]
Fig_Source='Fig. 3a'
columname='Animal'
roname='Pesticide'
res[['fig3a']]=linear_table(dat,Fig_Source,columname,roname)

### Figure 3b
dat=combined[combined$Taxonomic_group=='plants',]
Fig_Source='Fig. 3b'
columname='Plant'
roname='Pesticide'
res[['fig3b']]=linear_table(dat,Fig_Source,columname,roname)

### Figure 3c
dat=combined[combined$Taxonomic_group=='microorganisms',]
Fig_Source='Fig. 3c'
columname='Microorganism'
roname='Pesticide'
res[['fig3c']]=linear_table(dat,Fig_Source,columname,roname)

### Figure 3d
dat=subset(combined,Taxonomic_group=='animals'&pesticide_by_target_organisms=='insecticides')
Fig_Source='Fig. 3d'
columname='Animal'
roname='Insecticide'
res[['fig3d']]=linear_table(dat,Fig_Source,columname,roname)

### Figure 3e
dat=subset(combined,Taxonomic_group=='plants'&pesticide_by_target_organisms=='insecticides')
Fig_Source='Fig. 3e'
columname='Plant'
roname='Insecticide'
res[['fig3e']]=linear_table(dat,Fig_Source,columname,roname)

### Figure 3f
dat=subset(combined,Taxonomic_group=='microorganisms'&pesticide_by_target_organisms=='insecticides')
Fig_Source='Fig. 3f'
columname='Microorganism'
roname='Insecticide'
res[['fig3f']]=linear_table(dat,Fig_Source,columname,roname)

### Figure 3g
dat=subset(combined,Taxonomic_group=='animals'&pesticide_by_target_organisms=='fungicides')
Fig_Source='Fig. 3g'
columname='Animal'
roname='Fungicide'
res[['fig3g']]=linear_table(dat,Fig_Source,columname,roname)

### Figure 3h
dat=subset(combined,Taxonomic_group=='plants'&pesticide_by_target_organisms=='fungicides')
Fig_Source='Fig. 3h'
columname='Plant'
roname='Fungicide'
res[['fig3h']]=linear_table(dat,Fig_Source,columname,roname)

### Figure 3i
dat=subset(combined,Taxonomic_group=='microorganisms'&pesticide_by_target_organisms=='fungicides')
Fig_Source='Fig. 3i'
columname='Microorganism'
roname='Fungicide'
res[['fig3i']]=linear_table(dat,Fig_Source,columname,roname)

### Figure 3j
dat=subset(combined,Taxonomic_group=='animals'&pesticide_by_target_organisms=='herbicides')
Fig_Source='Fig. 3j'
columname='Animal'
roname='Herbicides'
res[['fig3j']]=linear_table(dat,Fig_Source,columname,roname)

### Figure 3k
dat=subset(combined,Taxonomic_group=='plants'&pesticide_by_target_organisms=='herbicides')
Fig_Source='Fig. 3k'
columname='Plant'
roname='Herbicides'
res[['fig3k']]=linear_table(dat,Fig_Source,columname,roname)

### Figure 3l
dat=subset(combined,Taxonomic_group=='microorganisms'&pesticide_by_target_organisms=='herbicides')
Fig_Source='Fig. 3l'
columname='Microorganism'
roname='Herbicides'
res[['fig3l']]=linear_table(dat,Fig_Source,columname,roname)

####### Supplementary Figure 1 ######
dat_F6=subset(combined,Experiment.type=='Laboratory')

### Supplementary Figure 1a
dat=dat_F6[dat_F6$Taxonomic_group=='animals',]
Fig_Source='Supplementary Fig. 1a'
columname='Animal'
roname='Pesticide'
res[['Supplementary Fig. 1a']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 1b
dat=dat_F6[dat_F6$Taxonomic_group=='plants',]
Fig_Source='Supplementary Fig. 1b'
columname='Plant'
roname='Pesticide'
res[['Supplementary Fig. 1b']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 1c
dat=dat_F6[dat_F6$Taxonomic_group=='microorganisms',]
Fig_Source='Supplementary Fig. 1c'
columname='Microorganism'
roname='Pesticide'
res[['Supplementary Fig. 1c']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 1d
dat=subset(dat_F6,Trophic_group=='animals'&pesticide_by_target_organisms=='insecticides')
Fig_Source='Supplementary Fig. 1d'
columname='Animal'
roname='Insecticide'
res[['Supplementary Fig. 1d']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 1e
dat=subset(dat_F6,Taxonomic_group=='plants'&pesticide_by_target_organisms=='insecticides')
Fig_Source='Supplementary Fig. 1e'
columname='Plant'
roname='Insecticide'
res[['Supplementary Fig. 1e']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 1f
dat=subset(dat_F6,Taxonomic_group=='microorganisms'&pesticide_by_target_organisms=='insecticides')
Fig_Source='Supplementary Fig. 1f'
columname='Microorganism'
roname='Insecticide'
res[['Supplementary Fig. 1f']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 1g
dat=subset(dat_F6,Taxonomic_group=='animals'&pesticide_by_target_organisms=='fungicides')
Fig_Source='Supplementary Fig. 1g'
columname='Animal'
roname='Fungicide'
res[['Supplementary Fig. 1g']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 1h
dat=subset(dat_F6,Taxonomic_group=='plants'&pesticide_by_target_organisms=='fungicides')
Fig_Source='Supplementary Fig. 1h'
columname='Plant'
roname='Fungicide'
res[['Supplementary Fig. 1h']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 1i
dat=subset(dat_F6,Taxonomic_group=='microorganisms'&pesticide_by_target_organisms=='fungicides')
Fig_Source='Supplementary Fig. 1i'
columname='Microorganism'
roname='Fungicide'
res[['Supplementary Fig. 1i']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 1j
dat=subset(dat_F6,Taxonomic_group=='animals'&pesticide_by_target_organisms=='herbicides')
Fig_Source='Supplementary Fig. 1j'
columname='Animal'
roname='Herbicides'
res[['Supplementary Fig. 1j']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 1k
dat=subset(dat_F6,Taxonomic_group=='plants'&pesticide_by_target_organisms=='herbicides')
Fig_Source='Supplementary Fig. 1k'
columname='Plant'
roname='Herbicides'
res[['Supplementary Fig. 1k']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 1l
dat=subset(dat_F6,Taxonomic_group=='microorganisms'&pesticide_by_target_organisms=='herbicides')
Fig_Source='Supplementary Fig. 1l'
columname='Microorganism'
roname='Herbicides'
res[['Supplementary Fig. 1l']]=linear_table(dat,Fig_Source,columname,roname)


####### Supplementary Figure 2 ######
dat_F7=subset(combined,Experiment.type=='Field')

### Supplementary Figure 2a
dat=dat_F7[dat_F7$Taxonomic_group=='animals',]
Fig_Source='Supplementary Fig. 2a'
columname='Animal'
roname='Pesticide'
res[['Supplementary Fig. 2a']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 2b
dat=dat_F7[dat_F7$Taxonomic_group=='plants',]
Fig_Source='Supplementary Fig. 2b'
columname='Plant'
roname='Pesticide'
res[['Supplementary Fig. 2b']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 2c
dat=dat_F7[dat_F7$Taxonomic_group=='microorganisms',]
Fig_Source='Supplementary Fig. 2c'
columname='Microorganism'
roname='Pesticide'
res[['Supplementary Fig. 2c']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 2d
dat=subset(dat_F7,Taxonomic_group=='animals'&pesticide_by_target_organisms=='insecticides')
Fig_Source='Supplementary Fig. 2d'
columname='Animal'
roname='Insecticide'
res[['Supplementary Fig. 2d']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 2e
dat=subset(dat_F7,Taxonomic_group=='plants'&pesticide_by_target_organisms=='insecticides')
Fig_Source='Supplementary Fig. 2e'
columname='Plant'
roname='Insecticide'
res[['Supplementary Fig. 2e']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 2f
dat=subset(dat_F7,Taxonomic_group=='microorganisms'&pesticide_by_target_organisms=='insecticides')
Fig_Source='Supplementary Fig. 2f'
columname='Microorganism'
roname='Insecticide'
res[['Supplementary Fig. 2f']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 2g
dat=subset(dat_F7,Taxonomic_group=='animals'&pesticide_by_target_organisms=='fungicides')
Fig_Source='Supplementary Fig. 2g'
columname='Animal'
roname='Fungicide'
res[['Supplementary Fig. 2g']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 2h
dat=subset(dat_F7,Taxonomic_group=='plants'&pesticide_by_target_organisms=='fungicides')
Fig_Source='Supplementary Fig. 2h'
columname='Plant'
roname='Fungicide'
res[['Supplementary Fig. 2h']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 2i
dat=subset(dat_F7,Taxonomic_group=='microorganisms'&pesticide_by_target_organisms=='fungicides')
Fig_Source='Supplementary Fig. 2i'
columname='Microorganism'
roname='Fungicide'
res[['Supplementary Fig. 2i']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 2j
dat=subset(dat_F7,Taxonomic_group=='animals'&pesticide_by_target_organisms=='herbicides')
Fig_Source='Supplementary Fig. 2j'
columname='Animal'
roname='Herbicides'
res[['Supplementary Fig. 2j']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 2k
dat=subset(dat_F7,Taxonomic_group=='plants'&pesticide_by_target_organisms=='herbicides')
Fig_Source='Supplementary Fig. 2k'
columname='Plant'
roname='Herbicides'
res[['Supplementary Fig. 2k']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 2l
dat=subset(dat_F7,Taxonomic_group=='microorganisms'&pesticide_by_target_organisms=='herbicides')
Fig_Source='Supplementary Fig. 2l'
columname='Microorganism'
roname='Herbicides'
res[['Supplementary Fig. 2l']]=linear_table(dat,Fig_Source,columname,roname)


####### Supplementary Figure 3 ######
dat_F8=subset(combined,Experiment.type=='Field'&Major.climatic.zones=='Temperate')

### Supplementary Figure 3a
dat=dat_F8[dat_F8$Taxonomic_group=='animals',]
Fig_Source='Supplementary Fig. 3a'
columname='Animal'
roname='Pesticide'
res[['Supplementary Fig. 3a']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 3b
dat=dat_F8[dat_F8$Taxonomic_group=='plants',]
Fig_Source='Supplementary Fig. 3b'
columname='Plant'
roname='Pesticide'
res[['Supplementary Fig. 3b']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 3c
dat=dat_F8[dat_F8$Taxonomic_group=='microorganisms',]
Fig_Source='Supplementary Fig. 3c'
columname='Microorganism'
roname='Pesticide'
res[['Supplementary Fig. 3c']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 3d
dat=subset(dat_F8,Taxonomic_group=='animals'&pesticide_by_target_organisms=='insecticides')
Fig_Source='Supplementary Fig. 3d'
columname='Animal'
roname='Insecticide'
res[['Supplementary Fig. 3d']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 3e
dat=subset(dat_F8,Taxonomic_group=='plants'&pesticide_by_target_organisms=='insecticides')
Fig_Source='Supplementary Fig. 3e'
columname='Plant'
roname='Insecticide'
res[['Supplementary Fig. 3e']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 3f
dat=subset(dat_F8,Taxonomic_group=='microorganisms'&pesticide_by_target_organisms=='insecticides')
Fig_Source='Supplementary Fig. 3f'
columname='Microorganism'
roname='Insecticide'
res[['Supplementary Fig. 3f']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 3g
dat=subset(dat_F8,Taxonomic_group=='animals'&pesticide_by_target_organisms=='fungicides')
Fig_Source='Supplementary Fig. 3g'
columname='Animal'
roname='Fungicide'
res[['Supplementary Fig. 3g']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 3h
dat=subset(dat_F8,Taxonomic_group=='plants'&pesticide_by_target_organisms=='fungicides')
Fig_Source='Supplementary Fig. 3h'
columname='Plant'
roname='Fungicide'
res[['Supplementary Fig. 3h']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 3i
dat=subset(dat_F8,Taxonomic_group=='microorganisms'&pesticide_by_target_organisms=='fungicides')
Fig_Source='Supplementary Fig. 3i'
columname='Microorganism'
roname='Fungicide'
res[['Supplementary Fig. 3i']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 3j
dat=subset(dat_F8,Taxonomic_group=='animals'&pesticide_by_target_organisms=='herbicides')
Fig_Source='Supplementary Fig. 3j'
columname='Animal'
roname='Herbicides'
res[['Supplementary Fig. 3j']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 3k
dat=subset(dat_F8,Taxonomic_group=='plants'&pesticide_by_target_organisms=='herbicides')
Fig_Source='Supplementary Fig. 3k'
columname='Plant'
roname='Herbicides'
res[['Supplementary Fig. 3k']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 3l
dat=subset(dat_F8,Taxonomic_group=='microorganisms'&pesticide_by_target_organisms=='herbicides')
Fig_Source='Supplementary Fig. 3l'
columname='Microorganism'
roname='Herbicides'
res[['Supplementary Fig. 3l']]=linear_table(dat,Fig_Source,columname,roname)

####### Supplementary Figure 4 ######
dat_F9=subset(combined,Experiment.type=='Field'&Major.climatic.zones=='Tropical')

### Supplementary Figure 4a
dat=dat_F9[dat_F9$Taxonomic_group=='animals',]
Fig_Source='Supplementary Fig. 4a'
columname='Animal'
roname='Pesticide'
res[['Supplementary Fig. 4a']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 4b
dat=dat_F9[dat_F9$Taxonomic_group=='plants',]
Fig_Source='Supplementary Fig. 4b'
columname='Plant'
roname='Pesticide'
res[['Supplementary Fig. 4b']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 4c
dat=dat_F9[dat_F9$Taxonomic_group=='microorganisms',]
Fig_Source='Supplementary Fig. 4c'
columname='Microorganism'
roname='Pesticide'
res[['Supplementary Fig. 4c']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 4d
dat=subset(dat_F9,Taxonomic_group=='animals'&pesticide_by_target_organisms=='insecticides')
Fig_Source='Supplementary Fig. 4d'
columname='Animal'
roname='Insecticide'
res[['Supplementary Fig. 4d']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 4e
dat=subset(dat_F9,Taxonomic_group=='plants'&pesticide_by_target_organisms=='insecticides')
Fig_Source='Supplementary Fig. 4e'
columname='Plant'
roname='Insecticide'
res[['Supplementary Fig. 4e']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 4f
dat=subset(dat_F9,Taxonomic_group=='microorganisms'&pesticide_by_target_organisms=='insecticides')
Fig_Source='Supplementary Fig. 4f'
columname='Microorganism'
roname='Insecticide'
res[['Supplementary Fig. 4f']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 4g
dat=subset(dat_F9,Taxonomic_group=='animals'&pesticide_by_target_organisms=='fungicides')
Fig_Source='Supplementary Fig. 4g'
columname='Animal'
roname='Fungicide'
res[['Supplementary Fig. 4g']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 4h
dat=subset(dat_F9,Taxonomic_group=='plants'&pesticide_by_target_organisms=='fungicides')
Fig_Source='Supplementary Fig. 4h'
columname='Plant'
roname='Fungicide'
res[['Supplementary Fig. 4h']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 4i
dat=subset(dat_F9,Taxonomic_group=='microorganisms'&pesticide_by_target_organisms=='fungicides')
Fig_Source='Supplementary Fig. 4i'
columname='Microorganism'
roname='Fungicide'
res[['Supplementary Fig. 4i']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 4j
dat=subset(dat_F9,Taxonomic_group=='animals'&pesticide_by_target_organisms=='herbicides')
Fig_Source='Supplementary Fig. 4j'
columname='Animal'
roname='Herbicides'
res[['Supplementary Fig. 4j']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 4k
dat=subset(dat_F9,Taxonomic_group=='plants'&pesticide_by_target_organisms=='herbicides')
Fig_Source='Supplementary Fig. 4k'
columname='Plant'
roname='Herbicides'
res[['Supplementary Fig. 4k']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 4l
dat=subset(dat_F9,Taxonomic_group=='microorganisms'&pesticide_by_target_organisms=='herbicides')
Fig_Source='Supplementary Fig. 4l'
columname='Microorganism'
roname='Herbicides'
res[['Supplementary Fig. 4l']]=linear_table(dat,Fig_Source,columname,roname)


####### Supplementary Figure 5 ######
dat_F10=subset(combined,Types.of.organism.exposure.to.pesticides=='Aquatic')

### Supplementary Figure 5a
dat=dat_F10[dat_F10$Taxonomic_group=='animals',]
Fig_Source='Supplementary Fig. 5a'
columname='Animal'
roname='Pesticide'
res[['Supplementary Fig. 5a']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 5b
dat=dat_F10[dat_F10$Taxonomic_group=='plants',]
Fig_Source='Supplementary Fig. 5b'
columname='Plant'
roname='Pesticide'
res[['Supplementary Fig. 5b']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 5c
dat=dat_F10[dat_F10$Taxonomic_group=='microorganisms',]
Fig_Source='Supplementary Fig. 5c'
columname='Microorganism'
roname='Pesticide'
res[['Supplementary Fig. 5c']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 5d
dat=subset(dat_F10,Taxonomic_group=='animals'&pesticide_by_target_organisms=='insecticides')
Fig_Source='Supplementary Fig. 5d'
columname='Animal'
roname='Insecticide'
res[['Supplementary Fig. 5d']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 5e
dat=subset(dat_F10,Taxonomic_group=='plants'&pesticide_by_target_organisms=='insecticides')
Fig_Source='Supplementary Fig. 5e'
columname='Plant'
roname='Insecticide'
res[['Supplementary Fig. 5e']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 5f
dat=subset(dat_F10,Taxonomic_group=='microorganisms'&pesticide_by_target_organisms=='insecticides')
Fig_Source='Supplementary Fig. 5f'
columname='Microorganism'
roname='Insecticide'
res[['Supplementary Fig. 5f']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 5g
dat=subset(dat_F10,Taxonomic_group=='animals'&pesticide_by_target_organisms=='fungicides')
Fig_Source='Supplementary Fig. 5g'
columname='Animal'
roname='Fungicide'
res[['Supplementary Fig. 5g']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 5h
dat=subset(dat_F10,Taxonomic_group=='plants'&pesticide_by_target_organisms=='fungicides')
Fig_Source='Supplementary Fig. 5h'
columname='Plant'
roname='Fungicide'
res[['Supplementary Fig. 5h']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 5i
dat=subset(dat_F10,Taxonomic_group=='microorganisms'&pesticide_by_target_organisms=='fungicides')
Fig_Source='Supplementary Fig. 5i'
columname='Microorganism'
roname='Fungicide'
res[['Supplementary Fig. 5i']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 5j
dat=subset(dat_F10,Taxonomic_group=='animals'&pesticide_by_target_organisms=='herbicides')
Fig_Source='Supplementary Fig. 5j'
columname='Animal'
roname='Herbicides'
res[['Supplementary Fig. 5j']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 5k
dat=subset(dat_F10,Taxonomic_group=='plants'&pesticide_by_target_organisms=='herbicides')
Fig_Source='Supplementary Fig. 5k'
columname='Plant'
roname='Herbicides'
res[['Supplementary Fig. 5k']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 5l
dat=subset(dat_F10,Taxonomic_group=='microorganisms'&pesticide_by_target_organisms=='herbicides')
Fig_Source='Supplementary Fig. 5l'
columname='Microorganism'
roname='Herbicides'
res[['Supplementary Fig. 5l']]=linear_table(dat,Fig_Source,columname,roname)


####### Supplementary Figure 6 ######
dat_F11=subset(combined,Types.of.organism.exposure.to.pesticides=='Terrestrial')

### Supplementary Figure 5a
dat=dat_F11[dat_F11$Taxonomic_group=='animals',]
Fig_Source='Supplementary Fig. 6a'
columname='Animal'
roname='Pesticide'
res[['Supplementary Fig. 6a']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 6b
dat=dat_F11[dat_F11$Taxonomic_group=='plants',]
Fig_Source='Supplementary Fig. 6b'
columname='Plant'
roname='Pesticide'
res[['Supplementary Fig. 6b']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 6c
dat=dat_F11[dat_F11$Taxonomic_group=='microorganisms',]
Fig_Source='Supplementary Fig. 6c'
columname='Microorganism'
roname='Pesticide'
res[['Supplementary Fig. 6c']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 6d
dat=subset(dat_F11,Taxonomic_group=='animals'&pesticide_by_target_organisms=='insecticides')
Fig_Source='Supplementary Fig. 6d'
columname='Animal'
roname='Insecticide'
res[['Supplementary Fig. 6d']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 6e
dat=subset(dat_F11,Taxonomic_group=='plants'&pesticide_by_target_organisms=='insecticides')
Fig_Source='Supplementary Fig. 6e'
columname='Plant'
roname='Insecticide'
res[['Supplementary Fig. 6e']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 6f
dat=subset(dat_F11,Taxonomic_group=='microorganisms'&pesticide_by_target_organisms=='insecticides')
Fig_Source='Supplementary Fig. 6f'
columname='Microorganism'
roname='Insecticide'
res[['Supplementary Fig. 6f']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 6g
dat=subset(dat_F11,Taxonomic_group=='animals'&pesticide_by_target_organisms=='fungicides')
Fig_Source='Supplementary Fig. 6g'
columname='Animal'
roname='Fungicide'
res[['Supplementary Fig. 6g']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 6h
dat=subset(dat_F11,Taxonomic_group=='plants'&pesticide_by_target_organisms=='fungicides')
Fig_Source='Supplementary Fig. 6h'
columname='Plant'
roname='Fungicide'
res[['Supplementary Fig. 6h']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 6i
dat=subset(dat_F11,Taxonomic_group=='microorganisms'&pesticide_by_target_organisms=='fungicides')
Fig_Source='Supplementary Fig. 6i'
columname='Microorganism'
roname='Fungicide'
res[['Supplementary Fig. 6i']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 6j
dat=subset(dat_F11,Taxonomic_group=='animals'&pesticide_by_target_organisms=='herbicides')
Fig_Source='Supplementary Fig. 6j'
columname='Animal'
roname='Herbicides'
res[['Supplementary Fig. 6j']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 6k
dat=subset(dat_F11,Taxonomic_group=='plants'&pesticide_by_target_organisms=='herbicides')
Fig_Source='Supplementary Fig. 6k'
columname='Plant'
roname='Herbicides'
res[['Supplementary Fig. 6k']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 6l
dat=subset(dat_F11,Taxonomic_group=='microorganisms'&pesticide_by_target_organisms=='herbicides')
Fig_Source='Supplementary Fig. 6l'
columname='Microorganism'
roname='Herbicides'
res[['Supplementary Fig. 6l']]=linear_table(dat,Fig_Source,columname,roname)

####### Supplementary Figure 7 ######
dat_F12=subset(combined,pesticide_by_source=='chemical')

### Supplementary Figure 7a
dat=dat_F12[dat_F12$Taxonomic_group=='animals',]
Fig_Source='Supplementary Fig. 7a'
columname='Animal'
roname='Chemical synthetic pesticide'
res[['Supplementary Fig. 7a']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 7b
dat=dat_F12[dat_F12$Taxonomic_group=='plants',]
Fig_Source='Supplementary Fig. 7b'
columname='Plant'
roname='Chemical synthetic pesticide'
res[['Supplementary Fig. 7b']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 7c
dat=dat_F12[dat_F12$Taxonomic_group=='microorganisms',]
Fig_Source='Supplementary Fig. 7c'
columname='Microorganism'
roname='Chemical synthetic pesticide'
res[['Supplementary Fig. 7c']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 7d
dat=subset(dat_F12,Taxonomic_group=='animals'&pesticide_by_target_organisms=='insecticides')
Fig_Source='Supplementary Fig. 7d'
columname='Animal'
roname='Chemical synthetic insecticide'
res[['Supplementary Fig. 7d']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 7e
dat=subset(dat_F12,Taxonomic_group=='plants'&pesticide_by_target_organisms=='insecticides')
Fig_Source='Supplementary Fig. 7e'
columname='Plant'
roname='Chemical synthetic insecticide'
res[['Supplementary Fig. 7e']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 7f
dat=subset(dat_F12,Taxonomic_group=='microorganisms'&pesticide_by_target_organisms=='insecticides')
Fig_Source='Supplementary Fig. 7f'
columname='Microorganism'
roname='Chemical synthetic insecticide'
res[['Supplementary Fig. 7f']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 7g
dat=subset(dat_F12,Taxonomic_group=='animals'&pesticide_by_target_organisms=='fungicides')
Fig_Source='Supplementary Fig. 7g'
columname='Animal'
roname='Chemical synthetic fungicide'
res[['Supplementary Fig. 7g']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 7h
dat=subset(dat_F12,Taxonomic_group=='plants'&pesticide_by_target_organisms=='fungicides')
Fig_Source='Supplementary Fig. 7h'
columname='Plant'
roname='Chemical synthetic fungicide'
res[['Supplementary Fig. 7h']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 7i
dat=subset(dat_F12,Taxonomic_group=='microorganisms'&pesticide_by_target_organisms=='fungicides')
Fig_Source='Supplementary Fig. 7i'
columname='Microorganism'
roname='Chemical synthetic fungicide'
res[['Supplementary Fig. 7i']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 7j
dat=subset(dat_F12,Taxonomic_group=='animals'&pesticide_by_target_organisms=='herbicides')
Fig_Source='Supplementary Fig. 7j'
columname='Animal'
roname='Chemical synthetic herbicides'
res[['Supplementary Fig. 7j']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 7k
dat=subset(dat_F12,Taxonomic_group=='plants'&pesticide_by_target_organisms=='herbicides')
Fig_Source='Supplementary Fig. 7k'
columname='Plant'
roname='Chemical synthetic herbicides'
res[['Supplementary Fig. 7k']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 7l
dat=subset(dat_F12,Taxonomic_group=='microorganisms'&pesticide_by_target_organisms=='herbicides')
Fig_Source='Supplementary Fig. 7l'
columname='Microorganism'
roname='Chemical synthetic herbicides'
res[['Supplementary Fig. 7l']]=linear_table(dat,Fig_Source,columname,roname)


####### Supplementary Figure 8 ######
dat_F13=subset(combined,pesticide_by_source=='mineral-based')

### Supplementary Figure 8a
dat=dat_F13[dat_F13$Taxonomic_group=='animals',]
Fig_Source='Supplementary Fig. 8a'
columname='Animal'
roname='Mineral-based pesticide'
res[['Supplementary Fig. 8a']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 8b
dat=dat_F13[dat_F13$Taxonomic_group=='plants',]
Fig_Source='Supplementary Fig. 8b'
columname='Plant'
roname='Mineral-based pesticide'
res[['Supplementary Fig. 8b']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 8c
dat=dat_F13[dat_F13$Taxonomic_group=='microorganisms',]
Fig_Source='Supplementary Fig. 8c'
columname='Microorganism'
roname='Mineral-based pesticide'
res[['Supplementary Fig. 8c']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 8d
dat=subset(dat_F13,Taxonomic_group=='animals'&pesticide_by_target_organisms=='insecticides')
Fig_Source='Supplementary Fig. 8d'
columname='Animal'
roname='Mineral-based insecticide'
res[['Supplementary Fig. 8d']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 8e
dat=subset(dat_F13,Taxonomic_group=='plants'&pesticide_by_target_organisms=='insecticides')
Fig_Source='Supplementary Fig. 8e'
columname='Plant'
roname='Mineral-based insecticide'
res[['Supplementary Fig. 8e']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 8f
dat=subset(dat_F13,Taxonomic_group=='microorganisms'&pesticide_by_target_organisms=='insecticides')
Fig_Source='Supplementary Fig. 8f'
columname='Microorganism'
roname='Mineral-based insecticide'
res[['Supplementary Fig. 8f']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 8g
dat=subset(dat_F13,Taxonomic_group=='animals'&pesticide_by_target_organisms=='fungicides')
Fig_Source='Supplementary Fig. 8g'
columname='Animal'
roname='Mineral-based fungicide'
res[['Supplementary Fig. 8g']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 8h
dat=subset(dat_F13,Taxonomic_group=='plants'&pesticide_by_target_organisms=='fungicides')
Fig_Source='Supplementary Fig. 8h'
columname='Plant'
roname='Mineral-based fungicide'
res[['Supplementary Fig. 8h']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 8i
dat=subset(dat_F13,Taxonomic_group=='microorganisms'&pesticide_by_target_organisms=='fungicides')
Fig_Source='Supplementary Fig. 8i'
columname='Microorganism'
roname='Mineral-based fungicide'
res[['Supplementary Fig. 8i']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 8j
dat=subset(dat_F13,Taxonomic_group=='animals'&pesticide_by_target_organisms=='herbicides')
Fig_Source='Supplementary Fig. 8j'
columname='Animal'
roname='Mineral-based herbicide'
res[['Supplementary Fig. 8j']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 8k
dat=subset(dat_F13,Taxonomic_group=='plants'&pesticide_by_target_organisms=='herbicides')
Fig_Source='Supplementary Fig. 8k'
columname='Plant'
roname='Mineral-based herbicide'
res[['Supplementary Fig. 8k']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 8l
dat=subset(dat_F13,Taxonomic_group=='microorganisms'&pesticide_by_target_organisms=='herbicides')
Fig_Source='Supplementary Fig. 8l'
columname='Microorganism'
roname='Mineral-based herbicide'
res[['Supplementary Fig. 8l']]=linear_table(dat,Fig_Source,columname,roname)

####### Supplementary Figure 9 ######
dat_F14=subset(combined,pesticide_by_source=='biogenic')

### Supplementary Figure 9a
dat=dat_F14[dat_F14$Taxonomic_group=='animals',]
Fig_Source='Supplementary Fig. 9a'
columname='Animal'
roname='Biogenic pesticide'
res[['Supplementary Fig. 9a']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 9b
dat=dat_F14[dat_F14$Taxonomic_group=='plants',]
Fig_Source='Supplementary Fig. 9b'
columname='Plant'
roname='Biogenic pesticide'
res[['Supplementary Fig. 9b']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 9c
dat=dat_F14[dat_F14$Taxonomic_group=='microorganisms',]
Fig_Source='Supplementary Fig. 9c'
columname='Microorganism'
roname='Biogenic pesticide'
res[['Supplementary Fig. 9c']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 9d
dat=subset(dat_F14,Taxonomic_group=='animals'&pesticide_by_target_organisms=='insecticides')
Fig_Source='Supplementary Fig. 9d'
columname='Animal'
roname='Biogenic insecticide'
res[['Supplementary Fig. 9d']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 9e
dat=subset(dat_F14,Taxonomic_group=='plants'&pesticide_by_target_organisms=='insecticides')
Fig_Source='Supplementary Fig. 9e'
columname='Plant'
roname='Biogenic insecticide'
res[['Supplementary Fig. 9e']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 9f
dat=subset(dat_F14,Taxonomic_group=='microorganisms'&pesticide_by_target_organisms=='insecticides')
Fig_Source='Supplementary Fig. 9f'
columname='Microorganism'
roname='Biogenic insecticide'
res[['Supplementary Fig. 9f']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 9g
dat=subset(dat_F14,Taxonomic_group=='animals'&pesticide_by_target_organisms=='fungicides')
Fig_Source='Supplementary Fig. 9g'
columname='Animal'
roname='Biogenic fungicide'
res[['Supplementary Fig. 9g']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 9h
dat=subset(dat_F14,Taxonomic_group=='plants'&pesticide_by_target_organisms=='fungicides')
Fig_Source='Supplementary Fig. 9h'
columname='Plant'
roname='Biogenic fungicide'
res[['Supplementary Fig. 9h']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 9i
dat=subset(dat_F14,Taxonomic_group=='microorganisms'&pesticide_by_target_organisms=='fungicides')
Fig_Source='Supplementary Fig. 9i'
columname='Microorganism'
roname='Biogenic fungicide'
res[['Supplementary Fig. 9i']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 9j
dat=subset(dat_F14,Taxonomic_group=='animals'&pesticide_by_target_organisms=='herbicides')
Fig_Source='Supplementary Fig. 9j'
columname='Animal'
roname='Biogenic herbicide'
res[['Supplementary Fig. 9j']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 9k
dat=subset(dat_F14,Taxonomic_group=='plants'&pesticide_by_target_organisms=='herbicides')
Fig_Source='Supplementary Fig. 9k'
columname='Plant'
roname='Biogenic herbicide'
res[['Supplementary Fig. 9k']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 9l
dat=subset(dat_F14,Taxonomic_group=='microorganisms'&pesticide_by_target_organisms=='herbicides')
Fig_Source='Supplementary Fig. 9l'
columname='Microorganism'
roname='Biogenic herbicide'
res[['Supplementary Fig. 9l']]=linear_table(dat,Fig_Source,columname,roname)


####### Supplementary Figure 10 ######
dat_F15_I=subset(combined,Primary.classification.of.animal......一级分类=='Invertebrate')
dat_F15_V=subset(combined,Primary.classification.of.animal......一级分类=='Vertebrate')

### Supplementary Figure 10a
dat=dat_F15_I
Fig_Source='Supplementary Fig. 10a'
columname='Invertebrate animal'
roname='Pesticide'
res[['Supplementary Fig. 10a']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 10b
dat=dat_F15_V
Fig_Source='Supplementary Fig. 10b'
columname='Vertebrate animal'
roname='Pesticide'
res[['Supplementary Fig. 10b']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 10c
dat=subset(dat_F15_I,pesticide_by_source=='chemical')
Fig_Source='Supplementary Fig. 10c'
columname='Invertebrate animal'
roname='Chemical synthetic pesticide'
res[['Supplementary Fig. 10c']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 10d
dat=subset(dat_F15_V,pesticide_by_source=='chemical')
Fig_Source='Supplementary Fig. 10d'
columname='Vertebrate animal'
roname='Chemical synthetic pesticide'
res[['Supplementary Fig. 10d']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 10e
dat=subset(dat_F15_I,pesticide_by_source=='mineral-based')
Fig_Source='Supplementary Fig. 10e'
columname='Invertebrate animal'
roname='Mineral-based pesticide'
res[['Supplementary Fig. 10e']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 10f
dat=subset(dat_F15_V,pesticide_by_source=='mineral-based')
Fig_Source='Supplementary Fig. 10f'
columname='Vertebrate animal'
roname='Mineral-based pesticide'
res[['Supplementary Fig. 10f']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 10g
dat=subset(dat_F15_I,pesticide_by_source=='biogenic')
Fig_Source='Supplementary Fig. 10g'
columname='Invertebrate animal'
roname='Biogenic pesticide'
res[['Supplementary Fig. 10g']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 10h
dat=subset(dat_F15_V,pesticide_by_source=='biogenic')
Fig_Source='Supplementary Fig. 10h'
columname='Vertebrate animal'
roname='Biogenic pesticide'
res[['Supplementary Fig. 10h']]=linear_table(dat,Fig_Source,columname,roname)


####### Supplementary Figure 11 ######
dat_F16_I=subset(combined,Primary.classification.of.animal......一级分类=='Invertebrate'&pesticide_by_target_organisms=='insecticides')
dat_F16_V=subset(combined,Primary.classification.of.animal......一级分类=='Vertebrate'&pesticide_by_target_organisms=='insecticides')

### Supplementary Figure 11a
dat=dat_F16_I
Fig_Source='Supplementary Fig. 11a'
columname='Invertebrate animal'
roname='insecticide'
res[['Supplementary Fig. 11a']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 11b
dat=dat_F16_V
Fig_Source='Supplementary Fig. 11b'
columname='Vertebrate animal'
roname='insecticide'
res[['Supplementary Fig. 11b']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 11c
dat=subset(dat_F16_I,pesticide_by_source=='chemical')
Fig_Source='Supplementary Fig. 11c'
columname='Invertebrate animal'
roname='Chemical synthetic insecticide'
res[['Supplementary Fig. 11c']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 11d
dat=subset(dat_F16_V,pesticide_by_source=='chemical')
Fig_Source='Supplementary Fig. 11d'
columname='Vertebrate animal'
roname='Chemical synthetic insecticide'
res[['Supplementary Fig. 11d']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 11e
dat=subset(dat_F16_I,pesticide_by_source=='mineral-based')
Fig_Source='Supplementary Fig. 11e'
columname='Invertebrate animal'
roname='Mineral-based insecticide'
res[['Supplementary Fig. 11e']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 11f
dat=subset(dat_F16_V,pesticide_by_source=='mineral-based')
Fig_Source='Supplementary Fig. 11f'
columname='Vertebrate animal'
roname='Mineral-based insecticide'
res[['Supplementary Fig. 11f']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 11g
dat=subset(dat_F16_I,pesticide_by_source=='biogenic')
Fig_Source='Supplementary Fig. 11g'
columname='Invertebrate animal'
roname='Biogenic insecticide'
res[['Supplementary Fig. 11g']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 11h
dat=subset(dat_F16_V,pesticide_by_source=='biogenic')
Fig_Source='Supplementary Fig. 11h'
columname='Vertebrate animal'
roname='Biogenic insecticide'
res[['Supplementary Fig. 11h']]=linear_table(dat,Fig_Source,columname,roname)


####### Supplementary Figure 12 ######
dat_F17_I=subset(combined,Primary.classification.of.animal......一级分类=='Invertebrate'&pesticide_by_target_organisms=='fungicides')
dat_F17_V=subset(combined,Primary.classification.of.animal......一级分类=='Vertebrate'&pesticide_by_target_organisms=='fungicides')

### Supplementary Figure 12a
dat=dat_F17_I
Fig_Source='Supplementary Fig. 12a'
columname='Invertebrate animal'
roname='fungicide'
res[['Supplementary Fig. 12a']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 12b
dat=dat_F17_V
Fig_Source='Supplementary Fig. 12b'
columname='Vertebrate animal'
roname='fungicide'
res[['Supplementary Fig. 12b']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 12c
dat=subset(dat_F17_I,pesticide_by_source=='chemical')
Fig_Source='Supplementary Fig. 12c'
columname='Invertebrate animal'
roname='Chemical synthetic fungicide'
res[['Supplementary Fig. 12c']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 12d
dat=subset(dat_F17_V,pesticide_by_source=='chemical')
Fig_Source='Supplementary Fig. 12d'
columname='Vertebrate animal'
roname='Chemical synthetic fungicide'
res[['Supplementary Fig. 12d']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 12e
dat=subset(dat_F17_I,pesticide_by_source=='mineral-based')
Fig_Source='Supplementary Fig. 12e'
columname='Invertebrate animal'
roname='Mineral-based fungicide'
res[['Supplementary Fig. 12e']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 12f
dat=subset(dat_F17_V,pesticide_by_source=='mineral-based')
Fig_Source='Supplementary Fig. 12f'
columname='Vertebrate animal'
roname='Mineral-based fungicide'
res[['Supplementary Fig. 12f']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 12g
dat=subset(dat_F17_I,pesticide_by_source=='biogenic')
Fig_Source='Supplementary Fig. 12g'
columname='Invertebrate animal'
roname='Biogenic fungicide'
res[['Supplementary Fig. 12g']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 12h
dat=subset(dat_F17_V,pesticide_by_source=='biogenic')
Fig_Source='Supplementary Fig. 12h'
columname='Vertebrate animal'
roname='Biogenic fungicide'
res[['Supplementary Fig. 12h']]=linear_table(dat,Fig_Source,columname,roname)


####### Supplementary Figure 13 ######
dat_F18_I=subset(combined,Primary.classification.of.animal......一级分类=='Invertebrate'&pesticide_by_target_organisms=='herbicides')
dat_F18_V=subset(combined,Primary.classification.of.animal......一级分类=='Vertebrate'&pesticide_by_target_organisms=='herbicides')

### Supplementary Figure 13a
dat=dat_F18_I
Fig_Source='Supplementary Fig. 13a'
columname='Invertebrate animal'
roname='herbicide'
res[['Supplementary Fig. 13a']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 13b
dat=dat_F18_V
Fig_Source='Supplementary Fig. 13b'
columname='Vertebrate animal'
roname='herbicide'
res[['Supplementary Fig. 13b']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 13c
dat=subset(dat_F18_I,pesticide_by_source=='chemical')
Fig_Source='Supplementary Fig. 13c'
columname='Invertebrate animal'
roname='Chemical synthetic herbicide'
res[['Supplementary Fig. 13c']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 13d
dat=subset(dat_F18_V,pesticide_by_source=='chemical')
Fig_Source='Supplementary Fig. 13d'
columname='Vertebrate animal'
roname='Chemical synthetic herbicide'
res[['Supplementary Fig. 13d']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 13e
dat=subset(dat_F18_I,pesticide_by_source=='mineral-based')
Fig_Source='Supplementary Fig. 13e'
columname='Invertebrate animal'
roname='Mineral-based herbicide'
res[['Supplementary Fig. 13e']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 13f
dat=subset(dat_F18_V,pesticide_by_source=='mineral-based')
Fig_Source='Supplementary Fig. 13f'
columname='Vertebrate animal'
roname='Mineral-based herbicide'
res[['Supplementary Fig. 13f']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 13g
dat=subset(dat_F18_I,pesticide_by_source=='biogenic')
Fig_Source='Supplementary Fig. 13g'
columname='Invertebrate animal'
roname='Biogenic herbicide'
res[['Supplementary Fig. 13g']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 13h
dat=subset(dat_F18_V,pesticide_by_source=='biogenic')
Fig_Source='Supplementary Fig. 13h'
columname='Vertebrate animal'
roname='Biogenic herbicide'
res[['Supplementary Fig. 13h']]=linear_table(dat,Fig_Source,columname,roname)


####### Supplementary Figure 14 ######
dat_F19_1=subset(combined,Primary.classification.of.animal......一级分类=='Seed plant')
dat_F19_2=subset(combined,Primary.classification.of.animal......一级分类=='Spore-producing  plant')

### Supplementary Figure 14a
dat=dat_F19_1
Fig_Source='Supplementary Fig. 14a'
columname='Seed plant'
roname='Pesticide'
res[['Supplementary Fig. 14a']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 14b
dat=dat_F19_2
Fig_Source='Supplementary Fig. 14b'
columname='Spore-producing plant'
roname='Pesticide'
res[['Supplementary Fig. 14b']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 14c
dat=subset(dat_F19_1,pesticide_by_source=='chemical')
Fig_Source='Supplementary Fig. 14c'
columname='Seed plant'
roname='Chemical synthetic pesticide'
res[['Supplementary Fig. 14c']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 14d
dat=subset(dat_F19_2,pesticide_by_source=='chemical')
Fig_Source='Supplementary Fig. 14d'
columname='Spore-producing plant'
roname='Chemical synthetic pesticide'
res[['Supplementary Fig. 14d']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 14e
dat=subset(dat_F19_1,pesticide_by_source=='mineral-based')
Fig_Source='Supplementary Fig. 14e'
columname='Seed plant'
roname='Mineral-based pesticide'
res[['Supplementary Fig. 14e']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 14f
dat=subset(dat_F19_2,pesticide_by_source=='mineral-based')
Fig_Source='Supplementary Fig. 14f'
columname='Spore-producing plant'
roname='Mineral-based pesticide'
res[['Supplementary Fig. 14f']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 14g
dat=subset(dat_F19_1,pesticide_by_source=='biogenic')
Fig_Source='Supplementary Fig. 14g'
columname='Seed plant'
roname='Biogenic pesticide'
res[['Supplementary Fig. 14g']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 14h
dat=subset(dat_F19_2,pesticide_by_source=='biogenic')
Fig_Source='Supplementary Fig. 14h'
columname='Spore-producing plant'
roname='Biogenic pesticide'
res[['Supplementary Fig. 14h']]=linear_table(dat,Fig_Source,columname,roname)


####### Supplementary Figure 15 ######
dat_F20_1=subset(combined,Primary.classification.of.animal......一级分类=='Seed plant' & pesticide_by_target_organisms=='insecticides')
dat_F20_2=subset(combined,Primary.classification.of.animal......一级分类=='Spore-producing  plant' & pesticide_by_target_organisms=='insecticides')

### Supplementary Figure 15a
dat=dat_F20_1
Fig_Source='Supplementary Fig. 15a'
columname='Seed plant'
roname='Insecticide'
res[['Supplementary Fig. 15a']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 15b
dat=dat_F20_2
Fig_Source='Supplementary Fig. 15b'
columname='Spore-producing plant'
roname='Insecticide'
res[['Supplementary Fig. 15b']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 15c
dat=subset(dat_F20_1,pesticide_by_source=='chemical')
Fig_Source='Supplementary Fig. 15c'
columname='Seed plant'
roname='Chemical synthetic insecticide'
res[['Supplementary Fig. 15c']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 15d
dat=subset(dat_F20_2,pesticide_by_source=='chemical')
Fig_Source='Supplementary Fig. 15d'
columname='Spore-producing plant'
roname='Chemical synthetic insecticide'
res[['Supplementary Fig. 15d']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 15e
dat=subset(dat_F20_1,pesticide_by_source=='mineral-based')
Fig_Source='Supplementary Fig. 15e'
columname='Seed plant'
roname='Mineral-based insecticide'
res[['Supplementary Fig. 15e']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 15f
dat=subset(dat_F20_2,pesticide_by_source=='mineral-based')
Fig_Source='Supplementary Fig. 15f'
columname='Spore-producing plant'
roname='Mineral-based insecticide'
res[['Supplementary Fig. 15f']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 15g
dat=subset(dat_F20_1,pesticide_by_source=='biogenic')
Fig_Source='Supplementary Fig. 15g'
columname='Seed plant'
roname='Biogenic insecticide'
res[['Supplementary Fig. 15g']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 15h
dat=subset(dat_F20_2,pesticide_by_source=='biogenic')
Fig_Source='Supplementary Fig. 15h'
columname='Spore-producing plant'
roname='Biogenic insecticide'
res[['Supplementary Fig. 15h']]=linear_table(dat,Fig_Source,columname,roname)


####### Supplementary Figure 16 ######
dat_F21_1=subset(combined,Primary.classification.of.animal......一级分类=='Seed plant' & pesticide_by_target_organisms=='fungicides')
dat_F21_2=subset(combined,Primary.classification.of.animal......一级分类=='Spore-producing  plant' & pesticide_by_target_organisms=='fungicides')

### Supplementary Figure 16a
dat=dat_F21_1
Fig_Source='Supplementary Fig. 16a'
columname='Seed plant'
roname='fungicide'
res[['Supplementary Fig. 16a']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 16b
dat=dat_F21_2
Fig_Source='Supplementary Fig. 16b'
columname='Spore-producing plant'
roname='fungicide'
res[['Supplementary Fig. 16b']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 16c
dat=subset(dat_F21_1,pesticide_by_source=='chemical')
Fig_Source='Supplementary Fig. 16c'
columname='Seed plant'
roname='Chemical synthetic fungicide'
res[['Supplementary Fig. 16c']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 16d
dat=subset(dat_F21_2,pesticide_by_source=='chemical')
Fig_Source='Supplementary Fig. 16d'
columname='Spore-producing plant'
roname='Chemical synthetic fungicide'
res[['Supplementary Fig. 16d']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 16e
dat=subset(dat_F21_1,pesticide_by_source=='mineral-based')
Fig_Source='Supplementary Fig. 16e'
columname='Seed plant'
roname='Mineral-based fungicide'
res[['Supplementary Fig. 16e']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 16f
dat=subset(dat_F21_2,pesticide_by_source=='mineral-based')
Fig_Source='Supplementary Fig. 16f'
columname='Spore-producing plant'
roname='Mineral-based fungicide'
res[['Supplementary Fig. 16f']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 16g
dat=subset(dat_F21_1,pesticide_by_source=='biogenic')
Fig_Source='Supplementary Fig. 16g'
columname='Seed plant'
roname='Biogenic fungicide'
res[['Supplementary Fig. 16g']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 16h
dat=subset(dat_F21_2,pesticide_by_source=='biogenic')
Fig_Source='Supplementary Fig. 16h'
columname='Spore-producing plant'
roname='Biogenic fungicide'
res[['Supplementary Fig. 16h']]=linear_table(dat,Fig_Source,columname,roname)


####### Supplementary Figure 17 ######
dat_F22_1=subset(combined,Primary.classification.of.animal......一级分类=='Seed plant' & pesticide_by_target_organisms=='herbicides')
dat_F22_2=subset(combined,Primary.classification.of.animal......一级分类=='Spore-producing  plant' & pesticide_by_target_organisms=='herbicides')

### Supplementary Figure 17a
dat=dat_F22_1
Fig_Source='Supplementary Fig. 17a'
columname='Seed plant'
roname='herbicide'
res[['Supplementary Fig. 17a']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 17b
dat=dat_F22_2
Fig_Source='Supplementary Fig. 17b'
columname='Spore-producing plant'
roname='herbicide'
res[['Supplementary Fig. 17b']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 17c
dat=subset(dat_F22_1,pesticide_by_source=='chemical')
Fig_Source='Supplementary Fig. 17c'
columname='Seed plant'
roname='Chemical synthetic herbicide'
res[['Supplementary Fig. 17c']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 17d
dat=subset(dat_F22_2,pesticide_by_source=='chemical')
Fig_Source='Supplementary Fig. 17d'
columname='Spore-producing plant'
roname='Chemical synthetic herbicide'
res[['Supplementary Fig. 17d']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 17e
dat=subset(dat_F22_1,pesticide_by_source=='mineral-based')
Fig_Source='Supplementary Fig. 17e'
columname='Seed plant'
roname='Mineral-based herbicide'
res[['Supplementary Fig. 17e']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 17f
dat=subset(dat_F22_2,pesticide_by_source=='mineral-based')
Fig_Source='Supplementary Fig. 17f'
columname='Spore-producing plant'
roname='Mineral-based herbicide'
res[['Supplementary Fig. 17f']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 17g
dat=subset(dat_F22_1,pesticide_by_source=='biogenic')
Fig_Source='Supplementary Fig. 17g'
columname='Seed plant'
roname='Biogenic herbicide'
res[['Supplementary Fig. 17g']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 17h
dat=subset(dat_F22_2,pesticide_by_source=='biogenic')
Fig_Source='Supplementary Fig. 17h'
columname='Spore-producing plant'
roname='Biogenic herbicide'
res[['Supplementary Fig. 17h']]=linear_table(dat,Fig_Source,columname,roname)


####### Supplementary Figure 18 ######
dat_F23_1=subset(combined,Primary.classification.of.animal......一级分类=='Bacteria')
dat_F23_2=subset(combined,Primary.classification.of.animal......一级分类=='Fungus')

### Supplementary Figure 18a
dat=dat_F23_1
Fig_Source='Supplementary Fig. 18a'
columname='Bacteria'
roname='Pesticide'
res[['Supplementary Fig. 18a']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 18b
dat=dat_F23_2
Fig_Source='Supplementary Fig. 18b'
columname='Fungus'
roname='Pesticide'
res[['Supplementary Fig. 18b']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 18c
dat=subset(dat_F23_1,pesticide_by_source=='chemical')
Fig_Source='Supplementary Fig. 18c'
columname='Bacteria'
roname='Chemical synthetic pesticide'
res[['Supplementary Fig. 18c']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 18d
dat=subset(dat_F23_2,pesticide_by_source=='chemical')
Fig_Source='Supplementary Fig. 18d'
columname='Fungus'
roname='Chemical synthetic pesticide'
res[['Supplementary Fig. 18d']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 18e
dat=subset(dat_F23_1,pesticide_by_source=='mineral-based')
Fig_Source='Supplementary Fig. 18e'
columname='Bacteria'
roname='Mineral-based pesticide'
res[['Supplementary Fig. 18e']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 18f
dat=subset(dat_F23_2,pesticide_by_source=='mineral-based')
Fig_Source='Supplementary Fig. 18f'
columname='Fungus'
roname='Mineral-based pesticide'
res[['Supplementary Fig. 18f']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 18g
dat=subset(dat_F23_1,pesticide_by_source=='biogenic')
Fig_Source='Supplementary Fig. 18g'
columname='Bacteria'
roname='Biogenic pesticide'
res[['Supplementary Fig. 18g']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 18h
dat=subset(dat_F23_2,pesticide_by_source=='biogenic')
Fig_Source='Supplementary Fig. 18h'
columname='Fungus'
roname='Biogenic pesticide'
res[['Supplementary Fig. 18h']]=linear_table(dat,Fig_Source,columname,roname)


####### Supplementary Figure 19 ######
dat_F24_1=subset(combined,Primary.classification.of.animal......一级分类=='Bacteria' & pesticide_by_target_organisms=='insecticides')
dat_F24_2=subset(combined,Primary.classification.of.animal......一级分类=='Fungus' & pesticide_by_target_organisms=='insecticides')

### Supplementary Figure 19a
dat=dat_F24_1
Fig_Source='Supplementary Fig. 19a'
columname='Bacteria'
roname='Insecticide'
res[['Supplementary Fig. 19a']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 19b
dat=dat_F24_2
Fig_Source='Supplementary Fig. 19b'
columname='Fungus'
roname='Insecticide'
res[['Supplementary Fig. 19b']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 19c
dat=subset(dat_F24_1,pesticide_by_source=='chemical')
Fig_Source='Supplementary Fig. 19c'
columname='Bacteria'
roname='Chemical synthetic insecticide'
res[['Supplementary Fig. 19c']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 19d
dat=subset(dat_F24_2,pesticide_by_source=='chemical')
Fig_Source='Supplementary Fig. 19d'
columname='Fungus'
roname='Chemical synthetic insecticide'
res[['Supplementary Fig. 19d']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 19e
dat=subset(dat_F24_1,pesticide_by_source=='mineral-based')
Fig_Source='Supplementary Fig. 19e'
columname='Bacteria'
roname='Mineral-based insecticide'
res[['Supplementary Fig. 19e']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 19f
dat=subset(dat_F24_2,pesticide_by_source=='mineral-based')
Fig_Source='Supplementary Fig. 19f'
columname='Fungus'
roname='Mineral-based insecticide'
res[['Supplementary Fig. 19f']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 19g
dat=subset(dat_F24_1,pesticide_by_source=='biogenic')
Fig_Source='Supplementary Fig. 19g'
columname='Bacteria'
roname='Biogenic insecticide'
res[['Supplementary Fig. 19g']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 19h
dat=subset(dat_F24_2,pesticide_by_source=='biogenic')
Fig_Source='Supplementary Fig. 19h'
columname='Fungus'
roname='Biogenic insecticide'
res[['Supplementary Fig. 19h']]=linear_table(dat,Fig_Source,columname,roname)


####### Supplementary Figure 20 ######
dat_F25_1=subset(combined,Primary.classification.of.animal......一级分类=='Bacteria' & pesticide_by_target_organisms=='fungicides')
dat_F25_2=subset(combined,Primary.classification.of.animal......一级分类=='Fungus' & pesticide_by_target_organisms=='fungicides')

### Supplementary Figure 20a
dat=dat_F25_1
Fig_Source='Supplementary Fig. 20a'
columname='Bacteria'
roname='Fungicide'
res[['Supplementary Fig. 20a']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 20b
dat=dat_F25_2
Fig_Source='Supplementary Fig. 20b'
columname='Fungus'
roname='Fungicide'
res[['Supplementary Fig. 20b']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 20c
dat=subset(dat_F25_1,pesticide_by_source=='chemical')
Fig_Source='Supplementary Fig. 20c'
columname='Bacteria'
roname='Chemical synthetic fungicide'
res[['Supplementary Fig. 20c']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 20d
dat=subset(dat_F25_2,pesticide_by_source=='chemical')
Fig_Source='Supplementary Fig. 20d'
columname='Fungus'
roname='Chemical synthetic fungicide'
res[['Supplementary Fig. 20d']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 20e
dat=subset(dat_F25_1,pesticide_by_source=='mineral-based')
Fig_Source='Supplementary Fig. 20e'
columname='Bacteria'
roname='Mineral-based fungicide'
res[['Supplementary Fig. 20e']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 20f
dat=subset(dat_F25_2,pesticide_by_source=='mineral-based')
Fig_Source='Supplementary Fig. 20f'
columname='Fungus'
roname='Mineral-based fungicide'
res[['Supplementary Fig. 20f']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 20g
dat=subset(dat_F25_1,pesticide_by_source=='biogenic')
Fig_Source='Supplementary Fig. 20g'
columname='Bacteria'
roname='Biogenic fungicide'
res[['Supplementary Fig. 20g']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 20h
dat=subset(dat_F25_2,pesticide_by_source=='biogenic')
Fig_Source='Supplementary Fig. 20h'
columname='Fungus'
roname='Biogenic fungicide'
res[['Supplementary Fig. 20h']]=linear_table(dat,Fig_Source,columname,roname)


####### Supplementary Figure 21 ######
dat_F26_1=subset(combined,Primary.classification.of.animal......一级分类=='Bacteria' & pesticide_by_target_organisms=='herbicides')
dat_F26_2=subset(combined,Primary.classification.of.animal......一级分类=='Fungus' & pesticide_by_target_organisms=='herbicides')

### Supplementary Figure 21a
dat=dat_F26_1
Fig_Source='Supplementary Fig. 21a'
columname='Bacteria'
roname='Herbicide'
res[['Supplementary Fig. 21a']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 21b
dat=dat_F26_2
Fig_Source='Supplementary Fig. 21b'
columname='Fungus'
roname='Herbicide'
res[['Supplementary Fig. 21b']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 21c
dat=subset(dat_F26_1,pesticide_by_source=='chemical')
Fig_Source='Supplementary Fig. 21c'
columname='Bacteria'
roname='Chemical synthetic herbicide'
res[['Supplementary Fig. 21c']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 21d
dat=subset(dat_F26_2,pesticide_by_source=='chemical')
Fig_Source='Supplementary Fig. 21d'
columname='Fungus'
roname='Chemical synthetic herbicide'
res[['Supplementary Fig. 21d']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 21e
dat=subset(dat_F26_1,pesticide_by_source=='mineral-based')
Fig_Source='Supplementary Fig. 21e'
columname='Bacteria'
roname='Mineral-based herbicide'
res[['Supplementary Fig. 21e']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 21f
dat=subset(dat_F26_2,pesticide_by_source=='mineral-based')
Fig_Source='Supplementary Fig. 21f'
columname='Fungus'
roname='Mineral-based herbicide'
res[['Supplementary Fig. 21f']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 21g
dat=subset(dat_F26_1,pesticide_by_source=='biogenic')
Fig_Source='Supplementary Fig. 21g'
columname='Bacteria'
roname='Biogenic herbicide'
res[['Supplementary Fig. 21g']]=linear_table(dat,Fig_Source,columname,roname)

### Supplementary Figure 21h
dat=subset(dat_F26_2,pesticide_by_source=='biogenic')
Fig_Source='Supplementary Fig. 21h'
columname='Fungus'
roname='Biogenic herbicide'
res[['Supplementary Fig. 21h']]=linear_table(dat,Fig_Source,columname,roname)

###########
S14=do.call('rbind',res)
colnames(S14)=c('Figure Source','Effect size of non-target organism category (Y)',
                'Pesticide, insecticide, fungicide or herbicide application rate added over the control (X)',
                'Number of observations','Number of Studies','Regression equation',
                '95%CL','95%CU','Std.Error','DF','T-value','P-value','Threshhold application rate over the control','R2')
write.xlsx(S14,'Supplementary Table 13.xls',row.names = F)

