library(metafor)
library(xlsx)
require(nlme)
require(effects)
setwd('...')
combined=read.csv('Meta combined data.csv',header=T)

Non_target_animal=read.xlsx('Latin names of specific animal species-2022-1021.xls',sheetName='specific animal species',header=F)

Non_target_plant=read.xlsx('Latin names of specific plant species-2022-1020.xls',sheetName='specific plant species',header=F)

Non_target_microorganism=read.xlsx('Latin names of microorganism species-2022-1023.xls',sheetName='Sheet1',header=F)

Pesticide_species=read.xlsx('Supplementary Table 15-effect sizes of each pesticide on each animal species-框架/农药种类汇总.xlsx',sheetName='Sheet1',header=F)

##### Supplementary Table 14 ######
dat_S15=subset(combined,Taxonomic_group=='animals')
dat_S15=c(list(dat_S15),lapply(unique(dat_S15$Taxonomic_group_response),
                               function(x){subset(dat_S15,Taxonomic_group_response==x)}))
S15=list()
for(i in 1:nrow(Pesticide_species)) {
  S15[[i]]=matrix(NA,nrow(Non_target_animal)*length(dat_S15),9)
  for(j in 1:nrow(Non_target_animal)) {
    for(k in 1:length(dat_S15)) {
      mydat=subset(dat_S15[[k]],Insecticide.name==Pesticide_species[i,2])
      if (nrow(mydat)<=1) {
        S15[[i]][(j-1)*length(dat_S15)+k,]=c(t(Pesticide_species[i,1:2]),
                                             Non_target_animal[j,],Nontarget_category[k],
                                             nrow(mydat),length(unique(mydat$Code)),rep(NA,3))
      } else {
        res=rma.mv(yi, vi,
                   data=mydat)
        S15[[i]][(j-1)*length(dat_S15)+k,]=c(t(Pesticide_species[i,1:2]),
                                             Non_target_animal[j,],Nontarget_category[k],
                                             nrow(mydat),length(unique(mydat$Code)),
                                             signif(res$b,4),signif(res$zval,4),signif(res$pval,4))
      }
    }
  }
}
S15=do.call('rbind',S15)
colnames(S15)<-c('Pesticide category','Pesticide species','Non-target animal species',
                 'Animal species response categories',
                 "Obs","Studies","Effect size","T-value","P-value")
S15=subset(as.data.frame(S15), Obs>0)
S15=na.omit(S15)
write.xlsx(S15,"Supplementary Table 14.xls",row.names = F)


##### Supplementary Table 15 ######
dat_S16=subset(combined,Taxonomic_group=='plants')
dat_S16=c(list(dat_S16),lapply(unique(dat_S16$Taxonomic_group_response),
                               function(x){subset(dat_S16,Taxonomic_group_response==x)}))
S16=list()
for(i in 1:nrow(Pesticide_species)) {
  S16[[i]]=matrix(NA,nrow(Non_target_plant)*length(dat_S16),9)
  for(j in 1:nrow(Non_target_plant)) {
    for(k in 1:length(dat_S16)) {
      mydat=subset(dat_S16[[k]],Insecticide.name==Pesticide_species[i,2])
      if (nrow(mydat)<=1) {
        S16[[i]][(j-1)*length(dat_S16)+k,]=c(t(Pesticide_species[i,1:2]),
                                             Non_target_plant[j,],Nontarget_category[k],
                                             nrow(mydat),length(unique(mydat$Code)),rep(NA,3))
      } else {
        res=rma.mv(yi, vi,
                   data=mydat)
        S16[[i]][(j-1)*length(dat_S16)+k,]=c(t(Pesticide_species[i,1:2]),
                                             Non_target_plant[j,],Nontarget_category[k],
                                             nrow(mydat),length(unique(mydat$Code)),
                                             signif(res$b,4),signif(res$zval,4),signif(res$pval,4))
      }
    }
  }
}
S16=do.call('rbind',S16)
colnames(S16)<-c('Pesticide category','Pesticide species','Non-target plant species',
                 'Plant species response categories',
                 "Obs","Studies","Effect size","T-value","P-value")
S16=subset(as.data.frame(S16), Obs>0)
S16=na.omit(S16)
write.xlsx(S16,"Supplementary Table 15.xls",row.names = F)


##### Supplementary Table 16 ######
dat_S17=subset(combined,Taxonomic_group=='microorganisms')
dat_S17=c(list(dat_S17),lapply(unique(dat_S17$Taxonomic_group_response),
                               function(x){subset(dat_S17,Taxonomic_group_response==x)}))
S17=list()
for(i in 1:nrow(Pesticide_species)) {
  S17[[i]]=matrix(NA,nrow(Non_target_microorganism)*length(dat_S17),9)
  for(j in 1:nrow(Non_target_microorganism)) {
    for(k in 1:length(dat_S17)) {
      mydat=subset(dat_S17[[k]],Insecticide.name==Pesticide_species[i,2])
      if (nrow(mydat)<=1) {
        S17[[i]][(j-1)*length(dat_S17)+k,]=c(t(Pesticide_species[i,1:2]),
                                             Non_target_microorganism[j,],Nontarget_category[k],
                                             nrow(mydat),length(unique(mydat$Code)),rep(NA,3))
      } else {
        res=rma.mv(yi, vi,
                   data=mydat)
        S17[[i]][(j-1)*length(dat_S17)+k,]=c(t(Pesticide_species[i,1:2]),
                                             Non_target_microorganism[j,],Nontarget_category[k],
                                             nrow(mydat),length(unique(mydat$Code)),
                                             signif(res$b,4),signif(res$zval,4),signif(res$pval,4))
      }
    }
  }
}
S17=do.call('rbind',S17)
colnames(S17)<-c('Pesticide category','Pesticide species','Non-target microorganism species',
                 'Microorganism species response categories',
                 "Obs","Studies","Effect size","T-value","P-value")
S17=subset(as.data.frame(S17), Obs>0)
S17=na.omit(S17)
write.xlsx(S17,"Supplementary Table 16.xls",row.names = F)


##### Supplementary Table 17 ######
category=read.xlsx('Animal family-2022-1020.xls',sheetName='Sheet1',header=T)

dat_S18=subset(combined,Taxonomic_group=='animals')
dat_S18=c(list(dat_S18),lapply(unique(dat_S18$Taxonomic_group_response),
                               function(x){subset(dat_S18,Taxonomic_group_response==x)}))
S18=list()
for(i in 1:nrow(Pesticide_species)) {
  S18[[i]]=matrix(NA,length(unique(category[,2]))*length(dat_S18),9)
  for(j in 1:length(unique(category$Animal.Family))) {
    for(k in 1:length(dat_S18)) {
      mydat=subset(dat_S18[[k]],Insecticide.name==Pesticide_species[i,2] & Animal.Latin.name %in% Name)
      if (nrow(mydat)<=1) {
        S18[[i]][(j-1)*length(dat_S18)+k,]=c(t(Pesticide_species[i,1:2]),
                                             unique(category[,2])[j],Nontarget_category[k],
                                             nrow(mydat),length(unique(mydat$Code)),rep(NA,3))
      } else {
        res=rma.mv(yi, vi,
                   data=mydat)
        S18[[i]][(j-1)*length(dat_S18)+k,]=c(t(Pesticide_species[i,1:2]),
                                             unique(category[,2])[j],Nontarget_category[k],
                                             nrow(mydat),length(unique(mydat$Code)),
                                             signif(res$b,4),signif(res$zval,4),signif(res$pval,4))
      }
    }
  }
}
S18=do.call('rbind',S18)
colnames(S18)<-c('Pesticide category','Pesticide species','Non-target animal family',
                 'Animal family response categories',
                 "Obs","Studies","Effect size","T-value","P-value")
S18=subset(as.data.frame(S18), Obs>0)
S18=na.omit(S18)
write.xlsx(S18,"Supplementary Table 17.xls",row.names = F)


##### Supplementary Table 18 ######
category=read.xlsx('plant family-2022-1021.xls',sheetName='Sheet1',header=T)

dat_S19=subset(combined,Taxonomic_group=='plants')
dat_S19=c(list(dat_S19),lapply(unique(dat_S19$Taxonomic_group_response),
                               function(x){subset(dat_S19,Taxonomic_group_response==x)}))
S19=list()
for(i in 1:nrow(Pesticide_species)) {
  S19[[i]]=matrix(NA,length(unique(category[,2]))*length(dat_S19),9)
  for(j in 1:length(unique(category[,2]))) {
    for(k in 1:length(dat_S19)) {
      mydat=subset(dat_S19[[k]],Insecticide.name==Pesticide_species[i,2] & Animal.Latin.name %in% Name)
      if (nrow(mydat)<=1) {
        S19[[i]][(j-1)*length(dat_S19)+k,]=c(t(Pesticide_species[i,1:2]),
                                             unique(category[,2])[j],Nontarget_category[k],
                                             nrow(mydat),length(unique(mydat$Code)),rep(NA,3))
      } else {
        res=rma.mv(yi, vi,
                   data=mydat)
        S19[[i]][(j-1)*length(dat_S19)+k,]=c(t(Pesticide_species[i,1:2]),
                                             unique(category[,2])[j],Nontarget_category[k],
                                             nrow(mydat),length(unique(mydat$Code)),
                                             signif(res$b,4),signif(res$zval,4),signif(res$pval,4))
      }
    }
  }
}
S19=do.call('rbind',S19)
colnames(S19)<-c('Pesticide category','Pesticide species','Non-target plant family',
                 'Plant family response categories',
                 "Obs","Studies","Effect size","T-value","P-value")
S19=subset(as.data.frame(S19), Obs>0)
S19=na.omit(S19)
write.xlsx(S19,"Supplementary Table 18.xls",row.names = F)


##### Supplementary Table 19 ######
category=read.xlsx('microorganism family-2022-1023.xls',sheetName='Sheet1',header=T)

dat_S20=subset(combined,Taxonomic_group=='microorganisms')
dat_S20=c(list(dat_S20),lapply(unique(dat_S20$Taxonomic_group_response),
                               function(x){subset(dat_S20,Taxonomic_group_response==x)}))
S20=list()
for(i in 1:nrow(Pesticide_species)) {
  S20[[i]]=matrix(NA,length(unique(category[,2]))*length(dat_S20),9)
  for(j in 1:length(unique(category[,2]))) {
    for(k in 1:length(dat_S20)) {
      mydat=subset(dat_S20[[k]],Insecticide.name==Pesticide_species[i,2] & Animal.Latin.name %in% Name)
      if (nrow(mydat)<=1) {
        S20[[i]][(j-1)*length(dat_S20)+k,]=c(t(Pesticide_species[i,1:2]),
                                             unique(category[,2])[j],Nontarget_category[k],
                                             nrow(mydat),length(unique(mydat$Code)),rep(NA,3))
      } else {
        res=rma.mv(yi, vi,
                   data=mydat)
        S20[[i]][(j-1)*length(dat_S20)+k,]=c(t(Pesticide_species[i,1:2]),
                                             unique(category[,2])[j],Nontarget_category[k],
                                             nrow(mydat),length(unique(mydat$Code)),
                                             signif(res$b,4),signif(res$zval,4),signif(res$pval,4))
      }
    }
  }
}
S20=do.call('rbind',S20)
colnames(S20)<-c('Pesticide category','Pesticide species','Non-target microorganism family',
                 'Microorganism family response categories',
                 "Obs","Studies","Effect size","T-value","P-value")
S20=subset(as.data.frame(S20), Obs>0)
S20=na.omit(S20)
write.xlsx(S20,"Supplementary Table 19.xls",row.names = F)


##### Supplementary Table 20 ######
category=read.xlsx('Animal order-2022-1023.xls',sheetName='Sheet1',header=T)

dat_S21=subset(combined,Taxonomic_group=='animals')
dat_S21=c(list(dat_S21),lapply(unique(dat_S21$Taxonomic_group_response),
                               function(x){subset(dat_S21,Taxonomic_group_response==x)}))
S21=list()
for(i in 1:nrow(Pesticide_species)) {
  S21[[i]]=matrix(NA,length(unique(category[,2]))*length(dat_S21),9)
  for(j in 1:length(unique(category[,2]))) {
    for(k in 1:length(dat_S21)) {
      mydat=subset(dat_S21[[k]],Insecticide.name==Pesticide_species[i,2] & Animal.Latin.name %in% Name)
      if (nrow(mydat)<=1) {
        S21[[i]][(j-1)*length(dat_S21)+k,]=c(t(Pesticide_species[i,1:2]),
                                             unique(category[,2])[j],Nontarget_category[k],
                                             nrow(mydat),length(unique(mydat$Code)),rep(NA,3))
      } else {
        res=rma.mv(yi, vi,
                   data=mydat)
        S21[[i]][(j-1)*length(dat_S21)+k,]=c(t(Pesticide_species[i,1:2]),
                                             unique(category[,2])[j],Nontarget_category[k],
                                             nrow(mydat),length(unique(mydat$Code)),
                                             signif(res$b,4),signif(res$zval,4),signif(res$pval,4))
      }
    }
  }
}
S21=do.call('rbind',S21)
colnames(S21)<-c('Pesticide category','Pesticide species','Non-target animal order',
                 'Animal order response categories',
                 "Obs","Studies","Effect size","T-value","P-value")
S21=subset(as.data.frame(S21), Obs>0)
S21=na.omit(S21)
write.xlsx(S21,"Supplementary Table 20.xls",row.names = F)


##### Supplementary Table 21 ######
category=read.xlsx('plant order-2022-1023.xls',sheetName='Sheet1',header=T)

dat_S22=subset(combined,Taxonomic_group=='plants')
dat_S22=c(list(dat_S22),lapply(unique(dat_S22$Taxonomic_group_response),
                               function(x){subset(dat_S22,Taxonomic_group_response==x)}))
S22=list()
for(i in 1:nrow(Pesticide_species)) {
  S22[[i]]=matrix(NA,length(unique(category[,2]))*length(dat_S22),9)
  for(j in 1:length(unique(category[,2]))) {
    for(k in 1:length(dat_S22)) {
      mydat=subset(dat_S22[[k]],Insecticide.name==Pesticide_species[i,2] & Animal.Latin.name %in% Name)
      if (nrow(mydat)<=1) {
        S22[[i]][(j-1)*length(dat_S22)+k,]=c(t(Pesticide_species[i,1:2]),
                                             unique(category[,2])[j],Nontarget_category[k],
                                             nrow(mydat),length(unique(mydat$Code)),rep(NA,3))
      } else {
        res=rma.mv(yi, vi,
                   data=mydat)
        S22[[i]][(j-1)*length(dat_S22)+k,]=c(t(Pesticide_species[i,1:2]),
                                             unique(category[,2])[j],Nontarget_category[k],
                                             nrow(mydat),length(unique(mydat$Code)),
                                             signif(res$b,4),signif(res$zval,4),signif(res$pval,4))
      }
    }
  }
}
S22=do.call('rbind',S22)
colnames(S22)<-c('Pesticide category','Pesticide species','Non-target plant order',
                 'Plant order response categories',
                 "Obs","Studies","Effect size","T-value","P-value")
S22=subset(as.data.frame(S22), Obs>0)
S22=na.omit(S22)
write.xlsx(S22,"Supplementary Table 21.xls",row.names = F)


##### Supplementary Table 22 ######
category=read.xlsx('microorganism order-2022-1023.xls',sheetName='Sheet1',header=T)

dat_S23=subset(combined,Taxonomic_group=='microorganisms')
dat_S23=c(list(dat_S23),lapply(unique(dat_S23$Taxonomic_group_response),
                               function(x){subset(dat_S23,Taxonomic_group_response==x)}))
S23=list()
for(i in 1:nrow(Pesticide_species)) {
  S23[[i]]=matrix(NA,length(unique(category[,2]))*length(dat_S23),9)
  for(j in 1:length(unique(category[,2]))) {
    for(k in 1:length(dat_S23)) {
      mydat=subset(dat_S23[[k]],Insecticide.name==Pesticide_species[i,2] & Animal.Latin.name %in% Name)
      if (nrow(mydat)<=1) {
        S23[[i]][(j-1)*length(dat_S23)+k,]=c(t(Pesticide_species[i,1:2]),
                                             unique(category[,2])[j],Nontarget_category[k],
                                             nrow(mydat),length(unique(mydat$Code)),rep(NA,3))
      } else {
        res=rma.mv(yi, vi,
                   data=mydat,)
        S23[[i]][(j-1)*length(dat_S23)+k,]=c(t(Pesticide_species[i,1:2]),
                                             unique(category[,2])[j],Nontarget_category[k],
                                             nrow(mydat),length(unique(mydat$Code)),
                                             signif(res$b,4),signif(res$zval,4),signif(res$pval,4))
      }
    }
  }
}
S23=do.call('rbind',S23)
colnames(S23)<-c('Pesticide category','Pesticide species','Non-target microorganism order',
                 'Microorganism order response categories',
                 "Obs","Studies","Effect size","T-value","P-value")
S23=subset(as.data.frame(S23), Obs>0)
S23=na.omit(S23)
write.xlsx(S23,"Supplementary Table 22.xls",row.names = F)


##### Supplementary Table 23 ######
category=read.xlsx('Animal class-2022-1023.xls',sheetName='Sheet1',header=T)

dat_S24=subset(combined,Taxonomic_group=='animals')
dat_S24=c(list(dat_S24),lapply(unique(dat_S24$Taxonomic_group_response),
                               function(x){subset(dat_S24,Taxonomic_group_response==x)}))
S24=list()
for(i in 1:nrow(Pesticide_species)) {
  S24[[i]]=matrix(NA,length(unique(category[,2]))*length(dat_S24),9)
  for(j in 1:length(unique(category[,2]))) {
    for(k in 1:length(dat_S24)) {
      mydat=subset(dat_S24[[k]],Insecticide.name==Pesticide_species[i,2] & Animal.Latin.name %in% Name)
      if (nrow(mydat)<=1) {
        S24[[i]][(j-1)*length(dat_S24)+k,]=c(t(Pesticide_species[i,1:2]),
                                             unique(category[,2])[j],Nontarget_category[k],
                                             nrow(mydat),length(unique(mydat$Code)),rep(NA,3))
      } else {
        res=rma.mv(yi, vi,
                   data=mydat)
        S24[[i]][(j-1)*length(dat_S24)+k,]=c(t(Pesticide_species[i,1:2]),
                                             unique(category[,2])[j],Nontarget_category[k],
                                             nrow(mydat),length(unique(mydat$Code)),
                                             signif(res$b,4),signif(res$zval,4),signif(res$pval,4))
      }
    }
  }
}
S24=do.call('rbind',S24)
colnames(S24)<-c('Pesticide category','Pesticide species','Non-target animal class',
                 'Animal class response categories',
                 "Obs","Studies","Effect size","T-value","P-value")
S24=subset(as.data.frame(S24), Obs>0)
S24=na.omit(S24)
write.xlsx(S24,"Supplementary Table 23.xls",row.names = F)


##### Supplementary Table 24 ######
category=read.xlsx('plant class-2022-1023.xls',sheetName='Sheet1',header=T)

dat_S25=subset(combined,Taxonomic_group=='plants')
dat_S25=c(list(dat_S25),lapply(unique(dat_S25$Taxonomic_group_response),
                               function(x){subset(dat_S25,Taxonomic_group_response==x)}))
S25=list()
for(i in 1:nrow(Pesticide_species)) {
  S25[[i]]=matrix(NA,length(unique(category[,2]))*length(dat_S25),9)
  for(j in 1:length(unique(category[,2]))) {
    for(k in 1:length(dat_S25)) {
      mydat=subset(dat_S25[[k]],Insecticide.name==Pesticide_species[i,2] & Animal.Latin.name %in% Name)
      if (nrow(mydat)<=1) {
        S25[[i]][(j-1)*length(dat_S25)+k,]=c(t(Pesticide_species[i,1:2]),
                                             unique(category[,2])[j],Nontarget_category[k],
                                             nrow(mydat),length(unique(mydat$Code)),rep(NA,3))
      } else {
        res=rma.mv(yi, vi,
                   data=mydat)
        S25[[i]][(j-1)*length(dat_S25)+k,]=c(t(Pesticide_species[i,1:2]),
                                             unique(category[,2])[j],Nontarget_category[k],
                                             nrow(mydat),length(unique(mydat$Code)),
                                             signif(res$b,4),signif(res$zval,4),signif(res$pval,4))
      }
    }
  }
}
S25=do.call('rbind',S25)
colnames(S25)<-c('Pesticide category','Pesticide species','Non-target plant class',
                 'Plant class response categories',
                 "Obs","Studies","Effect size","T-value","P-value")
S25=subset(as.data.frame(S25), Obs>0)
S25=na.omit(S25)
write.xlsx(S25,"Supplementary Table 24.xls",row.names = F)


##### Supplementary Table 25 ######
category=read.xlsx('microorganism class-2022-1023.xls',sheetName='Sheet1',header=T)

dat_S26=subset(combined,Taxonomic_group=='microorganisms')
dat_S26=c(list(dat_S26),lapply(unique(dat_S26$Taxonomic_group_response),
                               function(x){subset(dat_S26,Taxonomic_group_response==x)}))
S26=list()
for(i in 1:nrow(Pesticide_species)) {
  S26[[i]]=matrix(NA,length(unique(category[,2]))*length(dat_S26),9)
  for(j in 1:length(unique(category[,2]))) {
    for(k in 1:length(dat_S26)) {
      mydat=subset(dat_S26[[k]],Insecticide.name==Pesticide_species[i,2] & Animal.Latin.name %in% Name)
      if (nrow(mydat)<=1) {
        S26[[i]][(j-1)*length(dat_S26)+k,]=c(t(Pesticide_species[i,1:2]),
                                             unique(category[,2])[j],Nontarget_category[k],
                                             nrow(mydat),length(unique(mydat$Code)),rep(NA,3))
      } else {
        res=rma.mv(yi, vi,
                   data=mydat)
        S26[[i]][(j-1)*length(dat_S26)+k,]=c(t(Pesticide_species[i,1:2]),
                                             unique(category[,2])[j],Nontarget_category[k],
                                             nrow(mydat),length(unique(mydat$Code)),
                                             signif(res$b,4),signif(res$zval,4),signif(res$pval,4))
      }
    }
  }
}
S26=do.call('rbind',S26)
colnames(S26)<-c('Pesticide category','Pesticide species','Non-target microorganism class',
                 'Microorganism class response categories',
                 "Obs","Studies","Effect size","T-value","P-value")
S26=subset(as.data.frame(S26), Obs>0)
S26=na.omit(S26)
write.xlsx(S26,"Supplementary Table 25.xls",row.names = F)


##### Supplementary Table 26 ######
category=read.xlsx('Supplementary Table 27-effect sizes of each pesticide on animal phylum/Animal phylum-2022-1023.xls',sheetName='Sheet1',header=T)

dat_S27=subset(combined,Taxonomic_group=='animals')
dat_S27=c(list(dat_S27),lapply(unique(dat_S27$Taxonomic_group_response),
                               function(x){subset(dat_S27,Taxonomic_group_response==x)}))
S27=list()
for(i in 1:nrow(Pesticide_species)) {
  S27[[i]]=matrix(NA,length(unique(category[,2]))*length(dat_S27),9)
  for(j in 1:length(unique(category[,2]))) {
    for(k in 1:length(dat_S27)) {
      mydat=subset(dat_S27[[k]],Insecticide.name==Pesticide_species[i,2] & Animal.Latin.name %in% Name)
      if (nrow(mydat)<=1) {
        S27[[i]][(j-1)*length(dat_S27)+k,]=c(t(Pesticide_species[i,1:2]),
                                             unique(category[,2])[j],Nontarget_category[k],
                                             nrow(mydat),length(unique(mydat$Code)),rep(NA,3))
      } else {
        res=rma.mv(yi, vi,
                   data=mydat)
        S27[[i]][(j-1)*length(dat_S27)+k,]=c(t(Pesticide_species[i,1:2]),
                                             unique(category[,2])[j],Nontarget_category[k],
                                             nrow(mydat),length(unique(mydat$Code)),
                                             signif(res$b,4),signif(res$zval,4),signif(res$pval,4))
      }
    }
  }
}
S27=do.call('rbind',S27)
colnames(S27)<-c('Pesticide category','Pesticide species','Non-target animal phylum',
                 'Animal phylum response categories',
                 "Obs","Studies","Effect size","T-value","P-value")
S27=subset(as.data.frame(S27), Obs>0)
S27=na.omit(S27)
write.xlsx(S27,"Supplementary Table 26.xls",row.names = F)


##### Supplementary Table 27 ######
category=read.xlsx('Supplementary Table 28-effect sizes of each pesticide on plant phylum/plant phylum-2022-1023.xls',sheetName='Sheet1',header=T)

dat_S28=subset(combined,Taxonomic_group=='plants')
dat_S28=c(list(dat_S28),lapply(unique(dat_S28$Taxonomic_group_response),
                               function(x){subset(dat_S28,Taxonomic_group_response==x)}))
S28=list()
for(i in 1:nrow(Pesticide_species)) {
  S28[[i]]=matrix(NA,length(unique(category[,2]))*length(dat_S28),9)
  for(j in 1:length(unique(category[,2]))) {
    for(k in 1:length(dat_S28)) {
      mydat=subset(dat_S28[[k]],Insecticide.name==Pesticide_species[i,2] & Animal.Latin.name %in% Name)
      if (nrow(mydat)<=1) {
        S28[[i]][(j-1)*length(dat_S28)+k,]=c(t(Pesticide_species[i,1:2]),
                                             unique(category[,2])[j],Nontarget_category[k],
                                             nrow(mydat),length(unique(mydat$Code)),rep(NA,3))
      } else {
        res=rma.mv(yi, vi,
                   data=mydat)
        S28[[i]][(j-1)*length(dat_S28)+k,]=c(t(Pesticide_species[i,1:2]),
                                             unique(category[,2])[j],Nontarget_category[k],
                                             nrow(mydat),length(unique(mydat$Code)),
                                             signif(res$b,4),signif(res$zval,4),signif(res$pval,4))
      }
    }
  }
}
S28=do.call('rbind',S28)
colnames(S28)<-c('Pesticide category','Pesticide species','Non-target plant phylum',
                 'Plant phylum response categories',
                 "Obs","Studies","Effect size","T-value","P-value")
S28=subset(as.data.frame(S28), Obs>0)
S28=na.omit(S28)
write.xlsx(S28,"Supplementary Table 27.xls",row.names = F)


##### Supplementary Table 28 ######
category=read.xlsx('microorganism phylum-2022-1023.xls',sheetName='Sheet1',header=T)

dat_S29=subset(combined,Taxonomic_group=='microorganisms')
dat_S29=c(list(dat_S29),lapply(unique(dat_S29$Taxonomic_group_response),
                               function(x){subset(dat_S29,Taxonomic_group_response==x)}))
S29=list()
for(i in 1:nrow(Pesticide_species)) {
  S29[[i]]=matrix(NA,length(unique(category[,2]))*length(dat_S29),9)
  for(j in 1:length(unique(category[,2]))) {
    for(k in 1:length(dat_S29)) {
      mydat=subset(dat_S29[[k]],Insecticide.name==Pesticide_species[i,2] & Animal.Latin.name %in% Name)
      if (nrow(mydat)<=1) {
        S29[[i]][(j-1)*length(dat_S29)+k,]=c(t(Pesticide_species[i,1:2]),
                                             unique(category[,2])[j],Nontarget_category[k],
                                             nrow(mydat),length(unique(mydat$Code)),rep(NA,3))
      } else {
        res=rma.mv(yi, vi,
                   data=mydat)
        S29[[i]][(j-1)*length(dat_S29)+k,]=c(t(Pesticide_species[i,1:2]),
                                             unique(category[,2])[j],Nontarget_category[k],
                                             nrow(mydat),length(unique(mydat$Code)),
                                             signif(res$b,4),signif(res$zval,4),signif(res$pval,4))
      }
    }
  }
}
S29=do.call('rbind',S29)
colnames(S29)<-c('Pesticide category','Pesticide species','Non-target microorganism phylum',
                 'Microorganism phylum response categories',
                 "Obs","Studies","Effect size","T-value","P-value")
S29=subset(as.data.frame(S29), Obs>0)
S29=na.omit(S29)
write.xlsx(S29,"Supplementary Table 28.xls",row.names = F)


##### Supplementary Table 29 ######

S30=list()
for(i in 1:nrow(Pesticide_species)) {
  S30[[i]]=matrix(NA,nrow(Non_target_animal)*length(dat_S15),14)
  for(j in 1:nrow(Non_target_animal)) {
    for(k in 1:length(dat_S15)) {
      mydat=na.omit(mydat)
      if (length(unique(mydat)$xi)<=2) {
        S30[[i]][(j-1)*length(dat_S15)+k,]=c(t(Pesticide_species[i,1:2]),
                                             Non_target_animal[j,],Nontarget_category[k],
                                             nrow(mydat),length(unique(mydat$Code)),rep(NA,8))
      } else {
        x <- mydat[,'xi']
        y <- mydat[,'yi']
        mydata=data.frame(x=x,y=y,Code=mydat[,'Code'],
                          pesticide_by_target_organisms=mydat[,'pesticide_by_target_organisms'],
                          Insecticide.name=mydat[,'Insecticide.name'],vi=mydat[,'vi'])
        assign("mydata2", mydata, envir=.GlobalEnv) 
        fit<-try(lme(y~x,
                     random =list(~ 1|Code),
                     weights = varFixed(~ vi),data=mydata2,control = lmeControl(sigma = 1)))
          parameter=signif(summary(fit)$tTable[2,-1],4)
          Intercept_x=round(summary(fit)$tTable[,1],4)
          con95=round(intervals(fit, 'fixed', level=0.95)[[1]][2,c(1,3)],4)
          if (Intercept_x[1]<0) {
            equation=paste0('Y = ',Intercept_x[2],'X - ',abs(Intercept_x[1]))
          } else {
            equation=paste0('Y = ',Intercept_x[2],'X + ',abs(Intercept_x[1]))
          }
          Threshhold=round(-Intercept_x[1]/Intercept_x[2],4)
          S30[[i]][(j-1)*length(dat_S15)+k,]=c(t(Pesticide_species[i,1:2]),
                                               Non_target_animal[j,],Nontarget_category[k],
                                               nrow(mydat),length(unique(mydat$Code)),
                                               equation,con95,parameter,Threshhold)
      }
    }
  }
}
S30=do.call('rbind',S30)
colnames(S30)<-c('Pesticide category','Pesticide species','Non-target animal species',
                 'Effect size of non-target animal species category (Y)',
                 "Obs","Studies","Regression equation","95%CL","95%CU",
                 'Std. Error','DF','T-value','P-value','Threshhold application rate over the control')
S30=subset(as.data.frame(S30), Obs>0)
S30=na.omit(S30)
write.xlsx(S30,"Supplementary Table 29.xls",row.names = F)


##### Supplementary Table 30 ######
Non_target_plant=read.xlsx('Latin names of specific plant species-2022-1020.xls',sheetName='specific plant species',header=F)

dat_S31=subset(combined,Taxonomic_group=='plants')
dat_S31=c(list(dat_S31),lapply(unique(dat_S31$Taxonomic_group_response),
                               function(x){subset(dat_S31,Taxonomic_group_response==x)}))
S31=list()
for(i in 1:nrow(Pesticide_species)) {
  S31[[i]]=matrix(NA,nrow(Non_target_plant)*length(dat_S31),14)
  for(j in 1:nrow(Non_target_plant)) {
    for(k in 1:length(dat_S31)) {
      mydat=na.omit(mydat)
      if (length(unique(mydat)$xi)<=2) {
        S31[[i]][(j-1)*length(dat_S31)+k,]=c(t(Pesticide_species[i,1:2]),
                                             Non_target_plant[j,],Nontarget_category[k],
                                             nrow(mydat),length(unique(mydat$Code)),rep(NA,8))
      } else {
        x <- mydat[,'xi']
        y <- mydat[,'yi']
        mydata=data.frame(x=x,y=y,Code=mydat[,'Code'],
                          pesticide_by_target_organisms=mydat[,'pesticide_by_target_organisms'],
                          Insecticide.name=mydat[,'Insecticide.name'],vi=mydat[,'vi'])
        assign("mydata2", mydata, envir=.GlobalEnv) 
        fit<-try(lme(y~x,
                     random =list(~ 1|Code),
                     weights = varFixed(~ vi),data=mydata2,control = lmeControl(sigma = 1))
          parameter=signif(summary(fit)$tTable[2,-1],4)
          Intercept_x=round(summary(fit)$tTable[,1],4)
          con95=round(intervals(fit, 'fixed', level=0.95)[[1]][2,c(1,3)],4)
          if (Intercept_x[1]<0) {
            equation=paste0('Y = ',Intercept_x[2],'X - ',abs(Intercept_x[1]))
          } else {
            equation=paste0('Y = ',Intercept_x[2],'X + ',abs(Intercept_x[1]))
          }
          Threshhold=round(-Intercept_x[1]/Intercept_x[2],4)
          S31[[i]][(j-1)*length(dat_S31)+k,]=c(t(Pesticide_species[i,1:2]),
                                               Non_target_plant[j,],Nontarget_category[k],
                                               nrow(mydat),length(unique(mydat$Code)),
                                               equation,con95,parameter,Threshhold)
      }
    }
  }
}
S31=do.call('rbind',S31)
colnames(S31)<-c('Pesticide category','Pesticide species','Non-target animal species',
                 'Effect size of non-target animal species category (Y)',
                 "Obs","Studies","Regression equation","95%CL","95%CU",
                 'Std. Error','DF','T-value','P-value','Threshhold application rate over the control')
S31=subset(as.data.frame(S31), Obs>0)
S31=na.omit(S31)
write.xlsx(S31,"Supplementary Table 30.xls",row.names = F)


##### Supplementary Table 31 ######
Non_target_microorganism=read.xlsx('Latin names of microorganism species-2022-1023.xls',sheetName='Sheet1',header=F)

dat_S32=subset(combined,Taxonomic_group=='microorganisms')
dat_S32=c(list(dat_S32),lapply(unique(dat_S32$Taxonomic_group_response),
                               function(x){subset(dat_S32,Taxonomic_group_response==x)}))
S32=list()
for(i in 1:nrow(Pesticide_species)) {
  S32[[i]]=matrix(NA,nrow(Non_target_microorganism)*length(dat_S32),14)
  for(j in 1:nrow(Non_target_microorganism)) {
    for(k in 1:length(dat_S32)) {
      mydat=na.omit(mydat)
      if (length(unique(mydat)$xi)<=2) {
        S32[[i]][(j-1)*length(dat_S32)+k,]=c(t(Pesticide_species[i,1:2]),
                                             Non_target_microorganism[j,],Nontarget_category[k],
                                             nrow(mydat),length(unique(mydat$Code)),rep(NA,8))
      } else {
        x <- mydat[,'xi']
        y <- mydat[,'yi']
        mydata=data.frame(x=x,y=y,Code=mydat[,'Code'],
                          pesticide_by_target_organisms=mydat[,'pesticide_by_target_organisms'],
                          Insecticide.name=mydat[,'Insecticide.name'],vi=mydat[,'vi'])
        assign("mydata2", mydata, envir=.GlobalEnv) 
        fit<-try(lme(y~x,
                     random =list(~ 1|Code),
                     weights = varFixed(~ vi),data=mydata2,control = lmeControl(sigma = 1)))
          parameter=signif(summary(fit)$tTable[2,-1],4)
          Intercept_x=round(summary(fit)$tTable[,1],4)
          con95=round(intervals(fit, 'fixed', level=0.95)[[1]][2,c(1,3)],4)
          if (Intercept_x[1]<0) {
            equation=paste0('Y = ',Intercept_x[2],'X - ',abs(Intercept_x[1]))
          } else {
            equation=paste0('Y = ',Intercept_x[2],'X + ',abs(Intercept_x[1]))
          }
          Threshhold=round(-Intercept_x[1]/Intercept_x[2],4)
          S32[[i]][(j-1)*length(dat_S32)+k,]=c(t(Pesticide_species[i,1:2]),
                                               Non_target_microorganism[j,],Nontarget_category[k],
                                               nrow(mydat),length(unique(mydat$Code)),
                                               equation,con95,parameter,Threshhold)
      }
    }
  }
}
S32=do.call('rbind',S32)
colnames(S32)<-c('Pesticide category','Pesticide species','Non-target animal species',
                 'Effect size of non-target animal species category (Y)',
                 "Obs","Studies","Regression equation","95%CL","95%CU",
                 'Std. Error','DF','T-value','P-value','Threshhold application rate over the control')
S32=subset(as.data.frame(S32), Obs>0)
S32=na.omit(S32)
write.xlsx(S32,"Supplementary Table 31.xls",row.names = F)










