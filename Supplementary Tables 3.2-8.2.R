library(metafor)
library(xlsx)
setwd('...')
combined=read.csv('Meta combined data.csv',header=T)

library(rgbif)
library(V.PhyloMaker)
species=read.table("non-target organism species-2023-0515_combined.txt",h=T,sep="\t")
plants=subset(species,type=="plants" & species %in% unique(combined$Animal.Latin.name))
plants.taxa=name_backbone_checklist(plants$species)
a2=as.data.frame(plants.taxa)
a3=merge(plants,a2,by.x="species",by.y="verbatim_name")
mycor=phylo.maker(data.frame(species=a3$species,genus=a3$genus,
                                           family=a3$family))
mytree <- compute.brlen(mycor$scenario.3)
A <- vcv(mytree, corr=TRUE)


rm(species,plants,plants.taxa,a2,a3,mytree,A)

myphylo=function(Total.data2){
  species=subset(mycor[["species.list"]],species %in% unique(Total.data2$Animal.Latin.name))
  mycor1=phylo.maker(species)
  mytree <- compute.brlen(mycor1$scenario.3)
  A <- vcv(mytree, corr=TRUE)
  A2=round(A,1)
  levels(Total.data2$Plant.species.new)=sort(dimnames(A)[[1]])
  Total.data2$Plant.species.new.p=Total.data2$Plant.species.new
  return(list(Total.data2,A2))
}
# combined=na.omit(combined)
#  combined=combined[sample(nrow(combined))[1:500],]

########## Supplementary Table 3.2 ###########
S4=matrix(NA,12,8)
Pesticide_category=c('Pesticide',
                     unique(combined$pesticide_by_target_organisms))
dat_list1=lapply(unique(combined$Taxonomic_group_response),
                   function(x){subset(combined,Taxonomic_group_response==x)})
for (i in 1:length(Pesticide_category)) {
  for (j in 1:length(dat_list1)){
      mydat=subset(dat_list1[[j]],
                   pesticide_by_target_organisms==Pesticide_category[i])
      mydat=myphylo(mydat)
      res=rma.mv(yi, vi, 
                 random =list( ~ 1 | factor(Insecticide.name),
                               ~ 1 | Code,
                               ~ 1 | Publication_year,
                               ~ 1 | Plant.species.new,
                               ~ 1|Plant.species.new.p),
                 R=list(Plant.species.new.p=mydat[[2]]),
                 data=mydat[[1]])
    S4[(i-1)*length(dat_list1)+j,]=c(nrow(mydat[[1]]),length(unique(mydat[[1]]$Code)),
                                     res$b,res$zval,nrow(mydat[[1]])-1,
                                     res$pval,res$ci.lb,res$ci.ub)
  }
}
S4<-as.data.frame(S4,row.names<-Rnames)
colnames(S4)<-c("#Obs","#Studies","Effect size","T-value","Df",
                "P-value","CL.lb","CL.ub")
write.xlsx(S4,"Supplementary Table 3.2.xls")

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
        mydat=myphylo(mydat)
        res=rma.mv(yi, vi, 
                   random =list( ~ 1 | factor(Insecticide.name),
                                 ~ 1 | Code,
                                 ~ 1 | Publication_year,
                                 ~ 1 | Plant.species.new,
                                 ~ 1|Plant.species.new.p),
                   R=list(Plant.species.new.p=mydat[[2]]),
                   data=mydat[[1]])
      S5[(i-1)*length(dat_list3)+j,]=c(nrow(mydat[[1]]),length(unique(mydat[[1]]$Code)),
                                       res$b,res$zval,nrow(mydat[[1]])-1,
                                       res$pval,res$ci.lb,res$ci.ub)
    }
  }
}
colnames(S5)<-c("#Obs","#Studies","Effect size","T-value","Df",
                "P-value","CL.lb","CL.ub")
write.xlsx(S5,"Supplementary Table 4.2.xls")

########## Supplementary Table 5.2 ###########
dat_list1=lapply(unique(combined$Taxonomic_group_response),
                   function(x){subset(combined,Taxonomic_group_response==x)})
S6_category=c(unique(combined$pesticide_by_source),
              unique(combined$Primary.classification.of.insecticide........一级分类))
S6=matrix(NA,length(S6_category)*length(dat_list1),8)
for (i in 1:length(S6_category)) {
  for (j in 1:length(dat_list1)){
    if (i<=3) {
      mydat=subset(dat_list1[[j]],
                   pesticide_by_source==S6_category[i])
      if (nrow(mydat)<=1) {
        S6[(i-1)*length(dat_list1)+j,]=c(nrow(mydat),length(unique(mydat$Code)),
                                         rep(NA,6))
      } else {        
          mydat=myphylo(mydat)
          res=rma.mv(yi, vi, 
                     random =list( ~ 1 | factor(Insecticide.name),
                                   ~ 1 | Code,
                                   ~ 1 | Publication_year,
                                   ~ 1 | Plant.species.new,
                                   ~ 1|Plant.species.new.p),
                     R=list(Plant.species.new.p=mydat[[2]]),
                     data=mydat[[1]])
        S6[(i-1)*length(dat_list1)+j,]=c(nrow(mydat[[1]]),length(unique(mydat[[1]]$Code)),
                                         res$b,res$zval,nrow(mydat[[1]])-1,
                                         res$pval,res$ci.lb,res$ci.ub)
      }
    } else {
      mydat=subset(dat_list1[[j]],
                   Primary.classification.of.insecticide........一级分类==S6_category[i])
      if (nrow(mydat)<=1) {
        S6[(i-1)*length(dat_list1)+j,]=c(nrow(mydat),length(unique(mydat$Code)),
                                         rep(NA,6))
      } else {
          mydat=myphylo(mydat)
          res=rma.mv(yi, vi, 
                     random =list( ~ 1 | factor(Insecticide.name),
                                   ~ 1 | Code,
                                   ~ 1 | Publication_year,
                                   ~ 1 | Plant.species.new,
                                   ~ 1|Plant.species.new.p),
                     R=list(Plant.species.new.p=mydat[[2]]),
                     data=mydat[[1]])
        S6[(i-1)*length(dat_list1)+j,]=c(nrow(mydat[[1]]),length(unique(mydat[[1]]$Code)),
                                         res$b,res$zval,nrow(mydat[[1]])-1,
                                         res$pval,res$ci.lb,res$ci.ub)
      }
    }
  }
}
colnames(S6)<-c("#Obs","#Studies","Effect size","T-value","Df",
                "P-value","CL.lb","CL.ub")
write.xlsx(S6,"Supplementary Table 5.2.xls")

########## Supplementary Table 6.2 ###########
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
        S7[(k-1)*length(Pesticide_category)*length(dat_list1)+(i-1)*length(dat_list1)+j,]=c(nrow(mydat),length(unique(mydat$Code)),rep(NA,6))
      } else {
          mydat=myphylo(mydat)
          res=rma.mv(yi, vi, 
                     random =list( ~ 1 | factor(Insecticide.name),
                                   ~ 1 | Code,
                                   ~ 1 | Publication_year,
                                   ~ 1 | Plant.species.new,
                                   ~ 1|Plant.species.new.p),
                     R=list(Plant.species.new.p=mydat[[2]]),
                     data=mydat[[1]])
        }
        S7[(k-1)*length(Pesticide_category)*length(dat_list1)+(i-1)*length(dat_list1)+j,]=c(nrow(mydat[[1]]),length(unique(mydat[[1]]$Code)),
                                                                                            res$b,res$zval,nrow(mydat[[1]])-1,
                                                                                            res$pval,res$ci.lb,res$ci.ub)
    }
  }
}

colnames(S7)<-c("#Obs","#Studies","Effect size","T-value","Df",
                "P-value","CL.lb","CL.ub")
write.xlsx(S7,"Supplementary Table 6.2.xls")

########## Supplementary Table 7.2 ###########
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
        S8[(k-1)*length(Pesticide_category)*length(dat_list1)+(i-1)*length(dat_list1)+j,]=c(nrow(mydat),length(unique(mydat$Code)),
                                   rep(NA,6))
      } else {
          mydat=myphylo(mydat)
          res=rma.mv(yi, vi, 
                     random =list( ~ 1 | factor(Insecticide.name),
                                   ~ 1 | Code,
                                   ~ 1 | Publication_year,
                                   ~ 1 | Plant.species.new,
                                   ~ 1|Plant.species.new.p),
                     R=list(Plant.species.new.p=mydat[[2]]),
                     data=mydat[[1]])
        S8[(k-1)*length(Pesticide_category)*length(dat_list1)+(i-1)*length(dat_list1)+j,]=c(nrow(mydat[[1]]),length(unique(mydat[[1]]$Code)),
                                                                                            res$b,res$zval,nrow(mydat[[1]])-1,
                                                                                            res$pval,res$ci.lb,res$ci.ub)
      }
    }
  }
}

S8<-as.data.frame(S8,row.names<-Rnames)
colnames(S8)<-c("#Obs","#Studies","Effect size","T-value","Df",
                "P-value","CL.lb","CL.ub")
write.xlsx(S8,"Supplementary Table 7.2.xls")

########## Supplementary Table 8.2 ###########
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
        S9[(k-1)*length(Pesticide_category)*length(dat_list1)+(i-1)*length(dat_list1)+j,]=c(nrow(mydat),length(unique(mydat$Code)),
                                   rep(NA,6))
      } else {
          mydat=myphylo(mydat)
          res=rma.mv(yi, vi, 
                     random =list( ~ 1 | factor(Insecticide.name),
                                   ~ 1 | Code,
                                   ~ 1 | Publication_year,
                                   ~ 1 | Plant.species.new,
                                   ~ 1|Plant.species.new.p),
                     R=list(Plant.species.new.p=mydat[[2]]),
                     data=mydat[[1]])
        S9[(k-1)*length(Pesticide_category)*length(dat_list1)+(i-1)*length(dat_list1)+j,]=c(nrow(mydat[[1]]),length(unique(mydat[[1]]$Code)),
                                                                                            res$b,res$zval,nrow(mydat[[1]])-1,
                                                                                            res$pval,res$ci.lb,res$ci.ub)
      }
    }
  }
}
colnames(S9)<-c("#Obs","#Studies","Effect size","T-value","Df",
                "P-value","CL.lb","CL.ub")
write.xlsx(S9,"Supplementary Table 8.2.xls")

