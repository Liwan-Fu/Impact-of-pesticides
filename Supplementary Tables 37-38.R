library(metafor)
library(xlsx)
setwd('...')
combined=read.csv('Meta combined data.csv',header=T)

Old_new_pesticide=read.xlsx('Supplementary Data 2-Classification of old and new pesticides-2023-0505.xls',sheetName='classification',header=T)


Model_animal=read.xlsx('Supplementary Data 3-Classification of model and nonmodel species-2023-0505.xls',sheetName='animals',header=T)
Model_plant=read.xlsx('Supplementary Data 3-Classification of model and nonmodel species-2023-0505.xls',sheetName='plants',header=T)
Model_micro=read.xlsx('Supplementary Data 3-Classification of model and nonmodel species-2023-0505.xls',sheetName='microorganisms',header=T)
Model_comb=rbind(Model_animal[,c(1,6)],Model_plant[,c(1,6)],Model_micro[,c(1,6)])

combined[,'Model-non']=sapply(combined$Latin.name,function(x) {
  index=which(x==Model_comb$Species.or.other.groups)[1]
  Model_comb[index,'Model.or.non.model.species']
})

########## Supplementary Table 37 ###########
Pesticide_category=c('Pesticide',
                     unique(combined$pesticide_by_target_organisms))
dat_list1=c(lapply(unique(combined$Taxonmic_group),
                   function(x){subset(combined,Taxonomic_group==x)}),
            lapply(unique(combined$Taxonoic_group_response),
                   function(x){subset(combined,Taxnomic_group_response==x)}))
Experiment=unique(combined$`Old-new-pesticide`)
S37=matrix(NA,length(Experiment)*length(Pesticide_category)*length(dat_list1),8)
for (k in 1:length(Experiment)) {
  for (i in 1:length(Pesticide_category)) {
    for (j in 1:length(dat_list1)){
        mydat=subset(dat_list1[[j]],
                     `Old-new-pesticide`==Experiment[k]&pesticide_by_target_organisms==Pesticide_category[i])
      if (nrow(mydat)<=1) {
        S37[(k-1)*52+(i-1)*13+j,]=c(nrow(mydat),length(unique(mydat$Code)),
                                   rep(NA,6))
      } else {
        res=rma.mv(yi, vi, 
                   random = ~ 1 | factor(Insecticide.name),
                   data=mydat)
        S37[(k-1)*52+(i-1)*13+j,]=c(nrow(mydat),length(unique(mydat$Code)),
                                   res$b,res$zval,nrow(mydat)-1,res$pval,res$ci.lb,res$ci.ub)
      }
    }
  }
}
colnames(S37)<-c("#Obs","#Studies","Effect size","T-value","Df",
                "P-value","CL.lb","CL.ub")
write.xlsx(S37,"Supplementary Table 37.xls")

########## Supplementary Table 38 ###########
Pesticide_category=c('Pesticide',
                     unique(combined$pesticide_by_target_organisms))
dat_list1=c(lapply(unique(combined$Taxonomic_group),
                   function(x){subset(combind,Taonomic_group==x)}),
            lapply(unique(combined$Taxonomic_group_response),
                   function(x){subset(combined,Taxomic_group_response==x)}))
Experiment=unique(combined$`Model-non`)
S38=matrix(NA,length(Experiment)*length(Pesticide_category)*length(dat_list1),8)
for (k in 1:length(Experiment)) {
  for (i in 1:length(Pesticide_category)) {
    for (j in 1:length(dat_list1)){
        mydat=subset(dat_list1[[j]],
                     `Model-non`==Experiment[k]&pesticide_by_target_organisms==Pesticide_category[i])
      if (nrow(mydat)<=1) {
        S38[(k-1)*52+(i-1)*13+j,]=c(nrow(mydat),length(unique(mydat$Code)),
                                    rep(NA,6))
      } else {
        res=rma.mv(yi, vi, 
                   random = ~ 1 | factor(Insecticide.name),
                   data=mydat)
        S38[(k-1)*52+(i-1)*13+j,]=c(nrow(mydat),length(unique(mydat$Code)),
                                    res$b,res$zval,nrow(mydat)-1,res$pval,res$ci.lb,res$ci.ub)
      }
    }
  }
}
colnames(S38)<-c("#Obs","#Studies","Effect size","T-value","Df",
                 "P-value","CL.lb","CL.ub")
write.xlsx(S38,"Supplementary Table 38.xls")















