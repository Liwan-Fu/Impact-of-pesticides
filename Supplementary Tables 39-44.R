library(metafor)
library(xlsx)
setwd('...')
combined=read.csv('data for recomendation dose of aquatic.csv',header=T)


########## Supplementary Table 39 ###########
S32=matrix(NA,52,8)
Pesticide_category=c('Pesticide',
                     unique(combined$pestcide_by_target_organisms))
dat_list1=c(lapply(unique(combined$Taxonomic_group),
                   function(x){subset(combined,Taonomic_group==x)}),
            lapply(unique(combined$Taxonomic_group_response),
                   function(x){subset(combined,Taxonoic_group_response==x)}))
for (i in 1:length(Pesticide_category)) {
  for (j in 1:length(dat_list1)){
      mydat=subset(dat_list1[[j]],
                   pesticide_by_target_organisms==Pesticide_category[i])
    if (nrow(mydat)<=1) {
      S32[(i-1)*13+j,]=c(nrow(mydat),length(unique(mydat$Code)),
                         rep(NA,6))
    } else {
      res=rma.mv(yi, vi, 
                 random = ~ 1 | factor(Insecticide.name),
                 data=mydat)
      S32[(i-1)*13+j,]=c(nrow(mydat),length(unique(mydat$Code)),
                         res$b,res$zval,nrow(mydat)-1,res$pval,res$ci.lb,res$ci.ub)
    }
  }
}
colnames(S32)<-c("#Obs","#Studies","Effect size","T-value","Df",
                 "P-value","CL.lb","CL.ub")
write.xlsx(S32,"Supplementary Table 39.xls")

########## Supplementary Table 40 ###########
dat_list2=c(lapply(unique(combined$Primary.classification.of.animal......一级分类),
                   function(x){
                     subset(combined,Primary.classifcation.of.animal......一级分类==x)
                   }))
dat_list3=c(dat_list2[1],lapply(unique(dat_list2[[1]]$Taxonomic_group_response),
                                function(x) {
                                  subset(dat_list2[[1]],Taxonomic_group_response==x)
                                }),
            dat_list2[2],lapply(unique(dat_list2[[2]]$Taxnomic_group_response),
                                function(x) {
                                  subset(dat_list2[[2]],Taxonomic_group_response==x)
                                }),
            dat_list2[3],lapply(unique(dat_list2[[3]]$Taxonomic_group_response),
                                function(x) {
                                  subset(dat_list2[[3]],Taxoomic_group_response==x)
                                }),
            dat_list2[4],lapply(unique(dat_list2[[4]]$Taxonomic_group_response),
                                function(x) {
                                  subset(dat_list2[[4]],Taxonomic_group_response==x)
                                }),
            dat_list2[5],lapply(unique(dat_list2[[5]]$Taxonomc_group_response),
                                function(x) {
                                  subset(dat_list2[[5]],Taxonomic_group_response==x)
                                }),
            dat_list2[6],lapply(unique(dat_list2[[6]]$Taxonoic_group_response),
                                function(x) {
                                  subset(dat_list2[[6]],Taxnomic_group_response==x)
                                }))
S33=matrix(NA,96,8)
for (i in 1:length(Pesticide_category)) {
  for (j in 1:length(dat_list3)){
      mydat=subset(dat_list3[[j]],
                   pesticide_by_target_organisms==Pesticide_category[i])
    if (nrow(mydat)<=1) {
      S33[(i-1)*24+j,]=c(nrow(mydat),length(unique(mydat$Code)),
                        rep(NA,6))
    } else {
      res=rma.mv(yi, vi, 
                 random = ~ 1 | factor(Insecticide.name),
                 data=mydat)
      S33[(i-1)*24+j,]=c(nrow(mydat),length(unique(mydat$Code)),
                         res$b,res$zval,nrow(mydat)-1,res$pval,res$ci.lb,res$ci.ub)
    }
  }
}
S33<-as.data.frame(S33,row.names<-Rnames)
colnames(S33)<-c("#Obs","#Studies","Effect size","T-value","Df",
                 "P-value","CL.lb","CL.ub")
write.xlsx(S33,"Supplementary Table 40.xls")

########## Supplementary Table 41 ###########
Pesticide_category=c('Pesticide',
                     unique(combined$pesticie_by_target_organisms))
dat_list1=c(lapply(unique(combined$Taxonomic_group),
                   function(x){subset(combind,Taxonomi_group==x)}),
            lapply(unique(combined$Taxonomic_group_response),
                   function(x){subset(combined,Taxonomic_group_response==x)}))
S6_category=c(unique(combined$pesticide_by_source),
              unique(combined$Primary.classfication.of.insecticide........一级分类))[c(1,3,2,4,9,5,6,8,7)]
S34=matrix(NA,length(S6_category)*length(dat_list1),8)
for (i in 1:length(S6_category)) {
  for (j in 1:length(dat_list1)){
      mydat=subset(dat_list1[[j]],
                   pesticide_by_source==S6_category[i])
      if (nrow(mydat)<=1) {
        S34[(i-1)*13+j,]=c(nrow(mydat),length(unique(mydat$Code)),
                          rep(NA,6))
      } else {
        res=rma.mv(yi, vi, 
                   random = ~ 1 | factor(Insecticide.name),
                   data=mydat)
        S34[(i-1)*13+j,]=c(nrow(mydat),length(unique(mydat$Code)),
                           res$b,res$zval,nrow(mydat)-1,res$pval,res$ci.lb,res$ci.ub)
      }

    } else {
      mydat=subset(dat_list1[[j]],
                   Primary.classification.of.insecticide........一级分类==S6_category[i])
      if (nrow(mydat)<=1) {
        S34[(i-1)*13+j,]=c(nrow(mydat),length(unique(mydat$Code)),
                          rep(NA,6))
      } else {
        res=rma.mv(yi, vi, 
                   random = ~ 1 | factor(Insecticide.name),
                   data=mydat)
        S34[(i-1)*13+j,]=c(nrow(mydat),length(unique(mydat$Code)),
                           res$b,res$zval,nrow(mydat)-1,res$pval,res$ci.lb,res$ci.ub)
      }
  }
}
colnames(S34)<-c("#Obs","#Studies","Effect size","T-value","Df",
                 "P-value","CL.lb","CL.ub")
write.xlsx(S34,"Supplementary Table 41.xls")


########## Supplementary Table 42 ###########
Pesticide_category=c('Pesticide',
                     unique(combined$pesticide_by_target_organisms))
dat_list1=c(lapply(unique(combined$Taxonomic_group),
                   function(x){subset(combined,Taxonmic_group==x)}),
            lapply(unique(combined$Taxonomic_group_response)[c(4,1:3,5:7,10,8:9)],
                   function(x){subset(combined,Taxonomc_group_response==x)}))
Experiment=unique(combined$Major.climatic.zones)
S8=matrix(NA,length(Experiment)*length(Pesticide_category)*length(dat_list1),8)
for (k in 1:length(Experiment)) {
  for (i in 1:length(Pesticide_category)) {
    for (j in 1:length(dat_list1)){
        mydat=subset(dat_list1[[j]],
                     Major.climatic.zones==Experiment[k]&pesticide_by_target_organisms==Pesticide_category[i]&Experiment.type=='Field')
      if (nrow(mydat)<=1) {
        S8[(k-1)*52+(i-1)*13+j,]=c(nrow(mydat),length(unique(mydat$Code)),
                                   rep(NA,6))
      } else {
        res=rma.mv(yi, vi, 
                   random = ~ 1 | factor(Insecticide.name),
                   data=mydat)
        S8[(k-1)*52+(i-1)*13+j,]=c(nrow(mydat),length(unique(mydat$Code)),
                                   res$b,res$zval,nrow(mydat)-1,res$pval,res$ci.lb,res$ci.ub)
      }
    }
  }
}
colnames(S8)<-c("#Obs","#Studies","Effect size","T-value","Df",
                "P-value","CL.lb","CL.ub")
write.xlsx(S8,"Supplementary Table 42.xls")

########## Supplementary Table 43 ###########
Old_new_pesticide=read.xlsx('Supplementary Data 2-Classification of old and new pesticides-2023-0505.xls',sheetName='classification',header=T)
combined[,'Old-new-pesticide']=Old_new_pesticide[index,'Old.or.new.pesticides']

Pesticide_category=c('Pesticide',
                     unique(combined$pesticide_by_target_organisms))
dat_list1=c(lapply(unique(combined$Taxonomic_group),
                   function(x){subset(combined,Taxonomic_group==x)}),
            lapply(unique(combined$Taxonomic_group_response)[c(4,1:3,5:7,10,8:9)],
                   function(x){subset(combined,Taxonomic_group_response==x)}))

Experiment=unique(combined$`Old-new-pesticide`)
S9=matrix(NA,length(Experiment)*length(Pesticide_category)*length(dat_list1),8)
for (k in 1:length(Experiment)) {
  for (i in 1:length(Pesticide_category)) {
    for (j in 1:length(dat_list1)){
        mydat=subset(dat_list1[[j]],
                     `Old-new-pesticide`==Experiment[k]&pesticide_by_target_organisms==Pesticide_category[i])
      if (nrow(mydat)<=1) {
        S9[(k-1)*52+(i-1)*13+j,]=c(nrow(mydat),length(unique(mydat$Code)),
                                   rep(NA,6))
      } else {
        res=rma.mv(yi, vi, 
                   random = ~ 1 | factor(Insecticide.name),
                   data=mydat)
        S9[(k-1)*52+(i-1)*13+j,]=c(nrow(mydat),length(unique(mydat$Code)),
                                   res$b,res$zval,nrow(mydat)-1,res$pval,res$ci.lb,res$ci.ub)
      }
    }
  }
}
colnames(S9)<-c("#Obs","#Studies","Effect size","T-value","Df",
                "P-value","CL.lb","CL.ub")
write.xlsx(S9,"Supplementary Table 43.xls")


########## Supplementary Table 44 ###########
Model_animal=read.xlsx('Supplementary Data 3-Classification of model and nonmodel species-2023-0505.xls',sheetName='animals',header=T)
Model_plant=read.xlsx('Supplementary Data 3-Classification of model and nonmodel species-2023-0505.xls',sheetName='plants',header=T)
Model_micro=read.xlsx('Supplementary Data 3-Classification of model and nonmodel species-2023-0505.xls',sheetName='microorganisms',header=T)
Model_comb=rbind(Model_animal[,c(1,6)],Model_plant[,c(1,6)],Model_micro[,c(1,6)])

combined[,'Model-non']=sapply(combined$Latin.name,function(x) {
  index=which(x==Model_comb$Species.or.other.groups)[1]
  Model_comb[index,'Model.or.non.model.species']
})

Pesticide_category=c('Pesticide',
                     unique(combined$pesticide_by_target_organisms))
dat_list1=c(lapply(unique(combined$Taxonomic_group),
                   function(x){subset(combined,Taonomic_group==x)}),
            lapply(unique(combined$Taxonomic_group_response)[c(4,1:3,5:7,10,8:9)],
                   function(x){subset(combined,Taxonmic_group_response==x)}))

Experiment=unique(combined$`Model-non`)[c(2,1)]
S9=matrix(NA,length(Experiment)*length(Pesticide_category)*length(dat_list1),8)
for (k in 1:length(Experiment)) {
  for (i in 1:length(Pesticide_category)) {
    for (j in 1:length(dat_list1)){
        mydat=subset(dat_list1[[j]],
                     `Model-non`==Experiment[k]&pesticide_by_target_organisms==Pesticide_category[i])
      if (nrow(mydat)<=1) {
        S9[(k-1)*52+(i-1)*13+j,]=c(nrow(mydat),length(unique(mydat$Code)),
                                   rep(NA,6))
      } else {
        res=rma.mv(yi, vi, 
                   random = ~ 1 | factor(Insecticide.name),
                   data=mydat)
        S9[(k-1)*52+(i-1)*13+j,]=c(nrow(mydat),length(unique(mydat$Code)),
                                   res$b,res$zval,nrow(mydat)-1,res$pval,res$ci.lb,res$ci.ub)
      }
    }
  }
}

colnames(S9)<-c("#Obs","#Studies","Effect size","T-value","Df",
                "P-value","CL.lb","CL.ub")
write.xlsx(S9,"Supplementary Table 44.xls")

