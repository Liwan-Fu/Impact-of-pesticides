library(metafor)
library(xlsx)
library(optimParallel)
setwd('G:/博士1001/wannianfeng/3 农药多样性/Analysis files')
combined=read.csv('Meta combined data.csv',header=T)

pub_conf=read.xlsx('Publication year and Conflict of interest-2023-1009-version 2.xlsx',sheetIndex = 1,header = T)
combined$Publication_year=as.numeric(sapply(combined$Code,function(x){
  pub_conf[which(pub_conf$Item==x),'Publication.year']
}))
combined$Conflict_of_interst=as.numeric(sapply(combined$Code,function(x){
  pub_conf[which(pub_conf$Item==x),'Conflict.of.interst']
}))

### add variables for subsequent analyses
combined$Taxonomic_group_response=sapply(combined$category, function(x){
  if (x %in% grep('animal growth',unique(combined$category),value = T)){
    'animal growth'
  } else if (x %in% grep('animal reproduction',unique(combined$category),value = T)) {
    'animal reproduction'
  } else if (x %in% grep('animal behavior',unique(combined$category),value = T)) {
    'animal behavior'
  } else if (x %in% grep('animal biomarker',unique(combined$category),value = T)) {
    'animal biomarker'
  } else if (x %in% grep('plant growth',unique(combined$category),value = T)) {
    'plant growth'
  } else if (x %in% grep('plant reproduction',unique(combined$category),value = T)) {
    'plant reproduction'
  } else if (x %in% grep('plant biomarker',unique(combined$category),value = T)) {
    'plant biomarker'
  } else if (x %in% grep('microorgan growth',unique(combined$category),value = T)) {
    'microorganism growth'
  } else if (x %in% c('insecticide-microo reproduction',
                      'fungicide-microorg reproduction',
                      'herbicide-microorg reproduction')) {
    'microorganism reproduction'
  } else {'microorganism biomarker'}
})

combined$Taxonomic_group=sapply(combined$category, function(x){
  if (x %in% grep('animal',unique(combined$category),value = T)){
    'animals'
  } else if (x %in% grep('plant',unique(combined$category),value = T)) {
    'plants'
  } else {'microorganisms'}
})

combined$pesticide_by_target_organisms=sapply(combined$category, function(x){
  if (x %in% grep('insecticide',unique(combined$category),value = T)){
    'insecticides'
  } else if (x %in% grep('herbicide',unique(combined$category),value = T)) {
    'herbicides'
  } else {'fungicides'}
})

combined$pesticide_by_source=sapply(combined$Primary.classification.of.insecticide........一级分类,
                                    function(x) {
                                      if (x %in% grep('Mineral',unique(combined$Primary.classification.of.insecticide........一级分类),value = T)) {'mineral-based'} else if (x %in% grep('Biogenic',unique(combined$Primary.classification.of.insecticide........一级分类),value = T)) {'biogenic'} else {'chemical'}
                                    })

combined$xi=apply(combined,1,function(x) {
  dat=subset(combined,
             Code==x['Code'] & Insecticide.use.unit==x['Insecticide.use.unit'],
             select = Insecticide.dosage.treat)
  if (length(unique(dat$Insecticide.dosage.treat))==1) {
    NA
  } else {
    as.numeric(x['Insecticide.dosage.treat'])/min(dat)
  }
})
combined$xi=log2(round(combined$xi,3))

combined=escalc(measure="SMD",n1i=Treatment.n, n2i=Control.n, 
                m1i=Treatment..value., m2i=Control..value.,
                sd1i=Treatment.sd, sd2i=Control.sd,
                data=combined,append = TRUE)
for(i in 1:nrow(combined)) {
  if (grepl('biomarker',combined[i,'Taxonomic_group_response'])) {
    combined[i,'yi']=abs(combined[i,'yi'])
  }
}

combined0=subset(combined,Conflict_of_interst==0)
combined1=subset(combined,Conflict_of_interst==1)

########## Supplementary Table 58 ###########
Pesticide_category=c('Pesticide',
                     unique(combined$pesticide_by_target_organisms))
Nontarget_category=unique(combined$Taxonomic_group_response)
dat_list0=lapply(unique(combined$Taxonomic_group_response),
                 function(x){subset(combined0,Taxonomic_group_response==x)})
S0=matrix(NA,length(Pesticide_category)*length(dat_list0),8)
for (i in 1:length(Pesticide_category)) {
  for (j in 1:length(dat_list0)){
    if (i==1) {
      mydat=dat_list0[[j]]
    } else {
      mydat=subset(dat_list0[[j]],
                   pesticide_by_target_organisms==Pesticide_category[i])
    }
    if (nrow(mydat)<=1) {
      S0[(i-1)*length(dat_list0)+j,]=c(nrow(mydat),length(unique(mydat$Code)),
                                       rep(NA,6))
    } else {
      res=rma.mv(yi, vi, 
                 random = list(~ 1 | factor(Insecticide.name),~1|Code,
                               ~1|Publication_year),
                 data=mydat,
                 method = 'ML',
      S0[(i-1)*length(dat_list0)+j,]=c(nrow(mydat),length(unique(mydat$Code)),
                                       res$b,res$zval,nrow(mydat)-1,
                                       res$pval,res$ci.lb,res$ci.ub)
    }
  }
}

dat_list1=lapply(unique(combined$Taxonomic_group_response),
                 function(x){subset(combined1,Taxonomic_group_response==x)})
S1=matrix(NA,length(Pesticide_category)*length(dat_list1),8)
for (i in 1:length(Pesticide_category)) {
  for (j in 1:length(dat_list1)){
    if (i==1) {
      mydat=dat_list1[[j]]
    } else {
      mydat=subset(dat_list1[[j]],
                   pesticide_by_target_organisms==Pesticide_category[i])
    }
    if (nrow(mydat)<=1) {
      S1[(i-1)*length(dat_list1)+j,]=c(nrow(mydat),length(unique(mydat$Code)),
                                       rep(NA,6))
    } else {
      res=rma.mv(yi, vi, 
                 random = list(~ 1 | factor(Insecticide.name),~1|Code,
                               ~1|Publication_year),
                 data=mydat,
                 method = 'ML',
                 sparse=TRUE,verbose=TRUE)
      S1[(i-1)*length(dat_list1)+j,]=c(nrow(mydat),length(unique(mydat$Code)),
                                       res$b,res$zval,nrow(mydat)-1,
                                       res$pval,res$ci.lb,res$ci.ub)
    }
  }
}

S4=rbind(S0,S1)

Rnames=lapply(Pesticide_category,function(x){
  paste(x,Nontarget_category,sep = '_')
})
Rnames=do.call('c',Rnames)
Rnames=c(paste0('Non_',Rnames),paste0('Conf_',Rnames))
S4<-as.data.frame(S4,row.names<-Rnames)
colnames(S4)<-c("#Obs","#Studies","Effect size","T-value","Df",
                "P-value","CL.lb","CL.ub")
write.xlsx(S4,"Supplementary Table 58.xls")
