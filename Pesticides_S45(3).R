library(metafor)
library(xlsx)
library(optimParallel)
setwd('...')
combined=read.csv('Meta combined data.csv',header=T)

Old_new_pesticide=read.xlsx('Supplementary Data 2-Classification of old and new pesticides-2023-0505.xls',sheetName='classification',header=T)

combined[,'Old-new-pesticide']=Old_new_pesticide[index,'Old.or.new.pesticides']


Model_animal=read.xlsx('Supplementary Data 3-Classification of model and nonmodel species-2023-0505.xls',sheetName='animals',header=T)
Model_plant=read.xlsx('Supplementary Data 3-Classification of model and nonmodel species-2023-0505.xls',sheetName='plants',header=T)
Model_micro=read.xlsx('Supplementary Data 3-Classification of model and nonmodel species-2023-0505.xls',sheetName='microorganisms',header=T)
Model_comb=rbind(Model_animal[,c(1,6)],Model_plant[,c(1,6)],Model_micro[,c(1,6)])

combined[,'Model-non']=sapply(combined$Latin.name,function(x) {
  index=which(x==Model_comb$Species.or.other.groups)[1]
  Model_comb[index,'Model.or.non.model.species']
})

combined=combined[combined$yi>=-2&combined$yi<=2,]

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

## supplementary table 45

# supplementary table 45(3)
modelis=list()
myindex=c('Experiment.type','Major.climatic.zones',
          'pesticide_by_target_organisms','pesticide_by_source',
          'Types.of.organism.exposure.to.pesticides',
          'Old-new-pesticide',
          'Model-non','xi')
# combined=combined[sample(nrow(combined))[1:1000],]
combined=combined[,c("Code","Major.climatic.zones","Experiment.type",
                     "Insecticide.name",
                     "Types.of.organism.exposure.to.pesticides",
                     "No.insecticide..control",
                     "Insecticide.dosage.treat","category",
                     "Taxonomic_group_response","Taxonomic_group",
                     "pesticide_by_target_organisms","pesticide_by_source",
                     'Old-new-pesticide','Model-non',
                     "xi","yi","vi","Animal.Latin.name")]

mydat=myphylo(combined)
system.time(model1<-rma.mv(yi,vi,mods=~factor(Taxonomic_group)-1,
                           random =list( ~ 1 | factor(Insecticide.name),
                                         ~ 1 | Code,
                                         ~ 1 | Plant.species.new,
                                         ~ 1|Plant.species.new.p),
                           R=list(Plant.species.new.p=mydat[[2]]),
                           data=mydat[[1]])

S2_model1=matrix(NA,17,5)
for (i in 1:length(myindex)) {
  if(i<length(myindex)) {
        modelB<-rma.mv(yi,vi,mods=~factor(Taxonomic_group)+factor(mydat[[1]][,myindex[i]])+factor(Taxonomic_group):factor(mydat[[1]][,myindex[i]])-1,
                       random =list( ~ 1 | factor(Insecticide.name),
                                     ~ 1 | Code,
                                     ~ 1 | Plant.species.new,
                                     ~ 1|Plant.species.new.p),
                       R=list(Plant.species.new.p=mydat[[2]]),
                       data=mydat[[1]])
        modelA<-rma.mv(yi,vi,mods=~factor(Taxonomic_group)+factor(mydat[[1]][,myindex[i]])-1,
                       random =list( ~ 1 | factor(Insecticide.name),
                                     ~ 1 | Code,
                                     ~ 1 | Plant.species.new,
                                     ~ 1|Plant.species.new.p),
                       R=list(Plant.species.new.p=mydat[[2]]),
                       data=mydat[[1]])
        A1=anova(modelA,model1)
        BA=anova(modelB,modelA)
    }
  }
  else{
    modelB<-rma.mv(yi,vi,mods=~factor(Taxonomic_group)+xi+factor(Taxonomic_group):xi-1,
                   random =list( ~ 1 | factor(Insecticide.name),
                                 ~ 1 | Code,
                                 ~ 1 | Plant.species.new,
                                 ~ 1|Plant.species.new.p),
                   R=list(Plant.species.new.p=mydat[[2]]),
                   data=mydat[[1]])
    modelA<-rma.mv(yi,vi,mods=~factor(Taxonomic_group)+xi-1,
                   random =list( ~ 1 | factor(Insecticide.name),
                                 ~ 1 | Code,
                                 ~ 1 | Plant.species.new,
                                 ~ 1|Plant.species.new.p),
                   R=list(Plant.species.new.p=mydat[[2]]),
                   data=mydat[[1]])
    A1=anova(modelA,model1)
    BA=anova(modelB,modelA)
  }
  S2_model1[1,]=c(A1[["fit.stats.r"]][["AIC"]],A1[["fit.stats.r"]][["ll"]],'-',A1[["parms.r"]],'-')
  S2_model1[2*i,]=c(A1[["fit.stats.f"]][["AIC"]],A1[["fit.stats.f"]][["ll"]],
                    A1[["LRT"]],A1[["parms.f"]],A1[["pval"]])
  S2_model1[2*i+1,]=c(BA[["fit.stats.f"]][["AIC"]],BA[["fit.stats.f"]][["ll"]],
                      BA[["LRT"]],BA[["parms.f"]],BA[["pval"]])
}


model2<-rma.mv(yi,vi,mods=~factor(Taxonomic_group_response)-1,
               random =list( ~ 1 | factor(Insecticide.name),
                             ~ 1 | Code,
                             ~ 1 | Plant.species.new,
                             ~ 1|Plant.species.new.p),
               R=list(Plant.species.new.p=mydat[[2]]),
               data=mydat[[1]])
S2_model2=matrix(NA,17,5)
aa=anova(model2,model1)
S2_model2[1,]=c(aa[["fit.stats.f"]][["AIC"]],aa[["fit.stats.f"]][["ll"]],
                aa[["LRT"]],aa[["parms.f"]],aa[["pval"]])
for (i in 1:length(myindex)) {
  if(i<length(myindex)) {
        modelB<-rma.mv(yi,vi,mods=~factor(Taxonomic_group_response)+factor(mydat[[1]][,myindex[i]])+factor(Taxonomic_group_response):factor(mydat[[1]][,myindex[i]])-1,
                       random =list( ~ 1 | factor(Insecticide.name),
                                     ~ 1 | Code,
                                     ~ 1 | Plant.species.new,
                                     ~ 1|Plant.species.new.p),
                       R=list(Plant.species.new.p=mydat[[2]]),
                       data=mydat[[1]])
        modelA<-rma.mv(yi,vi,mods=~factor(Taxonomic_group_response)+factor(mydat[[1]][,myindex[i]])-1,
                       random =list( ~ 1 | factor(Insecticide.name),
                                     ~ 1 | Code,
                                     ~ 1 | Plant.species.new,
                                     ~ 1|Plant.species.new.p),
                       R=list(Plant.species.new.p=mydat[[2]]),
                       data=mydat[[1]])
        A1=anova(modelA,model2)
        BA=anova(modelB,modelA)
    }
  }
  else{
    modelB<-rma.mv(yi,vi,mods=~factor(Taxonomic_group_response)+xi+factor(Taxonomic_group_response):xi-1,
                   random =list( ~ 1 | factor(Insecticide.name),
                                 ~ 1 | Code,
                                 ~ 1 | Plant.species.new,
                                 ~ 1|Plant.species.new.p),
                   R=list(Plant.species.new.p=mydat[[2]]),
                   data=mydat[[1]])
    modelA<-rma.mv(yi,vi,mods=~factor(Taxonomic_group_response)+xi-1,
                   random =list( ~ 1 | factor(Insecticide.name),
                                 ~ 1 | Code,
                                 ~ 1 | Plant.species.new,
                                 ~ 1|Plant.species.new.p),
                   R=list(Plant.species.new.p=mydat[[2]]),
                   data=mydat[[1]])
    A1=anova(modelA,model2)
    BA=anova(modelB,modelA)
  }
  S2_model2[2*i,]=c(A1[["fit.stats.f"]][["AIC"]],A1[["fit.stats.f"]][["ll"]],
                    A1[["LRT"]],A1[["parms.f"]],A1[["pval"]])
  S2_model2[2*i+1,]=c(BA[["fit.stats.f"]][["AIC"]],BA[["fit.stats.f"]][["ll"]],
                      BA[["LRT"]],BA[["parms.f"]],BA[["pval"]])
}
Ref=c('-',rep(c(1,"A"),8),1,rep(c(2,"A"),8))
Predictor_category1=c('Taxonomic group + Type of experimental study',
                      'Taxonomic group + Type of experimental study + Taxonomic group × Type of experimental study',
                      'Taxonomic group + Climatic zone type',
                      'Taxonomic group + Climatic zone type + Taxonomic group × Climatic zone type',
                      'Taxonomic group + Types of pesticide classification by controlling target organisms',
                      'Taxonomic group + Types of pesticide classification by controlling target organisms + Taxonomic group × Types of pesticide classification by controlling target organisms',
                      'Taxonomic group + Types of pesticide classification by source',
                      'Taxonomic group + Types of pesticide classification by source + Taxonomic group × Types of pesticide classification by source',
                      'Taxonomic group + Types of organism exposure to pesticides',
                      'Taxonomic group + Types of organism exposure to pesticides + Taxonomic group × Types of organism exposure to pesticides',
                      'Taxonomic group + Types of pesticide use history',
                      'Taxonomic group + Types of pesticide use history + Taxonomic group × Types of pesticide use history',
                      'Taxonomic group + Types of experimental organisms',
                      'Taxonomic group + Types of experimental organisms + Taxonomic group × Types of experimental organisms',
                      'Taxonomic group + Log2 (added pesticide dosage over control)',
                      'Taxonomic group + Log2 (added pesticide dosage over control) + Taxonomic group × Log2 (added pesticide dosage over control)')
Predictor_category2=c('Taxonomic group response category + Type of experimental study',
                      'Taxonomic group response category + Type of experimental study + Taxonomic group response category × Type of experimental study',
                      'Taxonomic group response category + Climatic zone type',
                      'Taxonomic group response category + Climatic zone type + Taxonomic group response category × Climatic zone type',
                      'Taxonomic group response category + Types of pesticide classification by controlling target organisms',
                      'Taxonomic group response category + Types of pesticide classification by controlling target organisms + Taxonomic group response category × Types of pesticide classification by controlling target organisms',
                      'Taxonomic group response category + Types of pesticide classification by source',
                      'Taxonomic group response category + Types of pesticide classification by source + Taxonomic group response category × Types of pesticide classification by source',
                      'Taxonomic group response category + Types of organism exposure to pesticides',
                      'Taxonomic group response category + Types of organism exposure to pesticides + Taxonomic group response category × Types of organism exposure to pesticides',
                      'Taxonomic group response category + Types of pesticide use history',
                      'Taxonomic group response category + Types of pesticide use history + Taxonomic group response category × Types of pesticide use history',
                      'Taxonomic group response category + Types of experimental organisms',
                      'Taxonomic group response category + Types of experimental organisms + Taxonomic group response category × Types of experimental organisms',
                      'Taxonomic group response category + Log2 (added pesticide dosage over control)',
                      'Taxonomic group response category + Log2 (added pesticide dosage over control) + Taxonomic group response category × Log2 (added pesticide dosage over control)')
Predictor=c('Taxonomic group',Predictor_category1,
            'Taxonomic group response category',Predictor_category2)

n=c(model1[["k"]],unlist(lapply(modelis,function(x)x[["k"]]))[1:16],
    model2[["k"]],unlist(lapply(modelis,function(x)x[["k"]]))[17:32])
S2=cbind(Predictor,Ref,rbind(S2_model1,S2_model2),n)
colnames(S2)=c('Predictor','Ref','AIC','L-L','χ2','d.f.','P','n')
write.xlsx(S2,'Supplementary table 45(3).xls',row.names = F)
save.image(file = 'Supplementary Table 45(3).Rdata')





