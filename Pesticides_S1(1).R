library(metafor)
library(xlsx)
setwd('...')
combined=read.csv('Meta combined data.csv',header=T)

Old_new_pesticide=read.xlsx('Supplementary Data 2-Classification of old and new pesticides-2023-0505.xls',sheetName='classification',header=T)

combined[,'Old-new-pesticide']= Old_new_pesticide[index,'Old.or.new.pesticides']

combined[which(is.na(combined$`Old-new-pesticide`)),'Old-new-pesticide']="Old pesticide"

Model_animal=read.xlsx('Supplementary Data 3-Classification of model and nonmodel species-2023-0505.xls',sheetName='animals',header=T)
Model_plant=read.xlsx('Supplementary Data 3-Classification of model and nonmodel species-2023-0505.xls',sheetName='plants',header=T)
Model_micro=read.xlsx('Supplementary Data 3-Classification of model and nonmodel species-2023-0505.xls',sheetName='microorganisms',header=T)
Model_comb=rbind(Model_animal[,c(1,6)],Model_plant[,c(1,6)],Model_micro[,c(1,6)])

combined[,'Model-non']=sapply(combined$Latin.name,function(x) {
  index=which(x==Model_comb$Species.or.other.groups)
  Model_comb[index,'Model.or.non.model.species']
})

## supplementary table 1

# supplementary table 1(1)
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
                     "xi","yi","vi")]
system.time(model1<-rma.mv(yi,vi,mods=~factor(Taxonomic_group)-1,
                           random = list(~ 1 | factor(Insecticide.name),~1|Code),
                           data=combined)

S2_model1=matrix(NA,17,5)
for (i in 1:length(myindex)) {
  if(i<length(myindex)) {
        modelB<-rma.mv(yi,vi,mods=~factor(Taxonomic_group)+factor(combined[,myindex[i]])+factor(Taxonomic_group):factor(combined[,myindex[i]])-1,
                       random = list(~ 1 | factor(Insecticide.name),~1|Code),
                       data=combined)
        modelA<-rma.mv(yi,vi,mods=~factor(Taxonomic_group)+factor(combined[,myindex[i]])-1,
                       random = list(~ 1 | factor(Insecticide.name),~1|Code),
                       data=combined)
        A1=anova(modelA,model1)
        BA=anova(modelB,modelA)
    }
  }
  else{
    modelB<-rma.mv(yi,vi,mods=~factor(Taxonomic_group)+xi+factor(Taxonomic_group):xi-1,
                   random = list(~ 1 | factor(Insecticide.name),~1|Code),
                   data=combined)
    modelA<-rma.mv(yi,vi,mods=~factor(Taxonomic_group)+xi-1,
                   random = list(~ 1 | factor(Insecticide.name),~1|Code),
                   data=combined)
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
               random = list(~ 1 | factor(Insecticide.name),~1|Code),
               data=combined)
S2_model2=matrix(NA,17,5)
aa=anova(model2,model1)
S2_model2[1,]=c(aa[["fit.stats.f"]][["AIC"]],aa[["fit.stats.f"]][["ll"]],
                aa[["LRT"]],aa[["parms.f"]],aa[["pval"]])
for (i in 1:length(myindex)) {
  if(i<length(myindex)) {    
        modelB<-rma.mv(yi,vi,mods=~factor(Taxonomic_group_response)+factor(combined[,myindex[i]])+factor(Taxonomic_group_response):factor(combined[,myindex[i]])-1,
                       random = list(~ 1 | factor(Insecticide.name),~1|Code),
                       data=combined)
        modelA<-rma.mv(yi,vi,mods=~factor(Taxonomic_group_response)+factor(combined[,myindex[i]])-1,
                       random = list(~ 1 | factor(Insecticide.name),~1|Code),
                       data=combined)
        A1=anova(modelA,model1)
        BA=anova(modelB,modelA)
    }
  }
  else{
    modelB<-rma.mv(yi,vi,mods=~factor(Taxonomic_group_response)+xi+factor(Taxonomic_group_response):xi-1,
                   random = list(~ 1 | factor(Insecticide.name),~1|Code),
                   data=combined)
    modelA<-rma.mv(yi,vi,mods=~factor(Taxonomic_group_response)+xi-1,
                   random = list(~ 1 | factor(Insecticide.name),~1|Code),
                   data=combined)
    A1=anova(modelA,model1)
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
write.xlsx(S2,'Supplementary table 1(1).xls',row.names = F)
save.image(file = 'Supplementary Table 1(1).Rdata')












