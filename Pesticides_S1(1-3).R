rm(list = ls())
setwd("....")
DataPre_path <- "...."
library(readxl)
library(metafor)
library(ape)
library(rotl)

Animal_data <- read.table(paste(DataPre_path, "Animal_Combined_data.txt", sep = ""),
                          header = T)
Old_New_pesticide <- read_excel(paste("/lustre/huyueqing/ssy/Wan2023/DataPre/",
                                      "Supplementary Data 2-Classification of old and new pesticides-2023-0825.xls",
                                      sep = ""),sheet = 1)
Match_Odnp <- match(Animal_data$Hrbn, Old_New_pesticide$a.i.)
Animal_data[,"Old_New_pesticide"] <- Old_New_pesticide[Match_Odnp,"Old or new pesticides"]

Model_non <- read_excel(paste("/lustre/huyueqing/ssy/Wan2023/DataPre/",
                              "Supplementary Data 3-Classification of model and nonmodel species-2023-0825.xls",
                              sep = ""),sheet = 1)
Match_Mdnn <- match(Animal_data$AnLn, Model_non$`Species or other groups`)
Animal_data[,"Model_non"] <- Model_non[Match_Mdnn,
                                       which(colnames(Model_non) == "Model or non-model species")]


Micro_data <- read.table(paste(DataPre_path,"Microorgan_Combined_data.txt",sep = ""),
                         header = T)
Old_New_pesticide <- read_excel(paste("/lustre/huyueqing/ssy/Wan2023/DataPre/",
                                      "Supplementary Data 2-Classification of old and new pesticides-2023-0825.xls",
                                      sep = ""),sheet = 1)
Match_Odnp <- match(Micro_data$Hrbn, Old_New_pesticide$a.i.)
Micro_data[,"Old_New_pesticide"] <- Old_New_pesticide[Match_Odnp, "Old or new pesticides"]
Model_non <- read_excel(paste("/lustre/huyueqing/ssy/Wan2023/DataPre/",
                              "Supplementary Data 3-Classification of model and nonmodel species-2023-0825.xls",
                              sep = ""),sheet = 3)
Match_Mdnn <- match(Micro_data$Mcln, Model_non$`Species or other groups`)
Micro_data[,"Model_non"] <- Model_non[Match_Mdnn,"Model or non-model species"]

Plant_data <- read.table(paste(DataPre_path,"Plant_Combined_data.txt",sep = ""),
                         header = T)
Old_New_pesticide <- read_excel(paste("/lustre/huyueqing/ssy/Wan2023/DataPre/",
                                      "Supplementary Data 2-Classification of old and new pesticides-2023-0825.xls",
                                      sep = ""),sheet = 1)
Match_Odnp <- match(Plant_data$Hrbn, Old_New_pesticide$a.i.)
Plant_data[,"Old_New_pesticide"] <- Old_New_pesticide[Match_Odnp,"Old or new pesticides"]
Model_non <- read_excel(paste("/lustre/huyueqing/ssy/Wan2023/DataPre/",
                              "Supplementary Data 3-Classification of model and nonmodel species-2023-0825.xls",
                              sep = ""),sheet = 2)
Match_Mdnn <- match(Plant_data$PlLn,Model_non$`Species or other groups`)
Plant_data[,"Model_non"] <- Model_non[Match_Mdnn,"Model or non-model species"]


# "Spcn"---"动植物的名字,"SpLn"拉丁名
Colnam <- c("Code","Lctn","Lttd","Lngt","Mjcz","Expt","Expy","Hrbn",
            "Scoh.二","Pcoh.一","Spcn","SpLn","Tooetp","Scos.二","Pcos.一",
            "Scbi","Scbu","Hruu","Nhr.","Hds.","C.v.","Cntrl.s","Cntrl.n",
            "T.v.","Trtmnt.s","Trtmnt.n","LnR","yi","vi","pesticide_group",
            "response_category","Pbly","Cnoi","Old_New_pesticide","Model_non")
colnames(Animal_data) <- Colnam
colnames(Plant_data) <- Colnam
colnames(Micro_data) <- Colnam
Whole_data <- rbind(Animal_data, Plant_data, Micro_data)

Whole_data$Taxonomic_group <- c(rep("Animal", times = nrow(Animal_data)),
                                rep("Plant", times = nrow(Plant_data)),
                                rep("Microorganim", times = nrow(Micro_data)))
Whole_data$Taxonomic_group_response_category <- rep(1, times = nrow(Whole_data))
for (k in 1:nrow(Whole_data)){
  i <- Whole_data[k, "Taxonomic_group"]
  j <- Whole_data[k, "response_category"]
  response_category <- regmatches(j, regexec("growth|reproduction|biomarker|behavior",j))[[1]][1]
  flag <- paste(i, response_category, sep = "_")
  Whole_data[k, "Taxonomic_group_response_category"] <- flag
}

Whole_data$Tpcs<-rep(1,times=nrow(Whole_data))
for (k in 1:nrow(Whole_data)){
  i <- Whole_data[k,"Pcoh.一"]
  lable <- regmatches(i, regexec("Chemical synthetic|Biogenic|Mineral-based",i))[[1]][1]
  Whole_data[k,"Tpcs"] <- lable
}

Whole_data$Code <-factor(Whole_data$Code)
Whole_data$SpLn <- factor(Whole_data$SpLn)
Whole_data$Cnoi <- factor(Whole_data$Cnoi)
Whole_data$Hrbn <- factor(Whole_data$Hrbn)

Pre_vars <- c("Pbly","Cnoi")

model.null <- rma.mv(yi, vi, mods = ~Taxonomic_group-1,
                     random = list(~1|Code, ~1|Hrbn),
                     data = Whole_data,
                     method = "ML")
Results <- data.frame("Predictors variables",
                      NA,
                      AIC(model.null),
                      logLik(model.null),
                      NA,
                      summary(model.null)$p,
                      NA,
                      summary(model.null)$k,
                      fix.empty.names = F)

for (i in length(Pre_vars)){
  formula.A <- as.formula(paste("~",paste("Taxonomic_group",Pre_vars[i],sep="+"), "-1", sep = ""))
  formula.B <- as.formula(paste("~", paste("Taxonomic_group",Pre_vars[i], sep = "+"),
                                paste("+Taxonomic_group",Pre_vars[i],sep = ":"), "-1", sep = ""))
  Tem_data <- Whole_data
  if (length(table(is.na(Whole_data[,Pre_vars[i]])))==2){
    Tem_data <- Whole_data[!is.na(Whole_data[,Pre_vars[i]]),]
    model.null <- rma.mv(yi, vi, mods = ~Taxonomic_group-1,
                         random = list( ~1|Code, ~1|Hrbn),
                         data = Tem_data,
                         method = "ML")
  }
  model.A <- rma.mv(yi, vi, mods = formula.A,
                    random = list(~1|Code, ~1|Hrbn),
                    data = Tem_data,
                    method = "ML")
  model.B <- rma.mv(yi, vi, mods = formula.B,
                    random = list(~1|Code, ~1|Hrbn),
                    data = Tem_data,
                    method = "ML")
  A1 <- anova(model.A, model.null)
  Results <- rbind(Results, data.frame(as.character(formula.A)[-1],
                                       "1",
                                       A1[["fit.stats.f"]][["AIC"]],
                                       A1[["fit.stats.f"]][["ll"]],
                                       A1[["LRT"]],
                                       A1[["parms.f"]],
                                       A1[["pval"]],
                                       summary(model.A)$k,
                                       fix.empty.names = F))
  if (length(model.A$beta) != length(model.B$beta)){
    BA <- anova(model.B, model.A)
    Results <- rbind(Results, data.frame(as.character(formula.B)[-1],
                                         paste("1.", i, "A", sep = ""),
                                         BA[["fit.stats.f"]][["AIC"]],
                                         BA[["fit.stats.f"]][["ll"]],
                                         BA[["LRT"]],
                                         BA[["parms.f"]],
                                         BA[["pval"]],
                                         summary(model.B)$k,
                                         fix.empty.names = F))
  }
}

model.null2 <- rma.mv(yi, vi, mods = ~Taxonomic_group_response_category-1,
                      random = list(~1|Code, ~1|Hrbn),
                      data = Whole_data,
                      method = "ML")
null21 <- anova(model.null2, model.null)
Results <- rbind(Results, data.frame("Predictor variables",
                                     "2",
                                     null21[["fit.stats.f"]][["AIC"]],
                                     null21[["fit.stats.f"]][["ll"]],
                                     null21[["LRT"]],
                                     null21[["parms.f"]],
                                     null21[["pval"]],
                                     summary(model.null2)$k,
                                     fix.empty.names = F))

for (i in 1:length(Pre_vars)){
  formula.A <- as.formula(paste("~",paste("Taxonomic_group_response_category",Pre_vars[i],sep="+"), "-1", sep = ""))
  formula.B <- as.formula(paste("~", paste("Taxonomic_group_response_category",Pre_vars[i], sep = "+"),
                                paste("+Taxonomic_group_response_category",Pre_vars[i],sep = ":"), "-1", sep = ""))
  Tem_data <- Whole_data
  if (length(table(is.na(Whole_data[,Pre_vars[i]])))==2){
    Tem_data <- Whole_data[!is.na(Whole_data[,Pre_vars[i]]),]
    model.null2 <- rma.mv(yi, vi, mods = ~Taxonomic_group_response_category-1,
                         random = list( ~1|Code, ~1|Hrbn),
                         data = Tem_data,
                         method = "ML")
  }
  model.A <- rma.mv(yi, vi, mods = formula.A,
                    random = list(~1|Code, ~1|Hrbn),
                    data = Tem_data,
                    method = "ML")
  model.B <- rma.mv(yi, vi, mods = formula.B,
                    random = list(~1|Code, ~1|Hrbn),
                    data = Tem_data,
                    method = "ML")
  A1 <- anova(model.A, model.null2)

  Results <- rbind(Results, data.frame(as.character(formula.A)[-1],
                                       "2",
                                       A1[["fit.stats.f"]][["AIC"]],
                                       A1[["fit.stats.f"]][["ll"]],
                                       A1[["LRT"]],
                                       A1[["parms.f"]],
                                       A1[["pval"]],
                                       summary(model.A)$k,
                                       fix.empty.names = F))
  if (length(model.A$beta) != length(model.B$beta)){
    BA <- anova(model.B, model.A)
    Results <- rbind(Results, data.frame(as.character(formula.B)[-1],
                                         paste("2.", i, "A", sep = ""),
                                         BA[["fit.stats.f"]][["AIC"]],
                                         BA[["fit.stats.f"]][["ll"]],
                                         BA[["LRT"]],
                                         BA[["parms.f"]],
                                         BA[["pval"]],
                                         summary(model.B)$k,
                                         fix.empty.names = F))
  }

}


colnames(Results) <- c("Predictor variables",
                       "Ref.",
                       "AIC",
                       "L─L",
                       "Chisq",
                       "df",
                       "P",
                       "N")
write.csv(Results, "SuppTab1-1_1104.csv", row.names = F)


### SuppTab1-2
rm(list = ls())
setwd("/lustre/huyueqing/ssy/Wan2023/SuppTab1")
DataPre_path <- "/lustre/huyueqing/ssy/Wan2023/DataPre/"
library(readxl)
library(metafor)
library(ape)
library(rotl)

Animal_data <- read.table(paste(DataPre_path, "Animal_Combined_data.txt", sep = ""),
                          header = T)
Old_New_pesticide <- read_excel(paste(DataPre_path,
                                      "Supplementary Data 2-Classification of old and new pesticides-2023-0825.xls",
                                      sep = ""),sheet = 1)
Match_Odnp <- match(Animal_data$Hrbn, Old_New_pesticide$a.i.)
Animal_data[,"Old_New_pesticide"] <- Old_New_pesticide[Match_Odnp,"Old or new pesticides"]

Model_non <- read_excel(paste(DataPre_path,
                              "Supplementary Data 3-Classification of model and nonmodel species-2023-0825.xls",
                              sep = ""),sheet = 1)
Match_Mdnn <- match(Animal_data$AnLn, Model_non$`Species or other groups`)
Animal_data[,"Model_non"] <- Model_non[Match_Mdnn,
                                       which(colnames(Model_non) == "Model or non-model species")]


Micro_data <- read.table(paste(DataPre_path,"Microorgan_Combined_data.txt",sep = ""),
                         header = T)
Old_New_pesticide <- read_excel(paste(DataPre_path,
                                      "Supplementary Data 2-Classification of old and new pesticides-2023-0825.xls",
                                      sep = ""),sheet = 1)
Match_Odnp <- match(Micro_data$Hrbn, Old_New_pesticide$a.i.)
Micro_data[,"Old_New_pesticide"] <- Old_New_pesticide[Match_Odnp, "Old or new pesticides"]
Model_non <- read_excel(paste(DataPre_path,
                              "Supplementary Data 3-Classification of model and nonmodel species-2023-0825.xls",
                              sep = ""),sheet = 3)
Match_Mdnn <- match(Micro_data$Mcln, Model_non$`Species or other groups`)
Micro_data[,"Model_non"] <- Model_non[Match_Mdnn,"Model or non-model species"]

Plant_data <- read.table(paste(DataPre_path,"Plant_Combined_data.txt",sep = ""),
                         header = T)
Old_New_pesticide <- read_excel(paste(DataPre_path,
                                      "Supplementary Data 2-Classification of old and new pesticides-2023-0825.xls",
                                      sep = ""),sheet = 1)
Match_Odnp <- match(Plant_data$Hrbn, Old_New_pesticide$a.i.)
Plant_data[,"Old_New_pesticide"] <- Old_New_pesticide[Match_Odnp,"Old or new pesticides"]
Model_non <- read_excel(paste(DataPre_path,
                              "Supplementary Data 3-Classification of model and nonmodel species-2023-0825.xls",
                              sep = ""),sheet = 2)
Match_Mdnn <- match(Plant_data$PlLn,Model_non$`Species or other groups`)
Plant_data[,"Model_non"] <- Model_non[Match_Mdnn,"Model or non-model species"]

# "Spcn"---"动植物的名字,"SpLn"拉丁名
Colnam <- c("Code","Lctn","Lttd","Lngt","Mjcz","Expt","Expy","Hrbn",
            "Scoh.二","Pcoh.一","Spcn","SpLn","Tooetp","Scos.二","Pcos.一",
            "Scbi","Scbu","Hruu","Nhr.","Hds.","C.v.","Cntrl.s","Cntrl.n",
            "T.v.","Trtmnt.s","Trtmnt.n","LnR","yi","vi","pesticide_group",
            "response_category","Pbly","Cnoi","Old_New_pesticide","Model_non")
colnames(Animal_data) <- Colnam
colnames(Plant_data) <- Colnam
colnames(Micro_data) <- Colnam
Whole_data <- rbind(Animal_data, Plant_data, Micro_data)

Whole_data$Taxonomic_group <- Whole_data$Pcos.一
Whole_data$Taxonomic_group_response_category <- rep(1, times = nrow(Whole_data))
for (k in 1:nrow(Whole_data)){
  i <- Whole_data[k, "response_category"]
  response_category <- regmatches(i, regexec("growth|reproduction|biomarker|behavior",i))[[1]][1]
  j <- Whole_data[k, "Taxonomic_group"]
  taxonomic_group <- regmatches(j,regexec("Bacteria|Fungus|Invertebrate|Seed plant|Spore-producing plant|Vertebrate",j))[[1]][1]
  flag <- paste(taxonomic_group, response_category, sep = "_")
  Whole_data[k, "Taxonomic_group_response_category"] <- flag
}

Whole_data$Tpcs<-rep(1,times=nrow(Whole_data))
for (k in 1:nrow(Whole_data)){
  i <- Whole_data[k,"Pcoh.一"]
  lable <- regmatches(i, regexec("Chemical synthetic|Biogenic|Mineral-based",i))[[1]][1]
  Whole_data[k,"Tpcs"] <- lable
}


Whole_data$Code <- factor(Whole_data$Code)
Whole_data$SpLn <- factor(Whole_data$SpLn)
Whole_data$Cnoi <- factor(Whole_data$Cnoi)
Whole_data$Pbly <- factor(Whole_data$Pbly)
Whole_data$Hrbn <- factor(Whole_data$Hrbn)

Pre_vars <- c("Pbly","Cnoi")
model.null <- rma.mv(yi, vi, mods = ~Taxonomic_group-1,
                     random = list( ~1|Code, ~1|Hrbn),
                     data = Whole_data,
                     method = "ML")
Results <- data.frame("Predictors variables",
                      NA,
                      AIC(model.null),
                      logLik(model.null),
                      NA,
                      summary(model.null)$p,
                      NA,
                      summary(model.null)$k,
                      fix.empty.names = F) 

for(i in 1:length(Pre_vars)){
  formula.A <- as.formula(paste("~",paste("Taxonomic_group",Pre_vars[i],sep="+"), "-1", sep = ""))
  formula.B <- as.formula(paste("~", paste("Taxonomic_group",Pre_vars[i], sep = "+"), 
                                paste("+Taxonomic_group",Pre_vars[i],sep = ":"), "-1", sep = ""))
  Tem_data <- Whole_data
  if (length(table(is.na(Whole_data[,Pre_vars[i]])))==2){
    Tem_data <- Whole_data[!is.na(Whole_data[,Pre_vars[i]]),]
    model.null <- rma.mv(yi, vi, mods = ~Taxonomic_group-1,
                         random = list( ~1|Code, ~1|Hrbn),
                         data = Tem_data,
                         method = "ML")
  }
  model.A <- rma.mv(yi, vi, mods = formula.A,
                    random = list( ~1|Code, ~1|Hrbn),
                    data = Tem_data,
                    method = "ML")
  model.B <- rma.mv(yi, vi, mods = formula.B,
                    random = list(~1|Code, ~1|Hrbn),
                    data = Tem_data,
                    method = "ML")
  A1 <- anova(model.A, model.null)
  Results <- rbind(Results, data.frame(as.character(formula.A)[-1],
                                       "1",
                                       A1[["fit.stats.f"]][["AIC"]],
                                       A1[["fit.stats.f"]][["ll"]],
                                       A1[["LRT"]],
                                       A1[["parms.f"]],
                                       A1[["pval"]],
                                       summary(model.A)$k,
                                       fix.empty.names = F))
  if (length(model.A$beta) != length(model.B$beta)){
    BA <- anova(model.B, model.A)
    Results <- rbind(Results, data.frame(as.character(formula.B)[-1],
                                         paste("1.", i, "A", sep = ""),
                                         BA[["fit.stats.f"]][["AIC"]],
                                         BA[["fit.stats.f"]][["ll"]],
                                         BA[["LRT"]],
                                         BA[["parms.f"]],
                                         BA[["pval"]],
                                         summary(model.B)$k,
                                         fix.empty.names = F))
  }
  
}



model.null2 <- rma.mv(yi, vi, mods = ~Taxonomic_group_response_category-1,
                      random = list(~1|Code, ~1|Hrbn),
                      data = Whole_data,
                      method = "ML")
# null21 <- anova(model.null, model.null2)
# Results <-  rbind(Results, data.frame("Predictor variables",
#                                       "2",
#                                       null21[["fit.stats.f"]][["AIC"]],
#                                       null21[["fit.stats.f"]][["ll"]],
#                                       null21[["LRT"]],
#                                       null21[["parms.f"]],
#                                       null21[["pval"]],
#                                       summary(model.null2)$k,
#                                       fix.empty.names = F))

for(i in 1:length(Pre_vars)){
  formula.A <- as.formula(paste("~",paste("Taxonomic_group_response_category",Pre_vars[i],sep="+"), "-1", sep = ""))
  formula.B <- as.formula(paste("~", paste("Taxonomic_group_response_category",Pre_vars[i], sep = "+"), 
                                paste("+Taxonomic_group_response_category",Pre_vars[i],sep = ":"), "-1", sep = ""))
  Tem_data <- Whole_data
  if (length(table(is.na(Whole_data[,Pre_vars[i]])))==2){
    Tem_data <- Whole_data[!is.na(Whole_data[,Pre_vars[i]]),]
    model.null2 <- rma.mv(yi, vi, mods = ~Taxonomic_group_response_category-1,
                          random = list(~1|Code, ~1|Hrbn),
                          data = Tem_data,
                          method = "ML")
  }
  model.A <- rma.mv(yi, vi, mods = formula.A,
                    random = list(~1|Code, ~1|Hrbn),
                    data = Tem_data,
                    method = "ML")
  model.B <- rma.mv(yi, vi, mods = formula.B,
                    random = list( ~1|Code, ~1|Hrbn),
                    data = Tem_data,
                    method = "ML")
  A2 <- anova(model.A, model.null2)
  Results <- rbind(Results, data.frame(as.character(formula.A)[-1],
                                       "2",
                                       A2[["fit.stats.f"]][["AIC"]],
                                       A2[["fit.stats.f"]][["ll"]],
                                       A2[["LRT"]],
                                       A2[["parms.f"]],
                                       A2[["pval"]],
                                       summary(model.A)$k,
                                       fix.empty.names = F))
  if (length(model.A$beta) != length(model.B$beta)){
    BA <- anova(model.B, model.A)
    Results <- rbind(Results, data.frame(as.character(formula.B)[-1],
                                         paste("2.", i, "A", sep = ""),
                                         BA[["fit.stats.f"]][["AIC"]],
                                         BA[["fit.stats.f"]][["ll"]],
                                         BA[["LRT"]],
                                         BA[["parms.f"]],
                                         BA[["pval"]],
                                         summary(model.B)$k,
                                         fix.empty.names = F))
  }
  
}





colnames(Results) <- c("Predictor variables",
                       "Ref.",
                       "AIC",
                       "L─L",
                       "Chisq",
                       "df",
                       "P",
                       "N")
write.csv(Results, "SuppTab1-2_1104.csv", row.names = F)

### SuppTab1-3
rm(list = ls())
setwd("/lustre/huyueqing/ssy/Wan2023/SuppTab1")
DataPre_path <- "/lustre/huyueqing/ssy/Wan2023/DataPre/"
library(readxl)
library(metafor)
library(ape)
library(rotl)
library(stringr)

# data pre
Plant_data <- read.table(paste(DataPre_path,"Plant_Combined_data.txt",sep = ""),
                         header = T)
Old_New_pesticide <- read_excel(paste(DataPre_path,
                                      "Supplementary Data 2-Classification of old and new pesticides-2023-0825.xls",
                                      sep = ""),sheet = 1)
Match_Odnp <- match(Plant_data$Hrbn, Old_New_pesticide$a.i.)
Plant_data[,"Old_New_pesticide"] <- Old_New_pesticide[Match_Odnp,"Old or new pesticides"]
Model_non <- read_excel(paste(DataPre_path,
                              "Supplementary Data 3-Classification of model and nonmodel species-2023-0825.xls",
                              sep = ""),sheet = 2)
Match_Mdnn <- match(Plant_data$PlLn,Model_non$`Species or other groups`)
Plant_data[,"Model_non"] <- Model_non[Match_Mdnn,"Model or non-model species"]
Plant_data[,"Log_added_dosage"] <- log2(Plant_data$Hds. - Plant_data$Nhr.)
Plant_data$Taxonomic_group <- Plant_data$Scop.二
Plant_data$Taxonomic_group_response_category <- rep(1, times = nrow(Plant_data))
for (k in 1:nrow(Plant_data)){
  i <- Plant_data[k, "response_category"]
  response_category <- regmatches(i, regexec("growth|reproduction|biomarker",i))[[1]][1]
  if (Plant_data[k, "Taxonomic_group"] == "Dicotyledonous plant"){
    flag <- paste0("Dicotyledonous plant",sep = "_", response_category)
  }
  else if(Plant_data[k, "Taxonomic_group"] == "Algae plant"){
    flag <- paste0("Algae plant",sep = "_", response_category)
  }
  else{
    flag <- paste0("Monocotyledonous plant",sep = "_", response_category)
  }
  Plant_data[k, "Taxonomic_group_response_category"] <- flag
}

Plant_data$Tpcs <- rep(1,times=nrow(Plant_data))
for (k in 1:nrow(Plant_data)){
  i <- Plant_data[k,"Pcoh.一"]
  lable <- regmatches(i, regexec("Chemical synthetic|Biogenic|Mineral-based",i))[[1]][1]
  Plant_data[k,"Tpcs"] <- lable
}


## another method
library(V.PhyloMaker)
library(rgbif)


Plant_Cor <- read.table(paste(DataPre_path,"Plant_Cor.txt",sep = ""),header = T)

colnames(Plant_Cor) <- rownames(Plant_Cor)
Plant_data$PlLn <- str_to_lower(Plant_data$PlLn)
Matched_Plant_InCor <- match(Plant_data$PlLn,colnames(Plant_Cor))
Plant_data <- Plant_data[!is.na(Matched_Plant_InCor),]

Plant_Cor <- Plant_Cor[order(dimnames(Plant_Cor)[[1]]), order(dimnames(Plant_Cor)[[1]])]
Plant_data$Plant_Species <- factor(Plant_data$PlLn)
levels(Plant_data$Plant_Species) <- sort(dimnames(Plant_Cor)[[1]])


Plant_data$Code <- factor(Plant_data$Code)
Plant_data$PlLn <- factor(Plant_data$PlLn)
Plant_data$Expt <- factor(Plant_data$Expt)
Plant_data$Mjcz <- factor(Plant_data$Mjcz)
Plant_data$pesticide_group <- factor(Plant_data$pesticide_group)
Plant_data$Tpcs <- factor(Plant_data$Tpcs)
Plant_data$Tooept <- factor(Plant_data$Tooetp)
Plant_data$Old_New_pesticide <- factor(Plant_data$Old_New_pesticide)
Plant_data$Model_non <- factor(Plant_data$Model_non)
Plant_data$Cnoi <- factor(Plant_data$Cnoi)
Plant_data$Pbly <- factor(Plant_data$Pbly)
Plant_data$Hrbn <- factor(Plant_data$Hrbn)


# compute the model
Pre_vars <- c("Pbly","Cnoi")
t1 <- Sys.time()
model.null <- rma.mv(yi, vi, mods = ~Taxonomic_group-1,
                     random = list(~1|Plant_Species, ~1|Code, ~1|Hrbn),
                     R = list(Plant_Species = Plant_Cor[Plant_data$Plant_Species, Plant_data$Plant_Species]),
                     data = Plant_data,
                     method = "ML")
t2 <- Sys.time()
ti <- t2-t1

Results <- data.frame("Predictors variables",
                      NA,
                      AIC(model.null),
                      logLik(model.null),
                      NA,
                      summary(model.null)$p,
                      NA,
                      summary(model.null)$k,
                      fix.empty.names = F) 

for (i in 1:length(Pre_vars)){
  formula.A <- as.formula(paste("~",paste("Taxonomic_group",Pre_vars[i],sep="+"), "-1", sep = ""))
  formula.B <- as.formula(paste("~", paste("Taxonomic_group",Pre_vars[i], sep = "+"), 
                                paste("+Taxonomic_group",Pre_vars[i],sep = ":"), "-1", sep = ""))
  Tem_data <- Plant_data
  if (length(table(is.na(Plant_data[,Pre_vars[i]])))==2){
    Tem_data <- Plant_data[!is.na(Plant_data[,Pre_vars[i]]),]
    model.null <- rma.mv(yi, vi, mods = ~Taxonomic_group-1,
                         random = list(~1|Plant_Species, ~1|Code, ~1|Hrbn),
                         R = list(Plant_Species = Plant_Cor[Plant_data$Plant_Species, Plant_data$Plant_Species]),
                         data = Tem_data,
                         method = "ML")
  }
  model.A <- rma.mv(yi, vi, mods = formula.A,
                    random = list(~1|Plant_Species, ~1|Code, ~1|Hrbn),
                    R = list(Plant_Species = Plant_Cor[Tem_data$Plant_Species, Tem_data$Plant_Species]),
                    data = Tem_data,
                    method = "ML")
  model.B <- rma.mv(yi, vi, mods = formula.B,
                    random = list(~1|Plant_Species, ~1|Code, ~1|Hrbn),
                    R = list(Plant_Species = Plant_Cor[Tem_data$Plant_Species, Tem_data$Plant_Species]),
                    data = Tem_data,
                    method = "ML")
  A1 <- anova(model.A, model.null)
  Results <- rbind(Results, data.frame(as.character(formula.A)[-1],
                                       "1",
                                       A1[["fit.stats.f"]][["AIC"]],
                                       A1[["fit.stats.f"]][["ll"]],
                                       A1[["LRT"]],
                                       A1[["parms.f"]],
                                       A1[["pval"]],
                                       summary(model.A)$k,
                                       fix.empty.names = F))
  if (length(model.A$beta) != length(model.B$beta)){
    BA <- anova(model.B, model.A)
    Results <- rbind(Results, data.frame(as.character(formula.B)[-1],
                                         paste("1.", i, "A", sep = ""),
                                         BA[["fit.stats.f"]][["AIC"]],
                                         BA[["fit.stats.f"]][["ll"]],
                                         BA[["LRT"]],
                                         BA[["parms.f"]],
                                         BA[["pval"]],
                                         summary(model.B)$k,
                                         fix.empty.names = F))
  }
}


model.null2 <- rma.mv(yi, vi,mods = ~Taxonomic_group_response_category-1,
                      random = list(~1|Plant_Species, ~1|Code, ~1|Hrbn),
                      R = list(Plant_Species = Plant_Cor[Plant_data$Plant_Species, Plant_data$Plant_Species]),
                      data = Plant_data,
                      method = "ML")
# null21 <- anova(model.null2, model.null)
# Results <- rbind(Results, data.frame("Predictor variables",
#                                      "2",
#                                      null21[["fit.stats.f"]][["AIC"]],
#                                      null21[["fit.stats.f"]][["ll"]],
#                                      null21[["LRT"]],
#                                      null21[["parms.f"]],
#                                      null21[["pval"]],
#                                      summary(model.null2)$k,
#                                      fix.empty.names = F))

for (i in 1:length(Pre_vars)){
  formula.A <- as.formula(paste("~",paste("Taxonomic_group_response_category",Pre_vars[i],sep="+"), "-1", sep = ""))
  formula.B <- as.formula(paste("~", paste("Taxonomic_group_response_category",Pre_vars[i], sep = "+"), 
                                paste("+Taxonomic_group_response_category",Pre_vars[i],sep = ":"), "-1", sep = ""))
  
  Tem_data <- Plant_data
  if (length(table(is.na(Plant_data[,Pre_vars[i]])))==2){
    Tem_data <- Plant_data[!is.na(Plant_data[,Pre_vars[i]]),]
    model.null2 <- rma.mv(yi, vi,mods = ~Taxonomic_group_response_category-1,
                          random = list(~1|Plant_Species, ~1|Code, ~1|Hrbn),
                          R = list(Plant_Species = Plant_Cor[Plant_data$Plant_Species, Plant_data$Plant_Species]),
                          data = Tem_data,
                          method = "ML")
  }
  model.A <- rma.mv(yi, vi, mods = formula.A,
                    random = list(~1|Plant_Species, ~1|Code, ~1|Hrbn),
                    R = list(Plant_Species = Plant_Cor[Tem_data$Plant_Species, Tem_data$Plant_Species]),
                    data = Tem_data,
                    method = "ML")
  model.B <- rma.mv(yi, vi, mods = formula.B,
                    random = list(~1|Plant_Species, ~1|Code, ~1|Hrbn),
                    R = list(Plant_Species = Plant_Cor[Tem_data$Plant_Species, Tem_data$Plant_Species]),
                    data = Tem_data,
                    method = "ML")
  A2 <- anova(model.A, model.null2)
  
  Results <- rbind(Results, data.frame(as.character(formula.A)[-1],
                                       "2",
                                       A2[["fit.stats.f"]][["AIC"]],
                                       A2[["fit.stats.f"]][["ll"]],
                                       A2[["LRT"]],
                                       A2[["parms.f"]],
                                       A2[["pval"]],
                                       summary(model.A)$k,
                                       fix.empty.names = F))
  if (length(model.A$beta) != length(model.B$beta)){
    BA <- anova(model.B, model.A)
    Results <- rbind(Results, data.frame(as.character(formula.B)[-1],
                                         paste("2.", i, "A", sep = ""),
                                         BA[["fit.stats.f"]][["AIC"]],
                                         BA[["fit.stats.f"]][["ll"]],
                                         BA[["LRT"]],
                                         BA[["parms.f"]],
                                         BA[["pval"]],
                                         summary(model.B)$k,
                                         fix.empty.names = F))
  }
  
}

colnames(Results) <- c("Predictor variables",
                       "Ref.",
                       "AIC",
                       "L─L",
                       "Chisq",
                       "df",
                       "P",
                       "N")
write.csv(Results, "SuppTab1-3_1104.csv", row.names = F)


