rm(list = ls())
setwd("...")
DataPre_path <- "..."
library(readxl)
library(metafor)
library(ape)
library(stringr)
library(rotl)

# data preprocess
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
Animal_data$Log_added_dosage <- log2(Animal_data$Hds. - Animal_data$Nhr.)
Animal_data$Taxonomic_group <- Animal_data$Pcoa一
Animal_data$Taxonomic_group_response_category <- rep(1,times = nrow(Animal_data))
for (k in 1:nrow(Animal_data)){
  i <- Animal_data[k,"response_category"]
  response_category <- regmatches(i, regexec("growth|reproduction|behavior|biomarker",i))[[1]][1]
  if (Animal_data[k,"Taxonomic_group"] == "Vertebrate"){
    flag <- paste0("Vertebrate",sep = "_",response_category)
  }else{
    flag <- paste0("Invertebrate",sep = "_",response_category)
  }
  Animal_data[k,which(colnames(Animal_data) == "Taxonomic_group_response_category")] <- flag
}
Animal_data$Tpcs<-rep(1,times=nrow(Animal_data))
for (l in 1:nrow(Animal_data)){
  i <- Animal_data[k,"Pcoh.一"]
  lable <- regmatches(i, regexec("Chemical synthetic|Biogenic|Mineral-based",i))[[1]][1]
  Animal_data[k,"Tpcs"] <- lable
}



Animals_Cor <- read.table(paste(DataPre_path,
                                "Animal_Cor.txt",
                                sep = ""),sep = "\t")
colnames(Animals_Cor) <- rownames(Animals_Cor)
Animal_data$AnLn <- str_to_lower(Animal_data$AnLn)
Matched_Animal_InCor <- match(Animal_data$AnLn, rownames(Animals_Cor))
Animal_data <- Animal_data[!is.na(Matched_Animal_InCor),]

Animals_Cor <- Animals_Cor[order(dimnames(Animals_Cor)[[1]]), order(dimnames(Animals_Cor)[[1]])]
Animal_data$Animal_Species <- factor(Animal_data$AnLn)
levels(Animal_data$Animal_Species) <- sort(dimnames(Animals_Cor)[[1]])

# compute indicators
Animal_data$Code <- factor(Animal_data$Code)
Animal_data$AnLn <- factor(Animal_data$AnLn)
Animal_data$Taxonomic_group <- factor(Animal_data$Taxonomic_group)
Animal_data$Mjcz <- factor(Animal_data$Mjcz)
Animal_data$Expt <- factor(Animal_data$Expt)
Animal_data$pesticide_group <- factor(Animal_data$pesticide_group)
Animal_data$Tpcs <- factor(Animal_data$Tpcs)
Animal_data$Tooetp <- factor(Animal_data$Tooetp)
Animal_data$Old_New_pesticide <- factor(Animal_data$Old_New_pesticide)
Animal_data$Model_non <- factor(Animal_data$Model_non)
Animal_data$Cnoi <- factor(Animal_data$Cnoi)
Animal_data$Pbly <- factor(Animal_data$Pbly)

Pre_vars <- c("Expt","Mjcz","pesticide_group",
              "Tpcs","Old_New_pesticide",
              "Model_non","Pbly","Log_added_dosage") 


t1 <- Sys.time()
# Taxonomic group as covariant
model.null <- rma.mv(yi, vi, mods = ~Taxonomic_group-1, 
                     random = list(~1|Animal_Species, ~1|Code, ~1|Hrbn),
                     R = list(Animal_Species = Animals_Cor[Animal_data$Animal_Species, Animal_data$Animal_Species]),
                     data = Animal_data,
                     method = "ML")

t2 <- Sys.time()
t2-t1

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
  Tem_data <- Animal_data
  if(length(table(is.na(Animal_data[,Pre_vars[i]])))==2){
    Tem_data <- Animal_data[!is.na(Animal_data[,Pre_vars[i]]),]
  }
  if (Pre_vars[i] == "Mjcz"){
    Tem_data <- Tem_data[which(Tem_data$Expt == "Field"),]
    model.A <- rma.mv(yi, vi, mods = formula.A,
                      random = list(~1|Animal_Species, ~1|Code,~1|Hrbn),
                      R = list(Animal_Species = Animals_Cor[Tem_data$Animal_Species, Tem_data$Animal_Species]),
                      data = Tem_data,
                      method = "ML")
    model.B <- rma.mv(yi, vi, mods = formula.B,
                      random = list(~1|Animal_Species, ~1|Code,~1|Hrbn),
                      R = list(Animal_Species = Animals_Cor[Tem_data$Animal_Species, Tem_data$Animal_Species]),
                      data = Tem_data,
                      method = "ML")
    model.Mjcz <- rma.mv(yi,vi, mods = ~Taxonomic_group-1,
                         random = list(~1|Animal_Species, ~1|Code,~1|Hrbn),
                         R = list(Animal_Species = Animals_Cor[Tem_data$Animal_Species, Tem_data$Animal_Species]),
                         data = Tem_data,
                         method = "ML")
    A1 <- anova(model.A, model.Mjcz)
  }else{
    model.A <- rma.mv(yi, vi, mods = formula.A,
                      random = list(~1|Animal_Species, ~1|Code,~1|Hrbn),
                      R = list(Animal_Species = Animals_Cor[Tem_data$Animal_Species, Tem_data$Animal_Species]),
                      data = Tem_data,
                      method = "ML")
    model.B <- rma.mv(yi, vi, mods = formula.B,
                      random = list(~1|Animal_Species, ~1|Code,~1|Hrbn),
                      R = list(Animal_Species = Animals_Cor[Tem_data$Animal_Species, Tem_data$Animal_Species]),
                      data = Tem_data,
                      method = "ML")
    A1 <- anova(model.A, model.null)
    
  }
  Results <- rbind(Results, data.frame(as.character(formula.A)[-1],
                                       "1",
                                       A1[["fit.stats.f"]][["AIC"]],
                                       A1[["fit.stats.f"]][["ll"]],
                                       A1[["LRT"]],
                                       A1[["parms.f"]],
                                       A1[["pval"]],
                                       summary(model.A)$k,
                                       fix.empty.names = F))
  if(length(model.B$beta) != length(model.A$beta)){
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

# Response category as covariant
model.null2 <- rma.mv(yi, vi, mods = ~Taxonomic_group_response_category-1,
                      random = list(~1|Animal_Species, ~1|Code, ~1|Hrbn),
                      R = list(Animal_Species = Animals_Cor[Animal_data$Animal_Species, Animal_data$Animal_Species]),
                      data = Animal_data,
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
  formula.A <- as.formula(paste("~", paste("Taxonomic_group_response_category",Pre_vars[i],sep = "+"),"-1",sep = ""))
  formula.B <- as.formula(paste("~",paste("Taxonomic_group_response_category",Pre_vars[i],sep = "+"),
                                paste("+Taxonomic_group_response_category",Pre_vars[i], sep = ":"),"-1",sep = ""))
  Tem_data <- Animal_data
  if(length(table(is.na(Animal_data[,Pre_vars[i]])))==2){
    Tem_data <- Animal_data[!is.na(Animal_data[,Pre_vars[i]]),]
  }
  if (Pre_vars[i] == "Mjcz"){
    Tem_data <- Tem_data[which(Tem_data$Expt == "Field"),]
    model.A <- rma.mv(yi, vi, mods = formula.A,
                      random = list(~1|Animal_Species, ~1|Code,~1|Hrbn),
                      R = list(Animal_Species = Animals_Cor[Tem_data$Animal_Species, Tem_data$Animal_Species]),
                      data = Tem_data,
                      method = "ML")
    model.B <- rma.mv(yi, vi, mods = formula.B,
                      random = list(~1|Animal_Species, ~1|Code,~1|Hrbn),
                      R = list(Animal_Species = Animals_Cor[Tem_data$Animal_Species, Tem_data$Animal_Species]),
                      data = Tem_data,
                      method = "ML")
    model.Mjcz <- rma.mv(yi,vi, mods = ~Taxonomic_group-1,
                         random = list(~1|Animal_Species, ~1|Code,~1|Hrbn),
                         R = list(Animal_Species = Animals_Cor[Tem_data$Animal_Species, Tem_data$Animal_Species]),
                         data = Tem_data,
                         method = "ML")
    A1 <- anova(model.A, model.Mjcz)
  }else{
    model.A <- rma.mv(yi, vi, mods = formula.A,
                      random = list(~1|Animal_Species, ~1|Code,~1|Hrbn),
                      R = list(Animal_Species = Animals_Cor[Tem_data$Animal_Species, Tem_data$Animal_Species]),
                      data = Tem_data,
                      method = "ML")
    model.B <- rma.mv(yi, vi, mods = formula.B,
                      random = list(~1|Animal_Species, ~1|Code,~1|Hrbn),
                      R = list(Animal_Species = Animals_Cor[Tem_data$Animal_Species, Tem_data$Animal_Species]),
                      data = Tem_data,
                      method = "ML")
    A1 <- anova(model.A, model.null)
    
  }
  Results <- rbind(Results, data.frame(as.character(formula.A)[-1],
                                       "2",
                                       A1[["fit.stats.f"]][["AIC"]],
                                       A1[["fit.stats.f"]][["ll"]],
                                       A1[["LRT"]],
                                       A1[["parms.f"]],
                                       A1[["pval"]],
                                       summary(model.A)$k,
                                       fix.empty.names = F))
  if(length(model.B$beta) != length(model.A$beta)){
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

## Tooetp
Animal_data_Tooetp <- Animal_data[!is.na(Animal_data$Tooetp),]
formula.A <- as.formula(paste("~",paste("Taxonomic_group","Tooetp",sep="+"), "-1", sep = ""))
formula.B <- as.formula(paste("~", paste("Taxonomic_group","Tooetp", sep = "+"), 
                              paste("+Taxonomic_group","Tooetp",sep = ":"), "-1", sep = ""))
model.null <- rma.mv(yi, vi, mods = ~Taxonomic_group-1, 
                     random = list(~1|Animal_Species, ~1|Code,~1|Hrbn),
                     R = list(Animal_Species = Animals_Cor[Animal_data$Animal_Species, Animal_data$Animal_Species]),
                     data = Animal_data_Tooetp,method = "ML")
model.A <- rma.mv(yi, vi, mods = formula.A,
                  random = list(~1|Animal_Species, ~1|Code,~1|Hrbn),
                  R = list(Animal_Species = Animals_Cor[Animal_data$Animal_Species, Animal_data$Animal_Species]),
                  data = Animal_data_Tooetp,
                  method = "ML")
model.B <- rma.mv(yi, vi, mods = formula.B,
                  random = list(~1|Animal_Species, ~1|Code,~1|Hrbn),
                  R = list(Animal_Species = Animals_Cor[Animal_data$Animal_Species, Animal_data$Animal_Species]),
                  data = Animal_data_Tooetp,
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
if(length(model.B$beta) != length(model.A$beta)){
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


formula.A <- as.formula(paste("~",paste("Taxonomic_group_response_category","Tooetp",sep="+"), "-1", sep = ""))
formula.B <- as.formula(paste("~", paste("Taxonomic_group_response_category","Tooetp", sep = "+"), 
                              paste("+Taxonomic_group_response_category","Tooetp",sep = ":"), "-1", sep = ""))
model.null2 <- rma.mv(yi, vi, mods = ~Taxonomic_group_response_category-1, 
                      random = list(~1|Animal_Species, ~1|Code,~1|Hrbn),
                      R = list(Animal_Species = Animals_Cor[Animal_data$Animal_Species, Animal_data$Animal_Species]),
                      data = Animal_data_Tooetp,method = "ML")
model.A <- rma.mv(yi, vi, mods = formula.A,
                  random = list(~1|Animal_Species, ~1|Code,~1|Hrbn),
                  R = list(Animal_Species = Animals_Cor[Animal_data$Animal_Species, Animal_data$Animal_Species]),
                  data = Animal_data_Tooetp,
                  method = "ML")
model.B <- rma.mv(yi, vi, mods = formula.B,
                  random = list(~1|Animal_Species, ~1|Code,~1|Hrbn),
                  R = list(Animal_Species = Animals_Cor[Animal_data$Animal_Species, Animal_data$Animal_Species]),
                  data = Animal_data_Tooetp,
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
if(length(model.B$beta) != length(model.A$beta)){
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
#####
#####------------- Cnoi
#####-------------
Animal_data_Cnoi <- Animal_data[!is.na(Animal_data$Cnoi),]
formula.A <- as.formula(paste("~",paste("Taxonomic_group","Cnoi",sep="+"), "-1", sep = ""))
formula.B <- as.formula(paste("~", paste("Taxonomic_group","Cnoi", sep = "+"), 
                              paste("+Taxonomic_group","Cnoi",sep = ":"), "-1", sep = ""))
model.null <- rma.mv(yi, vi, mods = ~Taxonomic_group-1, 
                     random = list(~1|Animal_Species, ~1|Code,~1|Hrbn),
                     R = list(Animal_Species = Animals_Cor[Animal_data$Animal_Species, Animal_data$Animal_Species]),
                     data = Animal_data_Cnoi,method = "ML")
model.A <- rma.mv(yi, vi, mods = formula.A,
                  random = list(~1|Animal_Species, ~1|Code,~1|Hrbn),
                  R = list(Animal_Species = Animals_Cor[Animal_data$Animal_Species, Animal_data$Animal_Species]),
                  data = Animal_data_Cnoi,
                  method = "ML")
model.B <- rma.mv(yi, vi, mods = formula.B,
                  random = list(~1|Animal_Species, ~1|Code,~1|Hrbn),
                  R = list(Animal_Species = Animals_Cor[Animal_data$Animal_Species, Animal_data$Animal_Species]),
                  data = Animal_data_Cnoi,
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
if(length(model.B$beta) != length(model.A$beta)){
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


formula.A <- as.formula(paste("~",paste("Taxonomic_group_response_category","Cnoi",sep="+"), "-1", sep = ""))
formula.B <- as.formula(paste("~", paste("Taxonomic_group_response_category","Cnoi", sep = "+"), 
                              paste("+Taxonomic_group_response_category","Cnoi",sep = ":"), "-1", sep = ""))
model.null2 <- rma.mv(yi, vi, mods = ~Taxonomic_group_response_category-1, 
                      random = list(~1|Animal_Species, ~1|Code,~1|Hrbn),
                      R = list(Animal_Species = Animals_Cor[Animal_data$Animal_Species, Animal_data$Animal_Species]),
                      data = Animal_data_Cnoi,method = "ML")
model.A <- rma.mv(yi, vi, mods = formula.A,
                  random = list(~1|Animal_Species, ~1|Code,~1|Hrbn),
                  R = list(Animal_Species = Animals_Cor[Animal_data$Animal_Species, Animal_data$Animal_Species]),
                  data = Animal_data_Cnoi,
                  method = "ML")
model.B <- rma.mv(yi, vi, mods = formula.B,
                  random = list(~1|Animal_Species, ~1|Code,~1|Hrbn),
                  R = list(Animal_Species = Animals_Cor[Animal_data$Animal_Species, Animal_data$Animal_Species]),
                  data = Animal_data_Cnoi,
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
if(length(model.B$beta) != length(model.A$beta)){
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



colnames(Results) <- c("Predictor variables",
                       "Ref.",
                       "AIC",
                       "L─L",
                       "Chisq",
                       "df",
                       "P",
                       "N")
write.csv(Results, "SuppTab1-4_1104.csv", row.names = F)

### SUppTab1-5
rm(list = ls())
setwd("/lustre/huyueqing/ssy/Wan2023/SuppTab1")
DataPre_path <- "/lustre/huyueqing/ssy/Wan2023/DataPre/"
library(readxl)
library(metafor)
library(ape)
library(stringr)

# data preprocess
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
Match_Mdnn <- match(Micro_data$McLn, Model_non$`Species or other groups`)
Micro_data[,"Model_non"] <- Model_non[Match_Mdnn,"Model or non-model species"]
Micro_data$Log_added_dosage <- log2(Micro_data$Hds. - Micro_data$Nhr.)
Micro_data$Taxonomic_group <- Micro_data$Pcom一
Micro_data$Taxonomic_group_response_category<- rep(1,times = nrow(Micro_data))
for (k in 1:nrow(Micro_data)){
  i <- Micro_data[k,"response_category"]
  response_category <- regmatches(i, regexec("growth|reproduction|biomarker",i))[[1]][1]
  if (Micro_data[k,"Taxonomic_group"] == "Bacteria"){
    flag <- paste0("Bacteria",sep = "_",response_category)
  }else{
    flag <- paste0("Fungus",sep = "_",response_category)
  }
  Micro_data[k,"Taxonomic_group_response_category"] <- flag
}
Micro_data$Tpcs <- rep(1,times=nrow(Micro_data))
for (k in 1:nrow(Micro_data)){
  i <- Micro_data[k,"Pcoh.一"]
  lable <- regmatches(i, regexec("Chemical synthetic|Biogenic|Mineral-based",i))[[1]][1]
  Micro_data[k,"Tpcs"] <- lable
}



Micro_Cor <- read.table(paste(DataPre_path,
                              "Micro_Cor.txt",
                              sep = ""))
colnames(Micro_Cor) <- rownames(Micro_Cor)
Micro_data$McLn <- str_to_lower(Micro_data$McLn)
Matched_Micro_InCor <- match(Micro_data$McLn,rownames(Micro_Cor))
Micro_data <- Micro_data[!is.na(Matched_Micro_InCor),]

Micro_Cor <- Micro_Cor[order(dimnames(Micro_Cor)[[1]]), order(dimnames(Micro_Cor)[[1]])]
Micro_data$Micro_Species <- factor(Micro_data$McLn)
levels(Micro_data$Micro_Species) <- sort(dimnames(Micro_Cor)[[1]])




Micro_data$Code <- factor(Micro_data$Code)
Micro_data$McLn <- factor(Micro_data$McLn)
Micro_data$Expt <- factor(Micro_data$Expt)
Micro_data$Mjcz <- factor(Micro_data$Mjcz)
Micro_data$pesticide_group <- factor(Micro_data$pesticide_group)
Micro_data$Tpcs <- factor(Micro_data$Tpcs)
Micro_data$Tooetp <- factor(Micro_data$Tooetp)
Micro_data$Old_New_pesticide <- factor(Micro_data$Old_New_pesticide)
Micro_data$Model_non <- factor(Micro_data$Model_non)
Micro_data$Cnoi <- factor(Micro_data$Cnoi)
Micro_data$Pbly <- factor(Micro_data$Pbly)
Micro_data$Hrbn <- factor(Micro_data$Hrbn)

Pre_vars <- c("Expt","Mjcz","pesticide_group",
              "Tpcs","Old_New_pesticide",
              "Model_non","Pbly","Log_added_dosage")

model.null <- rma.mv(yi, vi, mods = ~Taxonomic_group-1,
                     random = list(~1|Micro_Species, ~1|Code, ~1|Hrbn),
                     R = list(Micro_Species=Micro_Cor[Micro_data$Micro_Species, Micro_data$Micro_Species]),
                     data = Micro_data,
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

for (i in 1:length(Pre_vars)){
  formula.A <- as.formula(paste("~",paste("Taxonomic_group",Pre_vars[i],sep="+"), "-1", sep = ""))
  formula.B <- as.formula(paste("~", paste("Taxonomic_group",Pre_vars[i], sep = "+"),
                                paste("+Taxonomic_group",Pre_vars[i],sep = ":"), "-1", sep = ""))
  Tem_data <- Micro_data
  if(length(table(is.na(Micro_data[,Pre_vars[i]])))==2){
    Tem_data <- Micro_data[!is.na(Micro_data[,Pre_vars[i]]),]
  }
  if (Pre_vars[i] == "Mjcz"){
    Tem_data <- Tem_data[which(Tem_data[,"Expt"]=="Field"),]
    model.A <- rma.mv(yi, vi,
                      mods = formula.A,
                      random = list(~1|Micro_Species, ~1|Code,~1|Hrbn),
                      R = list(Micro_Species = Micro_Cor[Tem_data$Micro_Species, Tem_data$Micro_Species]),
                      data = Tem_data,
                      method = "ML")
    model.B <- rma.mv(yi, vi,
                      mods = formula.B,
                      random = list(~1|Micro_Species, ~1|Code,~1|Hrbn),
                      R = list(Micro_Species = Micro_Cor[Tem_data$Micro_Species, Tem_data$Micro_Species]),
                      data = Tem_data,
                      method = "ML")
    model.Mjcz <- rma.mv(yi, vi,
                         mods = ~Taxonomic_group-1,
                         random = list(~1|Micro_Species, ~1|Code,~1|Hrbn),
                         R = list(Micro_Species = Micro_Cor[Tem_data$Micro_Species, Tem_data$Micro_Species]),
                         data = Tem_data,
                         method = "ML")
    A1 <- anova(model.A,model.Mjcz)
    if(length(model.B$beta) != length(model.A$beta)){
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
  }else{
    model.A <- rma.mv(yi, vi,
                      mods = formula.A,
                      random = list(~1|Micro_Species, ~1|Code ,~1|Hrbn),
                      R = list(Micro_Species = Micro_Cor[Tem_data$Micro_Species, Tem_data$Micro_Species]),
                      data = Tem_data,
                      method = "ML")
    model.B <- rma.mv(yi, vi,
                      mods = formula.B,
                      random = list(~1|Micro_Species, ~1|Code,~1|Hrbn),
                      R = list(Micro_Species = Micro_Cor[Tem_data$Micro_Species, Tem_data$Micro_Species]),
                      data = Tem_data,
                      method = "ML")
    A1 <- anova(model.A, model.null)
    if(length(model.B$beta) != length(model.A$beta)){
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
  Results <- rbind(Results, data.frame(as.character(formula.A)[-1],
                                       "1",
                                       A1[["fit.stats.f"]][["AIC"]],
                                       A1[["fit.stats.f"]][["ll"]],
                                       A1[["LRT"]],
                                       A1[["parms.f"]],
                                       A1[["pval"]],
                                       summary(model.A)$k,
                                       fix.empty.names = F))

}

# Response category as covariant
model.null2 <- rma.mv(yi, vi, mods = ~Taxonomic_group_response_category-1,
                      random = list(~1|Micro_Species, ~1|Code, ~1|Hrbn),
                      R = list(Micro_Species = Micro_Cor[Micro_data$Micro_Species, Micro_data$Micro_Species]),
                      data = Micro_data,
                      method = "ML")
null21 <- anova(model.null2, model.null)
Results <- rbind(Results, data.frame("Predictor variables",
                                     "1",
                                     null21[["fit.stats.f"]][["AIC"]],
                                     null21[["fit.stats.f"]][["ll"]],
                                     null21[["LRT"]],
                                     null21[["parms.f"]],
                                     null21[["pval"]],
                                     summary(model.null2)$k,
                                     fix.empty.names = F))
for (i in 1:length(Pre_vars)){
  formula.A <- as.formula(paste("~", paste("Taxonomic_group_response_category",Pre_vars[i],sep = "+"),"-1",sep = ""))
  formula.B <- as.formula(paste("~",paste("Taxonomic_group_response_category",Pre_vars[i],sep = "+"),
                                paste("+Taxonomic_group_response_category",Pre_vars[i], sep = ":"),"-1",sep = ""))
  Tem_data <- Micro_data
  if(length(table(is.na(Micro_data[,Pre_vars[i]])))==2){
    Tem_data <- Micro_data[!is.na(Micro_data[,Pre_vars[i]]),]
  }
  if (Pre_vars[i] == "Mjcz"){
    Tem_data <- Tem_data[which(Tem_data[,"Expt"]=="Field"),]
    model.A <- rma.mv(yi, vi,
                      mods = formula.A,
                      random = list(~1|Micro_Species, ~1|Code,~1|Hrbn),
                      R = list(Micro_Species = Micro_Cor[Tem_data$Micro_Species, Tem_data$Micro_Species]),
                      data = Tem_data,
                      method = "ML")
    model.B <- rma.mv(yi, vi,
                      mods = formula.B,
                      random = list(~1|Micro_Species, ~1|Code,~1|Hrbn),
                      R = list(Micro_Species = Micro_Cor[Tem_data$Micro_Species, Tem_data$Micro_Species]),
                      data = Tem_data,
                      method = "ML")
    model.Mjcz <- rma.mv(yi, vi,
                         mods = ~Taxonomic_group-1,
                         random = list(~1|Micro_Species, ~1|Code,~1|Hrbn),
                         R = list(Micro_Species = Micro_Cor[Tem_data$Micro_Species, Tem_data$Micro_Species]),
                         data = Tem_data,
                         method = "ML")
    A1 <- anova(model.A,model.Mjcz)
    if(length(model.B$beta) != length(model.A$beta)){
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
  }else{
    model.A <- rma.mv(yi, vi,
                      mods = formula.A,
                      random = list(~1|Micro_Species, ~1|Code ,~1|Hrbn),
                      R = list(Micro_Species = Micro_Cor[Tem_data$Micro_Species, Tem_data$Micro_Species]),
                      data = Tem_data,
                      method = "ML")
    model.B <- rma.mv(yi, vi,
                      mods = formula.B,
                      random = list(~1|Micro_Species, ~1|Code,~1|Hrbn),
                      R = list(Micro_Species = Micro_Cor[Tem_data$Micro_Species, Tem_data$Micro_Species]),
                      data = Tem_data,
                      method = "ML")
    A1 <- anova(model.A, model.null)
    if(length(model.B$beta) != length(model.A$beta)){
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
  Results <- rbind(Results, data.frame(as.character(formula.A)[-1],
                                       "2",
                                       A1[["fit.stats.f"]][["AIC"]],
                                       A1[["fit.stats.f"]][["ll"]],
                                       A1[["LRT"]],
                                       A1[["parms.f"]],
                                       A1[["pval"]],
                                       summary(model.A)$k,
                                       fix.empty.names = F))
}

## Tooetp
Animal_data_Tooetp <- Animal_data[!is.na(Animal_data$Tooetp),]
formula.A <- as.formula(paste("~",paste("Taxonomic_group","Tooetp",sep="+"), "-1", sep = ""))
formula.B <- as.formula(paste("~", paste("Taxonomic_group","Tooetp", sep = "+"), 
                              paste("+Taxonomic_group","Tooetp",sep = ":"), "-1", sep = ""))
model.null <- rma.mv(yi, vi, mods = ~Taxonomic_group-1, 
                     random = list(~1|Micro_Species, ~1|Code,~1|Hrbn),
                     R = list(Micro_Species = Micro_Cor[Micro_data$Micro_Species, Micro_data$Micro_Species]),
                     data = Micro_data_Tooetp,method = "ML")
model.A <- rma.mv(yi, vi, mods = formula.A,
                  random = list(~1|Micro_Species, ~1|Code,~1|Hrbn),
                  R = list(Micro_Species = Micro_Cor[Micro_data$Micro_Species, Micro_data$Micro_Species]),
                  data = Micro_data_Tooetp,
                  method = "ML")
model.B <- rma.mv(yi, vi, mods = formula.B,
                  random = list(~1|Micro_Species, ~1|Code,~1|Hrbn),
                  R = list(Micro_Species = Micro_Cor[Micro_data$Micro_Species, Micro_data$Micro_Species]),
                  data = Micro_data_Tooetp,
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
if(length(model.B$beta) != length(model.A$beta)){
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


formula.A <- as.formula(paste("~",paste("Taxonomic_group_response_category","Tooetp",sep="+"), "-1", sep = ""))
formula.B <- as.formula(paste("~", paste("Taxonomic_group_response_category","Tooetp", sep = "+"), 
                              paste("+Taxonomic_group_response_category","Tooetp",sep = ":"), "-1", sep = ""))
model.null2 <- rma.mv(yi, vi, mods = ~Taxonomic_group_response_category-1, 
                      random = list(~1|Micro_Species, ~1|Code,~1|Hrbn),
                      R = list(Micro_Species = Micro_Cor[Micro_data$Micro_Species, Micro_data$Micro_Species]),
                      data = Micro_data_Tooetp,method = "ML")
model.A <- rma.mv(yi, vi, mods = formula.A,
                  random = list(~1|Micro_Species, ~1|Code,~1|Hrbn),
                  R = list(Micro_Species = Micro_Cor[Micro_data$Micro_Species, Micro_data$Micro_Species]),
                  data = Micro_data_Tooetp,
                  method = "ML")
model.B <- rma.mv(yi, vi, mods = formula.B,
                  random = list(~1|Micro_Species, ~1|Code,~1|Hrbn),
                  R = list(Micro_Species = Micro_Cor[Micro_data$Micro_Species, Micro_data$Micro_Species]),
                  data = Micro_data_Tooetp,
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
if(length(model.B$beta) != length(model.A$beta)){
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
#####
#####------------- Cnoi
#####-------------
Animal_data_Cnoi <- Animal_data[!is.na(Animal_data$Cnoi),]
formula.A <- as.formula(paste("~",paste("Taxonomic_group","Cnoi",sep="+"), "-1", sep = ""))
formula.B <- as.formula(paste("~", paste("Taxonomic_group","Cnoi", sep = "+"), 
                              paste("+Taxonomic_group","Cnoi",sep = ":"), "-1", sep = ""))
model.null <- rma.mv(yi, vi, mods = ~Taxonomic_group-1, 
                     random = list(~1|Micro_Species, ~1|Code,~1|Hrbn),
                     R = list(Micro_Species = Micro_Cor[Micro_data$Micro_Species, Micro_data$Micro_Species]),
                     data = Micro_data_Cnoi,method = "ML")
model.A <- rma.mv(yi, vi, mods = formula.A,
                  random = list(~1|Micro_Species, ~1|Code,~1|Hrbn),
                  R = list(Micro_Species = Micro_Cor[Micro_data$Micro_Species, Micro_data$Micro_Species]),
                  data = Micro_data_Cnoi,
                  method = "ML")
model.B <- rma.mv(yi, vi, mods = formula.B,
                  random = list(~1|Micro_Species, ~1|Code,~1|Hrbn),
                  R = list(Micro_Species = Micro_Cor[Micro_data$Micro_Species, Micro_data$Micro_Species]),
                  data = Micro_data_Cnoi,
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
if(length(model.B$beta) != length(model.A$beta)){
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


formula.A <- as.formula(paste("~",paste("Taxonomic_group_response_category","Cnoi",sep="+"), "-1", sep = ""))
formula.B <- as.formula(paste("~", paste("Taxonomic_group_response_category","Cnoi", sep = "+"), 
                              paste("+Taxonomic_group_response_category","Cnoi",sep = ":"), "-1", sep = ""))
model.null2 <- rma.mv(yi, vi, mods = ~Taxonomic_group_response_category-1, 
                      random = list(~1|Micro_Species, ~1|Code,~1|Hrbn),
                      R = list(Micro_Species = Micro_Cor[Micro_data$Micro_Species, Micro_data$Micro_Species]),
                      data = Microl_data_Cnoi,method = "ML")
model.A <- rma.mv(yi, vi, mods = formula.A,
                  random = list(~1|Micro_Species, ~1|Code,~1|Hrbn),
                  R = list(Micro_Species = Micro_Cor[Micro_data$Micro_Species, Micro_data$Micro_Species]),
                  data = Micro_data_Cnoi,
                  method = "ML")
model.B <- rma.mv(yi, vi, mods = formula.B,
                  random = list(~1|Micro_Species, ~1|Code,~1|Hrbn),
                  R = list(Micro_Species = Micro_Cor[Micro_data$Micro_Species, Micro_data$Micro_Species]),
                  data = Micro_data_Cnoi,
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
if(length(model.B$beta) != length(model.A$beta)){
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




colnames(Results) <- c("Predictor variables",
                       "Ref.",
                       "AIC",
                       "L─L",
                       "Chisq",
                       "df",
                       "P",
                       "N")
write.csv(Results, "SuppTab1-5_1104.csv", row.names = T)
