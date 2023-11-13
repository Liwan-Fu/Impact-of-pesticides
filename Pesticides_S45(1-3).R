rm(list = ls())
setwd("...")
DataPre_path <- "..."
library(readxl)
library(metafor)
library(ape)
library(rotl)

Whole_data <- read.table(paste(DataPre_path,"AllData_Tab1.txt",sep=""),
                        header = T)
Whole_data <- Whole_data[which(Whole_data$yi>=-2 & Whole_data$yi<=2),]

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
write.csv(Results, "SuppTab45-1_1111.csv", row.names = F)


### SuppTab45-2
rm(list = ls())
setwd("/lustre/huyueqing/ssy/Wan2023/SuppTab1")
DataPre_path <- "/lustre/huyueqing/ssy/Wan2023/DataPre/"
library(readxl)
library(metafor)
library(ape)
library(rotl)

Whole_data <- read.table(paste(DataPre_path,"AllData_Tab2.txt",sep=""),
                        header = T)
Whole_data <- Whole_data[which(Whole_data$yi>=-2 & Whole_data$yi<=2),]


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
write.csv(Results, "SuppTab45-2-1111.csv", row.names = F)

### SuppTab45-3
rm(list = ls())
setwd("/lustre/huyueqing/ssy/Wan2023/SuppTab1")
DataPre_path <- "/lustre/huyueqing/ssy/Wan2023/DataPre/"
library(readxl)
library(metafor)
library(ape)
library(rotl)
library(stringr)

Plant_data <- read.table(paste(DataPre_path,"AllPlant.txt",sep = ""),
                        header = T)
Plant_Cor <- read.table(paste(DataPre_path,"Plant_Cor.txt",sep=""),
                        header = T)

colnames(Plant_Cor) <- rownames(Plant_Cor)
Plant_data$PlLn <- str_to_lower(Plant_data$PlLn)
Matched_Plant_InCor <- match(Plant_data$PlLn, rownames(Plant_Cor))
Plant_data <- Plant_data[!is.na(Matched_Plant_InCor),]

Plant_Cor <- Plant_Cor[order(dimnames(Plant_Cor)[[1]]), order(dimnames(Plant_Cor)[[1]])]
Plant_data$Plant_Species <- factor(Plant_data$PlLn)
levels(Plant_data$Plant_Species) <- sort(dimnames(Plant_Cor)[[1]])

plant <- read_excel(paste(DataPre_path,"matched plants-2023-0519.xls",sep=""))
plant <- rbind(colnames(plant),plant[,1])
colnames(plant) <- "plant_species"
plant$PlLn <- str_to_lower(plant$PlLn)

Plant_data <- Plant_data[which(Plant_data$yi >= -2 & Plant_data$yi <=2),]
Plant_data <- Plant_data[match(Plant_data$McLn, plant$PlLn),]



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
write.csv(Results, "SuppTab45-3_1111.csv", row.names = F)
