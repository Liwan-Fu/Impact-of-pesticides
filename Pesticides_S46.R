rm(list = ls())
setwd("...")
DataPre_path <- "..."
library(readxl)
library(metafor)
library(ape)
library(stringr)
library(rotl)

## SuppTab2-1
AllData_1 <- read.table(paste(DataPre_path,"AllData_Tab1.txt",sep=""),
                        header = T)
AllData_1 <- AllData_1[which(AllData_1$yi>=-2 & AllData_1$yi<=2),]

AllData_1$Pbly <- factor(AllData_1$Pbly)
AllData_1$Cnoi <- factor(AllData_1$Cnoi)
AllData_1$SpLn <- factor(AllData_1$SpLn)
AllData_1$Hrbn <- factor(AllData_1$Hrbn)

model.null <- rma.mv(yi,vi,mods = ~Taxonomic_group-1,
                     random = list(~1|Code, ~1|Hrbn),
                     data = AllData_1,
                     method = "ML")
egger_dt <- data.frame(r = residuals(model.null),
                       s = sign(residuals(model.null))*
                         sqrt(diag(vcov(model.null, type = "resid"))))
egg_test <- lm(r~s, egger_dt)
res <- data.frame("Taxonomic Group",
                  nrow(egg_test$model),
                  coef(summary(egg_test))[1, 3],
                  coef(summary(egg_test))[1, 4],
                  fix.empty.names = F)
Pre_vars <- c("Cnoi","Pbly")
for (i in 1:length(Pre_vars)){
  formula.A <- as.formula(paste("~",paste("Taxonomic_group",Pre_vars[i],sep="+"), "-1", sep = ""))
  formula.B <- as.formula(paste("~", paste("Taxonomic_group",Pre_vars[i], sep = "+"),
                                paste("+Taxonomic_group",Pre_vars[i],sep = ":"), "-1", sep = ""))
  TemData <- AllData_1
  if (length(table(is.na(AllData_1[,Pre_vars[i]])))==2){
    TemData <- AllData_1[!is.na(AllData_1[,Pre_vars[i]]),]
  }
  model.A <- rma.mv(yi,vi,mods = formula.A,
                    random = list(~1|Code, ~1|Hrbn),
                    data = TemData,
                    method = "ML")
  model.B <- rma.mv(yi, vi, mods = formula.B,
                    random = list(~ 1|Code, ~1|Hrbn),
                    data = TemData,
                    method = "ML")
  egger_dt <- data.frame(r = residuals(model.A),
                         s = sign(residuals(model.A))*
                           sqrt(diag(vcov(model.A, type = "resid"))))
  egg_test <- lm(r~s, egger_dt)
  res <- rbind(res, data.frame(as.character(formula.A)[-1],
                               nrow(egg_test$model),
                               coef(summary(egg_test))[1, 3],
                               coef(summary(egg_test))[1, 4],
                               fix.empty.names = F))
  egger_dt <- data.frame(r = residuals(model.B),
                         s = sign(residuals(model.B))*
                           sqrt(diag(vcov(model.B, type = "resid"))))
  egg_test <- lm(r~s, egger_dt)
  res <- rbind(res, data.frame(as.character(formula.B)[-1],
                               nrow(egg_test$model),
                               coef(summary(egg_test))[1, 3],
                               coef(summary(egg_test))[1, 4],
                               fix.empty.names = F))
}

for (i in 1:length(Pre_vars)){
  formula.A <- as.formula(paste("~",paste("Taxonomic_group_response_category",Pre_vars[i],sep="+"), "-1", sep = ""))
  formula.B <- as.formula(paste("~", paste("Taxonomic_group_response_category",Pre_vars[i], sep = "+"),
                                paste("+Taxonomic_group_response_category",Pre_vars[i],sep = ":"), "-1", sep = ""))
  TemData <- AllData_1
  if (length(table(is.na(AllData_1[,Pre_vars[i]])))==2){
    TemData <- AllData_1[!is.na(AllData_1[,Pre_vars[i]]),]
  }
  model.A <- rma.mv(yi,vi,mods = formula.A,
                    random = list(~1|Code, ~1|Hrbn),
                    data = TemData,
                    method = "ML")
  model.B <- rma.mv(yi, vi, mods = formula.B,
                    random = list(~ 1|Code, ~1|Hrbn),
                    data = TemData,
                    method = "ML")
  egger_dt <- data.frame(r = residuals(model.A),
                         s = sign(residuals(model.A))*
                           sqrt(diag(vcov(model.A, type = "resid"))))
  egg_test <- lm(r~s, egger_dt)
  res <- rbind(res, data.frame(as.character(formula.A)[-1],
                               nrow(egg_test$model),
                               coef(summary(egg_test))[1, 3],
                               coef(summary(egg_test))[1, 4],
                               fix.empty.names = F))
  egger_dt <- data.frame(r = residuals(model.B),
                         s = sign(residuals(model.B))*
                           sqrt(diag(vcov(model.B, type = "resid"))))
  egg_test <- lm(r~s, egger_dt)
  res <- rbind(res, data.frame(as.character(formula.B)[-1],
                               nrow(egg_test$model),
                               coef(summary(egg_test))[1, 3],
                               coef(summary(egg_test))[1, 4],
                               fix.empty.names = F))
}

colnames(res) <- c("Predictor variables",
                   "N",
                   "test-value",
                   "P-value")
write.csv(res, "SuppTab46-1_1111.csv", row.names = F)


## SuppTab2-2
rm(list = ls())
setwd("/lustre/huyueqing/ssy/Wan2023/SuppTab2")
DataPre_path <- "/lustre/huyueqing/ssy/Wan2023/DataPre/"
library(readxl)
library(metafor)
library(ape)
library(stringr)
library(rotl)

AllData_2 <- read.table(paste(DataPre_path,"AllData_Tab2.txt",sep=""),
                        header = T)
AllData_2 <- AllData_2[which(AllData_2$yi>=-2 & AllData_2$yi<=2),]

AllData_2$Pbly <- factor(AllData_2$Pbly)
AllData_2$Cnoi <- factor(AllData_2$Cnoi)
AllData_2$SpLn <- factor(AllData_2$SpLn)
AllData_2$Hrbn <- factor(AllData_2$Hrbn)

model.null <- rma.mv(yi, vi, mods = ~Taxonomic_group-1,
                     random = list(~1|Code, ~1|Hrbn),
                     data = AllData_2,
                     method = "ML")
egger_dt <- data.frame(r = residuals(model.null),
                       s = sign(residuals(model.null))*
                         sqrt(diag(vcov(model.null, type = "resid"))))
egg_test <- lm(r~s, egger_dt)
res <- data.frame("Taxonomic Group",
                  nrow(egg_test$model),
                  coef(summary(egg_test))[1, 3],
                  coef(summary(egg_test))[1, 4],
                  fix.empty.names = F)

Pre_vars <- c("Cnoi","Pbly")

for (i in 1:length(Pre_vars)){
  formula.A <- as.formula(paste("~",paste("Taxonomic_group",Pre_vars[i],sep="+"), "-1", sep = ""))
  formula.B <- as.formula(paste("~", paste("Taxonomic_group",Pre_vars[i], sep = "+"),
                                paste("+Taxonomic_group",Pre_vars[i],sep = ":"), "-1", sep = ""))
  TemData <- AllData_2
  if (length(table(is.na(AllData_2[,Pre_vars[i]])))==2){
    TemData <- AllData_2[!is.na(AllData_2[,Pre_vars[i]]),]
  }
  model.A <- rma.mv(yi,vi,mods = formula.A,
                    random = list(~1|Code, ~1|Hrbn),
                    data = TemData,
                    method = "ML")
  model.B <- rma.mv(yi, vi, mods = formula.B,
                    random = list(~ 1|Code, ~1|Hrbn),
                    data = TemData,
                    method = "ML")
  egger_dt <- data.frame(r = residuals(model.A),
                         s = sign(residuals(model.A))*
                           sqrt(diag(vcov(model.A, type = "resid"))))
  egg_test <- lm(r~s, egger_dt)
  res <- rbind(res, data.frame(as.character(formula.A)[-1],
                               nrow(egg_test$model),
                               coef(summary(egg_test))[1, 3],
                               coef(summary(egg_test))[1, 4],
                               fix.empty.names = F))
  egger_dt <- data.frame(r = residuals(model.B),
                         s = sign(residuals(model.B))*
                           sqrt(diag(vcov(model.B, type = "resid"))))
  egg_test <- lm(r~s, egger_dt)
  res <- rbind(res, data.frame(as.character(formula.B)[-1],
                               nrow(egg_test$model),
                               coef(summary(egg_test))[1, 3],
                               coef(summary(egg_test))[1, 4],
                               fix.empty.names = F))
}

for (i in 1:length(Pre_vars)){
  formula.A <- as.formula(paste("~",paste("Taxonomic_group_response_category",Pre_vars[i],sep="+"), "-1", sep = ""))
  formula.B <- as.formula(paste("~", paste("Taxonomic_group_response_category",Pre_vars[i], sep = "+"),
                                paste("+Taxonomic_group_response_category",Pre_vars[i],sep = ":"), "-1", sep = ""))
  TemData <- AllData_2
  if (length(table(is.na(AllData_2[,Pre_vars[i]])))==2){
    TemData <- AllData_2[!is.na(AllData_2[,Pre_vars[i]]),]
  }
  model.A <- rma.mv(yi,vi,mods = formula.A,
                    random = list(~1|Code, ~1|Hrbn),
                    data = TemData,
                    method = "ML")
  model.B <- rma.mv(yi, vi, mods = formula.B,
                    random = list(~ 1|Code, ~1|Hrbn),
                    data = TemData,
                    method = "ML")
  egger_dt <- data.frame(r = residuals(model.A),
                         s = sign(residuals(model.A))*
                           sqrt(diag(vcov(model.A, type = "resid"))))
  egg_test <- lm(r~s, egger_dt)
  res <- rbind(res, data.frame(as.character(formula.A)[-1],
                               nrow(egg_test$model),
                               coef(summary(egg_test))[1, 3],
                               coef(summary(egg_test))[1, 4],
                               fix.empty.names = F))
  egger_dt <- data.frame(r = residuals(model.B),
                         s = sign(residuals(model.B))*
                           sqrt(diag(vcov(model.B, type = "resid"))))
  egg_test <- lm(r~s, egger_dt)
  res <- rbind(res, data.frame(as.character(formula.B)[-1],
                               nrow(egg_test$model),
                               coef(summary(egg_test))[1, 3],
                               coef(summary(egg_test))[1, 4],
                               fix.empty.names = F))
}
colnames(res) <- c("Predictor variables",
                   "N",
                   "test-value",
                   "P-value")
write.csv(res, "SuppTab46-2_1111.csv", row.names = F)


##SuppTab2-3
rm(list = ls())
setwd("/lustre/huyueqing/ssy/Wan2023/SuppTab2")
DataPre_path <- "/lustre/huyueqing/ssy/Wan2023/DataPre/"
library(readxl)
library(metafor)
library(ape)
library(stringr)
library(rotl)

PlantData <- read.table(paste(DataPre_path,"AllPlant.txt",sep = ""),
                        header = T)
Plant_Cor <- read.table(paste(DataPre_path,"Plant_Cor.txt",sep=""),
                        header = T)

colnames(Plant_Cor) <- rownames(Plant_Cor)
PlantData$PlLn <- str_to_lower(PlantData$PlLn)
Matched_Plant_InCor <- match(PlantData$PlLn, rownames(Plant_Cor))
PlantData <- PlantData[!is.na(Matched_Plant_InCor),]

Plant_Cor <- Plant_Cor[order(dimnames(Plant_Cor)[[1]]), order(dimnames(Plant_Cor)[[1]])]
PlantData$Plant_Species <- factor(PlantData$PlLn)
levels(PlantData$Plant_Species) <- sort(dimnames(Plant_Cor)[[1]])

plant <- read_excel(paste(DataPre_path,"matched plants-2023-0519.xls",sep=""))
plant <- rbind(colnames(plant),plant[,1])
colnames(plant) <- "plant_species"
plant$PlLn <- str_to_lower(plant$PlLn)

PlantData <- PlantData[which(PlantData$yi >= -2 & PlantData$yi <=2),]
PlantData <- PlantData[match(PlantData$McLn, plant$PlLn),]



PlantData$Code <- factor(PlantData$Code)
PlantData$PlLn <- factor(PlantData$PlLn)
PlantData$Cnoi <- factor(PlantData$Cnoi)
PlantData$Pbly <- factor(PlantData$Pbly)
PlantData$Hrbn <- factor(PlantData$Hrbn)

model.null <- rma.mv(yi, vi, mods = ~Taxonomic_group-1,
                     random = list(~1|Plant_Species,~1|Code, ~1|Hrbn),
                     R = list(Plant_Species = Plant_Cor[PlantData$Plant_Species,PlantData$Plant_Species]),
                     data = PlantData,
                     method = "ML")
egger_dt <- data.frame(r = residuals(model.null),
                       s = sign(residuals(model.null))*
                         sqrt(diag(vcov(model.null, type = "resid"))))
egg_test <- lm(r~s, egger_dt)
res <- data.frame("Taxonomic Group",
                  nrow(egg_test$model),
                  coef(summary(egg_test))[1, 3],
                  coef(summary(egg_test))[1, 4],
                  fix.empty.names = F)

Pre_vars <- c("Cnoi","Pbly")

for (i in 1:length(Pre_vars)){
  formula.A <- as.formula(paste("~",paste("Taxonomic_group",Pre_vars[i],sep="+"), "-1", sep = ""))
  formula.B <- as.formula(paste("~", paste("Taxonomic_group",Pre_vars[i], sep = "+"),
                                paste("+Taxonomic_group",Pre_vars[i],sep = ":"), "-1", sep = ""))
  TemData <- PlantData
  if (length(table(is.na(PlantData[,Pre_vars[i]])))==2){
    TemData <- PlantData[!is.na(PlantData[,Pre_vars[i]]),]
  }
  model.A <- rma.mv(yi,vi,mods = formula.A,
                    random = list(~1|Plant_Species,~1|Code, ~1|Hrbn),
                    R = list(Plant_Species = Plant_Cor[TemData$Plant_Species,TemData$Plant_Species]),
                    data = TemData,
                    method = "ML")
  model.B <- rma.mv(yi, vi, mods = formula.B,
                    random = list(~1|Plant_Species,~1|Code, ~1|Hrbn),
                    R = list(Plant_Species = Plant_Cor[TemData$Plant_Species,TemData$Plant_Species]),
                    data = TemData,
                    method = "ML")
  egger_dt <- data.frame(r = residuals(model.A),
                         s = sign(residuals(model.A))*
                           sqrt(diag(vcov(model.A, type = "resid"))))
  egg_test <- lm(r~s, egger_dt)
  res <- rbind(res, data.frame(as.character(formula.A)[-1],
                               nrow(egg_test$model),
                               coef(summary(egg_test))[1, 3],
                               coef(summary(egg_test))[1, 4],
                               fix.empty.names = F))
  egger_dt <- data.frame(r = residuals(model.B),
                         s = sign(residuals(model.B))*
                           sqrt(diag(vcov(model.B, type = "resid"))))
  egg_test <- lm(r~s, egger_dt)
  res <- rbind(res, data.frame(as.character(formula.B)[-1],
                               nrow(egg_test$model),
                               coef(summary(egg_test))[1, 3],
                               coef(summary(egg_test))[1, 4],
                               fix.empty.names = F))
}

for (i in 1:length(Pre_vars)){
  formula.A <- as.formula(paste("~",paste("Taxonomic_group_response_category",Pre_vars[i],sep="+"), "-1", sep = ""))
  formula.B <- as.formula(paste("~", paste("Taxonomic_group_response_category",Pre_vars[i], sep = "+"),
                                paste("+Taxonomic_group_response_category",Pre_vars[i],sep = ":"), "-1", sep = ""))
  TemData <- PlantData
  if (length(table(is.na(PlantData[,Pre_vars[i]])))==2){
    TemData <- PlantData[!is.na(PlantData[,Pre_vars[i]]),]
  }
  model.A <- rma.mv(yi,vi,mods = formula.A,
                    random = list(~1|Plant_Species,~1|Code, ~1|Hrbn),
                    R = list(Plant_Species = Plant_Cor[TemData$Plant_Species,TemData$Plant_Species]),
                    data = TemData,
                    method = "ML")
  model.B <- rma.mv(yi, vi, mods = formula.B,
                    random = list(~1|Plant_Species,~1|Code, ~1|Hrbn),
                    R = list(Plant_Species = Plant_Cor[TemData$Plant_Species,TemData$Plant_Species]),
                    data = TemData,
                    method = "ML")
  egger_dt <- data.frame(r = residuals(model.A),
                         s = sign(residuals(model.A))*
                           sqrt(diag(vcov(model.A, type = "resid"))))
  egg_test <- lm(r~s, egger_dt)
  res <- rbind(res, data.frame(as.character(formula.A)[-1],
                               nrow(egg_test$model),
                               coef(summary(egg_test))[1, 3],
                               coef(summary(egg_test))[1, 4],
                               fix.empty.names = F))
  egger_dt <- data.frame(r = residuals(model.B),
                         s = sign(residuals(model.B))*
                           sqrt(diag(vcov(model.B, type = "resid"))))
  egg_test <- lm(r~s, egger_dt)
  res <- rbind(res, data.frame(as.character(formula.B)[-1],
                               nrow(egg_test$model),
                               coef(summary(egg_test))[1, 3],
                               coef(summary(egg_test))[1, 4],
                               fix.empty.names = F))
}

colnames(res) <- c("Predictor variables",
                   "N",
                   "test-value",
                   "P-value")
write.csv(res, "SuppTab46-3_1111.csv", row.names = F)


## SuppTab2-4
rm(list = ls())
setwd("/lustre/huyueqing/ssy/Wan2023/SuppTab2")
DataPre_path <- "/lustre/huyueqing/ssy/Wan2023/DataPre/"
library(readxl)
library(metafor)
library(ape)
library(stringr)
library(rotl)

AnimalData <- read.table(paste(DataPre_path,"AllAnimal.txt",sep = ""),
                         header = T)
Animal_Cor <- read.table(paste(DataPre_path,"Animal_Cor.txt",sep = ""),
                         header = T)

colnames(Animal_Cor) <- rownames(Animal_Cor)
AnimalData$AnLn <- str_to_lower(AnimalData$AnLn)
Matched_Plant_InCor <- match(AnimalData$AnLn, rownames(Animal_Cor))
AnimalData <- AnimalData[!is.na(Matched_Plant_InCor),]

Animal_Cor <- Animal_Cor[order(dimnames(Animal_Cor)[[1]]), order(dimnames(Animal_Cor)[[1]])]
AnimalData$Plant_Species <- factor(AnimalData$AnLn)
levels(AnimalData$Plant_Species) <- sort(dimnames(Animal_Cor)[[1]])

animal<-read_excel(paste(DataPre_path,"animal species -612种.xls",sep=""))
animal <- rbind(colnames(animal), animal[,1])
colnames(animal) <- "animal_species"

AnimalData <- AnimalData[which(AnimalData$yi >= -2 & AnimalData$yi <=2),]
AnimalData <- AnimalData[match(AnimalData$McLn, animal$animal_species),]


AnimalData$Code <- factor(AnimalData$Code)
AnimalData$AnLn <- factor(AnimalData$AnLn)
AnimalData$Taxonomic_group <- factor(AnimalData$Taxonomic_group)
AnimalData$Mjcz <- factor(AnimalData$Mjcz)
AnimalData$Expt <- factor(AnimalData$Expt)
AnimalData$pesticide_group <- factor(AnimalData$pesticide_group)
AnimalData$Tpcs <- factor(AnimalData$Tpcs)
AnimalData$Tooetp <- factor(AnimalData$Tooetp)
AnimalData$Old_New_pesticide <- factor(AnimalData$Old_New_pesticide)
AnimalData$Model_non <- factor(AnimalData$Model_non)
AnimalData$Cnoi <- factor(AnimalData$Cnoi)
AnimalData$Pbly <- factor(AnimalData$Pbly)
AnimalData$Hrbn <- factor(AnimalData$Hrbn)

Pre_vars <- c("Expt","Mjcz","pesticide_group",
              "Tpcs","Tooetp","Old_New_pesticide",
              "Model_non","Cnoi","Pbly","Log_added_dosage")

model.null <- rma.mv(yi, vi, mods = ~Taxonomic_group-1,
                     random = list(~1|Animal_Species,~1|Code, ~1|Hrbn),
                     R = list(Animal_Species = Animal_Cor[AnimalData$Animal_Species,AnimalData$Animal_Species]),
                     data = AnimalData,
                     method = "ML")
egger_dt <- data.frame(r = residuals(model.null),
                       s = sign(residuals(model.null))*
                         sqrt(diag(vcov(model.null, type = "resid"))))
egg_test <- lm(r~s, egger_dt)
res <- data.frame("Taxonomic Group",
                  nrow(egg_test$model),
                  coef(summary(egg_test))[1, 3],
                  coef(summary(egg_test))[1, 4],
                  fix.empty.names = F)

for (i in 1:length(Pre_vars)){
  formula.A <- as.formula(paste("~",paste("Taxonomic_group",Pre_vars[i],sep="+"), "-1", sep = ""))
  formula.B <- as.formula(paste("~", paste("Taxonomic_group",Pre_vars[i], sep = "+"),
                                paste("+Taxonomic_group",Pre_vars[i],sep = ":"), "-1", sep = ""))
  TemData <- AnimalData
  if (length(table(is.na(AnimalData[,Pre_vars[i]])))==2){
    TemData <- AnimalData[!is.na(AnimalData[,Pre_vars[i]]),]
  }
  if (Pre_vars[i]=="Mjcz"){
    TemData <- TemData[which(TemData$Expt=="Field"),]
  }
  model.A <- rma.mv(yi,vi,mods = formula.A,
                    random = list(~1|Animal_Species,~1|Code, ~1|Hrbn),
                    R = list(Animal_Species = Animal_Cor[TemData$Animal_Species,TemData$Animal_Species]),
                    data = TemData,
                    method = "ML")
  model.B <- rma.mv(yi, vi, mods = formula.B,
                    random = list(~1|Animal_Species,~1|Code, ~1|Hrbn),
                    R = list(Animal_Species = Animal_Cor[TemData$Animal_Species,TemData$Animal_Species]),
                    data = TemData,
                    method = "ML")
  egger_dt <- data.frame(r = residuals(model.A),
                         s = sign(residuals(model.A))*
                           sqrt(diag(vcov(model.A, type = "resid"))))
  egg_test <- lm(r~s, egger_dt)
  res <- rbind(res, data.frame(as.character(formula.A)[-1],
                               nrow(egg_test$model),
                               coef(summary(egg_test))[1, 3],
                               coef(summary(egg_test))[1, 4],
                               fix.empty.names = F))
  egger_dt <- data.frame(r = residuals(model.B),
                         s = sign(residuals(model.B))*
                           sqrt(diag(vcov(model.B, type = "resid"))))
  egg_test <- lm(r~s, egger_dt)
  res <- rbind(res, data.frame(as.character(formula.B)[-1],
                               nrow(egg_test$model),
                               coef(summary(egg_test))[1, 3],
                               coef(summary(egg_test))[1, 4],
                               fix.empty.names = F))
}

model.null2 <- rma.mv(yi, vi, mods = ~Taxonomic_group_response_category-1,
                      random = list(~1|Animal_Species,~1|Code, ~1|Hrbn),
                      R = list(Animal_Species = Animal_Cor[AnimalData$Animal_Species,AnimalData$Animal_Species]),
                      data = AnimalData,
                      method = "ML")
egger_dt <- data.frame(r = residuals(model.null),
                       s = sign(residuals(model.null))*
                         sqrt(diag(vcov(model.null, type = "resid"))))
egg_test <- lm(r~s, egger_dt)
res <- rbind(res,data.frame("Taxonomic Group Response Category",
                            nrow(egg_test$model),
                            coef(summary(egg_test))[1, 3],
                            coef(summary(egg_test))[1, 4],
                            fix.empty.names = F))

for (i in 1:length(Pre_vars)){
  formula.A <- as.formula(paste("~",paste("Taxonomic_group_response_category",Pre_vars[i],sep="+"), "-1", sep = ""))
  formula.B <- as.formula(paste("~", paste("Taxonomic_group_response_category",Pre_vars[i], sep = "+"),
                                paste("+Taxonomic_group_response_category",Pre_vars[i],sep = ":"), "-1", sep = ""))
  TemData <- AnimalData
  if (length(table(is.na(AnimalData[,Pre_vars[i]])))==2){
    TemData <- AnimalData[!is.na(AnimalData[,Pre_vars[i]]),]
  }
  if (Pre_vars[i]=="Mjcz"){
    TemData <- TemData[which(TemData$Expt=="Field"),]
  }
  model.A <- rma.mv(yi,vi,mods = formula.A,
                    random = list(~1|Animal_Species,~1|Code, ~1|Hrbn),
                    R = list(Animal_Species = Animal_Cor[TemData$Animal_Species,TemData$Animal_Species]),
                    data = TemData,
                    method = "ML")
  model.B <- rma.mv(yi, vi, mods = formula.B,
                    random = list(~1|Animal_Species,~1|Code, ~1|Hrbn),
                    R = list(Animal_Species = Animal_Cor[TemData$Animal_Species,TemData$Animal_Species]),
                    data = TemData,
                    method = "ML")
  egger_dt <- data.frame(r = residuals(model.A),
                         s = sign(residuals(model.A))*
                           sqrt(diag(vcov(model.A, type = "resid"))))
  egg_test <- lm(r~s, egger_dt)
  res <- rbind(res, data.frame(as.character(formula.A)[-1],
                               nrow(egg_test$model),
                               coef(summary(egg_test))[1, 3],
                               coef(summary(egg_test))[1, 4],
                               fix.empty.names = F))
  egger_dt <- data.frame(r = residuals(model.B),
                         s = sign(residuals(model.B))*
                           sqrt(diag(vcov(model.B, type = "resid"))))
  egg_test <- lm(r~s, egger_dt)
  res <- rbind(res, data.frame(as.character(formula.B)[-1],
                               nrow(egg_test$model),
                               coef(summary(egg_test))[1, 3],
                               coef(summary(egg_test))[1, 4],
                               fix.empty.names = F))
}

colnames(res) <- c("Predictor variables",
                   "N",
                   "test-value",
                   "P-value")
write.csv(res, "SuppTab46-4_1111.csv", row.names = F)

##SuppTab2-5
rm(list = ls())
setwd("/lustre/huyueqing/ssy/Wan2023/SuppTab2")
DataPre_path <- "/lustre/huyueqing/ssy/Wan2023/DataPre/"
library(readxl)
library(metafor)
library(ape)
library(stringr)
library(rotl)

MicroData <- read.table(paste(DataPre_path,"AllMicro.txt",sep = ""),
                        header = T)
Micro_Cor <- read.table(paste(DataPre_path,"Micro_Cor.txt",sep = ""),
                        header = T)

colnames(Micro_Cor) <- rownames(Micro_Cor)
MicroData$McLn <- str_to_lower(MicroData$McLn)
Matched_Micro_InCor <- match(MicroData$McLn,rownames(Micro_Cor))
MicroData <- MicroData[!is.na(Matched_Micro_InCor),]

Micro_Cor <- Micro_Cor[order(dimnames(Micro_Cor)[[1]]), order(dimnames(Micro_Cor)[[1]])]
MicroData$Micro_Species <- factor(MicroData$McLn)
levels(MicroData$Micro_Species) <- sort(dimnames(Micro_Cor)[[1]])

micro <- read_excel(paste(DataPre_path,"micro species - 96种.xls",sep=""))
micro <- rbind(colnames(micro), micro[,1])
colnames(micro) <- "micro_species"
micro$micro_species <- str_to_lower(micro$micro_species)

MicroData <- MicroData[which(MicroData$yi >= -2 & MicroData$yi <=2),]
MicroData <- MicroData[match(MicroData$McLn, micro$micro_species),]

MicroData$Code <- factor(MicroData$Code)
MicroData$McLn <- factor(MicroData$McLn)
MicroData$Expt <- factor(MicroData$Expt)
MicroData$Mjcz <- factor(MicroData$Mjcz)
MicroData$pesticide_group <- factor(MicroData$pesticide_group)
MicroData$Tpcs <- factor(MicroData$Tpcs)
MicroData$Tooetp <- factor(MicroData$Tooetp)
MicroData$Old_New_pesticide <- factor(MicroData$Old_New_pesticide)
MicroData$Model_non <- factor(MicroData$Model_non)
MicroData$Cnoi <- factor(MicroData$Cnoi)
MicroData$Pbly <- factor(MicroData$Pbly)
MicroData$Hrbn <- factor(MicroData$Hrbn)

Pre_vars <- c("Expt","Mjcz","pesticide_group",
              "Tpcs","Tooetp","Old_New_pesticide",
              "Model_non","Cnoi","Pbly","Log_added_dosage")

model.null <- rma.mv(yi, vi, mods = ~Taxonomic_group-1,
                     random = list(~1|Micro_Species,~1|Code, ~1|Hrbn),
                     R = list(Micro_Species = Micro_Cor[MicroData$Micro_Species, MicroData$Micro_Species]),
                     data = MicroData,
                     method = "ML")
egger_dt <- data.frame(r = residuals(model.null),
                       s = sign(residuals(model.null))*
                         sqrt(diag(vcov(model.null, type = "resid"))))
egg_test <- lm(r~s, egger_dt)
res <- data.frame("Taxonomic Group",
                  nrow(egg_test$model),
                  coef(summary(egg_test))[1, 3],
                  coef(summary(egg_test))[1, 4],
                  fix.empty.names = F)

for (i in 1:length(Pre_vars)){
  formula.A <- as.formula(paste("~",paste("Taxonomic_group",Pre_vars[i],sep="+"), "-1", sep = ""))
  formula.B <- as.formula(paste("~", paste("Taxonomic_group",Pre_vars[i], sep = "+"),
                                paste("+Taxonomic_group",Pre_vars[i],sep = ":"), "-1", sep = ""))
  TemData <- MicroData
  if (length(table(is.na(MicroData[,Pre_vars[i]])))==2){
    TemData <- MicroData[!is.na(MicroData[,Pre_vars[i]]),]
  }
  if (Pre_vars[i]=="Mjcz"){
    TemData <- TemData[which(TemData$Expt=="Field"),]
  }
  model.A <- rma.mv(yi,vi,mods = formula.A,
                    random = list(~1|Micro_Species,~1|Code, ~1|Hrbn),
                    R = list(Micro_Species = Micro_Cor[TemData$Micro_Species,TemData$Micro_Species]),
                    data = TemData,
                    method = "ML")
  model.B <- rma.mv(yi, vi, mods = formula.B,
                    random = list(~1|Micro_Species,~1|Code, ~1|Hrbn),
                    R = list(Micro_Species = Micro_Cor[TemData$Micro_Species,TemData$Micro_Species]),
                    data = TemData,
                    method = "ML")
  egger_dt <- data.frame(r = residuals(model.A),
                         s = sign(residuals(model.A))*
                           sqrt(diag(vcov(model.A, type = "resid"))))
  egg_test <- lm(r~s, egger_dt)
  res <- rbind(res, data.frame(as.character(formula.A)[-1],
                               nrow(egg_test$model),
                               coef(summary(egg_test))[1, 3],
                               coef(summary(egg_test))[1, 4],
                               fix.empty.names = F))
  egger_dt <- data.frame(r = residuals(model.B),
                         s = sign(residuals(model.B))*
                           sqrt(diag(vcov(model.B, type = "resid"))))
  egg_test <- lm(r~s, egger_dt)
  res <- rbind(res, data.frame(as.character(formula.B)[-1],
                               nrow(egg_test$model),
                               coef(summary(egg_test))[1, 3],
                               coef(summary(egg_test))[1, 4],
                               fix.empty.names = F))
}

model.null2 <- rma.mv(yi, vi, mods = ~Taxonomic_group_response_category-1,
                      random = list(~1|Micro_Species,~1|Code, ~1|Hrbn),
                      R = list(Micro_Species = Micro_Cor[MicroData$Micro_Species,MicroData$Micro_Species]),
                      data = MicroData,
                      method = "ML")
egger_dt <- data.frame(r = residuals(model.null),
                       s = sign(residuals(model.null))*
                         sqrt(diag(vcov(model.null, type = "resid"))))
egg_test <- lm(r~s, egger_dt)
res <- rbind(res,data.frame("Taxonomic Group Response Category",
                            nrow(egg_test$model),
                            coef(summary(egg_test))[1, 3],
                            coef(summary(egg_test))[1, 4],
                            fix.empty.names = F))

for (i in 1:length(Pre_vars)){
  formula.A <- as.formula(paste("~",paste("Taxonomic_group_response_category",Pre_vars[i],sep="+"), "-1", sep = ""))
  formula.B <- as.formula(paste("~", paste("Taxonomic_group_response_category",Pre_vars[i], sep = "+"),
                                paste("+Taxonomic_group_response_category",Pre_vars[i],sep = ":"), "-1", sep = ""))
  TemData <- MicroData
  if (length(table(is.na(MicroData[,Pre_vars[i]])))==2){
    TemData <- MicroData[!is.na(MicroData[,Pre_vars[i]]),]
  }
  if (Pre_vars[i]=="Mjcz"){
    TemData <- TemData[which(TemData$Expt=="Field"),]
  }
  model.A <- rma.mv(yi,vi,mods = formula.A,
                    random = list(~1|Micro_Species,~1|Code, ~1|Hrbn),
                    R = list(Micro_Species = Micro_Cor[TemData$Micro_Species,TemData$Micro_Species]),
                    data = TemData,
                    method = "ML")
  model.B <- rma.mv(yi, vi, mods = formula.B,
                    random = list(~1|Micro_Species,~1|Code, ~1|Hrbn),
                    R = list(Micro_Species = Micro_Cor[TemData$Micro_Species,TemData$Micro_Species]),
                    data = TemData,
                    method = "ML")
  egger_dt <- data.frame(r = residuals(model.A),
                         s = sign(residuals(model.A))*
                           sqrt(diag(vcov(model.A, type = "resid"))))
  egg_test <- lm(r~s, egger_dt)
  res <- rbind(res, data.frame(as.character(formula.A)[-1],
                               nrow(egg_test$model),
                               coef(summary(egg_test))[1, 3],
                               coef(summary(egg_test))[1, 4],
                               fix.empty.names = F))
  egger_dt <- data.frame(r = residuals(model.B),
                         s = sign(residuals(model.B))*
                           sqrt(diag(vcov(model.B, type = "resid"))))
  egg_test <- lm(r~s, egger_dt)
  res <- rbind(res, data.frame(as.character(formula.B)[-1],
                               nrow(egg_test$model),
                               coef(summary(egg_test))[1, 3],
                               coef(summary(egg_test))[1, 4],
                               fix.empty.names = F))
}

colnames(res) <- c("Predictor variables",
                   "N",
                   "test-value",
                   "P-value")
write.csv(res, "SuppTab46-5_1111.csv", row.names = F)


