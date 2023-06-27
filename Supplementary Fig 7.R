library(metafor)
library(xlsx)
library(Hmisc)
library(export)
setwd('...')
combined=read.csv('data for recomendation dose.csv',header=T)

########## Fig. S7 ###########
S32=matrix(NA,40,3)
Pesticide_category=c('Pesticide',
                     unique(combined$pesticide_by_target_organisms))
dat_list1=lapply(unique(combined$Taxonomic_group_response),
                   function(x){subset(combined,Taxonomic_group_response==x)})
for (i in 1:length(Pesticide_category)) {
  for (j in 1:length(dat_list1)){
    mydat=subset(dat_list1[[j]],
                 pesticide_by_target_organisms==Pesticide_category[i])
    if (nrow(mydat)<=1) {
      S32[(i-1)*10+j,]=rep(NA,3)
    } else {
      res=rma.mv(yi, vi, 
                 random = list(~ 1 | factor(Insecticide.name),~1|Code)
                 data=mydat)
      S32[(i-1)*10+j,]=c(res$b,res$ci.lb,res$ci.ub)
    }
  }
}
colnames(S32)<-c("Effect size","CL.lb","CL.ub")

### Fig. S7
errbar(1:nrow(Efig7),Efig7[,1], Efig7[,2], Efig7[,3],pch=19,
       ann = FALSE,lwd=9,cex=3,
       col=c(c('blue','green','black','red'),
             rep(c('blue','green','red'),2)),
       errbar.col=c(c('blue','green','black','red'),
                    rep(c('blue','green','red'),2)))
axis(3,cex.axis=6.5)
axis(1,at=1:nrow(Efig7a),labels=FALSE)
mtext("Effect size",side=4,line=3,cex=1.5,family="sans") #Arial font
abline(h=0,lty=2)
graph2jpg(file='Fig. S7 a-h',height=10,width=8.5,dpi=300)







