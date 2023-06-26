library(metafor)
library(xlsx)
library(Hmisc)
library(export)
library(reshape2)
setwd('...')
combined=read.csv('Meta combined data.csv',header=T)
F1b_dat=read.csv('Publication year-2023-0501.csv',header = T)

hist(F1b_dat$Publication.year,xlab = 'Publication year',
     xaxt="n",
     cex.lab=1.5, cex.axis=1.5,
     ylab = 'Number of papers',main = '')
axis(side=1,at=seq(1350,2925,25),cex.axis=1.5)
graph2jpg(file='Figure 1b',height=10,width=9,dpi=300)


########## Fig. S29a-d ###########
Pesticide_category=c('Pesticide',
                     unique(combined$pesticide_by_target_organisms))
dat_list1=lapply(unique(combined$Taxonomic_group_response),
                   function(x){subset(combined,Taxonomic_group_response==x)})
dat_fig2=matrix(NA,length(Pesticide_category)*length(dat_list1),2)
for (i in 1:length(Pesticide_category)) {
  for (j in 1:length(dat_list1)){
    mydat=subset(dat_list1[[j]],
                   pesticide_by_target_organisms==Pesticide_category[i])
    res=c(sum(mydat$yi>=0)/length(mydat$yi),sum(mydat$yi<0)/length(mydat$yi))
    dat_fig2[(i-1)*length(dat_list1)+j,]=res
  }
}
colnames(dat_fig2)<-c("ub","lb")

### Supplementary Figure 22
fig2a=melt(as.data.frame(dat_fig2))
ggplot(fig2a,aes(
  y=ifelse(variable=='ub',value,-value),
  fill=variable))+
  geom_bar(stat='identity')+
  coord_flip()+
  theme(axis.text.y = element_blank(),axis.ticks.y=element_blank(),
        axis.text.x = element_text(size = 15),
        legend.text = element_text(size = 15))+
  scale_y_continuous(limits=c(-1,1))+
  scale_fill_discrete(name="Percetage",
                      breaks=c("ub", "lb"),
                      labels=c("positive", "negative"))
graph2jpg(file='Fig. S29a-d',height=8,width=16,dpi=300)



########## Fig. S30a-i ########### 

combined=combined[combined$yi>=-2&combined$yi<=2,]


funnelplot=function(mydat,xlab,name) {
  funnel(mydat$yi,mydat$vi)
  graph2jpg(file=name,height=5,width=6,dpi=300)
}

########## Fig. S30a-i
mydat=subset(combined,
             pesticide_by_target_organisms=='insecticides'&Taxonomic_group=='animals')
funnelplot(mydat,'Insecticides: animal effect size','Fig. S30a-i')


