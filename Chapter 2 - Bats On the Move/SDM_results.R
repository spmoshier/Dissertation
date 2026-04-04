## Build boxplot to show all the different results
setwd("")

winter<-read.csv("winter_stats1.csv")
summer<-read.csv("summer_stats1.csv")

AUCnTSS<-merge(winter, summer, by="Species")
AUCnTSS2<-rbind(winter, summer)

winter1<-winter[,c(4:7,9:14,16:17)]
summer1<-summer[,c(4:7,9:14,16:17)]
#adjust the AUC values so they are all between 0.5 and 1
winter1[,1:6]<-winter1[,1:6]+0.5
summer1[,1:6]<-summer1[,1:6]+0.5

winter1[,1:6]<-winter1[,1:6]
summer1[,1:6]<-summer1[,1:6]

#format data for plotting
library(reshape)
winter2<-melt(as.data.frame(winter1))
winter2$Season<-"Winter"
summer2<-melt(as.data.frame(summer1))
summer2$Season<-"Summer"
alldata<-rbind(winter2,summer2)

write.csv(alldata, "alldata.csv")


## Plot of AUC & TSS results (Figure 3)
library(ggplot2)

alldata<-read.csv("alldata.csv")
levels<-unique(alldata$variable)

plot1<-ggplot(alldata, aes(x=factor(variable, levels=levels), y=value, fill=Season))+
  geom_boxplot(alpha=0.8) +
  geom_vline(xintercept = 6.5, linetype = "solid", color = "black", size = 1) +
  scale_fill_manual(values=c("#34909A","#A4548B")) +
  #ylim(0,1) +
  ggtitle("Evaluation of Species Distribution Models Using AUC and TSS") +
  labs(fill="") +
  xlab("AUC Evaluation                                                                                       TSS Evaluation") +
  ylab("") +
  theme(legend.position="right") +
  theme(axis.title.x = element_text(vjust = -1)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
#theme_classic(base_size = 12) # Increase base font size
plot1


## Variable importance for the different models
winter_weights<-read.csv("winter_stats2.csv")
#winter_weights<-melt(as.data.frame(winter_weights[,-1]))
winter_weights_rf<-winter_weights[,2:14]
winter_weights_glm<-winter_weights[,c(2:3,16:26)]
winter_weights_mx<-winter_weights[,c(2:3,27:37)]

summer2a<-read.csv("Summer_stats2a.csv")
summer2b<-read.csv("Summer_stats2b.csv")
summer2c<-read.csv("Summer_stats2c.csv")
summer2d<-read.csv("Summer_stats2d.csv")
summer2e<-read.csv("summer_stats2e.csv")
summer_weights<-rbind(summer2a, summer2b, summer2c, summer2d, summer2e)
#summer_weights<-melt(as.data.frame(summer_weights[,-1]))
summer_weights_rf<-summer_weights[,2:14]
summer_weights_glm<-summer_weights[,c(2:3,16:26)]
summer_weights_mx<-summer_weights[,c(2:3,27:37)]

all_weights_rf<-rbind(winter_weights_rf, summer_weights_rf)
all_weights_glm<-rbind(winter_weights_glm, summer_weights_glm)
all_weights_mx<-rbind(winter_weights_mx, summer_weights_mx)

all_weights_rf2<-melt(rbind(winter_weights_rf, summer_weights_rf))
all_weights_glm2<-melt(rbind(winter_weights_glm, summer_weights_glm))
all_weights_mx2<-melt(rbind(winter_weights_mx, summer_weights_mx))

## Test plots for the different models
library(ggplot2)

#random forest plot
plot2<-ggplot(all_weights_rf2, aes(x=factor(variable), y=value, fill=Season))+
  geom_boxplot(alpha=0.8) +
  #geom_vline(xintercept = 6.5, linetype = "solid", color = "black", size = 1) +
  scale_fill_manual(values=c("#34909A","#A4548B")) +
  ggtitle("Evaluation of Variable Rank Importance - Random Forest") +
  labs(fill="") +
  ylab("") +
  theme(legend.position="right") +
  theme(axis.title.x = element_text(vjust = -1)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
#theme_classic(base_size = 12) # Increase base font size
plot2

#glm plot
plot3<-ggplot(all_weights_glm2, aes(x=factor(variable), y=exp(value), fill=Season))+ #exp should fix the glm values?
  geom_boxplot(alpha=0.8) +
  #geom_vline(xintercept = 6.5, linetype = "solid", color = "black", size = 1) +
  scale_fill_manual(values=c("#34909A","#A4548B")) +
  ggtitle("Evaluation of Variable Rank Importance - Generalized Linear Model") +
  labs(fill="") +
  ylab("") +
  theme(legend.position="right") +
  theme(axis.title.x = element_text(vjust = -1)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
#theme_classic(base_size = 12) # Increase base font size
plot3

#maxent plot
plot4<-ggplot(all_weights_mx2, aes(x=factor(variable), y=value, fill=Season))+
  geom_boxplot(alpha=0.8) +
  #geom_vline(xintercept = 6.5, linetype = "solid", color = "black", size = 1) +
  scale_fill_manual(values=c("#34909A","#A4548B")) +
  ggtitle("Evaluation of Variable Rank Importance - MaxEnt Model") +
  labs(fill="") +
  ylab("") +
  theme(legend.position="right") +
  theme(axis.title.x = element_text(vjust = -1)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
#theme_classic(base_size = 12) # Increase base font size
plot4


## rank variables by average importance
#random forest
all_rf_rank<-t(as.data.frame(rank(-colMeans(all_weights_rf[,3:13], na.rm = TRUE))))
win_rf_rank<-t(as.data.frame(rank(-colMeans(winter_weights_rf[,3:13], na.rm = TRUE))))
sum_rf_rank<-t(as.data.frame(rank(-colMeans(summer_weights_rf[,3:13], na.rm = TRUE))))

#glm - exponentiate the coefficient values to help interpretation
all_glm_rank<-t(as.data.frame(rank(-colMeans(exp(all_weights_glm[,3:13]), na.rm = TRUE))))
win_glm_rank<-t(as.data.frame(rank(-colMeans(exp(winter_weights_glm[,3:13]), na.rm = TRUE))))
sum_glm_rank<-t(as.data.frame(rank(-colMeans(exp(summer_weights_glm[,3:13]), na.rm = TRUE))))

#maxent
all_mx_rank<-t(as.data.frame(rank(-colMeans(all_weights_mx[,3:13], na.rm = TRUE))))
win_mx_rank<-t(as.data.frame(rank(-colMeans(winter_weights_mx[,3:13], na.rm = TRUE))))
sum_mx_rank<-t(as.data.frame(rank(-colMeans(summer_weights_mx[,3:13], na.rm = TRUE))))


## build a table of variable importance
rf_importance<-as.data.frame(rbind(all_rf_rank, win_rf_rank, sum_rf_rank))
glm_importance<-as.data.frame(rbind(all_glm_rank, win_glm_rank, sum_glm_rank))
mx_importance<-as.data.frame(rbind(all_mx_rank, win_mx_rank, sum_mx_rank))

write.csv(rf_importance, file="rf_importance.csv")
write.csv(glm_importance, file="glm_importance.csv")
write.csv(mx_importance, file="mx_importance.csv")
