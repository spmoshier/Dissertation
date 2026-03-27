#### Plotting code for figures & tables
# load in initital datasets
library(dplyr)

setwd("./Final Files")

MO<-read.csv("migration_overlap_2-9-26.csv")
MO_filter<-read.csv("migration_overlap_filtered_2-9-26.csv")
MO_fig1<-read.csv("migration_overlap_figuredata_2-8-26.csv")
MO_fig2<-read.csv("migfig_2-8-26.csv")
NO<-read.csv("final_nicheoverlap.csv")

#create new dataframe using updated niche overlap & centroid calculations
#MO_new<-merge(MO_fig1[,1:9],NO, by="Latin.Name",all=T) #full join of the two datasets


### Figure 2: bat species sorted by literature categorization
library(ggplot2)
library(ggpattern)
library(ggbreak)

MO_dat<-MO_fig1$Fig_stat

ggplot(MO_fig1, aes(x = Colors, fill = Fig_stat))+
  geom_bar(position=position_stack(reverse=T))+
  geom_bar_pattern(position=position_stack(reverse=T),
                   aes(Colors, pattern_fill=Fig_stat),
                   pattern=c('none','crosshatch','none','crosshatch','none'),
                   pattern_density=0.02, pattern_alpha=0.3,
                   pattern_colour=c('#2E818A','#2E818A','#8A5672','#8A5672','#648B41'),
                   fill=c('#85CED6','#85CED6','#C09BAF','#C09BAF','#B6D19E'), #97BE74
                   colour=c('#2E818A','#2E818A','#8A5672','#8A5672','#648B41')) +
  geom_text(stat="count", position=position_stack(reverse=T), aes(label=Fig_stat), vjust=1.5, size=3)+ #color="white"
  labs(x = "Migratory Status", y = "Number of Species", title = "Migratory Classifications From the Literature") +
  theme(legend.position="none", legend.title=element_blank(), axis.text.x=element_blank())+
  scale_y_continuous(limits = c(0, 1275))+
  scale_y_cut(breaks=c(205), which=c(1,2), scales=c(0.5,2), expand=F)


### Figure 3: summary of AUC & TSS statistics for the modeling
library(ggplot2)

alldata<-read.csv("alldata.csv")
levels<-unique(alldata$variable)

plot<-ggplot(alldata, aes(x=factor(variable, levels=levels), y=value, fill=Season))+
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
plot


### Figure 4: boxplots of estimated migratory distances sorted by category (the bottom maps are from ArcGIS SDMToolbox)
library(ggplot2)
library(ggstatsplot)
library(NatParksPalettes)
library(dplyr)
library(gghighlight)

#sort MO_filter dataset
classes<-c("Long-distance Migrant","Short-distance Migrant","Nonmigratory")
MO_sort<-MO_filter[order(MO_filter$dist_km, decreasing=T),]
MO_sort$Index<-c(1:length(MO_sort$FE_Status))

#basic stats to add to boxplot
nLD<-sum(MO_sort$FE_Status=="Long-distance Migrant")
nSD<-sum(MO_sort$FE_Status=="Short-distance Migrant")
nNM<-sum(MO_sort$FE_Status=="Nonmigratory")
meanLD<-mean(MO_sort[MO_sort$FE_Status=="Long-distance Migrant", "dist_km"])
meanSD<-mean(MO_sort[MO_sort$FE_Status=="Short-distance Migrant", "dist_km"])
meanNM<-mean(MO_sort[MO_sort$FE_Status=="Nonmigratory", "dist_km"])

#colored rectangles in the background
colors=natparks.pals("Glacier",3)

#updated code for layering boxplot over top of the point distribution
ggplot(MO_sort, aes(x=Index, y=dist_km, fill=FE_Status))+
  annotate('rect', xmin=-5,xmax=46,ymin=0,ymax=3200, fill=colors[1], alpha=0.2)+
  annotate('rect', xmin=46,xmax=345,ymin=0,ymax=3200, fill=colors[2], alpha=0.2)+
  annotate('rect', xmin=345,xmax=397,ymin=0,ymax=3200, fill=colors[3], alpha=0.2)+
  geom_point(alpha=0.30, colour="gray20", size=2.5)+
  geom_point(data=MO_sort[c(1,61,380),],fill="black", size=4, pch=8, alpha=1)+
  geom_boxplot(alpha=0.30)+
  scale_fill_manual(values=natparks.pals("Glacier",3))+
  #scale_fill_manual(values=c("gray20","gray20","gray20"))+
  #stat_summary(fun.data=give.n, geom="text",vjust=-0.5)+
  scale_y_continuous(limits = c(0, 3200), breaks = c(0, 100, 500, 1000, 2000, 3000)) +
  labs(x="", y="")+ #Number of Species & Estimated Migration Distance (km)
  #xlim(0,400) + ylim(0,3200)+
  #ggplot2::scale_color_manual(values=natparks.pals("Glacier",3)) +
  #theme customizations
  theme(text=element_text(family="", size=0, color="black"),
        plot.title=element_text(family="",size=14,face="bold",color="black"),
        #stat annotations below main title    
        plot.subtitle=element_text(family="",size=0,face="bold",color="black"),
        plot.title.position="plot",
        axis.text=element_text(size=14, color="black"),
        axis.text.x=element_blank(),
        axis.title=element_text(size=14),
        axis.ticks = element_blank(),
        axis.line = element_line(colour = "grey50"),
        panel.grid = element_line(color = "#b4aea9"), #make color white for no lines
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(linetype = "dashed"),
        panel.background = element_rect(fill = "white", color = "white"),
        plot.background = element_rect(fill = "white", color = "white"),
        legend.position = "none") 
#geom_hline(yintercept = 0, linetype="solid", color = "darkblue") +
#geom_hline(yintercept = 100, linetype="solid", color = "darkblue") +
#geom_hline(yintercept = 1000, linetype="solid", color = "darkblue")

#figure out where species fit
ggplot(MO_sort, aes(x=Index, y=dist_km))+
  geom_point(alpha=0.30, colour="gray20", size=2.5)+
  geom_point(data=MO_sort[c(1,61,380),],fill="darkblue", size=6, pch=23)


### Figure 5: violin plots of our categories and Fleming2003 categories
library(ggplot2)
library(ggstatsplot)
library(NatParksPalettes)
library(dplyr)
library(gghighlight)

classes<-c("Long-distance Migrant","Short-distance Migrant","Nonmigratory")
MO_fig2$FE_status<-factor(MO_fig2$FE_status, levels=classes)
MO_filter$FE_Status<-factor(MO_filter$FE_Status, levels=classes)

#plot 1: D by literature categories
plot1<-ggbetweenstats(data=MO_fig2, x=Colors, y=D, p.adjust.method="bonferroni",
                      centrality.point.args = list(size=5, color="black"),
                      centrality.label.args = list(size = 5, nudge_x = 0.3, segment.linetype = 4, min.segment.length = 0),
                      point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.6), alpha = 0.6, size = 4, stroke = 0, na.rm = TRUE),
                      ggsignif.args = list(textsize = 5, tip_length = 0.01, na.rm = TRUE)) +
  labs(x="", y="Schoener's D",title="a) Migration Categories Based on Literature", axis.range.restrict=F) +
  #geom_text(data=MO, label=MO$Latin.Name, check_overlap=T) +
  scale_y_continuous(limits = c(0, 1.0), breaks = c(0, .2, .4, .6, .8, 1.0)) +
  #color palette
  ggplot2::scale_color_manual(values=natparks.pals("Glacier",10)) +
  #theme customizations
  theme(text=element_text(family="", size=16, color="black"),
        plot.title=element_text(family="",size=18,face="bold",color="black"),
        #stat annotations below main title    
        plot.subtitle=element_text(family="",size=16,face="bold",color="black"),
        plot.title.position="plot",
        axis.text=element_text(size=16, color="black"),
        axis.title=element_text(size=16),
        axis.ticks = element_blank(),
        axis.line = element_line(colour = "grey50"),
        panel.grid = element_line(color = "#b4aea9"), #make color white for no lines
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(linetype = "dashed"),
        panel.background = element_rect(fill = "white", color = "white"),
        plot.background = element_rect(fill = "white", color = "white"))

#plot 2: D by F&E categories
plot2<-ggbetweenstats(data=MO_filter, x=FE_Status, y=D, p.adjust.method="bonferroni",
                      centrality.point.args = list(size=5, color="black"),
                      centrality.label.args = list(size = 5, nudge_x = 0.3, segment.linetype = 4, min.segment.length = 0),
                      point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.6), alpha = 0.6, size = 4, stroke = 0, na.rm = TRUE),
                      ggsignif.args = list(textsize = 5, tip_length = 0.01, na.rm = TRUE)) +
  labs(x="", y="Schoener's D",title="b) Migration Categories Based on Fleming et al. (2003)", axis.range.restrict=F) +
  #geom_text(data=MO, label=MO$Latin.Name, check_overlap=T) +
  scale_y_continuous(limits = c(0, 1.2), breaks = c(0, .2, .4, .6, .8, 1.0)) +
  #color palette
  ggplot2::scale_color_manual(values=natparks.pals("Glacier",3)) +
  #theme customizations
  theme(text=element_text(family="", size=16, color="black"),
        plot.title=element_text(family="",size=18,face="bold",color="black"),
        #stat annotations below main title    
        plot.subtitle=element_text(family="",size=16,face="bold",color="black"),
        plot.title.position="plot",
        axis.text=element_text(size=16, color="black"),
        axis.title=element_text(size=16),
        axis.ticks = element_blank(),
        axis.line = element_line(colour = "grey50"),
        panel.grid = element_line(color = "#b4aea9"), #make color white for no lines
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(linetype = "dashed"),
        panel.background = element_rect(fill = "white", color = "white"),
        plot.background = element_rect(fill = "white", color = "white"))

grid.arrange(plot1, plot2, ncol=1)
```

### Figure 6: comparison of our species classification and the literature classification

library(forcats)
classes<-c("Long-distance Migrant","Short-distance Migrant","Nonmigratory")

ggplot(MO_fig2, aes(x=Colors, fill=factor(FE_status, levels=classes)))+
  geom_bar(stat="count", position=position_dodge(), width=0.75)+
  scale_fill_manual(values=c("#2E818A","#85CED6","#8A5672"))+
  labs(x = "Migratory Status From The Literature", y = "Number of Species", title = "")+
  theme(legend.position=c(0.18, 0.8),legend.background=element_blank())+
  guides(fill = guide_legend(title = "Our Classification"))
```

### all estimated mig distances plotted lowest to highest
NO_sort<-NO[order(NO$dist_km, decreasing=T),]
NO_sort$Index<-c(1:length(NO$Latin.Name))

geom_bar(stat="count", position=position_dodge(), width=0.75)+
  scale_fill_manual(values=c("#2E818A","#85CED6","#8A5672"))+
  labs(x = "Migratory Status From The Literature", y = "Number of Species", title = "")+
  theme(legend.position=c(0.18, 0.8),legend.background=element_blank())+
  guides(fill = guide_legend(title = "Our Classification"))


### Supplemental Figure 1: plot estimated migratory distances against some range variables
# midpoint latitude
# range size
# total lat extent
# total lon extent
library(ggplot2)
library(gridExtra)
library(ggpubr)

bat_stats<-read.csv("bat_range_stats.csv")
#merge with niche overlap dataset
bat_range_stats<-merge(bat_stats,MO_filter, by="Latin_Name",all=T) #full join of the two datasets

#set up linear regression analyses
lm_stats<-data.frame()
lm_stats[4,1:2]=NA
colnames(lm_stats)=c("intercept", "slope")

lm1<-lm(bat_range_stats$dist_km~bat_range_stats$total_area_m2)
lm2<-lm(bat_range_stats$dist_km~bat_range_stats$lat_extent_deg)
lm3<-lm(bat_range_stats$dist_km~bat_range_stats$mid_lat)
lm4<-lm(bat_range_stats$dist_km~bat_range_stats$mid_lon)
coeff1<-coefficients(lm1)
coeff2<-coefficients(lm2)
coeff3<-coefficients(lm3)
coeff4<-coefficients(lm4)
lm_stats[1,1:2]<-coeff1
lm_stats[2,1:2]<-coeff2
lm_stats[3,1:2]<-coeff3
lm_stats[4,1:2]<-coeff4

#set up linear regression plots
#plot0<-ggplot(bat_range_stats, aes(x=dist_km, y=total_area_m2))+
geom_point(alpha=0.30, colour="gray20", size=2.5)+
  labs(x = "Estimated Migration Distance (km)", y = expression("Total Area (m"^"2"*")"), title = "a) Total Range Area")+
  geom_smooth(method="lm", formula=y~x, color="#648B41")+
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`*`~")), label.y=20000000000000, label.x=2000) + # Add R^2 and p-value
  stat_regline_equation(label.y = 18000000000000, label.x=2000) # Add the equation, adjust label.y position as needed
#geom_abline(intercept=lm_stats[1,1],slope=lm_stats[1,2], color="red")

plot1<-ggplot(bat_range_stats, aes(x=dist_km, y=log10(total_area_m2)))+
  geom_point(alpha=0.30, colour="gray20", size=2.5)+
  labs(x = "Estimated Migration Distance (km)", y = expression("log"[10]*"(Total Area (m"^"2"*"))"), title = "a) Total Range Area")+
  geom_smooth(method="lm", formula=y~x, color="#648B41")+
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`*`~")), label.y=10.5, label.x=2000) + # Add R^2 and p-value
  stat_regline_equation(label.y = 10, label.x=2000) # Add the equation, adjust label.y position as needed
#geom_abline(intercept=lm_stats[1,1],slope=lm_stats[1,2], color="red")

plot2<-ggplot(bat_range_stats, aes(x=dist_km, y=lat_extent_deg))+
  geom_point(alpha=0.30, colour="gray20", size=2.5)+
  labs(x = "Estimated Migration Distance (km)", y = "Latitudinal Extent (degrees)", title = "b) Latitudinal Extent")+
  #geom_abline(intercept=lm_stats[2,1],slope=lm_stats[2,2], color="red")+
  geom_smooth(method="lm", formula=y~x, color="#648B41")+
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`*`~")), label.y=100, label.x=2000) + # Add R^2 and p-value
  stat_regline_equation(label.y = 90, label.x=2000) # Add the equation, adjust label.y position as needed

plot3<-ggplot(bat_range_stats, aes(x=dist_km, y=mid_lat))+
  geom_point(alpha=0.30, colour="gray20", size=2.5)+
  labs(x = "Estimated Migration Distance (km)", y = "Latitudinal Midpoint (degrees)", title = "c) Latitudinal Midpoint")+
  #geom_abline(intercept=lm_stats[3,1],slope=lm_stats[3,2], color="red")+
  geom_smooth(method="lm", formula=y~x, color="#648B41")+
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`*`~")), label.y=-20, label.x=2300) + # Add R^2 and p-value
  stat_regline_equation(label.y = -30, label.x=2300) # Add the equation, adjust label.y position as needed

plot4<-ggplot(bat_range_stats, aes(x=dist_km, y=mid_lon))+
  geom_point(alpha=0.30, colour="gray20", size=2.5)+
  labs(x = "Estimated Migration Distance (km)", y = "Longitudinal Midpoint (degrees)", title = "d) Longitudinal Midpoint")+
  #geom_abline(intercept=lm_stats[4,1],slope=lm_stats[4,2], color="red")+
  geom_smooth(method="lm", formula=y~x, color="#648B41")+
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`*`~")), label.y=130, label.x=1500) + # Add R^2 and p-value
  stat_regline_equation(label.y = 105, label.x=1500) # Add the equation, adjust label.y position as needed

grid.arrange(plot1, plot2, plot3, plot4, ncol=2)


#set up linear regression plots - try a different way
plot1<-ggplot(bat_range_stats, aes(x=total_area_m2, y=dist_km))+
  geom_point(alpha=0.30, colour="gray20", size=2.5)+
  labs(x = "Estimated Migration Distance (km)", y = expression("Total Area (m"^"2"*")"), title = "")

plot2<-ggplot(bat_range_stats, aes(x=lat_extent_deg, y=dist_km))+
  geom_point(alpha=0.30, colour="gray20", size=2.5)+
  labs(x = "Estimated Migration Distance (km)", y = "Latitudinal Extent (degrees)", title = "")

plot3<-ggplot(bat_range_stats, aes(x=mid_lat, y=dist_km))+
  geom_point(alpha=0.30, colour="gray20", size=2.5)+
  labs(x = "Estimated Migration Distance (km)", y = "Latitudinal Midpoint (degrees)", title = "")

plot4<-ggplot(bat_range_stats, aes(x=mid_lon, y=dist_km))+
  geom_point(alpha=0.30, colour="gray20", size=2.5)+
  labs(x = "Estimated Migration Distance (km)", y = "Longitudinal Midpoint (degrees)", title = "")

grid.arrange(plot1, plot2, plot3, plot4, ncol=2)
