#This script constructs linear models and mixed models for all data, conducts basic model tests relative to null models, and creates Figure S1

library(here) #set wd here
library(ggplot2)
library(tidyr)
library(dplyr)
library(data.table)
library(reshape2)
library(lubridate)
library(readxl)
library(grid)
library(gridExtra)
library(ggrepel)
library(lme4)
library(SIBER)
library(viridis)

#Load/QC Data
list.files(path="./") #list files in wd

scute<-read_excel("CR_Scute_lmm.xlsx") #load scute data
summary(scute) #check scute import

scute_factors<-c("TripCode", 
                 "FlipperTag1", 
                 "FlipperTag2", 
                 "PITTag", 
                 "Sex", 
                 "LifeStage", 
                 "TissueSampleUseable", 
                 "ScuteWeird", 
                 "TurtleID", 
                 "Depthum",
                 "Comment",
                 "IRMScomment",
                 "EAcomment") #select appropriate vars to convert into factors for scute data
scute[scute_factors]<-lapply(scute[scute_factors], as.factor) #convert selected vars to factors
scute$d13C<-as.numeric(scute$d13C)
scute$Cmicromol<-as.numeric(scute$Cmicromol)
scute$CNmass<-as.numeric(scute$CNmass)
scute$CNmolar<-as.numeric(scute$CNmolar)
scute<- scute %>% 
  filter(Layer!="1 2")
scute$Date <- as.Date(scute$Date , format = "%y/%m/%d") #convert date from char to date
scute$Month<- as.factor(month(scute$Date)) #extract month from date
scute$Year<-as.factor(year(scute$Date)) #extract year from date
scute$Layer<-as.numeric(scute$Layer)
scute$LayerTime<-scute$Date-days(scute$Layer*219)
scute$Layer<-factor(scute$Layer, levels = c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18')) #order scute layers
summary(scute) #check conversion

scute.id<-unique(scute$TurtleID) #create vector of turtle IDs for plotting later

skin<-read_xlsx("CR_Skin.xlsx") #load skin data
summary(skin) #check skin import

skin_factors<-c("TripCode", 
                "TurtleID",
                "FlipperTag1", 
                "FlipperTag2", 
                "PITTag", 
                "Sex", 
                "LifeStage",
                "ScuteSent",
                "ScuteDone", 
                "ScuteDataLocation", 
                "ScuteLayers",
                "TissueSent",
                "TissueDone",
                "TissueDataLocation",
                "TissueUse",
                "Comments") #select appropriate vars to convert into factors for skin data
skin[skin_factors]<-lapply(skin[skin_factors], as.factor) #convert selected vars to factors

colnames(skin)[which(names(skin) == 'Cnmass')] <- 'CNmass'
colnames(skin)[which(names(skin) == 'Cnmolar')] <- 'CNmolar'

skin_nums<-c("HeadWidth",
             "CCLn",
             "d15NEpi",
             "d13CEpi",
             "Nmg",
             "Cmg",
             "Npercent",
             "Cpercent",
             "Nmicromol",
             "Cmicromol",
             "CNmass",
             "CNmolar")
skin[skin_nums]<-lapply(skin[skin_nums], as.numeric)
skin$Date <- as.Date(skin$Date , format = "%m/%d/%y") #convert date from char to date
skin$Month<- as.factor(month(skin$Date)) #extract month from date
skin$Year<-as.factor(year(skin$Date)) #extract year from date
summary(skin) #check conversion

##Fix Life.Stage to Subadult<80<Adult
scute$LifeStage <-as.factor(if_else(scute$CCLn<80.0, "Subadult", "Adult"))
skin$LifeStage <-as.factor(if_else(skin$CCLn<80.0, "Subadult", "Adult"))

##Remove problematic points
scute_working<- scute %>%
  drop_na(CCLn) %>%
  filter(ScuteWeird == "No") %>% #removes points identified as having glue, or issues with QC
  filter(TurtleID != "263") %>% #this turtle only has 1 point
  filter(TurtleID != "147-2") %>% #this is a recapture
  filter(TurtleID != "212-2") %>% #this is a resample
  filter(TurtleID != "213-2") %>% #this is a resample
  filter(TurtleID != "137-2") %>%#this is a resample
  drop_na(d13C)#removes NA size data

skin_working<- skin %>%
  drop_na(d15NEpi) %>% #removes NAs, primarily from samples analyzed with lipids, plus one sample with bad QC
  drop_na(CCLn) #removes NA size data

##SCUTE LMM

lmmnullmodelC <- lmer(d13C ~ 1 + (1|TurtleID), data=scute_working)
lmmnullmodelN <- lmer(d15N ~ 1 + (1|TurtleID), data=scute_working)

lmmmodelC <- lmer(d13C ~ CCLn + (1|TurtleID), data=scute_working)
lmmmodelN <- lmer(d15N ~ CCLn + (1|TurtleID), data=scute_working)

plot(lmmmodelC)
hist(resid(lmmmodelC))
plot(order(resid(lmmmodelC)))

plot(lmmmodelN)
hist(resid(lmmmodelN))
plot(order(resid(lmmmodelN)))

anova(lmmmodelC, lmmnullmodelC, test="chisq")
anova(lmmmodelN, lmmnullmodelN, test="chisq")

summary(lmmmodelC)
summary(lmmmodelN)

#LMM Plots


lmmc<-scute_working %>%
  ggplot(aes(y=d13C, x=CCLn), group=TurtleID) +
  geom_point(alpha=0.6, size=2) +
  theme_classic () +
  scale_y_continuous(limits=c(-20,-10), breaks=c(-20, -18,-16,-14,-12, -10),
                     labels=c("20", "18", "16", "14", "12", "10"))+
  ylab(expression({delta}^13*C~'\u2030')) +
  theme(legend.position="none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank())+
  theme(plot.margin = unit(c(0,45,0,0), "pt")) +
  #scale_color_viridis_d(option="D") +
  geom_smooth(method="lm", color="black", size=1.5)+
  annotate(geom="text", x = 55, y = -10, label = "B", size=10, hjust=0) +
  annotate(geom="text", x = 55, y = -10.5, label = {delta}^13*C~'\u2030'~"=0.03 * CCLmin - 18.91", size=6, hjust=0)+
  annotate(geom="text", x = 55, y = -11, label = "p=0.02", size=6, hjust=0)

lmmn<-scute_working %>%
  ggplot(aes(y=d15N, x=CCLn), group=TurtleID)+
  geom_point(size=2, alpha=0.6)+
  theme_classic()+
  ylab(expression({delta}^15*N~'\u2030')) +
  scale_y_continuous(limits=c(3,13), breaks=c(4, 6, 8, 10, 12),
                     labels=c("4","6", "8", "10", "12"))+
  xlab ("Curved Carapace Length (cm)") +
  theme(legend.position="none",
        axis.title.x = element_text(size=20, hjust=-1),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=14),
        axis.text.y = element_blank()) +
  theme(plot.margin = unit(c(0,45,0,0), "pt")) +
  #scale_color_viridis_d(option="D") +
  geom_smooth(method="lm", color="black", size=1.5)+
  annotate(geom="text", x = 55, y = 13, label = "D", size=10, hjust=0) +
  annotate(geom="text", x = 55, y = 12.5, label = {delta}^13*N~'\u2030'~"=-0.04 * CCLmin + 9.79", size=6, hjust=0)+
  annotate(geom="text", x = 55, y = 12, label = "p=0.04", size=6, hjust=0)

#SKIN LM

nullmodelC <- lm(d13CEpi ~ 1, data=skin_working)
nullmodelN <- lm(d15NEpi ~ 1, data=skin_working)

modelC <- lm(d13CEpi ~ CCLn, data=skin_working)
modelN <- lm(d15NEpi ~ CCLn, data=skin_working)

hist(skin_working$d13CEpi)
hist(skin_working$d15NEpi)

plot(modelC)
hist(resid(modelC))
plot(modelN)
hist(resid(modelN))


anova(modelC, nullmodelC, test="Chisq")
anova(modelN, nullmodelN, test="Chisq")

summary(modelC)
summary(modelN)

#LM Plots

lmc<-skin_working%>%
  ggplot(aes(x=CCLn, y=d13CEpi))+
  geom_point(size=4, alpha=0.6)+
  theme_classic()+
  ylab(expression({delta}^13*C~'\u2030')) +
  xlab ("") +
  scale_y_continuous(limits=c(-20,-10), breaks=c(-20, -18,-16,-14,-12, -10),
                   labels=c("20", "18", "16", "14", "12", "10"))+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=20),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=14)) +
  theme(plot.margin = unit(c(0,0,0,0), "pt")) +
  geom_smooth(method="lm", color="black", size=1.5)+
  annotate(geom="text", x = 55, y = -10, label = "A", size=10, hjust=0) +
  annotate(geom="text", x = 55, y = -10.5, label = {delta}^13*C~'\u2030'~"=0.05 * CCLmin - 18.96", size=6, hjust=0)+
  annotate(geom="text", x = 55, y = -11, label = "p<0.001", size=6, hjust=0)

lmn<-skin_working %>%
  ggplot(aes(x=CCLn, y=d15NEpi))+
  geom_point(size=4, alpha=0.6)+
  theme_classic()+
  ylab(expression({delta}^15*N~'\u2030')) +
  scale_y_continuous(limits=c(3,13), breaks=c(4, 6, 8, 10, 12),
                     labels=c("4","6", "8", "10", "12"))+
  xlab ("") +
  theme(axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14)) +
  theme(plot.margin = unit(c(0,0,6,0), "pt")) +
  geom_smooth(method="lm", color="black", size=1.5)+
  annotate(geom="text", x = 55, y = 13, label = "C", size=10, hjust=0) +
  annotate(geom="text", x = 55, y = 12.5, label = {delta}^13*N~'\u2030'~"=-0.03 * CCLmin + 10.04", size=6, hjust=0)+
  annotate(geom="text", x = 55, y = 12, label = "p=0.017", size=6, hjust=0)

jpeg("lm_plots.jpg", width = 12, height = 14, units = "in", res = 350) 
grid.arrange(lmc, lmmc, lmn,  lmmn, ncol=2, nrow=2)
dev.off()
