#This script generates summary statistics for SIA data, as well as Figures 2 & 3

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
library(forcats)

#Load/QC Data
list.files(path="./") #list files in wd

scute<-read_excel("CR_Scute.xlsx") #load scute data
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
scute <- scute %>% 
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

skin.id<-unique(skin$TurtleID) #create vector of turtle IDs for plotting later

##Fix Life.Stage to Subadult<80<Adult
scute$LifeStage <-as.factor(if_else(scute$CCLn<80.0, "Subadult", "Adult"))
skin$LifeStage <-as.factor(if_else(skin$CCLn<80.0, "Subadult", "Adult"))

##Remove problematic points
scute_working<- scute %>%
  filter(ScuteWeird == "No") %>% #removes points identified as having glue, or issues with QC
  filter(TurtleID != "263") %>% #this turtle only has 1 point
  filter(TurtleID != "147-2") %>% #this is a recapture
  filter(TurtleID != "212-2") %>% #this is a resample
  filter(TurtleID != "213-2") %>% #this is a resample
  filter(TurtleID != "137-2") %>%#this is a resample
  drop_na(d13C)

skin_working<- skin %>%
  drop_na(d15NEpi) %>% #removes NAs, primarily from samples analyzed with lipids, plus one sample with bad QC
  drop_na(CCLn) #removes NA size data

##Skin Data
nskin<-length(skin_working$TurtleID)
uskinC<-mean(skin_working$d13CEpi)
sdskinC<-sd(skin_working$d13CEpi)
minskinC<-min(skin_working$d13CEpi)
maxskinC<-max(skin_working$d13CEpi)
uskinN<-mean(skin_working$d15NEpi)
sdskinN<-sd(skin_working$d15NEpi)
minskinN<-min(skin_working$d15NEpi)
maxskinN<-max(skin_working$d15NEpi)
adultskin<-length(skin_working$LifeStage[skin_working$LifeStage=="Adult"])
subadultskin<-length(skin_working$LifeStage[skin_working$LifeStage=="Subadult"])
maleskin<-length(skin_working$Sex[skin_working$Sex=="Male"])

##Skin Female
FEMuskinC<-mean(skin_working$d13CEpi[skin_working$Sex=="Female"])
FEMsdskinC<-sd(skin_working$d13CEpi[skin_working$Sex=="Female"])
FEMminskinC<-min(skin_working$d13CEpi[skin_working$Sex=="Female"])
FEMmaxskinC<-max(skin_working$d13CEpi[skin_working$Sex=="Female"])
FEMuskinN<-mean(skin_working$d15NEpi[skin_working$Sex=="Female"])
FEMsdskinN<-sd(skin_working$d15NEpi[skin_working$Sex=="Female"])
FEMminskinN<-min(skin_working$d15NEpi[skin_working$Sex=="Female"])
FEMmaxskinN<-max(skin_working$d15NEpi[skin_working$Sex=="Female"])

##Skin Male
MALEuskinC<-mean(skin_working$d13CEpi[skin_working$Sex=="Male"])
MALEsdskinC<-sd(skin_working$d13CEpi[skin_working$Sex=="Male"])
MALEminskinC<-min(skin_working$d13CEpi[skin_working$Sex=="Male"])
MALEmaxskinC<-max(skin_working$d13CEpi[skin_working$Sex=="Male"])
MALEuskinN<-mean(skin_working$d15NEpi[skin_working$Sex=="Male"])
MALEsdskinN<-sd(skin_working$d15NEpi[skin_working$Sex=="Male"])
MALEminskinN<-min(skin_working$d15NEpi[skin_working$Sex=="Male"])
MALEmaxskinN<-max(skin_working$d15NEpi[skin_working$Sex=="Male"])

#Skin Subadult
SUBuskinC<-mean(skin_working$d13CEpi[skin_working$Sex=="Subadult"])
SUBsdskinC<-sd(skin_working$d13CEpi[skin_working$Sex=="Subadult"])
SUBminskinC<-min(skin_working$d13CEpi[skin_working$Sex=="Subadult"])
SUBmaxskinC<-max(skin_working$d13CEpi[skin_working$Sex=="Subadult"])
SUBuskinN<-mean(skin_working$d15NEpi[skin_working$Sex=="Subadult"])
SUBsdskinN<-sd(skin_working$d15NEpi[skin_working$Sex=="Subadult"])
SUBminskinN<-min(skin_working$d15NEpi[skin_working$Sex=="Subadult"])
SUBmaxskinN<-max(skin_working$d15NEpi[skin_working$Sex=="Subadult"])

#Scute

nscute<-length(unique(scute_working$TurtleID))
uscuteC<-mean(scute_working$d13C)
sdscuteC<-sd(scute_working$d13C)
minscuteC<-min(scute_working$d13C)
maxscuteC<-max(scute_working$d13C)
uscuteN<-mean(scute_working$d15N)
sdscuteN<-sd(scute_working$d15N)
minscuteN<-min(scute_working$d15N)
maxscuteN<-max(scute_working$d15N)
maxlayers<-max(as.numeric(scute_working$Layer))
layers<- scute_working %>%
  group_by(TurtleID) %>%
  top_n(1, Layer)
avglayer<-mean(as.numeric(layers$Layer))
sdlayer<-sd(as.numeric(layers$Layer))
adultscute<-length(layers$LifeStage[layers$LifeStage=="Adult"])
subadultscute<-length(layers$LifeStage[layers$LifeStage=="Subadult"])
malescute<-length(layers$Sex[layers$Sex=="Male"])

#error bar calcs for skin
sbg <- skin_working %>% 
  group_by(Sex) %>% 
  summarise(count = n(),
            mC = mean(d13CEpi), 
            sdC = sd(d13CEpi), 
            mN = mean(d15NEpi), 
            sdN = sd(d15NEpi))

##Plots
#skin

jpeg(filename="skin_sex_plot.jpeg", width=8, height=8, units="in", res=350)
skin_working %>%
  mutate(Sex = as.character(Sex)) %>%
  mutate(Sex = recode(Sex, 
                          "Female" = "Female",
                          "Male" = "Male",
                          "NA" = "Subadult")) %>%
  mutate(Sex = as.factor(Sex)) %>%
  ggplot(aes(x=d13CEpi, y=d15NEpi, pch=Sex, col=Sex)) +
  geom_point(alpha=1, size = 4) +
  theme_classic()+
  scale_color_viridis_d()+
  theme(plot.title = element_blank(),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        legend.position = c(0.9, 0.9),
        legend.title.align = 0.5,
        legend.title=element_blank(), 
        legend.text=element_text(size=12))+
  xlab(expression({delta}^13*C~'\u2030')) +
  ylab(expression({delta}^15*N~'\u2030')) +
  geom_errorbar(data = sbg, 
                mapping = aes(x = mC, y = mN,
                              ymin = mN - 1.96*sdN, 
                              ymax = mN + 1.96*sdN), 
                width = 0, size=1.3, alpha=0.7) +
  geom_errorbarh(data = sbg, 
                 mapping = aes(x = mC, y = mN,
                               xmin = mC - 1.96*sdC,
                               xmax = mC + 1.96*sdC),
                 height = 0, size=1.3, alpha=0.7) + 
  geom_point(data = sbg, aes(x = mC, 
                             y = mN,
                             fill = Sex), 
             color = "black", shape = 22, size = 7,
             alpha = 1, show.legend = FALSE) +
  scale_fill_viridis_d()
dev.off()

#Scute
scuteNplot <-scute_working %>%
  ggplot(aes(x=Layer, y=d15N, group=TurtleID, pch=Sex, col=Sex)) +
  geom_line(alpha=0.8, size=1) +
  geom_point(alpha=0.8, size = 3) +
  #theme(legend.position="none") +
  scale_color_viridis_d()+
  theme_classic()+
  scale_y_continuous(limits=c(3.25,13), breaks=c(4, 6, 8, 10, 12),
                     labels=c("4","6", "8", "10", "12"))+
  theme(plot.margin = unit(c(1,1,1,1), "pt"),
        plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=20),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=14),
        legend.position = c(0.9, 0.9),
        legend.title.align = 0.5,
        legend.title=element_blank(), 
        legend.text=element_text(size=12))+
  ylab(expression({delta}^15*N~'\u2030')) +
  annotate(geom="text", x = 1, y = 13, label = "A", size=8, hjust=0) 

scuteCplot<-scute_working %>%
  ggplot(aes(x=Layer, y=d13C, group=TurtleID, pch=Sex, col=Sex)) +
  geom_line(alpha=0.8, size=1) +
  geom_point(alpha=0.8, size = 3) +
  scale_color_viridis_d()+
  theme_classic()+
  scale_y_continuous(limits=c(-20,-11.5), breaks=c(-20, -18,-16,-14,-12),
                     labels=c("20", "18", "16", "14", "12"))+
  theme(plot.margin = unit(c(1,1,1,1), "pt"),
        legend.position="none",
        plot.title = element_blank(),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14))+
  ylab(expression({delta}^13*C~'\u2030')) +
  annotate(geom="text", x = 1, y = -12, label = "B", size=8, hjust=0) 

jpeg("scute_sex_plots.jpg", width = 10, height = 10, units = "in", res = 350) 
grid.arrange(scuteNplot, scuteCplot, ncol=1, nrow=2)
dev.off()
