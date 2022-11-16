##Scute TNW Analyses##
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
  filter(ScuteWeird == "No") %>% #removes points identified as having glue, or issues with QC
  filter(TurtleID != "263") %>% #this turtle only has 1 point
  filter(TurtleID != "147-2") %>%  #this is a recapture
  filter(TurtleID != "212-2") %>% #this is a resample
  filter(TurtleID != "213-2") %>% #this is a resample
  filter(TurtleID != "137-2") %>%#this is a resample
  drop_na(d13C)
  
skin_working<- skin %>%
  filter(TurtleID != "147-2") %>% #this is a recapture
  drop_na(d15NEpi) #removes NAs, primarily from samples analyzed with lipids, plus one sample with bad QC

##TNW Metrics##
##SCUTE FIRST##
##AOV FIRST, LMM SECOND##
#Normality test + visualization. I think normality test not necessary, as we are primarily interested in residuals from model.
shapiro.test(scute_working$d13C)
shapiro.test(scute_working$d15N)
hist(scute_working$d13C)
hist(scute_working$d15N)

#Overall
model_all_C <- aov(d13C ~ Layer + TurtleID, data=scute_working) #same slope but different intercepts for all turtles
model_all_N <- aov(d15N ~ Layer + TurtleID, data=scute_working)

# Use the mean sum of squares to calculate WIC and BIC
WIC_all_C <- anova(model_all_C)["Residuals", "Mean Sq"] #within individual component
BIC_all_C <- anova(model_all_C)["TurtleID", "Mean Sq"] # between individual component
RIS_all_C <- BIC_all_C / WIC_all_C #relative index of specialization
TNW_all_C <- WIC_all_C + BIC_all_C #total niche width
WIC_all_C
BIC_all_C
RIS_all_C
TNW_all_C

WIC_all_N <- anova(model_all_N)["Residuals", "Mean Sq"] #within individual component
BIC_all_N <- anova(model_all_N)["TurtleID", "Mean Sq"] # between individual component
RIS_all_N <- BIC_all_N / WIC_all_N #relative index of specialization
TNW_all_N <- WIC_all_N + BIC_all_N #total niche width
WIC_all_N
BIC_all_N
RIS_all_N
TNW_all_N

#Adults
scute_adult <- scute_working %>%
  filter(LifeStage=="Adult")
  
model_adult_C <- aov(d13C ~ Layer + TurtleID, data=scute_adult)
model_adult_N <- aov(d15N ~ Layer + TurtleID, data=scute_adult)

# Use the mean sum of squares to calculate WIC and BIC
WIC_adult_C <- anova(model_adult_C)["Residuals", "Mean Sq"] #within individual component
BIC_adult_C <- anova(model_adult_C)["TurtleID", "Mean Sq"] # between individual component
RIS_adult_C <- BIC_adult_C / WIC_adult_C #relative index of specialization
TNW_adult_C <- WIC_adult_C + BIC_adult_C #total niche width
WIC_adult_C
BIC_adult_C
RIS_adult_C
TNW_adult_C

WIC_adult_N <- anova(model_adult_N)["Residuals", "Mean Sq"] #within individual component
BIC_adult_N <- anova(model_adult_N)["TurtleID", "Mean Sq"] # between individual component
RIS_adult_N <- BIC_adult_N / WIC_adult_N #relative index of specialization
TNW_adult_N <- WIC_adult_N + BIC_adult_N #total niche width
WIC_adult_N
BIC_adult_N
RIS_adult_N
TNW_adult_N

#Subadults
scute_subadult <- scute_working %>%
  filter(LifeStage=="Subadult")

model_subadult_C <- aov(d13C ~ Layer + TurtleID, data=scute_subadult)
model_subadult_N <- aov(d15N ~ Layer + TurtleID, data=scute_subadult)

# Use the mean sum of squares to calculate WIC and BIC
WIC_subadult_C <- anova(model_subadult_C)["Residuals", "Mean Sq"] #within individual component
BIC_subadult_C <- anova(model_subadult_C)["TurtleID", "Mean Sq"] # between individual component
RIS_subadult_C <- BIC_subadult_C / WIC_subadult_C #relative index of specialization
TNW_subadult_C <- WIC_subadult_C + BIC_subadult_C #total niche width
WIC_subadult_C
BIC_subadult_C
RIS_subadult_C
TNW_subadult_C

WIC_subadult_N <- anova(model_subadult_N)["Residuals", "Mean Sq"] #within individual component
BIC_subadult_N <- anova(model_subadult_N)["TurtleID", "Mean Sq"] # between individual component
RIS_subadult_N <- BIC_subadult_N / WIC_subadult_N #relative index of specialization
TNW_subadult_N <- WIC_subadult_N + BIC_subadult_N #total niche width
WIC_subadult_N
BIC_subadult_N
RIS_subadult_N
TNW_subadult_N

#LMM
scute_working$Layer<-as.numeric(scute_working$Layer)

memodel_all_C <- lmer(d13C ~ Layer + (1|TurtleID), data=scute_working)
memodel_all_N <- lmer(d15N ~ Layer + (1|TurtleID), data=scute_working)

memWIC_all_C <- as.data.frame(VarCorr(memodel_all_C))[2, "sdcor"] #within individual component
memBIC_all_C <- as.data.frame(VarCorr(memodel_all_C))[1, "sdcor"] # between individual component
memRIS_all_C <- memBIC_all_C / memWIC_all_C #relative index of specialization
memTNW_all_C <- memWIC_all_C + memBIC_all_C #total niche width

memWIC_all_N <- as.data.frame(VarCorr(memodel_all_N))[2, "sdcor"] #within individual component
memBIC_all_N <- as.data.frame(VarCorr(memodel_all_N))[1, "sdcor"] # between individual component
memRIS_all_N <- memBIC_all_N / memWIC_all_N #relative index of specialization
memTNW_all_N <- memWIC_all_N + memBIC_all_N #total niche width

scute_adult <- scute_working %>%
  filter(LifeStage=="Adult")

memodel_all_C_adult <- lmer(d13C ~ Layer + (1|TurtleID), data=scute_adult)
memodel_all_N_adult <- lmer(d15N ~ Layer + (1|TurtleID), data=scute_adult)

memWIC_all_C_adult <- as.data.frame(VarCorr(memodel_all_C_adult))[2, "sdcor"] #within individual component
memBIC_all_C_adult <- as.data.frame(VarCorr(memodel_all_C_adult))[1, "sdcor"] # between individual component
memRIS_all_C_adult <- memBIC_all_C_adult / memWIC_all_C_adult #relative index of specialization
memTNW_all_C_adult <- memWIC_all_C_adult + memBIC_all_C_adult #total niche width

memWIC_all_N_adult <- as.data.frame(VarCorr(memodel_all_N_adult))[2, "sdcor"] #within individual component
memBIC_all_N_adult <- as.data.frame(VarCorr(memodel_all_N_adult))[1, "sdcor"] # between individual component
memRIS_all_N_adult <- memBIC_all_N_adult / memWIC_all_N_adult #relative index of specialization
memTNW_all_N_adult <- memWIC_all_N_adult + memBIC_all_N_adult #total niche width

scute_subadult <- scute_working %>%
  filter(LifeStage=="Subadult")

memodel_all_C_sub <- lmer(d13C ~ Layer + (1|TurtleID), data=scute_subadult)
memodel_all_N_sub <- lmer(d15N ~ Layer + (1|TurtleID), data=scute_subadult)

memWIC_all_C_sub <- as.data.frame(VarCorr(memodel_all_C_sub))[2, "sdcor"] #within individual component
memBIC_all_C_sub <- as.data.frame(VarCorr(memodel_all_C_sub))[1, "sdcor"] # between individual component
memRIS_all_C_sub <- memBIC_all_C_sub / memWIC_all_C_sub #relative index of specialization
memTNW_all_C_sub <- memWIC_all_C_sub + memBIC_all_C_sub #total niche width

memWIC_all_N_sub <- as.data.frame(VarCorr(memodel_all_N_sub))[2, "sdcor"] #within individual component
memBIC_all_N_sub <- as.data.frame(VarCorr(memodel_all_N_sub))[1, "sdcor"] # between individual component
memRIS_all_N_sub <- memBIC_all_N_sub / memWIC_all_N_sub #relative index of specialization
memTNW_all_N_sub <- memWIC_all_N_sub + memBIC_all_N_sub #total niche width

#Export Metrics
niche_metrics<- list(A = c(WIC_all_C, WIC_adult_C,WIC_subadult_C, memWIC_all_C, memWIC_all_C_adult, memWIC_all_C_sub),
                     B = c(BIC_all_C, BIC_adult_C, BIC_subadult_C, memBIC_all_C, memBIC_all_C_adult, memBIC_all_C_sub),
                     C = c(RIS_all_C, RIS_adult_C, RIS_subadult_C, memRIS_all_C, memRIS_all_C_adult, memRIS_all_C_sub),
                     D = c(TNW_all_C, TNW_adult_C, TNW_subadult_C, memTNW_all_C, memTNW_all_C_adult, memTNW_all_C_sub),
                     E = c(WIC_all_N, WIC_adult_N, WIC_subadult_N, memWIC_all_N, memWIC_all_N_adult, memWIC_all_N_sub),
                     F = c(BIC_all_N, BIC_adult_N, BIC_subadult_N, memBIC_all_N, memBIC_all_N_adult, memBIC_all_N_sub),
                     G = c(RIS_all_N, RIS_adult_N, RIS_subadult_N, memRIS_all_N, memRIS_all_N_adult, memRIS_all_N_sub),
                     H = c(TNW_all_N, TNW_adult_N, TNW_subadult_N, memTNW_all_N, memTNW_all_N_adult, memTNW_all_N_sub))


niche_metrics<- as.data.frame(niche_metrics, col.names = c("WIC_C", "BIC_C", "RIS_C", "TNW_C", "WIC_N", "BIC_N", "RIS_N", "TNW_N"), row.names = c("All", "Adults", "Subadults", "LMM All", "LMM Adults", "LMM Subadults"))
write.csv(niche_metrics, file = "scute_niche_metrics.csv")

##SKIN##

#Overall
model_all_C <- aov(d13CEpi ~ TurtleID, data=skin_working) #same slope but different intercepts for all turtles
model_all_N <- aov(d15NEpi ~ TurtleID, data=skin_working)

# Use the mean sum of squares to calculate WIC and BIC
BIC_all_C <- anova(model_all_C)["TurtleID", "Mean Sq"] # between individual component
TNW_all_C <- BIC_all_C #total niche width
BIC_all_C
TNW_all_C

BIC_all_N <- anova(model_all_N)["TurtleID", "Mean Sq"] # between individual component
TNW_all_N <- BIC_all_N #total niche width
BIC_all_N
TNW_all_N

#Adults
skin_adult <- skin_working %>%
  filter(LifeStage=="Adult")

model_adult_C <- aov(d13CEpi ~ TurtleID, data=skin_adult)
model_adult_N <- aov(d15NEpi ~ TurtleID, data=skin_adult)

# Use the mean sum of squares to calculate WIC and BIC
BIC_adult_C <- anova(model_adult_C)["TurtleID", "Mean Sq"] # between individual component
TNW_adult_C <- BIC_adult_C #total niche width
BIC_adult_C
TNW_adult_C

BIC_adult_N <- anova(model_adult_N)["TurtleID", "Mean Sq"] # between individual component
TNW_adult_N <- BIC_adult_N #total niche width
BIC_adult_N
TNW_adult_N

#Subadults
skin_subadult <- skin_working %>%
  filter(LifeStage=="Subadult")

model_subadult_C <- aov(d13CEpi ~ TurtleID, data=skin_subadult)
model_subadult_N <- aov(d15NEpi ~ TurtleID, data=skin_subadult)

# Use the mean sum of squares to calculate WIC and BIC
BIC_subadult_C <- anova(model_subadult_C)["TurtleID", "Mean Sq"] # between individual component
TNW_subadult_C <- BIC_subadult_C #total niche width
BIC_subadult_C
TNW_subadult_C

BIC_subadult_N <- anova(model_subadult_N)["TurtleID", "Mean Sq"] # between individual component
TNW_subadult_N <- BIC_subadult_N #total niche width
BIC_subadult_N
TNW_subadult_N

#Export Metrics
niche_metrics<- list(A = c(BIC_all_C, BIC_adult_C, BIC_subadult_C),
                     B = c(TNW_all_C, TNW_adult_C, TNW_subadult_C),
                     C = c(BIC_all_N, BIC_adult_N, BIC_subadult_N),
                     D = c(TNW_all_N, TNW_adult_N, TNW_subadult_N))


niche_metrics<- as.data.frame(niche_metrics, col.names = c("BIC_C", "TNW_C", "BIC_N", "TNW_N"), row.names = c("All", "Adults", "Subadults"))
write.csv(niche_metrics, file = "skin_niche_metrics.csv")

