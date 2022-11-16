#This script transforms data into SIBER objects to examine niches as ellipses/convex hulls, and creates Figure 4

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
library(rjags)

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
  filter(TurtleID != "147-2") %>% #this is a recapture
  filter(TurtleID != "137-2") %>% #this is a resample and has only 3 layers
  drop_na(d13C)
  #filter(TurtleID != "181") %>% #sample has 4 layers
  #filter(TurtleID != "251") #sample has 5 layers but one NA

skin_working<- skin %>%
  filter(TurtleID != "147-2") %>%
  drop_na(d15NEpi) #removes NAs, primarily from samples analyzed with lipids, plus one sample with bad QC

##SIBER
#data("demo.siber.data") #check out formatting
#head(demo.siber.data) #check out formatting
#summary(demo.siber.data)

#subset and format data
scute_siber_df <- scute_working %>% 
  select(d13C, d15N, TurtleID, LifeStage) %>%
  mutate(LifeStage = recode(LifeStage, "Subadult" = 1, "Adult" = 2)) %>%
  mutate(TurtleID = as.numeric(TurtleID))
scute_working$SIBERID <- as.numeric(scute_working$TurtleID) #to link "group" metrics to individual turtles
names(scute_siber_df)<-c("iso1","iso2","group", "community")
summary(scute_siber_df)
scute_siber_df <- as.data.frame(scute_siber_df)
scute_siber_sub <- scute_siber_df %>%
  filter(community==1)
scute_siber_adult <- scute_siber_df %>%
  filter(community==2)

#make siber object
scute_siber_obj <- createSiberObject(scute_siber_df)
scute_siber_sub_obj <- createSiberObject(scute_siber_sub)
scute_siber_adult_obj <- createSiberObject(scute_siber_adult)

# Create lists of plotting arguments to be passed onwards to each 
# of the three plotting functions.
community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.95, lty = 1, lwd = 2)
group.hulls.args     <- list(lty = 2, col = "grey20")

#plot data, can change T/F to add/subtract ellipses/hulls
par(mfrow=c(1,1))
plotSiberObject(scute_siber_obj,
                ax.pad = 2, 
                hulls = T, community.hulls.args = community.hulls.args, 
                ellipses = F, group.ellipses.args = group.ellipses.args,
                group.hulls = F, group.hulls.args = group.hulls.args,
                bty = "L",
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'\u2030'),
                ylab = expression({delta}^15*N~'\u2030')
)

#calculate summary stats for each group (turtle)
group.ML <- groupMetricsML(scute_siber_obj)
print(group.ML)

#calculate layman metrics for communities - not bayesian
community.ML <- communityMetricsML(scute_siber_obj) 
print(community.ML)

#fit bayesian model to data
# options for running jags
parms <- list()
parms$n.iter <- 2 * 10^4   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains

# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

# fit the ellipses which uses an Inverse Wishart prior
# on the covariance matrix Sigma, and a vague normal prior on the 
# means. Fitting is via the JAGS method.
ellipses.posterior <- siberMVN(scute_siber_obj, parms, priors)

# The posterior estimates of the ellipses for each group can be used to
# calculate the SEA.B for each group.
SEA.B <- siberEllipses(ellipses.posterior)

#siberDensityPlot(SEA.B, xticklabels = colnames(group.ML), 
#                 xlab = c("Community | Group"),
#                ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
#                 bty = "L",
#                 las = 1,
#                 main = "SIBER ellipses on each group"
#)

# Add red x's for the ML estimated SEA-c
#points(1:ncol(SEA.B), group.ML[3,], col="red", pch = "x", lwd = 2)

# Calculate some credible intervals 
cr.p <- c(0.95, 0.99) # vector of quantiles

# call to hdrcde:hdr using lapply()
SEA.B.credibles <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = cr.p)

# do similar to get the modes, taking care to pick up multimodal posterior
# distributions if present
SEA.B.modes <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$mode},
  prob = cr.p, all.modes=T)

#extract SEA for individual turtles
SEA_group_modes <- as.data.frame(SEA.B.modes, row.names = "SEA", col.names = colnames(group.ML))
SEA_group_modes <- as.data.frame(t(SEA_group_modes))
IDS <- substr(row.names(SEA_group_modes), 4,5)
SEA_group_modes$SIBERID <- as.double(IDS)
SEA_group_modes <- SEA_group_modes[,c("SIBERID", "SEA")]
scute_working<-left_join(scute_working, SEA_group_modes)

SEA_size<- scute_working %>%
  select(CCLn, SEA)
SEA_size<-as.data.frame(unique(SEA_size))

SEA_size %>%
  ggplot(aes(x=CCLn, y=SEA)) +
  geom_point() +
  ggtitle("Mean SEA vs CCLn")+
  theme(plot.title = element_text(hjust = 0.5))

SEA_size %>%
  filter(SEA < 2) %>%
  ggplot(aes(x=CCLn, y=SEA)) +
  geom_point() +
  ggtitle("Mean SEA vs CCLn")+
  theme(plot.title = element_text(hjust = 0.5))

# extract the posterior means
mu.post <- extractPosteriorMeans(scute_siber_obj, ellipses.posterior)

# calculate the corresponding distribution of layman metrics
layman.B <- bayesianLayman(mu.post)


#plot hulls

with (scute_working, 
      plot(d13C, d15N, 
           col=alpha(c("#B8DE29", "#404788")[as.numeric(LifeStage)], 0.6),
           las=1,
           pch=16,
           cex=1.5,
           cex.axis=1.5,
           cex.lab=1.5,
           xlab=expression({delta}^13*C~'\u2030'),
           ylab=""))
title(ylab = expression({delta}^15*N~'\u2030'),
      line=2, 
      cex.lab=1.5)
legend(-13, 13, legend=c("Adult", "Subadult"), col=c("#B8DE29", "#404788"), pch=16, cex=1.2, bty="n")
plotCommunityHulls(scute_siber_sub_obj, plot.args=list(col="#404788", lty=1, lwd=3), iso.order=c(1,2))
plotCommunityHulls(scute_siber_adult_obj, plot.args=list(col="#B8DE29", lty=1, lwd=5), iso.order=c(1,2))

# --------------------------------------
# Visualise the first community: Subadults
# --------------------------------------
siberDensityPlot(layman.B[[1]], xticklabels = colnames(layman.B[[1]]), 
                 bty="L", ylim = c(0,30), xlab = "Subadults")

# add the ML estimates (if you want). Extract the correct means 
# from the appropriate array held within the overall array of means.
comm1.layman.ml <- laymanMetrics(scute_siber_obj$ML.mu[[1]][1,1,],
                                 scute_siber_obj$ML.mu[[1]][1,2,]
)
points(1:6, comm1.layman.ml$metrics, col = "red", pch = "x", lwd = 2)

# --------------------------------------
# Visualise the second community: Adults
# --------------------------------------
siberDensityPlot(layman.B[[2]], xticklabels = colnames(layman.B[[2]]), 
                 bty="L", ylim = c(0,30), xlab = "Adults")

# add the ML estimates. (if you want) Extract the correct means 
# from the appropriate array held within the overall array of means.
comm2.layman.ml <- laymanMetrics(scute_siber_obj$ML.mu[[2]][1,1,],
                                 scute_siber_obj$ML.mu[[2]][1,2,]
)
points(1:6, comm2.layman.ml$metrics, col = "red", pch = "x", lwd = 2)

# --------------------------------------
# Alternatively, pull out Layman metrics from both and aggregate them into a 
# single matrix using cbind() and plot them together on one graph.
# --------------------------------------

# go back to a 1x1 panel plot
par(mfrow=c(1,1))

siberDensityPlot(cbind(layman.B[[1]][,"TA"], layman.B[[2]][,"TA"]), ### bind both TAs from community 1 and 2
                 xticklabels = c("Subadults", "Adults"), 
                 bty="L", ylim = c(0,30),
                 las = 1,
                 ylab = "TA - Convex Hull Area",
                 xlab = "")

siberDensityPlot(cbind(layman.B[[1]][ , "NND"], ### bind both NNDs from community 1 and 2
                       layman.B[[2]][ , "NND"]),
                 xticklabels = c("Subadults", "Adults"), 
                 bty="L", ylim = c(0, 1.5),
                 las = 1,
                 ylab = "NND - nearest neighbor distance",
                 xlab = "")

siberDensityPlot(cbind(layman.B[[1]][ , "SDNND"], ### bind both SDNNDs from community 1 and 2
                       layman.B[[2]][ , "SDNND"]),
                 xticklabels = c("Subadults", "Adults"), 
                 bty="L", ylim = c(0, 1),
                 las = 1,
                 ylab = "SDNND - stdev of nearest neighbor distance",
                 xlab = "")

siberDensityPlot(cbind(layman.B[[1]][ , "CD"], ###bind both CDs fromcommunity 1 and 2
                       layman.B[[2]][ , "CD"]),
                 xticklabels = c("Subadults", "Adults"), 
                 bty="L", ylim = c(1, 3),
                 las = 1,
                 ylab = "CD - centroid dispersion",
                 xlab = "")

siberDensityPlot(cbind(layman.B[[1]][ , "dY_range"], ###bind both dY_ranges fromcommunity 1 and 2
                       layman.B[[2]][ , "dY_range"]),
                 xticklabels = c("Subadults", "Adults"), 
                 bty="L", ylim = c(0, 10),
                 las = 1,
                 ylab = "dY_range - range of mean d15N of groups within communities",
                 xlab = "")

siberDensityPlot(cbind(layman.B[[1]][ , "dX_range"], ###bind both dX_ranges fromcommunity 1 and 2
                       layman.B[[2]][ , "dX_range"]),
                 xticklabels = c("Subadults", "Adults"), 
                 bty="L", ylim = c(0, 10),
                 las = 1,
                 ylab = "dX_range - range of mean d13C of groups within communities",
                 xlab = "")

#Extract Means and CIs from Posteriors of Layman Metrics

SubadultTA <- mean(layman.B[[1]][,"TA"])
SubadultNND <- mean(layman.B[[1]][,"NND"])
SubadultSDNND <- mean(layman.B[[1]][,"SDNND"])
SubadultCD <- mean(layman.B[[1]][,"CD"])
SubadultdY <- mean(layman.B[[1]][,"dY_range"])
SubadultdX <- mean(layman.B[[1]][,"dX_range"])

AdultTA <- mean(layman.B[[2]][,"TA"])
AdultNND <- mean(layman.B[[2]][,"NND"])
AdultSDNND <- mean(layman.B[[2]][,"SDNND"])
AdultCD <- mean(layman.B[[2]][,"CD"])
AdultdY <- mean(layman.B[[2]][,"dY_range"])
AdultdX <- mean(layman.B[[2]][,"dX_range"])

#probabilities if we want, but not sure how useful?
dNr1.lt.dNr2 <- sum(layman.B[[1]][,"CD"] > 
                      layman.B[[2]][,"CD"]) / 
  length(layman.B[[1]][,"NND"])

print(dNr1.lt.dNr2)

#Export mean layman metrics
layman_metrics<- list(A = c(AdultdY, AdultdX, AdultTA, AdultNND, AdultSDNND, AdultCD),
                      B = c(SubadultdY, SubadultdX, SubadultTA, SubadultNND, SubadultSDNND, SubadultCD))
layman_metrics<-as.data.frame(t(as.data.frame(layman_metrics)))
colnames(layman_metrics)<-c("d15NRange", "d13CRange", "HullArea", "NearestNeighborDistance", "SDNND", "CentroidDispersion")
rownames(layman_metrics)<-c("Adults", "Subadults")
write.csv(layman_metrics, file = "layman_metrics.csv")


##########################################################################################################
#SKIN
##########################################################################################################

#subset and format data
skin_siber_df <- skin_working %>% 
  select(d13CEpi, d15NEpi, LifeStage) %>%
  mutate(LifeStage = recode(LifeStage, "Subadult" = 1, "Adult" = 2)) %>%
  drop_na(LifeStage)
skin_siber_df$community <- 1
names(skin_siber_df)<-c("iso1","iso2","group", "community")
summary(skin_siber_df)
skin_siber_df <- as.data.frame(skin_siber_df)
skin_siber_sub <- skin_siber_df %>%
  filter(group==1)
skin_siber_adult <- skin_siber_df %>%
  filter(group==2)

#make siber object
skin_siber_obj <- createSiberObject(skin_siber_df)
skin_siber_sub_obj <- createSiberObject(skin_siber_sub)
skin_siber_adult_obj <- createSiberObject(skin_siber_adult)

# Create lists of plotting arguments to be passed onwards to each 
# of the three plotting functions.
community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.95, lty = 1, lwd = 2)
group.hulls.args     <- list(lty = 2, col = "grey20")

#plot data, can change T/F to add/subtract ellipses/hulls
par(mfrow=c(1,1))
plotSiberObject(skin_siber_obj,
                ax.pad = 2, 
                hulls = F, community.hulls.args = community.hulls.args, 
                ellipses = T, group.ellipses.args = group.ellipses.args,
                group.hulls = F, group.hulls.args = group.hulls.args,
                bty = "L",
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'\u2030'),
                ylab = expression({delta}^15*N~'\u2030')
)


#calculate summary stats for each group (lifestage)
group.ML <- groupMetricsML(skin_siber_obj)
print(group.ML)

#calculate layman metrics for communities - not bayesian
community.ML <- communityMetricsML(skin_siber_obj) 
print(community.ML)

#fit bayesian model to data
# options for running jags
parms <- list()
parms$n.iter <- 2 * 10^4   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains

# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

# fit the ellipses which uses an Inverse Wishart prior
# on the covariance matrix Sigma, and a vague normal prior on the 
# means. Fitting is via the JAGS method.
ellipses.posterior <- siberMVN(skin_siber_obj, parms, priors)

# The posterior estimates of the ellipses for each group can be used to
# calculate the SEA.B for each group.
SEA.B <- siberEllipses(ellipses.posterior)

siberDensityPlot(SEA.B, xticklabels = colnames(group.ML), 
                 xlab = c("Community | Group"),
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses on each group"
)

# Add red x's for the ML estimated SEA-c
points(1:ncol(SEA.B), group.ML[3,], col="red", pch = "x", lwd = 2)

# Calculate some credible intervals 
cr.p <- c(0.95, 0.99) # vector of quantiles

# call to hdrcde:hdr using lapply()
SEA.B.credibles <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = cr.p)

# do similar to get the modes, taking care to pick up multimodal posterior
# distributions if present
SEA.B.modes <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$mode},
  prob = cr.p, all.modes=T)

#extract SEA for life stages
SEA_group_modes <- as.data.frame(SEA.B.modes, row.names = "SEA", col.names = c("Subadults", "Adults"))
SEA_group_modes <- as.data.frame(t(SEA_group_modes))
SEA_group_modes

#Prob that Subadults greater than adults
Pg1.lt.g2 <- sum( SEA.B[,1] > SEA.B[,2] ) / nrow(SEA.B)
print(Pg1.lt.g2)###comparing 1 with 2

#Overlap in ppm sq
overlap.G1.G2 <- maxLikOverlap("1.1", "1.2", skin_siber_obj, p = 0.95, n =)
overlap.G1.G2

#Percent overlap
prop.of.both <- as.numeric(overlap.G1.G2["overlap"] / (overlap.G1.G2["area.1"] + overlap.G1.G2["area.2"]))
print(prop.of.both) ####10%of the total area of both elipse areas are shared

##Export Metrics
skin_siber_metrics<-as.data.frame(t(group.ML))
rownames(skin_siber_metrics)<-c("Subadults", "Adults")
skin_siber_metrics$SEAb<-SEA_group_modes[]
names(skin_siber_metrics[,4])<-"SEAb"
write.csv(skin_siber_metrics, file="skin_siber_metrics.csv")

#Plot
jpeg(filename="siberplot.jpeg", width=8, height=9, units="in", res=350)

par(mfrow=c(2,1))
par(mar=c(1,4,0,1))
par(oma=c(3,2,2,2))

with (skin_working, 
      plot(d13CEpi, d15NEpi,
           xlim=c(-20, -11),
           ylim=c(3, 13),
           col=alpha(c("#440154", "#21918c","#fde725")[as.numeric(Sex)], 0.6),
           las=1,
           pch=16,
           cex=1.5,
           cex.axis=1.5,
           cex.lab=1.5,
           xlab="",
           ylab="",
           xaxt="n"))
title(ylab = expression({delta}^15*N~'\u2030'),
      line=2.25, 
      cex.lab=1.5)
legend(-12.5, 13.5, legend=c("Female","Male","Subadult"), col=c("#440154", "#21918c","#fde725"), pch=16, cex=1.2, bty="n")
palette(c("#fde725", "#440154", "#21918c")) #set viridis yellow for subadults (first color)
plotGroupEllipses(skin_siber_sub_obj, plot.args=list(lty=1, lwd=5), iso.order=c(1,2)) #plot subadult ellipse
palette(c("#21918c","#440154")) #set viridis blue for males (first color)
plotGroupEllipses(skin_siber_adult_obj, plot.args=list(lty=1, lwd=3), iso.order=c(1,2)) #plot adult (male ellipse)
palette(c("#440154", "#21918c")) #set viridis blue for females (first color)
plotGroupEllipses(skin_siber_adult_obj, plot.args=list(lty=3, lwd=3), iso.order=c(1,2)) #plot adult (female ellipse)
mtext("A", cex=1.5, side = 3, adj = 0.01,line = -1.5)

with (scute_working, 
      plot(d13C, d15N,
           xlim=c(-20, -11),
           ylim=c(3, 13),
           col=alpha(c("#440154", "#21918c","#fde725")[as.numeric(Sex)], 0.6),
           las=1,
           pch=16,
           cex=1.5,
           cex.axis=1.5,
           xlab="",
           ylab=""))
title(ylab = expression({delta}^15*N~'\u2030'),
      line=2.25, 
      cex.lab=1.5)
plotCommunityHulls(scute_siber_sub_obj, plot.args=list(col="#fde725", lty=1, lwd=5), iso.order=c(1,2))
plotCommunityHulls(scute_siber_adult_obj, plot.args=list(col="#21918c", lty=1, lwd=3), iso.order=c(1,2))
plotCommunityHulls(scute_siber_adult_obj, plot.args=list(col="#440154", lty=3, lwd=3), iso.order=c(1,2))
mtext("B", cex=1.5, side = 3, adj = 0.01,line = -1.5)
mtext(expression({delta}^13*C~'\u2030'), side=1, line=1.5, cex=1.5, outer=TRUE)

dev.off()
