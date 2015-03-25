#Code to import, organize, and explore the SEM-EBSD data obtained for Mytilus edulis (Me5-1)
#2014-02-07
#C.B. Salling

#Set and check current working directory, list current files
setwd("//Users//cbscientist//Desktop//shell_paleothermometer//Mc5-1 Angle//7-Analysis//R")
library(gplots)
library(Hmisc)

#Import angle spread data
angle_spread <- as.data.frame(read.table("angleSpread.csv",header=TRUE,sep=","))

#Import layer thickness data
lay_th_mic <- as.data.frame(read.table("Me5-1_layerThickness.csv",header=TRUE,sep=","))

#Subset data and calculate mean layer thicknesses for each time point
#----------------------------------------------------------------------------------
tp1_lt  <- subset(lay_th_mic,lay_th_mic$tp == 1)
tp3_lt  <- subset(lay_th_mic,lay_th_mic$tp == 3)
tp4_lt  <- subset(lay_th_mic,lay_th_mic$tp == 4)
tp5_lt  <- subset(lay_th_mic,lay_th_mic$tp == 5)
tp6_lt  <- subset(lay_th_mic,lay_th_mic$tp == 6)
tp7_lt  <- subset(lay_th_mic,lay_th_mic$tp == 7)
tp8_lt  <- subset(lay_th_mic,lay_th_mic$tp == 8)
tp9_lt  <- subset(lay_th_mic,lay_th_mic$tp == 9)
tp10_lt <- subset(lay_th_mic,lay_th_mic$tp ==10)
tp11_lt <- subset(lay_th_mic,lay_th_mic$tp ==11)

tp1_meanlt  <-  mean(tp1_lt$layer_thick_mic)
tp3_meanlt  <-  mean(tp3_lt$layer_thick_mic)
tp4_meanlt  <-  mean(tp4_lt$layer_thick_mic)
tp5_meanlt  <-  mean(tp5_lt$layer_thick_mic)
tp6_meanlt  <-  mean(tp6_lt$layer_thick_mic)
tp7_meanlt  <-  mean(tp7_lt$layer_thick_mic)
tp8_meanlt  <-  mean(tp8_lt$layer_thick_mic)
tp9_meanlt  <-  mean(tp9_lt$layer_thick_mic)
tp10_meanlt <- mean(tp10_lt$layer_thick_mic)
tp11_meanlt <- mean(tp11_lt$layer_thick_mic)

tp1_sdlt  <-  sd(tp1_lt$layer_thick_mic)
tp3_sdlt  <-  sd(tp3_lt$layer_thick_mic)
tp4_sdlt  <-  sd(tp4_lt$layer_thick_mic)
tp5_sdlt  <-  sd(tp5_lt$layer_thick_mic)
tp6_sdlt  <-  sd(tp6_lt$layer_thick_mic)
tp7_sdlt  <-  sd(tp7_lt$layer_thick_mic)
tp8_sdlt  <-  sd(tp8_lt$layer_thick_mic)
tp9_sdlt  <-  sd(tp9_lt$layer_thick_mic)
tp10_sdlt <- sd(tp10_lt$layer_thick_mic)
tp11_sdlt <- sd(tp11_lt$layer_thick_mic)

mean_lts <- c(tp1_meanlt,tp3_meanlt,tp4_meanlt,tp5_meanlt,tp6_meanlt,
           tp7_meanlt,tp8_meanlt,tp9_meanlt,tp10_meanlt,tp11_meanlt)
sd_lts <- c(tp1_sdlt,tp3_sdlt,tp4_sdlt,tp5_sdlt,tp6_sdlt,tp7_sdlt,
            tp8_sdlt,tp9_sdlt,tp10_sdlt,tp11_sdlt)
weights_lts <- 1/sd_lts  # create weights using the inverse standard deviation

#----------------------------------------------------------------------------------

#Continuity - established visually (see powerpoint _bycontinuity)
cont <- c(0.6,0.2,0.5,0.9,0.8,1,0.7,0.4,0.1,0.3)

#Merge layer thickness and angle spread data and subset resulting data frame
angle_spread <- data.frame(angle_spread,lt_mic=mean_lts)
angle_spread <- data.frame(angle_spread,conty=cont)
angle_spread <- data.frame(angle_spread,lt_mic_sd=sd_lts)
angle_spread$lt_mic <- round(angle_spread$lt_mic,2)
angle_spread$lt_mic_sd <- round(angle_spread$lt_mic_sd,2)
low_temp <- subset(angle_spread,angle_spread$temp_cel>4)

plotCI(angle_spread$conty,angle_spread$angle_deg,
       uiw=.05,liw=.05,
       err="x",
       ylim=c(10,45),xlim=c(0,1.2),
       xlab="Continuity of Pole Figures(AU)",
       ylab="Angle Spread Footprint(°)",
       main="Angle Spread Footprint vs. Continuity\nof Pole Figures in Me",
       pch=19)
par(new=TRUE)
errbar(angle_spread$conty,angle_spread$angle_deg,
       angle_spread$angle_deg+3,angle_spread$angle_deg-3,
       ylim=c(10,45),xlim=c(0,1.2),
       yaxt='n',
       xaxt='n',
       xlab='',
       ylab='',
       pch='.')
#plot(angle_spread$conty,angle_spread$lt_mic)
#plot(angle_spread$conty,angle_spread$temp_cel)


#Create linear fit of high temp angle spread data and isolate R^2 value
lm.ang <- lm(angle_spread$angle_deg ~ angle_spread$temp_cel)
rs.ang <- summary(lm.ang)$r.squared
slope.ang <- summary(lm.ang)$coefficients[2, 1]
summary(lm.ang)
slope.angtext <- as.character(round(slope.ang,3))
corco.ang <- cor(angle_spread$angle_deg,angle_spread$temp_cel) # Isolates Pearson product correlation
#plot(lm.ang) #Checks fitting results

#Plot angle spread vs. environmental temperature along with high-temp data fit
plot(angle_spread$angle_deg~angle_spread$temp_cel,
     main="Angle Spread Footprint vs. Environmental\nTemperature in Me",
     ylim=c(0,90),xlim=c(0,22),
     xlab="Environmental Temperature(°C)",
     ylab="Angle Spread Footprint(°)",
     pch=19)
abline(lm.ang)
text(13,45, bquote(R == .(round(corco.ang, 2))))
text(13,53,paste("Linear Fit,\nSlope = ",slope.angtext))
#write.table(angle_spread,file="Me5_1_angleSpread.csv",append=FALSE,quote=FALSE,sep=",",
#            eol="\n",na="NA",dec=".",row.names=FALSE)

# plotCI(angle_spread$temp_cel,angle_spread$angle_deg,
#        uiw=1.5,liw=1.5,
#        err="x",
#        ylim=c(0,90),xlim=c(0,22),
#        xlab="Environmental Temperature(°C)",
#        ylab="Angle Spread Footprint(°)",
#        main="Angle Spread Footprint vs. Environmental\nTemperature in Me",
#        pch=19)
# par(new=TRUE)
# errbar(angle_spread$temp_cel,angle_spread$angle_deg,
#        angle_spread$angle_deg+3,angle_spread$angle_deg-3,
#        ylim=c(0,90),xlim=c(0,22),
#        yaxt='n',
#        xaxt='n',
#        xlab='',
#        ylab='',
#        pch='.')
# abline(lm.ang)
# text(13,45, bquote(R^2 == .(round(rs.ang, 2))))
# text(13,53,paste("Linear Fit,\nSlope = ",slope.angtext))

#There exists no clear correlation between distance from first formed nacre layer
#and angle spread footprint, even after removing the low temperature data
plot(angle_spread$dist1st_mic,angle_spread$angle_deg,
     ylim=c(0,90),
     xlab="Distance from First Formed Nacre Layer (µm)",
     ylab="Angle Spread Footprint(°)",
     main="Angle Spread Footprint vs. Distance from\nFirst Formed Nacre Layer")
lm.dist1 <- lm(angle_spread$angle_deg ~ angle_spread$dist1st_mic)
lm.dist2 <- lm(low_temp$angle_deg ~ low_temp$dist1st_mic)
r2dist1 <- summary(lm.dist1)$r.squared
r2dist2 <- summary(lm.dist2)$r.squared

#Create linear fit of layer thickness data vs. environmental temperature and isolate R^2 value
#Poor R^2 value at 0.46
lm.lt <- lm(angle_spread$lt_mic ~ angle_spread$temp_cel,weights=weights_lts)
rs.lt <- summary(lm.lt)$r.squared
slope.lt <- summary(lm.lt)$coefficients[2, 1]
slope.text <- as.character(round(slope.lt,3))
corco.lt <- cor(angle_spread$lt_mic,angle_spread$temp_cel)
#plot(lm.lt) #Checks fitting results
no_out <- subset(angle_spread,angle_spread$TP !=8)
lm.lt2 <- lm(no_out$lt_mic ~ no_out$temp_cel)
rs.lt2 <- summary(lm.lt2)$r.squared

#Plot layer thickness vs. environmental temperature
#plotCI(angle_spread$temp_cel,angle_spread$lt_mic,
#       uiw=1.5,liw=1.5,
#       err="x",
#       ylim=c(0.4,1),xlim=c(0,22),
#       xlab="Environmental Temperature(°C)",
#       ylab="Mean Layer Thickness(µm)",
#       pch=19)
#par(new=TRUE)
plotCI(angle_spread$temp_cel,angle_spread$lt_mic,
       uiw=angle_spread$lt_mic_sd,liw=angle_spread$lt_mic_sd,
       err="y",
       ylim=c(0.4,1),xlim=c(0,22),
       xlab="Environmental Temperature(°C)",
       ylab="Mean Layer Thickness(µm)",
       pch=19)
abline(lm.lt)
title(main="Mean Layer Thickness(µm) vs. Environmental\nTemperature in Me")
text(5,.85, bquote(R == .(round(corco.lt, 2))))
text(5,.93,paste("Weighted Linear fit,\nSlope = ",slope.text))

write.table(angle_spread,file="Me5-1_plottedData.csv",append=FALSE,quote=FALSE,sep=",",
            eol="\n",na="NA",dec=".",row.names=FALSE)
