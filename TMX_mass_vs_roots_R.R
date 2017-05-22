rm(list = ls())
library(dplyr)
library(ggplot2)
library(Rmisc)
library(gridExtra)
library(grid)
library(multcompView)
library(lawstat)
library(car)


MvR <- read.csv("C:/Users/Jesse/Desktop/Neonicotinoids_15/Data/Column-Study/Raw-Data-Outputs/TMX_mass_vs_roots-R.csv")

attach(MvR)

#install.packages("dplyr")
#install.packages("ggplot2")
#install.packages("Rmisc")
#install.packages("grid")
#install.packages("gridExtra")
#install.packages("multcompView")
#install.packages("lawstat")
#install.packages("car")




MvR

MvR <- MvR[-c(24:41),] ##remove extra junk

MvR <- MvR[, -8] ##remove extra junk


MvR


###transform + calcualtions

MvR$log.TMX <- log10(MvR$TMX.ng.30.60.cm)

MvR$col.length.cm <- 60

MvR$col.volume.cc <- (pi * (10)^2 * MvR$col.length.cm) ##root length in cm

MvR

MvR$rt.length.dens.cm_cc <- MvR$Root.depth.cm/MvR$col.volume.cc ##length density (root length-cm)/(volume of soil-cm^3)

MvR$biomass.per.cm.root <- MvR$Dry.weight.g/MvR$Root.depth.cm ### root biomass/root depth

MvR$TMX.micrg.30.60.cm <- MvR$TMX.ng.30.60.cm/1000 ###convert to micrograms

##now subset v3-v5

late.stage <- subset(MvR, Stage != "V1")

late.stage

##plot

late.stage


####try log of mass below 30 cm as a function of root depth
plot(MvR$Root.depth.cm, MvR$log.TMX) ###just trying 

abline(lm(MvR$log.TMX ~ MvR$Root.depth.cm), col="red")

summary(lm(MvR$log.TMX ~ MvR$Root.depth.cm))

###Adjusted R-squared:  0.8221 

MvR




####now log mass vs root length density (root length/volume soil)
plot(MvR$rt.length.dens.cm_cc, MvR$log.TMX) ###just trying 

abline(lm(MvR$log.TMX ~ MvR$rt.length.dens.cm_cc), col="red")

summary(lm(MvR$log.TMX ~ MvR$rt.length.dens.cm_cc))

###Adjusted R-squared:  0.8221 same result




###try TMX mass vs root dry mass

qqnorm(MvR$TMX.ng.30.60.cm) ##not normal

qqline(MvR$TMX.ng.30.60.cm) ##not normal

qqnorm(MvR$Dry.weight.g) ##not normal

qqline(MvR$Dry.weight.g) ##not normal.. but can still plot the relationship

plot(MvR$Dry.weight.g, MvR$TMX.ng.30.60.cm) ###just trying 

abline(MvR$TMX.ng.30.60.cm ~ MvR$Dry.weight.g, col="red")

summary(lm(MvR$TMX.ng.30.60.cm ~ MvR$Dry.weight.g))

#####deeeeyum! :Adjusted R-square = 0.8837 *************************** probably becasue this includes v1

###logtMX vs dry root mass?

plot(MvR$Dry.weight.g, MvR$log.TMX) ###just trying 

abline(MvR$log.TMX ~ MvR$Dry.weight.g, col="red")

summary(lm(MvR$log.TMX ~ MvR$Dry.weight.g))

####weaker relationship Adjusted R-squared:  0.7461





###try root dry mass/cm root

plot(MvR$biomass.per.cm.root, MvR$TMX.ng.30.60.cm) ###just trying 

abline(MvR$TMX.ng.30.60.cm ~ MvR$biomass.per.cm.root, col="red")

summary(lm(MvR$TMX.ng.30.60.cm ~ MvR$biomass.per.cm.root))

#### Adjusted R-squared:  0.8294 ---pretty good

###adjust this with log TMX

plot(MvR$biomass.per.cm.root, MvR$log.TMX) ###just trying 

abline(MvR$log.TMX ~ MvR$biomass.per.cm.root, col="red")

summary(lm(MvR$log.TMX ~ MvR$biomass.per.cm.root))

##not as good  Adjusted R-squared:  0.6953 



##focus on biomass as predictor of TMX mass below 30 cm, and root length as a predictor of log(TMX)**********************

MvR



####log of mass below 30 cm as a function of root depth

plot(MvR$Root.depth.cm, MvR$log.TMX) ###

abline(lm(MvR$log.TMX ~ MvR$Root.depth.cm), col="black")

summary(lm(MvR$log.TMX ~ MvR$Root.depth.cm))

###Adjusted R-squared:  0.8221 
dev.off()



pdf("C:/Users/Jesse/Desktop/Neonicotinoids_15/Data/Column-Study/Raw-Data-Outputs/TMX_mass_vs_roots.plot1.pdf")

par(mar = c(5.1,6,4.1,4)) ##in base its c(bottom, left, top, right) in ggplot its c(top, right, bottom, left)


 plot(x = NULL, y = NULL, axes = FALSE, xlab = "Rooting Depth (cm)",
      
      ylab = "log[TMX (ng)]", ###blank plot with title and range specified
      
      xlim = c(0, 70), 
     
      ylim = c(0, 6), cex.lab = 2)

points(MvR$Root.depth.cm[Stage == "V3"], MvR$log.TMX[Stage == "V3"],
       
       pch = 21, col = "black", bg = "red",  cex = 2, lwd = 2)#add points

points(MvR$Root.depth.cm[Stage == "V5"], MvR$log.TMX[Stage == "V5"],
       
       pch = 21, col = "black", bg = "yellow", cex = 2, lwd = 2)

abline(lm(MvR$log.TMX ~ MvR$Root.depth.cm), col="black", lwd = 2) ###add line to plot

axis(1, lwd = 2, cex.axis = 1.5)

axis(2, lwd = 2, cex.axis = 1.5)

box(col = "black", lwd = 2)

title( main = "Mass of TMX in Bulk Soil (0-60 cm) vs Rooting Depth", cex.main = 1.5)

legend("bottomright", inset = 0.05, c("V3", "V5"), pch = c(21, 21), col = c("black", "black"),
       
       pt.bg = c("red", "yellow"), cex = c(2, 2), pt.lwd = c(2, 2), horiz = T, bty = "n")

label1 <- bquote(italic(R)^2 == 0.822)

text(x = 11, y = 5.75, label1, cex = 2, font = 2)

dev.off()







###now plot TMX mass vs root dry mass

plot(MvR$Dry.weight.g, MvR$TMX.ng.30.60.cm) ###just trying 

abline(MvR$TMX.ng.30.60.cm ~ MvR$Dry.weight.g, col="red")

summary(lm(MvR$TMX.ng.30.60.cm ~ MvR$Dry.weight.g))

#####deeeeyum! :Adjusted R-square = 0.8837 *************************** probably becasue this includes v1
dev.off()


pdf("C:/Users/Jesse/Desktop/Neonicotinoids_15/Data/Column-Study/Raw-Data-Outputs/TMX_mass_vs_roots.plot2.pdf")

par(mar = c(5.1,6,4.1,4)) ##in base its c(bottom, left, top, right) in ggplot its c(top, right, bottom, left)


plot(x = NULL, y = NULL, axes = FALSE, xlab = "Root Dry Mass (g)",
     
     ylab = expression(TMX ~(mu * g)), ###blank plot with title and range specified
     
     xlim = c(0, 5.5), 
     
     ylim = c(0, 85), cex.lab = 2)

points(MvR$Dry.weight.g[Stage == "V1"], MvR$TMX.micrg.30.60.cm[Stage == "V1"],
       
       pch = 21, col = "black", bg = "blue",  cex = 2, lwd = 2)#add points

points(MvR$Dry.weight.g[Stage == "V3"], MvR$TMX.micrg.30.60.cm[Stage == "V3"],
       
       pch = 21, col = "black", bg = "red", cex = 2, lwd = 2)

points(MvR$Dry.weight.g[Stage == "V5"], MvR$TMX.micrg.30.60.cm[Stage == "V5"],
       
       pch = 21, col = "black", bg = "yellow", cex = 2, lwd = 2)

abline(lm(MvR$TMX.micrg.30.60.cm ~ MvR$Dry.weight.g), col="black", lwd = 2) ###add line to plot

axis(1, lwd = 2, cex.axis = 1.5)

axis(2, lwd = 2, cex.axis = 1.5)

box(col = "black", lwd = 2)

title( main = "Mass of TMX in Bulk Soil (0-60 cm) vs Root Dry mass", cex.main = 1.5)

legend("bottomright", inset = 0.05, c("V1", "V3", "V5"), pch = c(21, 21, 21), col = c("black", "black", "black"),
       
       pt.bg = c("blue", "red", "yellow"), cex = c(2, 2, 2), pt.lwd = c(2, 2, 2), horiz = T, bty = "n")

label1 <- bquote(italic(R)^2 == 0.884)

text(x = 1, y = 80, label1, cex = 2, font = 2)

dev.off()



















#######now arrange in a two panel figure
dev.off()

pdf("C:/Users/Jesse/Desktop/Neonicotinoids_15/Data/Column-Study/Raw-Data-Outputs/TMX_mass_vs_roots.combo.pdf",
    
    paper = "USr")

par(mfrow = c(1,2), oma = c(4, 1, 1, 1))
#layout(matrix(c(1,2,3,3), ncol=2, byrow=TRUE), heights=c(3, 0.5))



##plot1
par(mar = c(15,4.5,1,0.5)) ##in base its c(bottom, left, top, right) in ggplot its c(top, right, bottom, left)


plot(x = NULL, y = NULL, axes = FALSE, xlab = "Rooting Depth (cm)",
     
     ylab = "log[TMX (ng)] Below 30 cm", ###blank plot with title and range specified
     
     xlim = c(0, 70), 
     
     ylim = c(0, 6), cex.lab = 1.25)

points(MvR$Root.depth.cm[Stage == "V3"], MvR$log.TMX[Stage == "V3"],
       
       pch = 21, col = "black", bg = "red",  cex = 2, lwd = 2)#add points

points(MvR$Root.depth.cm[Stage == "V5"], MvR$log.TMX[Stage == "V5"],
       
       pch = 21, col = "black", bg = "green", cex = 2, lwd = 2)

abline(lm(MvR$log.TMX ~ MvR$Root.depth.cm), col="black", lwd = 2) ###add line to plot

axis(1, lwd = 2, cex.axis = 1.25)

axis(2, lwd = 2, cex.axis = 1.25)

box(col = "black", lwd = 2)

#title( main = "Mass of TMX in Bulk Soil (0-60 cm) vs Rooting Depth", cex.main = 1.5)

#legend("bottomright", inset = 0.05, c("V3", "V5"), pch = c(21, 21), col = c("black", "black"),
       
#       pt.bg = c("red", "yellow"), cex = c(1.5, 1.5), pt.lwd = c(2, 2), horiz = T, bty = "n")

label1 <- bquote(italic(R)^2 == 0.822)

text(x = 17.5, y = 5.75, label1, cex = 1.25, font = 2)




###plot2

par(mar = c(15,4.5,1,0.5)) ##in base its c(bottom, left, top, right) in ggplot its c(top, right, bottom, left)


plot(x = NULL, y = NULL, axes = FALSE, xlab = "Root Dry Mass (g)",
     
     ylab = expression(TMX ~(mu * g) ~ 
                         
                         Below ~ 30 ~ cm), ###blank plot with title and range specified
     
     xlim = c(0, 5.5), 
     
     ylim = c(0, 85), cex.lab = 1.25)

points(MvR$Dry.weight.g[Stage == "V1"], MvR$TMX.micrg.30.60.cm[Stage == "V1"],
       
       pch = 21, col = "black", bg = "blue",  cex = 2, lwd = 2)#add points

points(MvR$Dry.weight.g[Stage == "V3"], MvR$TMX.micrg.30.60.cm[Stage == "V3"],
       
       pch = 21, col = "black", bg = "red", cex = 2, lwd = 2)

points(MvR$Dry.weight.g[Stage == "V5"], MvR$TMX.micrg.30.60.cm[Stage == "V5"],
       
       pch = 21, col = "black", bg = "green", cex = 2, lwd = 2)

abline(lm(MvR$TMX.micrg.30.60.cm ~ MvR$Dry.weight.g), col="black", lwd = 2) ###add line to plot

axis(1, lwd = 2, cex.axis = 1.25)

axis(2, lwd = 2, cex.axis = 1.25)

box(col = "black", lwd = 2)

#title( main = "Mass of TMX in Bulk Soil (0-60 cm) vs Root Dry mass", cex.main = 1.5)

#legend("bottomright", inset = 0.05, c("V1", "V3", "V5"), pch = c(21, 21, 21), col = c("black", "black", "black"),
       
#       pt.bg = c("blue", "red", "yellow"), cex = c(1.5, 1.5, 1.5), pt.lwd = c(2, 2, 2), horiz = T, bty = "n")

label1 <- bquote(italic(R)^2 == 0.884)

text(x = 1.5, y = 80, label1, cex = 1.25, font = 2)






###3

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 4, 10, 0), new = TRUE)

plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

legend("center", inset = 0.05, c("V1", "V3", "V5"), pch = c(21, 21, 21),
       
       col = c("black", "black", "black"),
       
       pt.bg = c("blue", "red", "green"), cex = c(1.75, 1.75, 1.75),
       
       pt.lwd = c(2, 2, 2), horiz = T,
       
       yjust = 0.25, border = "black", box.lwd = 2)

dev.off()

?legend()



