rm(list = ls())

library(xlsx)
library(dplyr)
library(ggplot2)
library(Rmisc)
library(gridExtra)
library(grid)
library(multcompView)
library(lawstat)

RR <- read.csv("C:/Users/Jesse/Desktop/Neonicotinoids_15/Data/Column-Study/Raw-Data-Outputs/TMX_root_rhiz_soil_R.csv")
attach(RR)


##JBR
##5/31/17

###Root and rhizosphere soil analysis

##for "no plant" volumns seed soil .. ~ seed soil = rhizosphere soil

###zeros were also added to "Ctr" TMX mass and conc columns to make compatible df lenths in r
##in reality no root soil or rhizosphere soil existed for Ctr treatemts




###get ng to micrograms

RR$TMX.microg <- RR$TMX.ng/1000



##start subsetting

#root soil
root.soil <- subset(RR, Section == "root soil")

#rhizosphere soil

rhiz.soil <- subset(RR, Section == "rhiz")



##make a lumped average: "plant-associated soil" for TMX conc

plant_ass <- root.soil[, c(1, 3:5)] ##new df with these columns



plant_ass$TMX.ppb.lumped <- ifelse(root.soil$Stage == "V5_Ctr", rhiz.soil$TMX.ppb,
                                   
                                (root.soil$TMX.ppb + rhiz.soil$TMX.ppb)/2) ##just keep "seed soil" conc


##can't do this for mass and include "ctr" treatment

# plant_ass$TMX.microg.lumped <- ifelse(root.soil$Stage == "V5-Ctr", rhiz.soil$TMX.microg,
#                                       
#                                       (root.soil$TMX.microg + rhiz.soil$TMX.microg)/2) ##



##fix with ifelse statements





  
##make a lumped average: "plant-associated soil" for TMX conc
  
#root.soil$TMX.ppb.lumped <- (root.soil$TMX.ppb + root.soil$TMX.ppb)/2
  
  
# rhiz.soil$TMX.ppb.lumped <- (rhiz.soil$TMX.ppb + rhiz.soil$TMX.ppb)/2
# 
# 
# 
# ##make a lumped average: "plant-associated soil" for TMX mass
# 
# root.soil$TMX.microg.lumped <- (root.soil$TMX.microg + root.soil$TMX.microg)/2
# 
# 
# rhiz.soil$TMX.microg.lumped <- (rhiz.soil$TMX.microg + rhiz.soil$TMX.microg)/2

  
  
  
  

###get R to treat these as characters:

#  RR$Texture <- as.character(RR$Texture)
# # 
#  RR$Stage <- as.character(RR$Stage)
# # 
#  RR$Combo <- as.character(RR$Combo)
# 
# 
# 
# 
# 
#  plant_ass <- data.frame("Sample.ID" = levels(RR$Sample.ID),
#                        "Texture" = NA,
#                        "Stage" = NA,
#                        "Combo" = NA,
#                        "TMX.ppb.lumped" = NA,
#                        "TMX.microg.lumped" = NA)
# 
# keep_cols <- c("Sample.ID", "Texture", "Stage", "Combo")
# 
# for (i in levels(RR$Sample.ID)) {
# 
#   temp_df <- RR[RR$Sample.ID == i,]
# 
#           plant_ass[plant_ass$Sample.ID == i, c(5)] <- (root.soil$TMX.ppb[i] + rhiz.soil$TMX.ppb[i])/2
# 
#           plant_ass[plant_ass$Sample.ID == i, c(6)] <- (root.soil$TMX.microg[i] + rhiz.soil$TMX.microg[i])/2
# 
# 
#           plant_ass[plant_ass$Sample.ID == i, c(1:4)] <- temp_df[1, c(1,3,4,5)]
# 
#         }
# 
#         plant_ass

#    
#         
        
        


#####################################################################################################

#### TMX concentration  stats on lumped TMX conc**************************************************


plant_ass

boxplot(plant_ass$TMX.ppb.lumped)#non-normal.. skewed left

hist(plant_ass$TMX.ppb.lumped) # non-normal.. skewed left

qqnorm(plant_ass$TMX.ppb.lumped) ## non-normal.. skewed left

qqline(plant_ass$TMX.ppb.lumped)## non-normal.. skewed left

plant_ass$log.TMX.ppb.lumped <- log(plant_ass$TMX.ppb.lumped) ##looks good

hist(plant_ass$log.TMX.ppb.lumped)

qqnorm(plant_ass$log.TMX.ppb.lumped) ## 

qqline(plant_ass$log.TMX.ppb.lumped)



bartlett.test(plant_ass$log.TMX.ppb.lumped ~ interaction(plant_ass$Texture, plant_ass$Stage))  ## 

fligner.test(plant_ass$log.TMX.ppb.lumped ~ interaction(plant_ass$Texture, plant_ass$Stage)) ## HOV met

levene.test(plant_ass$log.TMX.ppb.lumped ~ plant_ass$Texture * plant_ass$Stage) ##
# 
# boxplot(BS_top$log10.TMX ~ BS_top$Texture * BS_top$Stage, ylab = "TMX (ng g^-1)", xlab = "Trt") ##
# 
# interaction.plot(BS_top$Texture, BS_top$Stage, BS_top$log10.TMX) ## 

ANOVA.plant.ass <- aov(plant_ass$log.TMX.ppb.lumped ~ plant_ass$Texture * plant_ass$Stage)

summary.plant.ass <- summary(ANOVA.plant.ass)

summary.plant.ass ## stage and texture sig

Tukey.plant.ass <- TukeyHSD(ANOVA.plant.ass)

Tukey.plant.ass

letters.plant.ass <- multcompLetters4(ANOVA.plant.ass, Tukey.plant.ass) ## 

letters.plant.ass

# $`plant_ass$Texture`
# Sand Clay 
# "a"  "b" 
# 
# $`plant_ass$Stage`
# V1     V3 V5_Ctr     V5 
# "a"    "b"    "c"    "d" 
# 
# $`plant_ass$Texture:plant_ass$Stage`
# Sand:V1     Clay:V1     Sand:V3     Clay:V3 Sand:V5_Ctr Clay:V5_Ctr     Clay:V5     Sand:V5 
# "a"        "ab"         "b"        "bc"        "cd"        "de"         "e"         "e" 

#####################################################################################################







#########do stats without lumped averages


rhiz.soil

###starting with rhizosphere... stats on TMX conc


boxplot(rhiz.soil$TMX.ppb) #

hist(rhiz.soil$TMX.ppb) # 

qqnorm(rhiz.soil$TMX.ppb) ## 

qqline(rhiz.soil$TMX.ppb) ## 

rhiz.soil$log.TMX.ppb <- log(rhiz.soil$TMX.ppb) ##

hist(rhiz.soil$log.TMX.ppb)

qqnorm(rhiz.soil$log.TMX.ppb) ## 

qqline(rhiz.soil$log.TMX.ppb) ## normalish

bartlett.test(rhiz.soil$log.TMX.ppb ~ interaction(rhiz.soil$Texture, plant_ass$Stage))  ## 

fligner.test(rhiz.soil$log.TMX.ppb ~ interaction(rhiz.soil$Texture, plant_ass$Stage)) ## HOV met

levene.test(hiz.soil$log.TMX.ppb ~ rhiz.soil$Texture * rhiz.soil$Stage) ##
# 



ANOVA.rhiz.soil.TMX.ppb <- aov(rhiz.soil$log.TMX.ppb ~ rhiz.soil$Texture * rhiz.soil$Stage)

summary.rhiz.soil.TMX.ppb <- summary(ANOVA.rhiz.soil.TMX.ppb)

summary.rhiz.soil.TMX.ppb ## Stage signficant

Tukey.rhiz.soil.TMX.ppb <- TukeyHSD(ANOVA.rhiz.soil.TMX.ppb)

Tukey.rhiz.soil.TMX.ppb

letters.rhiz.soil.TMX.ppb <- multcompLetters4(ANOVA.rhiz.soil.TMX.ppb, Tukey.rhiz.soil.TMX.ppb) ## 

letters.rhiz.soil.TMX.ppb 

#`rhiz.soil$Stage`
# V1     V3 V5_Ctr     V5 
# "a"    "b"    "c"    "d" 
# 
# $`rhiz.soil$Texture:rhiz.soil$Stage`
# Sand:V1     Clay:V1     Sand:V3     Clay:V3 Sand:V5_Ctr Clay:V5_Ctr     Clay:V5     Sand:V5 
# "a"        "ab"        "ab"        "bc"        "cd"        "de"        "de"         "e" 



####now work on rhizpshphere TMX mass stats


boxplot(rhiz.soil$TMX.microg) #

hist(rhiz.soil$TMX.microg) # left skewed

qqnorm(rhiz.soil$TMX.microg) ## 

qqline(rhiz.soil$TMX.microg) ## 

rhiz.soil$log.TMX.microg <- log(rhiz.soil$TMX.microg) ##

hist(rhiz.soil$log.TMX.microg) ##

qqnorm(rhiz.soil$log.TMX.microg) ## 

qqline(rhiz.soil$log.TMX.microg) ## very skewed


##nah, try rank transforming

rhiz.soil$rank.TMX.microg <- rank(rhiz.soil$TMX.microg)


bartlett.test(rhiz.soil$rank.TMX.microg  ~ interaction(rhiz.soil$Texture, plant_ass$Stage))  ## 

fligner.test(rhiz.soil$rank.TMX.microg  ~ interaction(rhiz.soil$Texture, plant_ass$Stage)) ## HOV met.. barely

levene.test(rhiz.soil$rank.TMX.microg  ~ rhiz.soil$Texture * rhiz.soil$Stage) ##
# 



ANOVA.rhiz.soil.TMX.microg  <- aov(rhiz.soil$rank.TMX.microg ~ rhiz.soil$Texture * rhiz.soil$Stage)

summary.rhiz.soil.TMX.microg <- summary(ANOVA.rhiz.soil.TMX.microg)

summary.rhiz.soil.TMX.microg ## Stage signficant

Tukey.rhiz.soil.TMX.microg <- TukeyHSD(ANOVA.rhiz.soil.TMX.microg)

Tukey.rhiz.soil.TMX.microg

letters.rhiz.soil.TMX.microg <- multcompLetters4(ANOVA.rhiz.soil.TMX.microg, Tukey.rhiz.soil.TMX.microg) ## 

letters.rhiz.soil.TMX.microg



# 
# $`rhiz.soil$Texture`$LetterMatrix
# a
# Sand TRUE
# Clay TRUE
# 
# 
# $`rhiz.soil$Stage`
# V3     V1     V5 V5_Ctr 
# "a"   "ab"   "ab"    "b" 
# 
# $`rhiz.soil$Texture:rhiz.soil$Stage`
# Clay:V3     Sand:V1     Sand:V3     Clay:V1     Clay:V5     Sand:V5 Sand:V5_Ctr Clay:V5_Ctr 
# "a"        "ab"        "ab"        "ab"        "ab"        "ab"        "ab"         "b" 

#####################################################################################################




#####################################################################################################
root.soil

root.soil <- root.soil[root.soil$Stage != "V5_Ctr",]

root.soil

####now root soil concentration (without the Ctr trt)

boxplot(root.soil$TMX.ppb) #

hist(root.soil$TMX.ppb) # skewed left

qqnorm(root.soil$TMX.ppb) ## 

qqline(root.soil$TMX.ppb) ## 

root.soil$log.TMX.ppb <- log(root.soil$TMX.ppb) ##

hist(root.soil$log.TMX.ppb)

qqnorm(root.soil$log.TMX.ppb) ## 

qqline(root.soil$log.TMX.ppb) ## normal

bartlett.test(root.soil$log.TMX.ppb ~ interaction(root.soil$Texture, root.soil$Stage))  ## 

fligner.test(root.soil$log.TMX.ppb ~ interaction(root.soil$Texture, root.soil$Stage)) ## HOV met

levene.test(root.soil$log.TMX.ppb ~ root.soil$Texture * root.soil$Stage) ##
# 



ANOVA.root.soil.TMX.ppb <- aov(root.soil$log.TMX.ppb ~ root.soil$Texture * root.soil$Stage)

summary.root.soil.TMX.ppb <- summary(ANOVA.root.soil.TMX.ppb)

summary.root.soil.TMX.ppb ## Stage signficant and texture

Tukey.root.soil.TMX.ppb <- TukeyHSD(ANOVA.root.soil.TMX.ppb)

Tukey.root.soil.TMX.ppb

letters.root.soil.TMX.ppb <- multcompLetters4(ANOVA.root.soil.TMX.ppb, Tukey.root.soil.TMX.ppb) ## 

letters.root.soil.TMX.ppb 


# $`root.soil$Texture`
# Sand Clay 
# "a"  "b" 
# 
# $`root.soil$Stage`
# V1  V3  V5 
# "a" "b" "c" 
# 
# $`root.soil$Texture:root.soil$Stage`
# Sand:V1 Clay:V1 Sand:V3 Clay:V3 Sand:V5 Clay:V5 
# "a"    "ab"    "bc"     "c"     "d"     "d" 
# 


#######now mass of TMX in root soil


boxplot(root.soil$TMX.microg) #

hist(root.soil$TMX.microg) # probably alright but a little skewed left

qqnorm(root.soil$TMX.microg) ## 

qqline(root.soil$TMX.microg) ## 

root.soil$log.TMX.microg <- log(root.soil$TMX.microg) ##

hist(root.soil$log.TMX.microg)## normal

qqnorm(root.soil$log.TMX.microg) ## 

qqline(root.soil$log.TMX.microg) ## normal

bartlett.test(root.soil$log.TMX.microg ~ interaction(root.soil$Texture, root.soil$Stage))  ## 

fligner.test(root.soil$log.TMX.microg ~ interaction(root.soil$Texture, root.soil$Stage)) ## HOV met

levene.test(root.soil$log.TMX.microg ~ root.soil$Texture * root.soil$Stage) ##
# 



ANOVA.root.soil.TMX.microg <- aov(root.soil$log.TMX.microg ~ root.soil$Texture * root.soil$Stage)

summary.root.soil.TMX.microg <- summary(ANOVA.root.soil.TMX.microg)

summary.root.soil.TMX.microg ## Stage signficant and texture

Tukey.root.soil.TMX.microg <- TukeyHSD(ANOVA.root.soil.TMX.ppb)

Tukey.root.soil.TMX.microg

letters.root.soil.TMX.microg <- multcompLetters4(ANOVA.root.soil.TMX.microg, Tukey.root.soil.TMX.microg) ## 

letters.root.soil.TMX.microg 

# 
# $`root.soil$Texture`
# Sand Clay 
# "a"  "b" 
# 
# $`root.soil$Stage`
# V1  V3  V5 
# "a" "b" "c" 
# 
# $`root.soil$Texture:root.soil$Stage`
# Sand:V1 Sand:V3 Sand:V5 Clay:V3 Clay:V1 Clay:V5 
# "a"    "bc"     "d"     "b"    "ac"     "d" 
# 

#####################################################################################################







#####################################################################################################




###########now prep for plots******

###rhiz TMX.ppb

rhiz.soil

rhiz.soil.TMX.ppb <- summarySE(rhiz.soil, measurevar = "TMX.ppb", groupvars = c("Texture", "Stage"))
  
rhiz.soil.TMX.ppb 


letters.rhiz.soil.TMX.ppb ##now add letters

rhiz.soil.TMX.ppb$letters <- c("ab", "bc", "de", "de", "a", "ab", "e", "cd")


##rhiz TMX mass

rhiz.soil

rhiz.soil.TMX.microg <- summarySE(rhiz.soil, measurevar = "TMX.microg", groupvars = c("Texture", "Stage"))

rhiz.soil.TMX.microg 

letters.rhiz.soil.TMX.microg ##now add letters

rhiz.soil.TMX.microg 

rhiz.soil.TMX.microg$letters <- c("ab", "a", "ab", "b", "ab", "ab", "ab", "ab")


####root soil TMX conc

root.soil

root.soil.TMX.ppb <- summarySE(root.soil, measurevar = "TMX.ppb", groupvars = c("Texture", "Stage"))

root.soil.TMX.ppb

letters.root.soil.TMX.ppb

root.soil.TMX.ppb

root.soil.TMX.ppb$letters <- c("ab", "c", "d", "a", "bc", "d")

root.soil.TMX.ppb


####root soil TMX mass

root.soil

root.soil.TMX.microg <- summarySE(root.soil, measurevar = "TMX.microg", groupvars = c("Texture", "Stage"))

root.soil.TMX.microg

letters.root.soil.TMX.microg

root.soil.TMX.microg

root.soil.TMX.microg$letters <- c("ac", "b", "d", "a", "bc", "d")

root.soil.TMX.microg



###now on plant associated soil

plant_ass

plant_ass.TMX.ppb.lumped <- summarySE(plant_ass, measurevar = "TMX.ppb.lumped", groupvars = c("Texture", "Stage"))

plant_ass.TMX.ppb.lumped

letters.plant.ass

plant_ass.TMX.ppb.lumped

plant_ass.TMX.ppb.lumped$letters <- c("ab", "bc", "e", "de", "a", "b", "e", "cd")

plant_ass.TMX.ppb.lumped
  
##################################################################################################### 
  

##now plot
#####################################################################################################  


cleanup <- theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background =  element_blank(),
                 axis.line = element_line(color = "black"))



plant.ass <- ggplot(plant_ass.TMX.ppb.lumped, aes(Texture, TMX.ppb.lumped, fill = Stage)) +
  
  geom_bar(position = position_dodge(), stat = "identity") +
  
  geom_bar(position = "dodge", stat = "identity", color="black", show.legend =FALSE, lwd = 0.75) + ##add outline to bars
  
  geom_errorbar(aes(ymin = TMX.ppb.lumped, ymax = TMX.ppb.lumped + se),      ##error bars
                
                width = 0.2, position = position_dodge(0.9), lwd = 0.75) +
  
  cleanup +                                          #cleanup
  
  scale_fill_manual(name = "",
                    
                    labels = c("V1", "V3", "V5", "V5, No Plant"), 
                    
                    values = c("blue", "red", "green", "yellow")) +
  
  xlab("") +
  
  ylab(expression(TMX ~ (mu * g ~ kg^{-1}))) + ### pain in the ass way to express superscripts
  
  theme(axis.text.y = element_text(size = 16),  ## change font size of y axis
        
        axis.title = element_text(size = 16, face = "bold")) +
  
  #geom_label(aes(label = letters)) +
  
  #geom_text(aes(label = letters)) +
  
  geom_text(aes(label = letters, y = TMX.ppb.lumped + se, fontface = "italic"), position = position_dodge(width=0.9), vjust= -0.5) +  #####add in signficance letter--adjusted for errror bars
  
  theme(axis.text.x = element_text(size = 16),  ## change font size of y axis
        
        axis.title = element_text(size = 16, face = "bold")) +
  
  theme(axis.line.x = element_line(size = 1), ##change axis line size
        
        axis.line.y = element_line(size = 1),
        
        panel.border = element_rect(color = "black", fill = NA, size = 1.75), ##add border
        
        legend.position = c(0.2, 0.8), ## move legend
        
        legend.margin = margin(0, 0, 0, 0, "mm"), ##remove margin on legend
        
        legend.text = element_text( size = 16), ##legend text
        
        #legend.background = element_rect(color = "black", size = 1), ## add legend border
        
        #legend.key.height=unit(1,"line"), 
        
        #legend.key.width=unit(1,"line")
        
        legend.key = element_rect(color = "black", size = 0.75), ## put line around legend key
        
        legend.direction = "vertical", 
        
        axis.ticks.y = element_line(size = 1, color = "black"),  ##enlargen ticks--problem at op right corner
        
        legend.key.size = unit(7, 'mm'),
        
        legend.key.height = unit(7, "mm"), 
        
        plot.margin = unit(c(1.5, 1,1.5,-1.5),"mm"),
        
        axis.title.y = element_text(vjust = -1)) +
  
  #plot.title = element_text( Size = 18)) +
  
  #guides(guide_legend(nrow=4)) +
  
  #legend.spacing.y = unit(1, "mm") +
  
  #legend.key.size = unit(1, "line")) 
  
  #scale_fill_manual(values=values, labels=setNames(paste(labels, " "), entries)) +
  
ggtitle("TMX in Plant-Associated Soil") +
  
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5, vjust = 1)) + ## added title and adjusted
  
  #scale_y_continuous(expand = c(0,0), limit = c(0, 110)) + scale_x_discrete(expand = c(0.03,0)) ### removed excess space btw bars and axis
  
  scale_y_continuous(expand = c(0,0), limit = c(0, 1800)) +

  scale_x_discrete(expand = c(0.03,0)) ### removed excess space btw bars and axis



plant.ass


pdf("C:/Users/Jesse/Desktop/Neonicotinoids_15/Data/Column-Study/Raw-Data-Outputs/plant.associated.soil.pdf")
plant.ass
dev.off()










