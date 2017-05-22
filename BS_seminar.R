rm(list = ls())

BS <- read.csv("C:/Users/Jesse/Desktop/Neonicotinoids_15/Data/Column-Study/Raw-Data-Outputs/Bulk_soil_R.csv")
attach(BS)

library(car)

library(multcompView)

BS

BS <- BS[, c("Sample.ID", "Interval", "Texture", "Stage", "Combo", "TMX")] ## remove missing values

na.omit(BS) ## remove missing values

BS_top <- subset(BS, Interval == "0-30 cm")

BS_top ## 0-30 cm

BS_middle <- subset(BS, Interval == "30-45 cm")

BS_middle ## 30-45 cm

BS_bottom <- subset(BS, Interval == "45-60 cm")

BS_bottom ## 45-60 cm




########## BS_top Stats---2-1ay interaction + trt combo---"V5-Ctr"

BS_top

boxplot(BS_top$TMX)

hist(BS_top$TMX) # distribution in structure study

qqnorm(BS_top$TMX) ## not normal

qqline(BS_top$TMX)

BS_top$log.TMX <- log(BS_top$TMX)

hist(BS_top$log.TMX)

qqnorm(BS_top$log.TMX) ## not normal after natural log transform

qqline(BS_top$log.TMX)


BS_top$log10.TMX <- log10(BS_top$TMX) 

hist(BS_top$log10.TMX) ###normal**** after log base 10 transform

bartlett.test(BS_top$log10.TMX ~ interaction(BS_top$Texture, BS_top$Stage))  ## no HOV

fligner.test(BS_top$log10.TMX ~ interaction(BS_top$Texture, BS_top$Stage)) ## barely meets HOV here P=0.05

levene.test(BS_top$log10.TMX ~ BS_top$Texture * BS_top$Stage) ##

boxplot(BS_top$log10.TMX ~ BS_top$Texture * BS_top$Stage, ylab = "TMX (ng g^-1)", xlab = "Trt") ##

interaction.plot(BS_top$Texture, BS_top$Stage, BS_top$log10.TMX) ## looks like interactions btw  with plant and w/out plant at sand

ANOVA.top <- aov(BS_top$log10.TMX ~ BS_top$Texture * BS_top$Stage)

summary.top <- summary(ANOVA.top)

summary.top ## 

Tukey.top <- TukeyHSD(ANOVA.top)

Tukey.top

letters.top <- multcompLetters4(ANOVA.top, Tukey.top) ## Texture and Stage significant---- V1 V5.Ctr     V5     V3 
##                                                                          "a"    "b"    "b"    "b"

###Try rank Transforming

BS_top$Rank.TMX <- rank(BS_top$TMX)

hist(BS_top$Rank.TMX) ###normal**** 

qqnorm(BS_top$Rank.TMX) ## 

qqline(BS_top$Rank.TMX)

bartlett.test(BS_top$Rank.TMX ~ interaction(BS_top$Texture, BS_top$Stage))  ## 

fligner.test(BS_top$Rank.TMX ~ interaction(BS_top$Texture, BS_top$Stage)) ## meets HOV
levene.test(BS_top$Rank.TMX ~ BS_top$Texture * BS_top$Stage) ##

boxplot(BS_top$Rank.TMX ~ BS_top$Texture * BS_top$Stage, ylab = "TMX (ng g^-1)", xlab = "Trt") ##

interaction.plot(BS_top$Texture, BS_top$Stage, BS_top$Rank.TMX) ## looks like interactions btw  with plant and w/out plant at sand

ANOVA.top1 <- aov(BS_top$Rank.TMX ~ BS_top$Texture * BS_top$Stage)

summary.top1 <- summary(ANOVA.top1)

summary.top1 ## 

Tukey.top1 <- TukeyHSD(ANOVA.top1)

Tukey.top1

letters.top <- multcompLetters4(ANOVA.top1, Tukey.top1) ## Texture and Stage significant   V1 V5.Ctr     V5     V3 
##                                                                          "a"   "ab"    "b"    "b" 

## here V1 and V5 ctr not different

BS_top


######################BS_middle stats (30-45 cm)

BS_middle

boxplot(BS_middle$TMX)

boxplot(BS_middle$TMX ~ BS_middle$Texture * BS_middle$Stage)

hist(BS_middle$TMX) ## very skewed left

qqnorm(BS_middle$TMX)

qqline(BS_middle$TMX)

BS_middle$log10.TMX <- log10(BS_middle$TMX) ## log 10 transform

hist(BS_middle$log10.TMX)

qqnorm(BS_middle$log10.TMX)

qqline(BS_middle$log10.TMX)## looks fine

BS_middle$lsqrt.TMX <- sqrt(BS_middle$TMX) ## sqrt transform

hist(BS_middle$lsqrt.TMX)

qqnorm(BS_middle$lsqrt.TMX)

qqline(BS_middle$lsqrt.TMX)

BS_middle$Rank.TMX <- rank(BS_middle$TMX)

hist(BS_middle$Rank.TMX)

qqnorm(BS_middle$Rank.TMX)

qqline(BS_middle$Rank.TMX)

###stick with log10 TMX ??

bartlett.test(BS_middle$log10.TMX ~ interaction(BS_middle$Texture, BS_middle$Stage))  ## 

fligner.test(BS_middle$log10.TMX ~ interaction(BS_middle$Texture, BS_middle$Stage))## meets HOV

levene.test(BS_middle$log10.TMX ~ BS_middle$Texture * BS_middle$Stage) ##

boxplot(BS_middle$log10.TMX ~ BS_middle$Texture * BS_middle$Stage, ylab = "TMX (ng g^-1)", xlab = "Trt") ##

interaction.plot(BS_middle$Texture, BS_middle$Stage, BS_middle$log10.TMX) ## 

ANOVA.middle <- aov(BS_middle$log10.TMX ~ BS_middle$Texture * BS_middle$Stage)

summary.middle <- summary(ANOVA.middle)

summary.middle ## 

Tukey.middle <- TukeyHSD(ANOVA.middle)

Tukey.middle

multcompLetters4(ANOVA.middle, Tukey.middle)

BS_middle


##remove outlier "nk-v5-ctr-5" ---had large preferential flow pathway in this region attach image...

outlier.V5.sand <- subset(BS_middle, Combo == "Sandy-V5-Ctr")

outlier.V5.sand

boxplot(outlier.V5.sand$TMX) ##visualize the outlier--could be skewing 2 way interaction tests?

BS_middle1 <- subset(BS_middle, Sample.ID != "nk-v5-ctr-5")

BS_middle1

## rerun with log10.TMX

hist(BS_middle1$log10.TMX)

qqnorm(BS_middle1$log10.TMX)

qqline(BS_middle1$log10.TMX)

bartlett.test(BS_middle1$log10.TMX ~ interaction(BS_middle1$Texture, BS_middle$Stage))  ## 

fligner.test(BS_middle1$log10.TMX ~ interaction(BS_middle1$Texture, BS_middle1$Stage))## meets HOV

levene.test(BS_middle1$log10.TMX ~ BS_middle1$Texture * BS_middle1$Stage) ##

boxplot(BS_middle1$log10.TMX ~ BS_middle1$Texture * BS_middle1$Stage, ylab = "TMX (ng g^-1)", xlab = "Trt") ##

interaction.plot(BS_middle1$Texture, BS_middle1$Stage, BS_middle1$log10.TMX) ## 

ANOVA.middle1 <- aov(BS_middle1$log10.TMX ~ BS_middle1$Texture * BS_middle1$Stage)

summary.middle1 <- summary(ANOVA.middle1)

summary.middle1 ## Stage (including no plant-v5 combo) significant---interaction signficant

Tukey.middle1 <- TukeyHSD(ANOVA.middle1)

Tukey.middle1

letters.middle <- multcompLetters4(ANOVA.middle1, Tukey.middle1)###### Good shit!! differences plant effect in both textures! 

##so use BS_middle1


########################################  BS_bottom stats

BS_bottom

hist(BS_bottom$TMX)

BS_bottom$log10.TMX <- log10(BS_bottom$TMX)

hist(BS_bottom$log10.TMX) ###great!

qqnorm(BS_bottom$log10.TMX)

qqline(BS_bottom$log10.TMX) ###great!

bartlett.test(BS_bottom$log10.TMX ~ interaction(BS_bottom$Texture, BS_bottom$Stage))  ## 

fligner.test(BS_bottom$log10.TMX ~ interaction(BS_bottom$Texture, BS_bottom$Stage))## No HOV

levene.test(BS_bottom$log10.TMX ~ BS_bottom$Texture * BS_bottom$Stage) ##

boxplot(BS_bottom$log10.TMX ~ BS_bottom$Texture * BS_bottom$Stage, ylab = "TMX (ng g^-1)", xlab = "Trt") ##

boxplot(BS_bottom$log10.TMX)

interaction.plot(BS_bottom$Texture, BS_bottom$Stage, BS_bottom$log10.TMX) ## 

ANOVA.bottom <- aov(BS_bottom$log10.TMX ~ BS_bottom$Texture * BS_bottom$Stage)

summary.bottom <- summary(ANOVA.bottom)

summary.bottom ## 

Tukey.bottom <- TukeyHSD(ANOVA.bottom)

Tukey.bottom

multcompLetters4(ANOVA.bottom, Tukey.bottom)



## try natural log

BS_bottom$log.TMX <- log(BS_bottom$TMX)

hist(BS_bottom$log.TMX) ###great!

qqnorm(BS_bottom$log.TMX)

qqline(BS_bottom$log.TMX) ###great!

bartlett.test(BS_bottom$log.TMX ~ interaction(BS_bottom$Texture, BS_bottom$Stage))  ## 

fligner.test(BS_bottom$log.TMX ~ interaction(BS_bottom$Texture, BS_bottom$Stage))## No HOV again

levene.test(BS_bottom$log.TMX ~ BS_bottom$Texture * BS_bottom$Stage) ##


### sqrt?

BS_bottom$sqrt.TMX <- sqrt(BS_bottom$TMX)

hist(BS_bottom$sqrt.TMX) ###not normal

qqnorm(BS_bottom$sqrt.TMX)

qqline(BS_bottom$sqrt.TMX) ###not normal

bartlett.test(BS_bottom$sqrt.TMX ~ interaction(BS_bottom$Texture, BS_bottom$Stage))  ## 

fligner.test(BS_bottom$sqrt.TMX ~ interaction(BS_bottom$Texture, BS_bottom$Stage))## No HOV again

levene.test(BS_bottom$sqrt.TMX ~ BS_bottom$Texture * BS_bottom$Stage) ##


###rank transform?

BS_bottom$Rank.TMX <- rank(BS_bottom$TMX)

hist(BS_bottom$Rank.TMX) ###

qqnorm(BS_bottom$Rank.TMX)

qqline(BS_bottom$Rank.TMX) ###

bartlett.test(BS_bottom$Rank.TMX ~ interaction(BS_bottom$Texture, BS_bottom$Stage))  ## 

fligner.test(BS_bottom$Rank.TMX ~ interaction(BS_bottom$Texture, BS_bottom$Stage))## meets HOV

levene.test(BS_bottom$Rank.TMX ~ BS_bottom$Texture * BS_bottom$Stage) ##

boxplot(BS_bottom$Rank.TMX ~ BS_bottom$Texture * BS_bottom$Stage, ylab = "TMX (ng g^-1)", xlab = "Trt") ##

boxplot(BS_bottom$Rank.TMX)

interaction.plot(BS_bottom$Texture, BS_bottom$Stage, BS_bottom$Rank.TMX) ## 

ANOVA.bottom1 <- aov(BS_bottom$Rank.TMX ~ BS_bottom$Texture * BS_bottom$Stage)

summary.bottom1 <- summary(ANOVA.bottom1)

summary.bottom1 ## Stage significant--- interaction signficant

Tukey.bottom1 <- TukeyHSD(ANOVA.bottom1)

Tukey.bottom1

letters.bottom <- multcompLetters4(ANOVA.bottom1, Tukey.bottom1) #####Clay:V5 Clay:V5.Ctr     Sand:V1     Sand:V5 Sand:V5.Ctr     Sand:V3     Clay:V1     Clay:V3 
#####################################################"a"         "a"        "ab"        "ab"        "ab"         "b"         "b"         "b" 

####have to use rank transformation due to no HOV

###############End Stats on TMX conc***********************************************














##############TMX Mass in BS***************************************************************************************

BS_TMX.mass <- read.csv("C:/Users/Jesse/Desktop/Neonicotinoids_15/Data/Column-Study/Raw-Data-Outputs/Bulk_soil_mass_TMX_R.csv")

attach(BS_TMX.mass)

BS_TMX.mass

##change ng to micrograms

BS_TMX.mass$TMX.micrg <- BS_TMX.mass$TMX.ng/1000

BS_TMX.mass

BS_TMX.mass.top <- subset(BS_TMX.mass, Interval == "0-30 cm")

BS_TMX.mass.top  ## 0-30 cm

BS_TMX.mass.middle <- subset(BS_TMX.mass, Interval == "30-45 cm")

BS_TMX.mass.middle 

##remove outlier "nk-v5-ctr-5" ---had large preferential flow pathway in 30-45 cm 

BS_TMX.mass.middle <- subset(BS_TMX.mass.middle, Sample.ID != "nk-v5-ctr-5")

BS_TMX.mass.middle ## 30-45 cm

BS_TMX.mass.bottom <- subset(BS_TMX.mass, Interval == "45-60 cm")

BS_TMX.mass.bottom ## 45-60 cm




################################## BS_top Stats---2-way interaction + trt combo---"V5-Ctr" for TMX MASS

BS_TMX.mass.top

boxplot(BS_TMX.mass.top$TMX.ng)

hist(BS_TMX.mass.top$TMX.ng) # 

qqnorm(BS_TMX.mass.top$TMX.ng) ## not normal

qqline(BS_TMX.mass.top$TMX.ng)

BS_TMX.mass.top$log.TMX <- log(BS_TMX.mass.top$TMX.ng)

hist(BS_TMX.mass.top$log.TMX) ####perfect !

qqnorm(BS_TMX.mass.top$log.TMX) ## not normal after natural log transform

qqline(BS_TMX.mass.top$log.TMX) ###ok 





BS_TMX.mass.top$log10.TMX <- log10(BS_TMX.mass.top$TMX.ng) 

hist(BS_TMX.mass.top$log10.TMX) ###normal**** after log base 10 transform

bartlett.test(BS_TMX.mass.top$log10.TMX ~ interaction(BS_TMX.mass.top$Texture, BS_TMX.mass.top$Stage))  ## no HOV

fligner.test(BS_TMX.mass.top$log10.TMX ~ interaction(BS_TMX.mass.top$Texture, BS_TMX.mass.top$Stage)) ##  meets HOV 

levene.test(BS_TMX.mass.top$log10.TMX ~ BS_TMX.mass.top$Texture * BS_TMX.mass.top$Stage) ##

boxplot(BS_TMX.mass.top$log10.TMX ~ BS_TMX.mass.top$Texture * BS_TMX.mass.top$Stage,
        
        ylab = expression(TMX ~ (mu * g)), xlab = "Trt") ## 

interaction.plot(BS_TMX.mass.top$Texture,
                 
                 BS_TMX.mass.top$Stage, 
                 
                 BS_TMX.mass.top$log10.TMX) ## looks like texture and stage effect

ANOVA.top.TMX.mass <- aov(BS_TMX.mass.top$log10.TMX ~ BS_TMX.mass.top$Texture * BS_TMX.mass.top$Stage)

summary.top.TMX.mass <- summary(ANOVA.top.TMX.mass)

summary.top.TMX.mass ## Texture and stage significant, no interaction

Tukey.top.TMX.mass <- TukeyHSD(ANOVA.top.TMX.mass)

Tukey.top.TMX.mass

letters.top.TMX.mass <- multcompLetters4(ANOVA.top.TMX.mass,
                                         
                                         Tukey.top.TMX.mass) ##

letters.top.TMX.mass ##clear differences btw texture


########Try .TMX

BS_TMX.mass.top$Rank.TMX <- rank(BS_TMX.mass.top$TMX.micrg)

hist(BS_TMX.mass.top$Rank.TMX) ###

qqnorm(BS_TMX.mass.top$Rank.TMX)

qqline(BS_TMX.mass.top$Rank.TMX) ###

bartlett.test(BS_TMX.mass.top$Rank.TMX ~ interaction(BS_TMX.mass.top$Texture,
                                                        
                                                     BS_TMX.mass.top$Stage))  ## 

fligner.test(BS_TMX.mass.top$Rank.TMX ~ interaction(BS_TMX.mass.top$Texture,
                                                       
                                                    BS_TMX.mass.top$Stage))## meets HOV

levene.test(BS_TMX.mass.top$Rank.TMX ~ BS_TMX.mass.top$Texture * BS_TMX.mass.top$Stage) ##

boxplot(BS_TMX.mass.top$Rank.TMX ~ BS_TMX.mass.top$Texture * BS_TMX.mass.top$Stage,
        
        ylab = expression(TMX ~ (mu * g)), xlab = "Trt") ##

boxplot(BS_TMX.mass.top$Rank.TMX)

interaction.plot(BS_TMX.mass.top$Texture, BS_TMX.mass.top$Stage, BS_TMX.mass.top$Rank.TMX) ## 

ANOVA.top.TMX.mass1 <- aov(BS_TMX.mass.top$Rank.TMX ~ BS_TMX.mass.top$Texture * BS_TMX.mass.top$Stage)

summary.top.TMX.mass1 <- summary(ANOVA.top.TMX.mass1)

summary.top.TMX.mass1 ## Stage significant--- interaction signficant

Tukey.top.TMX.mass1 <- TukeyHSD(ANOVA.bottom.TMX.mass)

Tukey.top.TMX.mass1

letters.top.TMX.mass1 <- multcompLetters4(ANOVA.top.TMX.mass, Tukey.top.TMX.mass)

letters.top.TMX.mass1


####rank transform not needed








######################## Middle stats TMX MASS---2-way interaction + trt combo---"V5-Ctr" 

BS_TMX.mass.middle

BS_TMX.mass.middle$log10.TMX <- log10(BS_TMX.mass.middle$TMX.micrg)

BS_TMX.mass.middle


hist(BS_TMX.mass.middle$log10.TMX)

qqnorm(BS_TMX.mass.middle$log10.TMX) ### ok

qqline(BS_TMX.mass.middle$log10.TMX) ##ok

bartlett.test(BS_TMX.mass.middle$log10.TMX ~ interaction(BS_TMX.mass.middle$Texture,
                                                         
                                                         BS_TMX.mass.middle$Stage))  ## 


fligner.test(BS_TMX.mass.middle$log10.TMX ~ interaction(BS_TMX.mass.middle$Texture,
                                                        
                                                        BS_TMX.mass.middle$Stage))## meets HOV

levene.test(BS_TMX.mass.middle$log10.TMX ~ BS_TMX.mass.middle$Texture * BS_TMX.mass.middle$Stage) ## 

boxplot(BS_TMX.mass.middle$log10.TMX ~ BS_TMX.mass.middle$Texture * BS_TMX.mass.middle$Stage,
        
        ylab = expression(TMX ~ (mu * g)), xlab = "Trt") ##

interaction.plot(BS_TMX.mass.middle$Texture, BS_TMX.mass.middle$Stage, BS_TMX.mass.middle$log10.TMX) ##  no crossing

ANOVA.middle.TMX.mass <- aov(BS_TMX.mass.middle$log10.TMX ~ BS_TMX.mass.middle$Texture * BS_TMX.mass.middle$Stage)

summary.middle.TMX.mass <- summary(ANOVA.middle.TMX.mass)

summary.middle.TMX.mass ## Stage and interaction effect sifnificant, not texture

Tukey.middle.TMX.mass <- TukeyHSD(ANOVA.middle.TMX.mass)

Tukey.middle.TMX.mass

letters.middle.TMX.mass <- multcompLetters4(ANOVA.middle.TMX.mass, Tukey.middle.TMX.mass)###### Good shit!! differences plant effect in both textures! 

letters.middle.TMX.mass 

### plant effect in sand, not sifnificant in clay

## Sand:V5     Clay:V5 Clay:V5.Ctr Sand:V5.Ctr     Sand:V1     Clay:V1     Clay:V3     Sand:V3 
######"a"        "ab"        "bc"        "cd"        "de"        "ef"        "fg"         "g" 


######################## Bottom stats----- TMX MASS---2-way interaction + trt combo---"V5-Ctr" 

BS_TMX.mass.bottom  ## something wrong.... too much mass?  No checks out fine



hist(BS_TMX.mass.bottom$TMX.micrg)

BS_TMX.mass.bottom


BS_TMX.mass.bottom$log10.TMX <- log10(BS_TMX.mass.bottom$TMX.micrg)

hist(BS_TMX.mass.bottom$log10.TMX) ###great!

qqnorm(BS_TMX.mass.bottom$log10.TMX)

qqline(BS_TMX.mass.bottom$log10.TMX) ###great!

bartlett.test(BS_TMX.mass.bottom$log10.TMX ~ interaction(BS_TMX.mass.bottom$Texture,
                                                         
                                                         BS_TMX.mass.bottom$Stage))  ## 

fligner.test(BS_TMX.mass.bottom$log10.TMX ~ interaction(BS_TMX.mass.bottom$Texture,
                                                       
                                                       BS_TMX.mass.bottom$Stage))## No HOV

levene.test(BS_TMX.mass.bottom$log10.TMX ~ BS_TMX.mass.bottom$Texture * BS_TMX.mass.bottom$Stage)

## damn, need to rank transform just as in BS TMX- conc

boxplot(BS_TMX.mass.bottom$log10.TMX ~ BS_TMX.mass.bottom$Texture * BS_TMX.mass.bottom$Stage,
        
        ylab = expression(TMX ~ (mu * g)), xlab = "Trt") ##



#boxplot(BS_bottom$log10.TMX)

#interaction.plot(BS_bottom$Texture, BS_bottom$Stage, BS_bottom$log10.TMX) ## 

#ANOVA.bottom <- aov(BS_bottom$log10.TMX ~ BS_bottom$Texture * BS_bottom$Stage)

#summary.bottom <- summary(ANOVA.bottom)

#summary.bottom ## 

#Tukey.bottom <- TukeyHSD(ANOVA.bottom)

#Tukey.bottom

#multcompLetters4(ANOVA.bottom, Tukey.bottom)

##############Rank transforming

BS_TMX.mass.bottom$Rank.TMX <- rank(BS_TMX.mass.bottom$TMX.micrg)

hist(BS_TMX.mass.bottom$Rank.TMX) ###

qqnorm(BS_TMX.mass.bottom$Rank.TMX)

qqline(BS_TMX.mass.bottom$Rank.TMX) ###

bartlett.test(BS_TMX.mass.bottom$Rank.TMX ~ interaction(BS_TMX.mass.bottom$Texture,
                                                        
                                                        BS_TMX.mass.bottom$Stage))  ## 

fligner.test(BS_TMX.mass.bottom$Rank.TMX ~ interaction(BS_TMX.mass.bottom$Texture,
                                                       
                                                       BS_TMX.mass.bottom$Stage))## meets HOV

levene.test(BS_TMX.mass.bottom$Rank.TMX ~ BS_TMX.mass.bottom$Texture * BS_TMX.mass.bottom$Stage) ##

boxplot(BS_TMX.mass.bottom$Rank.TMX ~ BS_TMX.mass.bottom$Texture * BS_TMX.mass.bottom$Stage,
        
        ylab = expression(TMX ~ (mu * g)), xlab = "Trt") ##

boxplot(BS_TMX.mass.bottom$Rank.TMX)

interaction.plot(BS_bottom$Texture, BS_bottom$Stage, BS_TMX.mass.bottom$Rank.TMX) ## 

ANOVA.bottom.TMX.mass <- aov(BS_TMX.mass.bottom$Rank.TMX ~ BS_TMX.mass.bottom$Texture * BS_TMX.mass.bottom$Stage)

summary.bottom.TMX.mass <- summary(ANOVA.bottom.TMX.mass)

summary.bottom.TMX.mass ## Stage significant--- interaction signficant

Tukey.bottom.TMX.mass <- TukeyHSD(ANOVA.bottom.TMX.mass)

Tukey.bottom.TMX.mass

letters.bottom.TMX.mass <- multcompLetters4(ANOVA.bottom.TMX.mass, Tukey.bottom.TMX.mass)

letters.bottom.TMX.mass

##Clay:V5  Clay:V5.Ctr     Sandy:V1 Sandy:V5.Ctr     Sandy:V5      Clay:V1     Sandy:V3      Clay:V3 
####"a"          "a"         "ab"         "ab"         "ab"          "b"          "b"          "b" 



##########################***********************************************************************************



###################prepping TMX conc for plots

###make 3 bar graphs

library(dplyr)
library(ggplot2)
library(Rmisc)
library(gridExtra)
library(grid)
library(lattice)

###Data top

BS_top

BS_top.bar <- summarySE(BS_top, measurevar="TMX", groupvars=c("Texture","Stage"))

BS_top.bar



letters.top

BS_top.bar$letters <- c("abc", "cd", "d", "bcd", "a", "bc", "ab", "abc") #### add letters vector to df

BS_top.bar

###Data Middle

BS_middle.bar <- summarySE(BS_middle1, measurevar="TMX", groupvars=c("Texture","Stage"))

BS_middle.bar

letters.middle

BS_middle.bar$letters <-  c("cd", "d", "a", "b", "c", "d", "a", "b")  ###### add letters vector to df

BS_middle.bar


###Data Bottom
BS_bottom.bar <- summarySE(BS_bottom, measurevar="TMX", groupvars=c("Texture","Stage"))

BS_bottom.bar

letters.bottom

BS_bottom.bar$letters <- c("b", "b", "a", "a", "ab", "b", "ab", "ab")

BS_bottom.bar


##################### prepping BS TMX Mass for plots

###Data top

BS_TMX.mass.top

BS_TMX.mass.top.bar <- summarySE(BS_TMX.mass.top, measurevar="TMX.micrg", groupvars=c("Texture","Stage"))

BS_TMX.mass.top.bar



letters.top.TMX.mass

##Sand:V1     Clay:V1     Sand:V5     Sand:V3 Sand:V5.Ctr Clay:V5.Ctr     Clay:V3     Clay:V5 
##"a"        "ab"        "ab"        "ab"         "b"         "b"         "b"         "b" 



BS_TMX.mass.top.bar$letters <- c("ab", "b", "b", "b", "a", "ab", "ab", "b") #### add letters vector to df

BS_TMX.mass.top.bar

###Data Middle

BS_TMX.mass.middle

BS_TMX.mass.middle.bar <- summarySE(BS_TMX.mass.middle, measurevar="TMX.micrg", groupvars=c("Texture","Stage"))

BS_TMX.mass.middle.bar

letters.middle.TMX.mass

##Sand:V5     Clay:V5 Clay:V5.Ctr Sand:V5.Ctr     Sand:V1     Clay:V1     Clay:V3     Sand:V3 
####"a"        "ab"        "bc"        "cd"        "de"        "ef"        "fg"         "g" 
##
BS_TMX.mass.middle.bar$letters <-  c("ef", "fg", "ab", "bc", "de", "g", "a", "cd")  ###### add letters vector to df

BS_TMX.mass.middle.bar


###Data Bottom

BS_TMX.mass.bottom

BS_TMX.mass.bottom.bar <- summarySE(BS_TMX.mass.bottom, measurevar="TMX.micrg", groupvars=c("Texture","Stage"))

BS_TMX.mass.bottom.bar

letters.bottom.TMX.mass

##Clay:V5  Clay:V5.Ctr     Sandy:V1 Sandy:V5.Ctr     Sandy:V5      Clay:V1     Sandy:V3      Clay:V3 
####"a"          "a"         "ab"         "ab"         "ab"          "b"          "b"          "b" 

BS_TMX.mass.bottom.bar$letters <- c("b", "b", "a", "a", "ab", "b", "ab", "ab")

BS_TMX.mass.bottom.bar


######################################*********************************************************************
















##position_dodge--prevents bar overlap

###clean up function-cleans up ggplot defaults

###cleanup <- theme(panel.grid.major = element_blank(),
#                   panel.grid.minor = element_blank(),
#                   panel.background =  element_blank(),
#                   axis.line = element_line(color - "black"))

#scale_x_discrete() usesd tio change continuous variables on x
#fill controls the color of bar, color controls the outline color of bar

#postion_dodge() or psotion = "dodge" prevents a stacked bar and arranges as side-by-side columns

##to change legend and color specs:  scale_fill_manual( name = "name of legend"
#                                                     labels = 
#                                                      values = "specify color here" )










cleanup <- theme(panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(),
                                  panel.background =  element_blank(),
                                  axis.line = element_line(color = "black"))
 


########making plots             



bar.top <- ggplot(BS_top.bar, aes(Texture, TMX, fill = Stage)) +
  
         geom_bar(position = position_dodge(), stat = "identity") +

         geom_bar(position = "dodge", stat = "identity", color="black", show.legend =FALSE, lwd = 0.75) + ##add outline to bars
           
         geom_errorbar(aes(ymin = TMX, ymax = TMX + se),      ##error bars
         
         width = 0.2, position = position_dodge(0.9), lwd = 0.75) +
  
         cleanup +                                          #cleanup
  
         scale_fill_manual(name = "",
                           
                           labels = c("V1", "V3", "V5", "V5, No Plant"), 
                                     
                           values = c("blue", "red", "green", "yellow")) +
                                                                                                  
         xlab("") +
        
         ylab(expression(TMX ~ (ng ~ g^{-1}))) + ### pain in the ass way to express superscripts

         theme(axis.text.y = element_text(size = 16),  ## change font size of y axis
               
         axis.title = element_text(size = 16, face = "bold")) +
  
         #geom_label(aes(label = letters)) +

         #geom_text(aes(label = letters)) +
      
         geom_text(aes(label = letters, y = TMX + se, fontface = "italic"), position = position_dodge(width=0.9), vjust= -0.5) +  #####add in signficance letter--adjusted for errror bars
  
         theme(axis.text.x = element_text(size = 16),  ## change font size of y axis
               
               axis.title = element_text(size = 16, face = "bold")) +
  
         theme(axis.line.x = element_line(size = 1), ##change axis line size
               
               axis.line.y = element_line(size = 1),
               
               panel.border = element_rect(color = "black", fill = NA, size = 1.75), ##add border
               
               legend.position = c(0.16,0.87), ## move legend
               
               legend.margin = margin(0, 0, 0, 0, "mm"), ##remove margin on legend
               
               legend.text = element_text( size = 16), ##legend text
               
               #legend.background = element_rect(color = "black", size = 1), ## add legend border
               
               #legend.key.height=unit(1,"line"), 

               #legend.key.width=unit(1,"line")
               
               legend.key = element_rect(color = "black", size = 0.75), ## put line around legend key
               
               legend.direction = "vertical", 

               axis.ticks.y = element_line(size = 1, color = "black"),  ##enlargen ticks--problem at op right corner
                
               legend.key.size = unit(7, 'mm'),
               
               legend.key.height = unit(6, "mm")) +
               
               #plot.title = element_text( Size = 18)) +
  
               #guides(guide_legend(nrow=4)) +

               #legend.spacing.y = unit(1, "mm") +
               
               #legend.key.size = unit(1, "line")) 
  
          #scale_fill_manual(values=values, labels=setNames(paste(labels, " "), entries)) +
      
          ggtitle("Bulk Soil, 0-30 cm") +
 
          theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5, vjust = 1)) + ## added title and adjusted
  
          scale_y_continuous(expand = c(0,0), limit = c(0, 100)) + scale_x_discrete(expand = c(0.03,0)) ### removed excess space btw bars and axis


bar.middle <- ggplot(BS_middle.bar, aes(Texture, TMX, fill = Stage)) +
  
  geom_bar(position = position_dodge(), stat = "identity") +
  
  geom_bar(position = "dodge", stat = "identity", color="black", show.legend =FALSE, lwd = 0.75) + ##add outline to bars
  
  geom_errorbar(aes(ymin = TMX, ymax = TMX + se),      ##error bars
                
                width = 0.2, position = position_dodge(0.9), lwd = 0.75) +
  
  cleanup +                                          #cleanup
  
  scale_fill_manual(name = "",
                    
                    labels = c("V1", "V3", "V5", "V5, No Plant"), 
                    
                    values = c("blue", "red", "green", "yellow")) +
  
  xlab("") +
  
  ylab(expression(TMX ~ (ng ~ g^{-1}))) + ### pain in the ass way to express superscripts
  
  theme(axis.text.y = element_text(size = 16),  ## change font size of y axis
        
        axis.title = element_text(size = 16, face = "bold")) +
  
  #geom_label(aes(label = letters)) +
  
  #geom_text(aes(label = letters)) +
  
  geom_text(aes(label = letters, y = TMX + se, fontface = "italic"), position = position_dodge(width=0.9), vjust= -0.5) +  #####add in signficance letter--adjusted for errror bars
  
  theme(axis.text.x = element_text(size = 16),  ## change font size of y axis
        
        axis.title = element_text(size = 16, face = "bold")) +
  
  theme(axis.line.x = element_line(size = 1), ##change axis line size
        
        axis.line.y = element_line(size = 1),
        
        panel.border = element_rect(color = "black", fill = NA, size = 1.75), ##add border
        
        legend.position = c(0.16,0.87), ## move legend
        
        legend.margin = margin(0, 0, 0, 0, "mm"), ##remove margin on legend
        
        legend.text = element_text( size = 16), ##legend text
        
        #legend.background = element_rect(color = "black", size = 1), ## add legend border
        
        #legend.key.height=unit(1,"line"), 
        
        #legend.key.width=unit(1,"line")
        
        legend.key = element_rect(color = "black", size = 0.75), ## put line around legend key
        
        legend.direction = "vertical", 
        
        axis.ticks.y = element_line(size = 1, color = "black"),  ##enlargen ticks--problem at op right corner
        
        legend.key.size = unit(7, 'mm'),
        
        legend.key.height = unit(6, "mm")) +
  
  #plot.title = element_text( Size = 18)) +
  
  #guides(guide_legend(nrow=4)) +
  
  #legend.spacing.y = unit(1, "mm") +
  
  #legend.key.size = unit(1, "line")) 
  
  #scale_fill_manual(values=values, labels=setNames(paste(labels, " "), entries)) +
  
ggtitle("Bulk Soil, 30-45 cm") +
  
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5, vjust = 1)) + ## added title and adjusted
  
  scale_y_continuous(expand = c(0,0), limit = c(0, 100)) + scale_x_discrete(expand = c(0.03,0)) ### removed excess space btw bars and axis






bar.bottom <- ggplot(BS_bottom.bar, aes(Texture, TMX, fill = Stage)) +
  
  geom_bar(position = position_dodge(), stat = "identity") +
  
  geom_bar(position = "dodge", stat = "identity", color="black", show.legend =FALSE, lwd = 0.75) + ##add outline to bars
  
  geom_errorbar(aes(ymin = TMX, ymax = TMX + se),      ##error bars
                
                width = 0.2, position = position_dodge(0.9), lwd = 0.75) +
  
  cleanup +                                          #cleanup
  
  scale_fill_manual(name = "",
                    
                    labels = c("V1", "V3", "V5", "V5, No Plant"), 
                    
                    values = c("blue", "red", "green", "yellow")) +
  
  xlab("") +
  
  ylab(expression(TMX ~ (ng ~ g^{-1}))) + ### pain in the ass way to express superscripts
  
  theme(axis.text.y = element_text(size = 16),  ## change font size of y axis
        
        axis.title = element_text(size = 16, face = "bold")) +
  
  #geom_label(aes(label = letters)) +
  
  #geom_text(aes(label = letters)) +
  
  geom_text(aes(label = letters, y = TMX + se, fontface = "italic"), position = position_dodge(width=0.9), vjust= -0.5) +  #####add in signficance letter--adjusted for errror bars
  
  theme(axis.text.x = element_text(size = 16),  ## change font size of y axis
        
        axis.title = element_text(size = 16, face = "bold")) +
  
  theme(axis.line.x = element_line(size = 1), ##change axis line size
        
        axis.line.y = element_line(size = 1),
        
        panel.border = element_rect(color = "black", fill = NA, size = 1.75), ##add border
        
        legend.position = c(0.16,0.87), ## move legend
        
        legend.margin = margin(0, 0, 0, 0, "mm"), ##remove margin on legend
        
        legend.text = element_text( size = 16), ##legend text
        
        #legend.background = element_rect(color = "black", size = 1), ## add legend border
        
        #legend.key.height=unit(1,"line"), 
        
        #legend.key.width=unit(1,"line")
        
        legend.key = element_rect(color = "black", size = 0.75), ## put line around legend key
        
        legend.direction = "vertical", 
        
        axis.ticks.y = element_line(size = 1, color = "black"),  ##enlargen ticks--problem at op right corner
        
        legend.key.size = unit(7, 'mm'),
        
        legend.key.height = unit(6, "mm")) +
  
  #plot.title = element_text( Size = 18)) +
  
  #guides(guide_legend(nrow=4)) +
  
  #legend.spacing.y = unit(1, "mm") +
  
  #legend.key.size = unit(1, "line")) 
  
  #scale_fill_manual(values=values, labels=setNames(paste(labels, " "), entries)) +
  
ggtitle("Bulk Soil, 45-60 cm") +
  
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5, vjust = 1)) + ## added title and adjusted
  
  scale_y_continuous(expand = c(0,0), limit = c(0, 100)) + scale_x_discrete(expand = c(0.03,0)) ### removed excess space btw bars and axis



bar.top

bar.middle

bar.bottom





?citation()






###################Make TMX mass plots*******************************************************

bar.top.mass <- ggplot(BS_TMX.mass.top.bar, aes(Texture, TMX.micrg, fill = Stage)) +
  
  geom_bar(position = position_dodge(), stat = "identity") +
  
  geom_bar(position = "dodge", stat = "identity", color="black", show.legend =FALSE, lwd = 0.75) + ##add outline to bars
  
  geom_errorbar(aes(ymin = TMX.micrg, ymax = TMX.micrg + se),      ##error bars
                
                width = 0.2, position = position_dodge(0.9), lwd = 0.75) +
  
  cleanup +                                          #cleanup
  
  scale_fill_manual(name = "",
                    
                    labels = c("V1", "V3", "V5", "V5, No Plant"), 
                    
                    values = c("blue", "red", "green", "yellow")) +
  
  xlab("") +
  
  ylab(expression(TMX ~ (mu * g))) + ### pain in the ass way to express superscripts
  
  theme(axis.text.y = element_text(size = 16),  ## change font size of y axis
        
        axis.title = element_text(size = 16, face = "bold")) +
  
  #geom_label(aes(label = letters)) +
  
  #geom_text(aes(label = letters)) +
  
  geom_text(aes(label = letters, y = TMX.micrg + se, fontface = "italic"), position = position_dodge(width=0.9), vjust= -0.5) +  #####add in signficance letter--adjusted for errror bars
  
  theme(axis.text.x = element_text(size = 16),  ## change font size of y axis
        
        axis.title = element_text(size = 16, face = "bold")) +
  
  theme(axis.line.x = element_line(size = 1), ##change axis line size
        
        axis.line.y = element_line(size = 1),
        
        panel.border = element_rect(color = "black", fill = NA, size = 1.75), ##add border
        
        legend.position = c(0.16,0.87), ## move legend
        
        legend.margin = margin(0, 0, 0, 0, "mm"), ##remove margin on legend
        
        legend.text = element_text( size = 16), ##legend text
        
        #legend.background = element_rect(color = "black", size = 1), ## add legend border
        
        #legend.key.height=unit(1,"line"), 
        
        #legend.key.width=unit(1,"line")
        
        legend.key = element_rect(color = "black", size = 0.75), ## put line around legend key
        
        legend.direction = "vertical", 
        
        axis.ticks.y = element_line(size = 1, color = "black"),  ##enlargen ticks--problem at op right corner
        
        legend.key.size = unit(7, 'mm'),
        
        legend.key.height = unit(6, "mm")) +
  
  #plot.title = element_text( Size = 18)) +
  
  #guides(guide_legend(nrow=4)) +
  
  #legend.spacing.y = unit(1, "mm") +
  
  #legend.key.size = unit(1, "line")) 
  
  #scale_fill_manual(values=values, labels=setNames(paste(labels, " "), entries)) +
  
ggtitle("Bulk Soil, 0-30 cm") +
  
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5, vjust = 1)) + ## added title and adjusted
  
  scale_y_continuous(expand = c(0,0), limit = c(0, 1550)) + scale_x_discrete(expand = c(0.03,0)) ### removed excess space btw bars and axis

bar.top.mass



bar.middle.mass <- ggplot(BS_TMX.mass.middle.bar, aes(Texture, TMX.micrg, fill = Stage)) +
  
  geom_bar(position = position_dodge(), stat = "identity") +
  
  geom_bar(position = "dodge", stat = "identity", color="black", show.legend =FALSE, lwd = 0.75) + ##add outline to bars
  
  geom_errorbar(aes(ymin = TMX.micrg, ymax = TMX.micrg + se),      ##error bars
                
                width = 0.2, position = position_dodge(0.9), lwd = 0.75) +
  
  cleanup +                                          #cleanup
  
  scale_fill_manual(name = "",
                    
                    labels = c("V1", "V3", "V5", "V5, No Plant"), 
                    
                    values = c("blue", "red", "green", "yellow")) +
  
  xlab("") +
  
  ylab(expression(TMX ~ (mu * g))) + ### pain in the ass way to express superscripts
  
  theme(axis.text.y = element_text(size = 16),  ## change font size of y axis
        
        axis.title = element_text(size = 16, face = "bold")) +
  
  #geom_label(aes(label = letters)) +
  
  #geom_text(aes(label = letters)) +
  
  geom_text(aes(label = letters, y = TMX.micrg + se, fontface = "italic"), position = position_dodge(width=0.9), vjust= -0.5) +  #####add in signficance letter--adjusted for errror bars
  
  theme(axis.text.x = element_text(size = 16),  ## change font size of y axis
        
        axis.title = element_text(size = 16, face = "bold")) +
  
  theme(axis.line.x = element_line(size = 1), ##change axis line size
        
        axis.line.y = element_line(size = 1),
        
        panel.border = element_rect(color = "black", fill = NA, size = 1.75), ##add border
        
        legend.position = c(0.16,0.87), ## move legend
        
        legend.margin = margin(0, 0, 0, 0, "mm"), ##remove margin on legend
        
        legend.text = element_text( size = 16), ##legend text
        
        #legend.background = element_rect(color = "black", size = 1), ## add legend border
        
        #legend.key.height=unit(1,"line"), 
        
        #legend.key.width=unit(1,"line")
        
        legend.key = element_rect(color = "black", size = 0.75), ## put line around legend key
        
        legend.direction = "vertical", 
        
        axis.ticks.y = element_line(size = 1, color = "black"),  ##enlargen ticks--problem at op right corner
        
        legend.key.size = unit(7, 'mm'),
        
        legend.key.height = unit(6, "mm")) +
  
  #plot.title = element_text( Size = 18)) +
  
  #guides(guide_legend(nrow=4)) +
  
  #legend.spacing.y = unit(1, "mm") +
  
  #legend.key.size = unit(1, "line")) 
  
  #scale_fill_manual(values=values, labels=setNames(paste(labels, " "), entries)) +
  
ggtitle("Bulk Soil, 30-45 cm") +
  
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5, vjust = 1)) + ## added title and adjusted
  
  scale_y_continuous(expand = c(0,0), limit = c(0, 1500)) + scale_x_discrete(expand = c(0.03,0)) ### removed excess space btw bars and axis

bar.middle.mass ####plant effect only sig in sand




bar.bottom.mass <- ggplot(BS_TMX.mass.bottom.bar, aes(Texture, TMX.micrg, fill = Stage)) +
  
  geom_bar(position = position_dodge(), stat = "identity") +
  
  geom_bar(position = "dodge", stat = "identity", color="black", show.legend =FALSE, lwd = 0.75) + ##add outline to bars
  
  geom_errorbar(aes(ymin = TMX.micrg, ymax = TMX.micrg + se),      ##error bars
                
                width = 0.2, position = position_dodge(0.9), lwd = 0.75) +
  
  cleanup +                                          #cleanup
  
  scale_fill_manual(name = "",
                    
                    labels = c("V1", "V3", "V5", "V5, No Plant"), 
                    
                    values = c("blue", "red", "green", "yellow")) +
  
  xlab("") +
  
  ylab(expression(TMX ~ (mu * g))) + ### pain in the ass way to express superscripts
  
  theme(axis.text.y = element_text(size = 16),  ## change font size of y axis
        
        axis.title = element_text(size = 16, face = "bold")) +
  
  #geom_label(aes(label = letters)) +
  
  #geom_text(aes(label = letters)) +
  
  geom_text(aes(label = letters, y = TMX.micrg + se, fontface = "italic"), position = position_dodge(width=0.9), vjust= -0.5) +  #####add in signficance letter--adjusted for errror bars
  
  theme(axis.text.x = element_text(size = 16),  ## change font size of y axis
        
        axis.title = element_text(size = 16, face = "bold")) +
  
  theme(axis.line.x = element_line(size = 1), ##change axis line size
        
        axis.line.y = element_line(size = 1),
        
        panel.border = element_rect(color = "black", fill = NA, size = 1.75), ##add border
        
        legend.position = c(0.16,0.87), ## move legend
        
        legend.margin = margin(0, 0, 0, 0, "mm"), ##remove margin on legend
        
        legend.text = element_text( size = 16), ##legend text
        
        #legend.background = element_rect(color = "black", size = 1), ## add legend border
        
        #legend.key.height=unit(1,"line"), 
        
        #legend.key.width=unit(1,"line")
        
        legend.key = element_rect(color = "black", size = 0.75), ## put line around legend key
        
        legend.direction = "vertical", 
        
        axis.ticks.y = element_line(size = 1, color = "black"),  ##enlargen ticks--problem at op right corner
        
        legend.key.size = unit(7, 'mm'),
        
        legend.key.height = unit(6, "mm")) +
  
  #plot.title = element_text( Size = 18)) +
  
  #guides(guide_legend(nrow=4)) +
  
  #legend.spacing.y = unit(1, "mm") +
  
  #legend.key.size = unit(1, "line")) 
  
  #scale_fill_manual(values=values, labels=setNames(paste(labels, " "), entries)) +
  
ggtitle("Bulk Soil, 45-60 cm") +
  
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5, vjust = 1)) + ## added title and adjusted
  
  scale_y_continuous(expand = c(0,0), limit = c(0, 1500)) + scale_x_discrete(expand = c(0.03,0)) ### removed excess space btw bars and axis



bar.bottom.mass



####Try to arrange everything on one graph? ******************************************

#####make seprate objects for combination/comparison

bar.top1 <- ggplot(BS_top.bar, aes(Texture, TMX, fill = Stage)) +
  
  geom_bar(position = position_dodge(), stat = "identity") +
  
  geom_bar(position = "dodge", stat = "identity", color="black", show.legend =FALSE, lwd = 0.75) + ##add outline to bars
  
  geom_errorbar(aes(ymin = TMX, ymax = TMX + se),      ##error bars
                
                width = 0.2, position = position_dodge(0.9), lwd = 0.75) +
  
  cleanup +                                          #cleanup
  
  scale_fill_manual(name = "",
                    
                    labels = c("V1", "V3", "V5", "V5, No Plant"), 
                    
                    values = c("blue", "red", "green", "yellow")) +
  
  xlab("") +
  
  ylab(expression(TMX ~ (ng ~ g^{-1}))) + ### pain in the ass way to express superscripts
  
  theme(axis.text.y = element_text(size = 16),  ## change font size of y axis
        
        axis.title = element_text(size = 16, face = "bold")) +
  
  #geom_label(aes(label = letters)) +
  
  #geom_text(aes(label = letters)) +
  
  geom_text(aes(label = letters, y = TMX + se, fontface = "italic"), position = position_dodge(width=0.9), vjust= -0.5) +  #####add in signficance letter--adjusted for errror bars
  
  theme(axis.text.x = element_text(size = 16),  ## change font size of y axis
        
        axis.title = element_text(size = 16, face = "bold")) +
  
  theme(axis.line.x = element_line(size = 1), ##change axis line size
        
        axis.line.y = element_line(size = 1),
        
        panel.border = element_rect(color = "black", fill = NA, size = 1.75), ##add border
        
        legend.position = "none", ## move legend
        
        legend.margin = margin(0, 0, 0, 0, "mm"), ##remove margin on legend
        
        legend.text = element_text( size = 16), ##legend text
        
        #legend.background = element_rect(color = "black", size = 1), ## add legend border
        
        #legend.key.height=unit(1,"line"), 
        
        #legend.key.width=unit(1,"line")
        
        legend.key = element_rect(color = "black", size = 0.75), ## put line around legend key
        
        legend.direction = "vertical", 
        
        axis.ticks.y = element_line(size = 1, color = "black"),  ##enlargen ticks--problem at op right corner
        
        legend.key.size = unit(7, 'mm'),
        
        legend.key.height = unit(6, "mm")) +
  
  #plot.title = element_text( Size = 18)) +
  
  #guides(guide_legend(nrow=4)) +
  
  #legend.spacing.y = unit(1, "mm") +
  
  #legend.key.size = unit(1, "line")) 
  
  #scale_fill_manual(values=values, labels=setNames(paste(labels, " "), entries)) +
  
ggtitle("Bulk Soil, 0-30 cm") +
  
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5, vjust = 1),
        
        plot.margin = unit( c(2, 0, 0, 2), "mm")) + ## added title and adjusted
  
  scale_y_continuous(expand = c(0,0), limit = c(0, 110)) + scale_x_discrete(expand = c(0.03,0)) ### removed excess space btw bars and axis





bar.middle1 <- ggplot(BS_middle.bar, aes(Texture, TMX, fill = Stage)) +
  
  geom_bar(position = position_dodge(), stat = "identity") +
  
  geom_bar(position = "dodge", stat = "identity", color="black", show.legend =FALSE, lwd = 0.75) + ##add outline to bars
  
  geom_errorbar(aes(ymin = TMX, ymax = TMX + se),      ##error bars
                
                width = 0.2, position = position_dodge(0.9), lwd = 0.75) +
  
  cleanup +                                          #cleanup
  
  scale_fill_manual(name = "",
                    
                    labels = c("V1", "V3", "V5", "V5, No Plant"), 
                    
                    values = c("blue", "red", "green", "yellow")) +
  
  xlab("") +
  
  ylab("") + ### pain in the ass way to express superscripts
  
  theme(axis.text.y = element_text(size = 16),  ## change font size of y axis
        
        axis.title = element_text(size = 16, face = "bold")) +
  
  #geom_label(aes(label = letters)) +
  
  #geom_text(aes(label = letters)) +
  
  geom_text(aes(label = letters, y = TMX + se, fontface = "italic"), position = position_dodge(width=0.9), hjust = -0.5) +  #####add in signficance letter--adjusted for errror bars
  
  theme(axis.text.x = element_text(size = 16),  ## change font size of y axis
        
        axis.title = element_text(size = 16, face = "bold")) +
  
  theme(axis.line.x = element_line(size = 1), ##change axis line size
        
        axis.line.y = element_line(size = 1),
        
        panel.border = element_rect(color = "black", fill = NA, size = 1.75), ##add border
        
        legend.position = c(0.37,0.84), ## move legend
        
        legend.margin = margin(0, 0, 0, 0, "mm"), ##remove margin on legend
        
        legend.text = element_text( size = 16), ##legend text
        
        #legend.background = element_rect(color = "black", size = 1), ## add legend border
        
        #legend.key.height=unit(1,"line"), 
        
        #legend.key.width=unit(1,"line")
        
        legend.key = element_rect(color = "black", size = 0.75), ## put line around legend key
        
        legend.direction = "vertical", 
        
        axis.ticks.y = element_line(size = 1, color = "black"),  ##enlargen ticks--problem at op right corner
        
        legend.key.size = unit(8, 'mm'),
        
        legend.key.height = unit(8, "mm")) +
  
  #plot.title = element_text( Size = 18)) +
  
  #guides(guide_legend(nrow=4)) +
  
  #legend.spacing.y = unit(1, "mm") +
  
  #legend.key.size = unit(1, "line")) 
  
  #scale_fill_manual(values=values, labels=setNames(paste(labels, " "), entries)) +
  
ggtitle("Bulk Soil, 30-45 cm") +
  
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5, vjust = 1),
        
        plot.margin = unit( c(2, 2, 0, 0), "mm")) + ## added title and adjusted
  
  scale_y_continuous(expand = c(0,0), limit = c(0, 100)) + scale_x_discrete(expand = c(0.03,0)) +

 coord_flip() ### removed excess space btw bars and axis + swapped axes

bar.middle1




bar.bottom1 <- ggplot(BS_bottom.bar, aes(Texture, TMX, fill = Stage)) +
  
  geom_bar(position = position_dodge(), stat = "identity") +
  
  geom_bar(position = "dodge", stat = "identity", color="black", show.legend =FALSE, lwd = 0.75) + ##add outline to bars
  
  geom_errorbar(aes(ymin = TMX, ymax = TMX + se),      ##error bars
                
                width = 0.2, position = position_dodge(0.9), lwd = 0.75) +
  
  cleanup +                                          #cleanup
  
  scale_fill_manual(name = "",
                    
                    labels = c("V1", "V3", "V5", "V5, No Plant"), 
                    
                    values = c("blue", "red", "green", "yellow")) +
  
  xlab("") +
  
  ylab("") + ### pain in the ass way to express superscripts
  
  theme(axis.text.y = element_text(size = 16),  ## change font size of y axis
        
        axis.title = element_text(size = 16, face = "bold")) +
  
  #geom_label(aes(label = letters)) +
  
  #geom_text(aes(label = letters)) +
  
  geom_text(aes(label = letters, y = TMX + se, fontface = "italic"), position = position_dodge(width=0.9), hjust= -0.5) +  #####add in signficance letter--adjusted for errror bars
  
  theme(axis.text.x = element_text(size = 16),  ## change font size of y axis
        
        axis.title = element_text(size = 16, face = "bold")) +
  
  theme(axis.line.x = element_line(size = 1), ##change axis line size
        
        axis.line.y = element_line(size = 1),
        
        panel.border = element_rect(color = "black", fill = NA, size = 1.75), ##add border
        
        legend.position = "none", ## move legend
        
        legend.margin = margin(0, 0, 0, 0, "mm"), ##remove margin on legend
        
        legend.text = element_text( size = 16), ##legend text
        
        #legend.background = element_rect(color = "black", size = 1), ## add legend border
        
        #legend.key.height=unit(1,"line"), 
        
        #legend.key.width=unit(1,"line")
        
        legend.key = element_rect(color = "black", size = 0.75), ## put line around legend key
        
        legend.direction = "horizontal", 
        
        axis.ticks.y = element_line(size = 1, color = "black"),  ##enlargen ticks--problem at op right corner
        
        legend.key.size = unit(7, 'mm'),
        
        legend.key.height = unit(6, "mm")) +
  
  #plot.title = element_text( Size = 18)) +
  
  #guides(guide_legend(nrow=4)) +
  
  #legend.spacing.y = unit(1, "mm") +
  
  #legend.key.size = unit(1, "line")) 
  
  #scale_fill_manual(values=values, labels=setNames(paste(labels, " "), entries)) +
  
ggtitle("Bulk Soil, 45-60 cm") +
  
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5, vjust = 1),
        
        plot.margin = unit( c(2, 2, 0, 0), "mm")) + ## added title and adjusted
  
  scale_y_continuous(expand = c(0,0), limit = c(0, 100)) + scale_x_discrete(expand = c(0.03,0)) +
  
  coord_flip() ### removed excess space btw bars and axis and swapped axes

bar.bottom1


#################Mass plots***************************************************************************

bar.top.mass1 <- ggplot(BS_TMX.mass.top.bar, aes(Texture, TMX.micrg, fill = Stage)) +
  
  geom_bar(position = position_dodge(), stat = "identity") +
  
  geom_bar(position = "dodge", stat = "identity", color="black", show.legend =FALSE, lwd = 0.75) + ##add outline to bars
  
  geom_errorbar(aes(ymin = TMX.micrg, ymax = TMX.micrg + se),      ##error bars
                
                width = 0.2, position = position_dodge(0.9), lwd = 0.75) +
  
  cleanup +                                          #cleanup
  
  scale_fill_manual(name = "",
                    
                    labels = c("V1", "V3", "V5", "V5, No Plant"), 
                    
                    values = c("blue", "red", "green", "yellow")) +
  
  xlab("") +
  
  ylab(expression(TMX ~ (mu * g))) + ### pain in the ass way to express superscripts
  
  theme(axis.text.y = element_text(size = 16),  ## change font size of y axis
        
        axis.title = element_text(size = 16, face = "bold")) +
  
  #geom_label(aes(label = letters)) +
  
  #geom_text(aes(label = letters)) +
  
  geom_text(aes(label = letters, y = TMX.micrg + se, fontface = "italic"), position = position_dodge(width=0.9), vjust= -0.5) +  #####add in signficance letter--adjusted for errror bars
  
  theme(axis.text.x = element_text(size = 16),  ## change font size of y axis
        
        axis.title = element_text(size = 16, face = "bold"),
        
        axis.title.y = element_text(vjust = -0.5)) +
  
  theme(axis.line.x = element_line(size = 1), ##change axis line size
        
        axis.line.y = element_line(size = 1),
        
        panel.border = element_rect(color = "black", fill = NA, size = 1.75), ##add border
        
        legend.position = "none", ## move legend
        
        legend.margin = margin(0, 0, 0, 0, "mm"), ##remove margin on legend
        
        legend.text = element_text( size = 16), ##legend text
        
        #legend.background = element_rect(color = "black", size = 1), ## add legend border
        
        #legend.key.height=unit(1,"line"), 
        
        #legend.key.width=unit(1,"line")
        
        legend.key = element_rect(color = "black", size = 0.75), ## put line around legend key
        
        legend.direction = "vertical", 
        
        axis.ticks.y = element_line(size = 1, color = "black"),  ##enlargen ticks--problem at op right corner
        
        legend.key.size = unit(7, 'mm'),
        
        legend.key.height = unit(6, "mm"),
        
        plot.margin = unit(c(0,0,0,0.5),"mm")) +
  
  #plot.title = element_text( Size = 18)) +
  
  #guides(guide_legend(nrow=4)) +
  
  #legend.spacing.y = unit(1, "mm") +
  
  #legend.key.size = unit(1, "line")) 
  
  #scale_fill_manual(values=values, labels=setNames(paste(labels, " "), entries)) +
  
#ggtitle("Bulk Soil, 0-30 cm") +
  
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5, vjust = 1)) + ## added title and adjusted
  
  scale_y_continuous(expand = c(0,0), limit = c(0, 1600)) + scale_x_discrete(expand = c(0.03,0)) ### removed excess space btw bars and axis

bar.top.mass1




bar.middle.mass1 <- ggplot(BS_TMX.mass.middle.bar, aes(Texture, TMX.micrg, fill = Stage)) +
  
  geom_bar(position = position_dodge(), stat = "identity") +
  
  geom_bar(position = "dodge", stat = "identity", color="black", show.legend =FALSE, lwd = 0.75) + ##add outline to bars
  
  geom_errorbar(aes(ymin = TMX.micrg, ymax = TMX.micrg + se),      ##error bars
                
                width = 0.2, position = position_dodge(0.9), lwd = 0.75) +
  
  cleanup +                                          #cleanup
  
  scale_fill_manual(name = "",
                    
                    labels = c("V1", "V3", "V5", "V5, No Plant"), 
                    
                    values = c("blue", "red", "green", "yellow")) +
  
  xlab("") +
  
  ylab(expression(TMX ~ (mu * g))) + ### pain in the ass way to express superscripts
  
  theme(axis.text.y = element_text(size = 16),  ## change font size of y axis
        
        axis.title = element_text(size = 16, face = "bold")) +
  
  #geom_label(aes(label = letters)) +
  
  #geom_text(aes(label = letters)) +
  
  geom_text(aes(label = letters, y = TMX.micrg + se, fontface = "italic"), position = position_dodge(width=0.9), hjust= -0.5) +  #####add in signficance letter--adjusted for errror bars
  
  theme(axis.text.x = element_text(size = 16),  ## change font size of y axis
        
        axis.title = element_text(size = 16, face = "bold")) +
  
  theme(axis.line.x = element_line(size = 1), ##change axis line size
        
        axis.line.y = element_line(size = 1),
        
        panel.border = element_rect(color = "black", fill = NA, size = 1.75), ##add border
        
        legend.position = c(0.16,0.87), ## move legend
        
        legend.margin = margin(0, 0, 0, 0, "mm"), ##remove margin on legend
        
        legend.text = element_text( size = 16), ##legend text
        
        #legend.background = element_rect(color = "black", size = 1), ## add legend border
        
        #legend.key.height=unit(1,"line"), 
        
        #legend.key.width=unit(1,"line")
        
        legend.key = element_rect(color = "black", size = 0.75), ## put line around legend key
        
        legend.direction = "vertical", 
        
        axis.ticks.y = element_line(size = 1, color = "black"),  ##enlargen ticks--problem at op right corner
        
        legend.key.size = unit(7, 'mm'),
        
        legend.key.height = unit(6, "mm")) +
  
  #plot.title = element_text( Size = 18)) +
  
  #guides(guide_legend(nrow=4)) +
  
  #legend.spacing.y = unit(1, "mm") +
  
  #legend.key.size = unit(1, "line")) 
  
  #scale_fill_manual(values=values, labels=setNames(paste(labels, " "), entries)) +
  
#ggtitle("Bulk Soil, 30-45 cm") +
  
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5, vjust = 1)) + ## added title and adjusted
  
  scale_y_continuous(expand = c(0,0), limit = c(0, 100)) + scale_x_discrete(expand = c(0.03,0)) +
                                                                               
                                  coord_flip() ### removed excess space btw bars and axis and swapped axes

bar.middle.mass1 ####plant effect only sig in sand




bar.bottom.mass1 <- ggplot(BS_TMX.mass.bottom.bar, aes(Texture, TMX.micrg, fill = Stage)) +
  
  geom_bar(position = position_dodge(), stat = "identity") +
  
  geom_bar(position = "dodge", stat = "identity", color="black", show.legend =FALSE, lwd = 0.75) + ##add outline to bars
  
  geom_errorbar(aes(ymin = TMX.micrg, ymax = TMX.micrg + se),      ##error bars
                
                width = 0.2, position = position_dodge(0.9), lwd = 0.75) +
  
  cleanup +                                          #cleanup
  
  scale_fill_manual(name = "",
                    
                    labels = c("V1", "V3", "V5", "V5, No Plant"), 
                    
                    values = c("blue", "red", "green", "yellow")) +
  
  xlab("") +
  
  ylab(expression(TMX ~ (mu * g))) + ### pain in the ass way to express superscripts
  
  theme(axis.text.y = element_text(size = 16),  ## change font size of y axis
        
        axis.title = element_text(size = 16, face = "bold")) +
  
  #geom_label(aes(label = letters)) +
  
  #geom_text(aes(label = letters)) +
  
  geom_text(aes(label = letters, y = TMX.micrg + se, fontface = "italic"), position = position_dodge(width=0.9), hjust= -0.5) +  #####add in signficance letter--adjusted for errror bars
  
  theme(axis.text.x = element_text(size = 16),  ## change font size of y axis
        
        axis.title = element_text(size = 16, face = "bold")) +
  
  theme(axis.line.x = element_line(size = 1), ##change axis line size
        
        axis.line.y = element_line(size = 1),
        
        panel.border = element_rect(color = "black", fill = NA, size = 1.75), ##add border
        
        legend.position = "none", ## move legend
        
        legend.margin = margin(0, 0, 0, 0, "mm"), ##remove margin on legend
        
        legend.text = element_text( size = 16), ##legend text
        
        #legend.background = element_rect(color = "black", size = 1), ## add legend border
        
        #legend.key.height=unit(1,"line"), 
        
        #legend.key.width=unit(1,"line")
        
        legend.key = element_rect(color = "black", size = 0.75), ## put line around legend key
        
        legend.direction = "vertical", 
        
        axis.ticks.y = element_line(size = 1, color = "black"),  ##enlargen ticks--problem at op right corner
        
        legend.key.size = unit(7, 'mm'),
        
        legend.key.height = unit(6, "mm")) +
  
  #plot.title = element_text( Size = 18)) +
  
  #guides(guide_legend(nrow=4)) +
  
  #legend.spacing.y = unit(1, "mm") +
  
  #legend.key.size = unit(1, "line")) 
  
  #scale_fill_manual(values=values, labels=setNames(paste(labels, " "), entries)) +
  
#ggtitle("Bulk Soil, 45-60 cm") +
  
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5, vjust = 1)) + ## added title and adjusted
  
  scale_y_continuous(expand = c(0,0), limit = c(0, 100)) + scale_x_discrete(expand = c(0.03,0)) +

 coord_flip()### removed excess space btw bars and axis and swapped x and y


bar.bottom.mass1


#multiplot(bar.top, bar.middle, bar.bottom, cols = 3)



?lapply

?grid.arrange


grid.arrange(bar.top1, bar.middle1, bar.bottom1, bar.top.mass1,
                         
            bar.middle.mass1, bar.bottom.mass1,
             
             layout_matrix = rbind( c(1, 2, 2, NA ),
                                    c(1, 3, 3, NA),
                                    c(4, 5, 5, NA),
                                    c(4, 6, 6, NA))) 
                                    
                          
                                                                                    



#grid.arrange(bar.top1, bar.middle1, bar.bottom1, bar.top.mass1,
             
            #bar.middle.mass1, bar.bottom.mass1, ncol = 3)

?grid.arrange()


########### retry multiplotting


bar.top2 <- ggplot(BS_top.bar, aes(Texture, TMX, fill = Stage)) +
  
  geom_bar(position = position_dodge(), stat = "identity") +
  
  geom_bar(position = "dodge", stat = "identity", color="black", show.legend =FALSE, lwd = 0.75) + ##add outline to bars
  
  geom_errorbar(aes(ymin = TMX, ymax = TMX + se),      ##error bars
                
                width = 0.2, position = position_dodge(0.9), lwd = 0.75) +
  
  cleanup +                                          #cleanup
  
  scale_fill_manual(name = "",
                    
                    labels = c("V1", "V3", "V5", "V5, No Plant"), 
                    
                    values = c("blue", "red", "green", "yellow")) +
  
  xlab("") +
  
  ylab(expression(TMX ~ (ng ~ g^{-1}))) + ### pain in the ass way to express superscripts
  
  theme(axis.text.y = element_text(size = 16),  ## change font size of y axis
        
        axis.title = element_text(size = 16, face = "bold")) +
  
  #geom_label(aes(label = letters)) +
  
  #geom_text(aes(label = letters)) +
  
  geom_text(aes(label = letters, y = TMX + se, fontface = "italic"), position = position_dodge(width=0.9), vjust= -0.5) +  #####add in signficance letter--adjusted for errror bars
  
  theme(axis.text.x = element_text(size = 16),  ## change font size of y axis
        
        axis.title = element_text(size = 16, face = "bold")) +
  
  theme(axis.line.x = element_line(size = 1), ##change axis line size
        
        axis.line.y = element_line(size = 1),
        
        panel.border = element_rect(color = "black", fill = NA, size = 1.75), ##add border
        
        legend.position = "none", ## move legend
        
        legend.margin = margin(0, 0, 0, 0, "mm"), ##remove margin on legend
        
        legend.text = element_text( size = 16), ##legend text
        
        #legend.background = element_rect(color = "black", size = 1), ## add legend border
        
        #legend.key.height=unit(1,"line"), 
        
        #legend.key.width=unit(1,"line")
        
        legend.key = element_rect(color = "black", size = 0.75), ## put line around legend key
        
        legend.direction = "vertical", 
        
        axis.ticks.y = element_line(size = 1, color = "black"),  ##enlargen ticks--problem at op right corner
        
        legend.key.size = unit(7, 'mm'),
        
        legend.key.height = unit(6, "mm")) +
  
  #plot.title = element_text( Size = 18)) +
  
  #guides(guide_legend(nrow=4)) +
  
  #legend.spacing.y = unit(1, "mm") +
  
  #legend.key.size = unit(1, "line")) 
  
  #scale_fill_manual(values=values, labels=setNames(paste(labels, " "), entries)) +
  
ggtitle("a) Bulk Soil, 0-30 cm") +
  
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.7, vjust = 1)) + ## added title and adjusted
  
  scale_y_continuous(expand = c(0,0), limit = c(0, 110)) + scale_x_discrete(expand = c(0.03,0)) ### removed excess space btw bars and axis




bar.middle2 <- ggplot(BS_middle.bar, aes(Texture, TMX, fill = Stage)) +
  
  geom_bar(position = position_dodge(), stat = "identity") +
  
  geom_bar(position = "dodge", stat = "identity", color="black", show.legend =FALSE, lwd = 0.75) + ##add outline to bars
  
  geom_errorbar(aes(ymin = TMX, ymax = TMX + se),      ##error bars
                
                width = 0.2, position = position_dodge(0.9), lwd = 0.75) +
  
  cleanup +                                          #cleanup
  
  scale_fill_manual(name = "",
                    
                    labels = c("V1", "V3", "V5", "V5, No Plant"), 
                    
                    values = c("blue", "red", "green", "yellow")) +
  
  xlab("") +
  
  ylab(expression(TMX ~ (ng ~ g^{-1}))) + ### pain in the ass way to express superscripts
  
  theme(axis.text.y = element_text(size = 16),  ## change font size of y axis
        
        axis.title = element_text(size = 16, face = "bold")) +
  
  #geom_label(aes(label = letters)) +
  
  #geom_text(aes(label = letters)) +
  
  geom_text(aes(label = letters, y = TMX + se, fontface = "italic"), position = position_dodge(width=0.9), vjust= -0.5) +  #####add in signficance letter--adjusted for errror bars
  
  theme(axis.text.x = element_text(size = 16),  ## change font size of y axis
        
        axis.title = element_text(size = 16, face = "bold")) +
  
  theme(axis.line.x = element_line(size = 1), ##change axis line size
        
        axis.line.y = element_line(size = 1),
        
        panel.border = element_rect(color = "black", fill = NA, size = 1.75), ##add border
        
        legend.position = c(0.35,0.6), ## move legend
        
        legend.margin = margin(0, 0, 0, 0, "mm"), ##remove margin on legend
        
        legend.text = element_text( size = 20), ##legend text
        
        #legend.background = element_rect(color = "black", size = 1), ## add legend border
        
        #legend.key.height=unit(1,"line"), 
        
        #legend.key.width=unit(1,"line")
        
        legend.key = element_rect(color = "black", size = 0.75), ## put line around legend key
        
        legend.direction = "vertical", 
        
        axis.ticks.y = element_line(size = 1, color = "black"),  ##enlargen ticks--problem at op right corner
        
        legend.key.size = unit(7, 'mm'),
        
        legend.key.height = unit(6, "mm")) +
  
  #plot.title = element_text( Size = 18)) +
  
  #guides(guide_legend(nrow=4)) +
  
  #legend.spacing.y = unit(1, "mm") +
  
  #legend.key.size = unit(1, "line")) 
  
  #scale_fill_manual(values=values, labels=setNames(paste(labels, " "), entries)) +
  
ggtitle("b)  Bulk Soil, 30-45 cm") + annotate("text", x = 1.5, y= 25, label = "*** 1/4 scale in 0-30 cm section") +
  
  theme(plot.title = element_text(size = 16, face = "bold", vjust = 1)) + ## added title and adjusted
  
  scale_y_continuous(expand = c(0,0), limit = c(0, 27.5)) + scale_x_discrete(expand = c(0.03,0)) ### removed excess space btw bars and axis


bar.middle2



bar.bottom2 <- ggplot(BS_bottom.bar, aes(Texture, TMX, fill = Stage)) +
  
  geom_bar(position = position_dodge(), stat = "identity") +
  
  geom_bar(position = "dodge", stat = "identity", color="black", show.legend =FALSE, lwd = 0.75) + ##add outline to bars
  
  geom_errorbar(aes(ymin = TMX, ymax = TMX + se),      ##error bars
                
                width = 0.2, position = position_dodge(0.9), lwd = 0.75) +
  
  cleanup +                                          #cleanup
  
  scale_fill_manual(name = "",
                    
                    labels = c("V1", "V3", "V5", "V5, No Plant"), 
                    
                    values = c("blue", "red", "green", "yellow")) +
  
  xlab("") +
  
  ylab(expression(TMX ~ (ng ~ g^{-1}))) + ### pain in the ass way to express superscripts
  
  theme(axis.text.y = element_text(size = 16),  ## change font size of y axis
        
        axis.title = element_text(size = 16, face = "bold")) +
  
  #geom_label(aes(label = letters)) +
  
  #geom_text(aes(label = letters)) +
  
  geom_text(aes(label = letters, y = TMX + se, fontface = "italic"), position = position_dodge(width=0.9), vjust= -0.5) +  #####add in signficance letter--adjusted for errror bars
  
  theme(axis.text.x = element_text(size = 16),  ## change font size of y axis
        
        axis.title = element_text(size = 16, face = "bold")) +
  
  theme(axis.line.x = element_line(size = 1), ##change axis line size
        
        axis.line.y = element_line(size = 1),
        
        panel.border = element_rect(color = "black", fill = NA, size = 1.75), ##add border
        
        legend.position = "none", ## move legend
        
        legend.margin = margin(0, 0, 0, 0, "mm"), ##remove margin on legend
        
        legend.text = element_text( size = 16), ##legend text
        
        #legend.background = element_rect(color = "black", size = 1), ## add legend border
        
        #legend.key.height=unit(1,"line"), 
        
        #legend.key.width=unit(1,"line")
        
        legend.key = element_rect(color = "black", size = 0.75), ## put line around legend key
        
        legend.direction = "vertical", 
        
        axis.ticks.y = element_line(size = 1, color = "black"),  ##enlargen ticks--problem at op right corner
        
        legend.key.size = unit(7, 'mm'),
        
        legend.key.height = unit(6, "mm")) +
  
  #plot.title = element_text( Size = 18)) +
  
  #guides(guide_legend(nrow=4)) +
  
  #legend.spacing.y = unit(1, "mm") +
  
  #legend.key.size = unit(1, "line")) 
  
  #scale_fill_manual(values=values, labels=setNames(paste(labels, " "), entries)) +
  
ggtitle("c)  Bulk Soil, 45-60 cm") + annotate("text", x = 1.5, y= 25, label = "*** 1/4 scale in 0-30 cm section") +
  
  theme(plot.title = element_text(size = 16, face = "bold", vjust = 1)) + ## added title and adjusted
  
  scale_y_continuous(expand = c(0,0), limit = c(0, 27.5)) + scale_x_discrete(expand = c(0.03,0)) ### removed excess space btw bars and axis


pdf("C:/Users/Jesse/Desktop/Neonicotinoids_15/Data/Column-Study/Raw-Data-Outputs/bar.top2.pdf")
bar.top2
dev.off()

pdf("C:/Users/Jesse/Desktop/Neonicotinoids_15/Data/Column-Study/Raw-Data-Outputs/bar.middle2.pdf")
bar.middle2
dev.off()

pdf("C:/Users/Jesse/Desktop/Neonicotinoids_15/Data/Column-Study/Raw-Data-Outputs/bar.bottom2.pdf")
bar.bottom2
dev.off()


################################################MASS

bar.top.mass2 <- ggplot(BS_TMX.mass.top.bar, aes(Texture, TMX.micrg, fill = Stage)) +
  
  geom_bar(position = position_dodge(), stat = "identity") +
  
  geom_bar(position = "dodge", stat = "identity", color="black", show.legend =FALSE, lwd = 0.75) + ##add outline to bars
  
  geom_errorbar(aes(ymin = TMX.micrg, ymax = TMX.micrg + se),      ##error bars
                
                width = 0.2, position = position_dodge(0.9), lwd = 0.75) +
  
  cleanup +                                          #cleanup
  
  scale_fill_manual(name = "",
                    
                    labels = c("V1", "V3", "V5", "V5, No Plant"), 
                    
                    values = c("blue", "red", "green", "yellow")) +
  
  xlab("") +
  
  ylab(expression(TMX ~ (mu * g))) + ### pain in the ass way to express superscripts
  
  theme(axis.text.y = element_text(size = 16),  ## change font size of y axis
        
        axis.title = element_text(size = 16, face = "bold")) +
  
  #geom_label(aes(label = letters)) +
  
  #geom_text(aes(label = letters)) +
  
  geom_text(aes(label = letters, y = TMX.micrg + se, fontface = "italic"), position = position_dodge(width=0.9), vjust= -0.5) +  #####add in signficance letter--adjusted for errror bars
  
  theme(axis.text.x = element_text(size = 16),  ## change font size of y axis
        
        axis.title = element_text(size = 16, face = "bold"),
        
        axis.title.y = element_text(vjust = -0.7)) +
  
  theme(axis.line.x = element_line(size = 1), ##change axis line size
        
        axis.line.y = element_line(size = 1),
        
        panel.border = element_rect(color = "black", fill = NA, size = 1.75), ##add border
        
        legend.position = "none", ## move legend
        
        legend.margin = margin(0, 0, 0, 0, "mm"), ##remove margin on legend
        
        legend.text = element_text( size = 16), ##legend text
        
        #legend.background = element_rect(color = "black", size = 1), ## add legend border
        
        #legend.key.height=unit(1,"line"), 
        
        #legend.key.width=unit(1,"line")
        
        legend.key = element_rect(color = "black", size = 0.75), ## put line around legend key
        
        legend.direction = "vertical", 
        
        axis.ticks.y = element_line(size = 1, color = "black"),  ##enlargen ticks--problem at op right corner
        
        legend.key.size = unit(7, 'mm'),
        
        legend.key.height = unit(6, "mm"),
        
        plot.margin = unit(c(0,1,0,0),"mm")) +
  
  #plot.title = element_text( Size = 18)) +
  
  #guides(guide_legend(nrow=4)) +
  
  #legend.spacing.y = unit(1, "mm") +
  
  #legend.key.size = unit(1, "line")) 
  
  #scale_fill_manual(values=values, labels=setNames(paste(labels, " "), entries)) +
  
ggtitle("d)") +
  
  theme(plot.title = element_text(size = 16, face = "bold", vjust = 1)) + ## added title and adjusted
  
  scale_y_continuous(expand = c(0,0), limit = c(0, 1620)) + scale_x_discrete(expand = c(0.03,0)) ### removed excess space btw bars and axis

bar.top.mass2



bar.middle.mass2 <- ggplot(BS_TMX.mass.middle.bar, aes(Texture, TMX.micrg, fill = Stage)) +
  
  geom_bar(position = position_dodge(), stat = "identity") +
  
  geom_bar(position = "dodge", stat = "identity", color="black", show.legend =FALSE, lwd = 0.75) + ##add outline to bars
  
  geom_errorbar(aes(ymin = TMX.micrg, ymax = TMX.micrg + se),      ##error bars
                
                width = 0.2, position = position_dodge(0.9), lwd = 0.75) +
  
  cleanup +                                          #cleanup
  
  scale_fill_manual(name = "",
                    
                    labels = c("V1", "V3", "V5", "V5, No Plant"), 
                    
                    values = c("blue", "red", "green", "yellow")) +
  
  xlab("") +
  
  ylab("") + ### pain in the ass way to express superscripts
  
  theme(axis.text.y = element_text(size = 16),  ## change font size of y axis
        
        axis.title = element_text(size = 16, face = "bold")) +
  
  #geom_label(aes(label = letters)) +
  
  #geom_text(aes(label = letters)) +
  
  geom_text(aes(label = letters, y = TMX.micrg + se, fontface = "italic"), position = position_dodge(width=0.9), vjust= -0.5) +  #####add in signficance letter--adjusted for errror bars
  
  theme(axis.text.x = element_text(size = 16),  ## change font size of y axis
        
        axis.title = element_text(size = 16, face = "bold")) +
  
  theme(axis.line.x = element_line(size = 1), ##change axis line size
        
        axis.line.y = element_line(size = 1),
        
        panel.border = element_rect(color = "black", fill = NA, size = 1.75), ##add border
        
        legend.position = c(0.35,0.6), ## move legend
        
        legend.margin = margin(0, 0, 0, 0, "mm"), ##remove margin on legend
        
        legend.text = element_text( size = 16), ##legend text
        
        #legend.background = element_rect(color = "black", size = 1), ## add legend border
        
        #legend.key.height=unit(1,"line"), 
        
        #legend.key.width=unit(1,"line")
        
        legend.key = element_rect(color = "black", size = 0.75), ## put line around legend key
        
        legend.direction = "vertical", 
        
        axis.ticks.y = element_line(size = 1, color = "black"),  ##enlargen ticks--problem at op right corner
        
        legend.key.size = unit(7, 'mm'),
        
        legend.key.height = unit(6, "mm"),
        
        plot.margin = unit(c(0,1,0,0),"mm")) +
  
  #plot.title = element_text( Size = 18)) +
  
  #guides(guide_legend(nrow=4)) +
  
  #legend.spacing.y = unit(1, "mm") +
  
  #legend.key.size = unit(1, "line")) 
  
  #scale_fill_manual(values=values, labels=setNames(paste(labels, " "), entries)) +
  
ggtitle("e)") + annotate("text", x = 1.5, y= 380, label = "*** 1/4 scale in 0-30 cm section") +
  
  theme(plot.title = element_text(size = 16, face = "bold", vjust = 1)) + ## added title and adjusted
  
  scale_y_continuous(expand = c(0,0), limit = c(0, 420)) + scale_x_discrete(expand = c(0.03,0)) ### removed excess space btw bars and axis

bar.middle.mass2 ####plant effect only sig in sand




bar.bottom.mass2 <- ggplot(BS_TMX.mass.bottom.bar, aes(Texture, TMX.micrg, fill = Stage)) +
  
  geom_bar(position = position_dodge(), stat = "identity") +
  
  geom_bar(position = "dodge", stat = "identity", color="black", show.legend =FALSE, lwd = 0.75) + ##add outline to bars
  
  geom_errorbar(aes(ymin = TMX.micrg, ymax = TMX.micrg + se),      ##error bars
                
                width = 0.2, position = position_dodge(0.9), lwd = 0.75) +
  
  cleanup +                                          #cleanup
  
  scale_fill_manual(name = "",
                    
                    labels = c("V1", "V3", "V5", "V5, No Plant"), 
                    
                    values = c("blue", "red", "green", "yellow")) +
  
  xlab("") +
  
  ylab("") + ### pain in the ass way to express superscripts
  
  theme(axis.text.y = element_text(size = 16),  ## change font size of y axis
        
        axis.title = element_text(size = 16, face = "bold")) +
  
  #geom_label(aes(label = letters)) +
  
  #geom_text(aes(label = letters)) +
  
  geom_text(aes(label = letters, y = TMX.micrg + se, fontface = "italic"), position = position_dodge(width=0.9), vjust= -0.5) +  #####add in signficance letter--adjusted for errror bars
  
  theme(axis.text.x = element_text(size = 16),  ## change font size of y axis
        
        axis.title = element_text(size = 16, face = "bold")) +
  
  theme(axis.line.x = element_line(size = 1), ##change axis line size
        
        axis.line.y = element_line(size = 1),
        
        panel.border = element_rect(color = "black", fill = NA, size = 1.75), ##add border
        
        legend.position = "none", ## move legend
        
        legend.margin = margin(0, 0, 0, 0, "mm"), ##remove margin on legend
        
        legend.text = element_text( size = 16), ##legend text
        
        #legend.background = element_rect(color = "black", size = 1), ## add legend border
        
        #legend.key.height=unit(1,"line"), 
        
        #legend.key.width=unit(1,"line")
        
        legend.key = element_rect(color = "black", size = 0.75), ## put line around legend key
        
        legend.direction = "vertical", 
        
        axis.ticks.y = element_line(size = 1, color = "black"),  ##enlargen ticks--problem at op right corner
        
        legend.key.size = unit(7, 'mm'),
        
        legend.key.height = unit(6, "mm"),
        
        plot.margin = unit(c(0,1,0,0),"mm")) +
  
  #plot.title = element_text( Size = 18)) +
  
  #guides(guide_legend(nrow=4)) +
  
  #legend.spacing.y = unit(1, "mm") +
  
  #legend.key.size = unit(1, "line")) 
  
  #scale_fill_manual(values=values, labels=setNames(paste(labels, " "), entries)) +
  
ggtitle("f)") + annotate("text", x = 1.5, y= 380, label = "*** 1/4 scale in 0-30 cm section") +
  
  theme(plot.title = element_text(size = 16, face = "bold", vjust = 1)) + ## added title and adjusted
  
  scale_y_continuous(expand = c(0,0), limit = c(0, 420)) + scale_x_discrete(expand = c(0.03,0)) ### removed excess space btw bars and axis

?annotate

bar.bottom.mass2




BS.combo <- grid.arrange(bar.top2, bar.middle2, bar.bottom2, ncol = 1)
BS.combo

ggsave(file = "C:/Users/Jesse/Desktop/Neonicotinoids_15/Data/Column-Study/Raw-Data-Outputs/BS.combo.seminar.pdf",
       
       BS.combo)
dev.off()



