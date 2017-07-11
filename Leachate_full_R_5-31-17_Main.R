rm(list = ls())
##test

library(xlsx)
library(dplyr)
library(ggplot2)
library(Rmisc)
library(gridExtra)
library(grid)
library(multcompView)
library(lawstat)
#library(dplyr)


##JBR
##5/16/17
###Estimating ET for Column Experiment... and analyzing leachate/TMX transport



###Assumptions:
#             - no change in storage over time after initial drainage following
#               saturation event, or theta was ~ equivalent to field capacity
#             - no evaporation during rain events... (probably realistic since, the simulator "chamber" was closed
#               and rain event was 7 mins in duration)
#             - so I = P
#             -  I - Leachate = ET


#######start Data cleaning and initial calculations (using code from Leachate_full)

###includes TMX mass and conc analysis from Leachate_full...keep same organization/structure

###note that 300 mL will be maxiumum "ET" from events 1-11 and 3000 mL Will be max ET for event 12,
#whereas in reality torage was probably changing before final event and ET > Inf


Leachate_full <- read.csv("C:/Users/Jesse/Desktop/Neonicotinoids_15/Data/Column-Study/Raw-Data-Outputs/Leachate_full_R.csv")
attach(Leachate_full)


Leachate_full

str(Leachate_full)

head(Leachate_full)



###first clean data source

###leachate volumes <0.3 = 0

#if(Leachate_full$Leach.Vol.mL <= 0.3) {
#  Leachate_full$Leach.Vol.mL <- 0
#}

Leachate_full$Leach.Vol.mL[Leach.Vol.mL <= 0.3] <- 0 #### this works!

#Leachate_full$Leach.Vol.mL[Leach.Vol.mL == 0] <- NA   ###replace 0 with NA in volume leachate

#Leachate_full[Leachate_full$Leach.Vol.mL == 0, Leachate_full$TMX_PPB] <- 0

#Leachate_full$TMX_PPB[TMX_PPB == 0] <- NA   ###replace 0 with NA in leachate conc

Leachate_full$TMX_PPB <- ifelse(Leachate_full$Leach.Vol.mL == 0, NA, Leachate_full$TMX_PPB)

View(Leachate_full)




cum.na <- function(x) { 
  
  x[which(is.na(x))] <- 0 
  
  return(cumsum(x))  ####### use this to treat na as zero in cumulative calculation. IMPORTANT!!!!!
} 

Leachate_full$cum.Leach.vol.mL <- ave(Leachate_full$Leach.Vol.mL,
                                      
                                      Leachate_full$Sample.ID, FUN=cum.na) 

### cumulative function but NAs not removed in calculation, "KF-60 cm-LS/SS-5" leached less than 400 mL total!

### ave() is a shorthand way of treating a function by group (for stacked data)----just adjust function within body^^



####TMX_mass first

Leachate_full$TMX.ng <- (Leachate_full$TMX_PPB * Leachate_full$Leach.Vol.mL)

Leachate_full$TMX.micrg <- Leachate_full$TMX.ng/1000



####now calculate the cumulative TMX mass transported 

Leachate_full$cum.TMX.ng <- ave(Leachate_full$TMX.ng,
                                
                                Leachate_full$Sample.ID, FUN=cum.na) ####good!

Leachate_full$cum.TMX.micrg <- Leachate_full$cum.TMX.ng/1000



##now remove outliers in final leachate event not present in stat analysis: "KF-60 cm-LS-Exp-4",

#"KF-60 cm-LS/SS-4", "NK-20 cm-LS-Exp-4"

Leachate_full[Leachate_full$Sample.ID == "KF-60 cm-LS-Exp-4",] ##row 436

Leachate_full[Leachate_full$Sample.ID == "KF-60 cm-LS/SS-4",]  ##row 467

Leachate_full[Leachate_full$Sample.ID == "NK-20 cm-LS-Exp-4",]  ## row 463

#Leachate_full <- subset(Leachate_full, Sample.ID != "KF-60 cm-LS-Exp-4" & Event != 12)

#Leachate_full[Leachate_full$Sample.ID == "KF-60 cm-LS-Exp-4",]

Leachate_full <- Leachate_full[-c(436, 467, 463),] ###this is a useful code---- from dplyr? removed selected!

Leachate_full


tail(Leachate_full)




###ET estimates##########################################################################

###need ET(rate) and ET(cumualtive)

Leachate_full$Leach.Vol.mL[is.na(Leachate_full$Leach.Vol.mL)] <- 0 ##change NA's back to 0 for calculation





### ET(cumulative) first (mL)

Leachate_full$ET.cum.mL <- Leachate_full$cum.inf.mL - Leachate_full$cum.Leach.vol.mL 

Leachate_full

###a couple of negative values: KF-20 cm-LS-Ctr-5 and KF-20 cm-LS-Exp-3

### Leachate exceeded input... these columns may have been slightly above field capacity initailly?

###should be fine when averaged





###convert to mm ... more conventional to use length
##column area

radius <- 10 ## radius in cm 
area <- pi * (radius)^2

Leachate_full$ET.cum.mm <- (Leachate_full$ET.cum.mL / area) * 10 ## cumulative ET in mm

tail(Leachate_full)

Leachate_full[Leachate_full$Event == 11,]



##ET(t) mL
Leachate_full$ET.mL <- Leachate_full$Inf.mL - Leachate_full$Leach.Vol.mL
##negative values towards the end... why???

Leachate_full[Leachate_full$Event == 1, ]

Leachate_full[Leachate_full$ET.mL <= 0, ] ## where Vleach > Vinf
##seems highest (most negative in Clay-20cm-no plant)
##wait and see if it affects the average but... should be expressed anyway






##now ET(rate) (mm/day)

Leachate_full$ET.rate.mm.d <- ((Leachate_full$ET.mL/area)* 10)/Leachate_full$dt

Leachate_full$ET.rate.cm.d <- ((Leachate_full$ET.mL/area))/Leachate_full$dt ##in cm/day


Leachate_full[Leachate_full$Event == 10,]
Leachate_full[Leachate_full$Event == 11,]
Leachate_full[Leachate_full$Event == 12,]


######try p rate - et rate (in cm day^-1)


Leachate_full$p-et.cm.day <- 




















######prepare for plot



######all 60 cm columns

tall.columns <- Leachate_full[Leachate_full$Size_cm == "60",]

tall.columns

tall.columns.stats <- summarySE(tall.columns, measurevar="TMX_PPB",
                                
                                groupvars=c("Texture","Plant.Influence",
                                            
                                            "Time.day"), na.rm =TRUE) ## use na.rm = TRUE to omit NA from calculation

#tall.columns.stats

#Letters.tall.leach

#tall.columns.stats$Letters <- ifelse(tall.columns.stats$TMX_PPB == 10.25525743, 6, 

#                                    1)
tall.columns.stats






#####tall columns without structureless trt*********************************

#tall.columns[tall.columns$Structure == "NS",]


tall.columns.no.SS <- subset(tall.columns, Structure == "S" | Texture == "Sand")  ####good, works

tall.columns.no.SS

tall.columns.no.SS.stats <- summarySE(tall.columns.no.SS, measurevar="TMX_PPB",
                                      
                                      groupvars=c("Texture","Plant.Influence",
                                                  
                                                  "Time.day"), na.rm =TRUE) ## 
tall.columns.no.SS.stats

### need to add column with cum.inf.cm

tall.columns.no.SS.stats$Cum.inf.cm <- ifelse(tall.columns.no.SS.stats$Time.day == 3, 0.9,
                                              ifelse(tall.columns.no.SS.stats$Time.day == 6, 1.8,
                                                     ifelse(tall.columns.no.SS.stats$Time.day == 9, 2.7,
                                                            ifelse(tall.columns.no.SS.stats$Time.day == 12, 3.6,
                                                                   ifelse(tall.columns.no.SS.stats$Time.day == 15, 4.5,
                                                                          ifelse(tall.columns.no.SS.stats$Time.day == 18, 5.4,
                                                                                 ifelse(tall.columns.no.SS.stats$Time.day == 21, 6.3,
                                                                                        ifelse(tall.columns.no.SS.stats$Time.day == 24, 7.2,
                                                                                               ifelse(tall.columns.no.SS.stats$Time.day == 27, 8.1,
                                                                                                      ifelse(tall.columns.no.SS.stats$Time.day == 30, 9,
                                                                                                             ifelse(tall.columns.no.SS.stats$Time.day == 31, 9.9,
                                                                                                                    ifelse(tall.columns.no.SS.stats$Time.day == 33, 18.9,
                                                                                                                           NA
                                                                                                                    ))))))))))))




tall.columns.no.SS.stats

Letters.tall.leach


row.names(tall.columns.no.SS.stats) = "12" 

######## add in sig letters and make all exclusions = NA

tall.columns.no.SS.stats$Letters <-  ifelse(tall.columns.no.SS.stats$Texture == "Clay" & 
                                              
                                              tall.columns.no.SS.stats$Time.day == "33", "a",
                                            
                                            ifelse (tall.columns.no.SS.stats$Texture == "Sand" & 
                                                      
                                                      tall.columns.no.SS.stats$Time.day == "33", "b",
                                                    
                                                    NA)
)


######Good 



tall.columns.no.SS.stats

##now make factorial combo (eg. clay No Plant)

tall.columns.no.SS.stats$Combo <- paste(tall.columns.no.SS.stats$Texture,
                                        
                                        tall.columns.no.SS.stats$Plant.Influence)

tall.columns.no.SS.stats ##good



###########Structure Study*******************

structure.study <- subset(Leachate_full, Size_cm == "60" &
                            
                            Structure == "S" |
                            
                            Structure == "NS")

structure.study

structure.study <- subset(structure.study, structure.study$Plant.Influence == "Plant") ## exclude "No Plant"

tail(structure.study)

structure.study

structure.study.stats <- summarySE(structure.study, measurevar="TMX_PPB",
                                   
                                   groupvars=c("Structure",
                                               
                                               "Time.day"), na.rm =TRUE) ## stat summary

structure.study.stats

structure.study.stats$Letters <-  ifelse(structure.study.stats$Structure == "NS" & 
                                           
                                           structure.study.stats$Time.day == "33", "a",
                                         
                                         ifelse (structure.study.stats$Structure == "S" & 
                                                   
                                                   structure.study.stats$Time.day == "33", "b",
                                                 
                                                 NA)
)

structure.study.stats

structure.study.stats$Cum.inf.cm <- ifelse(structure.study.stats$Time.day == 3, 0.9,
                                           ifelse(structure.study.stats$Time.day == 6, 1.8,
                                                  ifelse(structure.study.stats$Time.day == 9, 2.7,
                                                         ifelse(structure.study.stats$Time.day == 12, 3.6,
                                                                ifelse(structure.study.stats$Time.day == 15, 4.5,
                                                                       ifelse(structure.study.stats$Time.day == 18, 5.4,
                                                                              ifelse(structure.study.stats$Time.day == 21, 6.3,
                                                                                     ifelse(structure.study.stats$Time.day == 24, 7.2,
                                                                                            ifelse(structure.study.stats$Time.day == 27, 8.1,
                                                                                                   ifelse(structure.study.stats$Time.day == 30, 9,
                                                                                                          ifelse(structure.study.stats$Time.day == 31, 9.9,
                                                                                                                 ifelse(structure.study.stats$Time.day == 33, 18.9,
                                                                                                                        NA
                                                                                                                 ))))))))))))

structure.study.stats

##now make factorial combo (eg. clay No Plant)

#structure.study.stats$Combo <- paste(structure.study.stats$Texture,

#                                     structure.study.stats$Plant.Influence)

#structure.study.stats ##good






###good


###small columns

small.columns <- Leachate_full[Leachate_full$Size_cm == "20",]

small.columns

tail(small.columns)

small.columns.stats <- summarySE(small.columns, measurevar="TMX_PPB",
                                 
                                 groupvars=c("Texture", "Plant.Influence",
                                             
                                             "Time.day"), na.rm =TRUE)


Letters.small

small.columns.stats$Letters <- ifelse(small.columns.stats$Texture == "Clay" & 
                                        
                                        small.columns.stats$Time.day == "33" &
                                        
                                        small.columns.stats$Plant.Influence == "No Plant", "ab",
                                      
                                      ifelse (small.columns.stats$Texture == "Clay" & 
                                                
                                                small.columns.stats$Time.day == "33" &
                                                
                                                small.columns.stats$Plant.Influence == "Plant",
                                              
                                              "b", 
                                              
                                              ifelse(small.columns.stats$Texture == "Sand" & 
                                                       
                                                       small.columns.stats$Time.day == "33" &
                                                       
                                                       small.columns.stats$Plant.Influence == "No Plant",
                                                     
                                                     "ab",
                                                     
                                                     ifelse(small.columns.stats$Texture == "Sand" & 
                                                              
                                                              small.columns.stats$Time.day == "33" &
                                                              
                                                              small.columns.stats$Plant.Influence == "Plant",
                                                            
                                                            "a", NA)
                                              )
                                      )
)


small.columns.stats

Letters.small

small.columns.stats$Cum.inf.cm <- ifelse(small.columns.stats$Time.day == 3, 0.9,
                                         ifelse(small.columns.stats$Time.day == 6, 1.8,
                                                ifelse(small.columns.stats$Time.day == 9, 2.7,
                                                       ifelse(small.columns.stats$Time.day == 12, 3.6,
                                                              ifelse(small.columns.stats$Time.day == 15, 4.5,
                                                                     ifelse(small.columns.stats$Time.day == 18, 5.4,
                                                                            ifelse(small.columns.stats$Time.day == 21, 6.3,
                                                                                   ifelse(small.columns.stats$Time.day == 24, 7.2,
                                                                                          ifelse(small.columns.stats$Time.day == 27, 8.1,
                                                                                                 ifelse(small.columns.stats$Time.day == 30, 9,
                                                                                                        ifelse(small.columns.stats$Time.day == 31, 9.9,
                                                                                                               ifelse(small.columns.stats$Time.day == 33, 18.9,
                                                                                                                      NA
                                                                                                               ))))))))))))

small.columns.stats

##now make factorial combo (eg. clay No Plant)

small.columns.stats$Combo <- paste(small.columns.stats$Texture,
                                   
                                   small.columns.stats$Plant.Influence)

small.columns.stats ##good








#############stats on mass TMX transported

#### tall columns no SS first

tall.columns.no.SS

tall.columns.no.SS.final <- subset(tall.columns.no.SS, Event == 12) ###need to isolate final event


tall.columns.no.SS.final
###***

tall.columns.no.SS.final

tall.columns.no.SS.final$log10.TMX.micrg <- log10(tall.columns.no.SS.final$cum.TMX.micrg)## log 10

tall.columns.no.SS.final$rank.TMX.micrg <- rank(tall.columns.no.SS.final$cum.TMX.micrg)## rank transform data to normality

boxplot(tall.columns.no.SS.final$rank.TMX.micrg ~
          
          tall.columns.no.SS.final$Texture * tall.columns.no.SS.final$Plant.Influence,
        
        ylab = expression(TMX ~ mu*g), xlab = "Trt")## 

bartlett.test(tall.columns.no.SS.final$rank.TMX.micrg
              
              ~ interaction(tall.columns.no.SS.final$Texture,
                            
                            tall.columns.no.SS.final$Plant.Influence)) ## 

fligner.test(tall.columns.no.SS.final$rank.TMX.micrg
             
             ~ interaction(tall.columns.no.SS.final$Texture,
                           
                           tall.columns.no.SS.final$Plant.Influence)) ## 

levene.test(tall.columns.no.SS.final$rank.TMX.micrg ~
              
              tall.columns.no.SS.final$Texture * tall.columns.no.SS.final$Plant.Influence) ## no HOV but i guess it doesn't mater since npar




ANOVA.tall.mass <- aov(tall.columns.no.SS.final$rank.TMX.micrg ~
                         
                         tall.columns.no.SS.final$Texture * tall.columns.no.SS.final$Plant.Influence)

summary.tall.mass <- summary(ANOVA.tall.mass)

summary.tall.mass ## texture and plant influence. No interaction

Tukey.tall.mass <- TukeyHSD(ANOVA.tall.mass)

Tukey.tall.mass

Letters.tall.leach.mass <- multcompLetters4(ANOVA.tall.mass, Tukey.tall.mass)

Letters.tall.leach.mass

##Clay:Plant Clay:No Plant    Sand:Plant Sand:No Plant 
#####"a"           "a"           "b"           "b" 



################Stats on structure TMX mass transport

structure.study.final <- subset(structure.study, Event == 12)

boxplot(structure.study.final$cum.TMX.micrg ~
          
          structure.study.final$Structure, ylab = expression(TMX ~ mu*g),
        
        xlab = "TRT") ## 

hist(structure.study.final$cum.TMX.micrg) ## not normal 

qqnorm(structure.study.final$cum.TMX.micrg) ## not normal

qqline(structure.study.final$cum.TMX.micrg) ## not normal

bartlett.test(structure.study.final$cum.TMX.micrg ~ structure.study.final$Structure) ## meets HOV

fligner.test(structure.study.final$cum.TMX.micrg ~ structure.study.final$Structure) ## meets HOV

levene.test(structure.study.final$cum.TMX.micrg, structure.study.final$Structure) ## meets HOV



##########################just use KW test

kruskal.test(structure.study.final$cum.TMX.micrg ~ structure.study.final$Structure) ## difference!





################Stats on short columns ---tmx mass transport

small.columns ##somehow there are a few 60 cm columns

small.columns.final <- subset(small.columns, Event == 12)

small.columns.final

small.columns.final$log10.TMX.mcrg <- log10(small.columns.final$cum.TMX.micrg) ##log transform

fligner.test(small.columns.final$log10.TMX.mcrg ~
               
               interaction(small.columns.final$Texture, small.columns.final$Plant.Influence))##HOV met

qqnorm(small.columns.final$log10.TMX.mcrg) ## normal

qqline(small.columns.final$log10.TMX.mcrg) ## normal 

hist(small.columns.final$log10.TMX.mcrg)###stick with log10 trans


ANOVA.small.mass <- aov(small.columns.final$log10.TMX.mcrg ~
                          
                          small.columns.final$Texture *
                          
                          small.columns.final$Plant.Influence)

summary.small.columns.mass <- summary(ANOVA.small.mass)

summary.small.columns.mass ## no sig


Tukey.small.columns.mass <- TukeyHSD(ANOVA.small.mass)

Tukey.small.columns.mass


Letters.small.mass <- multcompLetters4(ANOVA.small.mass,Tukey.small.columns.mass)

Letters.small.mass ###all "a"

small.columns.final


### no need to rank transform

#small.columns.final$Rank.TMX.mcrg <- rank(small.columns.final$cum.TMX.micrg)

#hist(small.columns.final$Rank.TMX.mcrg) ## normal

#qqnorm(small.columns.final$Rank.TMX.mcrg) ## 

#qqline(small.columns.final$Rank.TMX.mcrg) ##normal


#fligner.test(small.columns.final$Rank.TMX.mcrg ~

#               interaction(small.columns.final$Texture,

#                           small.columns.final$Plant.Influence))## meets HOV

#ANOVA.small.mass <- aov(small.columns.final$Rank.TMX.mcrg ~

#                         small.columns.final$Texture *

#                        small.columns.final$Plant.Influence)

#summary.small.columns.mass <- summary(ANOVA.small.mass)

#summary.small.columns.mass ## only interaction signficant


#Tukey.small.columns.mass <- TukeyHSD(ANOVA.small.mass)

#Tukey.small.columns.mass


#Letters.small.mass <- multcompLetters4(ANOVA.small.mass,Tukey.small.columns.mass)

#Letters.small.mass ###all "a"

#(small.columns.final)


###################Stats on cumulative leachate volume (mL)

##tall columns no structure first

tall.columns.no.SS.final

tall.columns.no.SS.final$log10.cum.Leach.vol.mL <- log10(tall.columns.no.SS.final$cum.Leach.vol.mL)## log 10

hist(tall.columns.no.SS.final$log10.cum.Leach.vol.mL) ## normal 

qqnorm(tall.columns.no.SS.final$log10.cum.Leach.vol.mL) ##  normal

qqline(tall.columns.no.SS.final$log10.cum.Leach.vol.mL) ## normal



tall.columns.no.SS.final$rank.cum.Leach.vol.mL <- rank(tall.columns.no.SS.final$cum.Leach.vol.mL)## rank transform data to normality

####use log 10

boxplot(tall.columns.no.SS.final$log10.cum.Leach.vol.mL ~
          
          tall.columns.no.SS.final$Texture * tall.columns.no.SS.final$Plant.Influence,
        
        ylab = "Cumualtive Leachate (mL)", xlab = "Trt")## 

fligner.test(tall.columns.no.SS.final$log10.cum.Leach.vol.mL ~
               
               interaction(tall.columns.no.SS.final$Texture, tall.columns.no.SS.final$Plant.Influence)) 
###HOV met


ANOVA.tall.Leach.vol <- aov(tall.columns.no.SS.final$log10.cum.Leach.vol.mL ~
                              
                              tall.columns.no.SS.final$Texture * tall.columns.no.SS.final$Plant.Influence)

summary.tall.Leach.vol <- summary(ANOVA.tall.Leach.vol)

summary.tall.Leach.vol ## texture and plant influence. No interaction 

Tukey.tall.Leach.vol <- TukeyHSD(ANOVA.tall.Leach.vol)

Tukey.tall.Leach.vol

Letters.tall.Leach.vol <- multcompLetters4(ANOVA.tall.Leach.vol, Tukey.tall.Leach.vol)

Letters.tall.Leach.vol

##Clay:No Plant    Clay:Plant Sand:No Plant    Sand:Plant 
#####"a"          "ab"          "ab"           "b" 


########now the structure study***********************************

structure.study.final

boxplot(structure.study.final$cum.Leach.vol.mL ~
          
          structure.study.final$Structure, ylab = "Cumulative Leachate (mL)",
        
        xlab = "TRT") ## 

hist(structure.study.final$cum.Leach.vol.mL) ## not normal 

qqnorm(structure.study.final$cum.Leach.vol.mL) ## not normal

qqline(structure.study.final$cum.Leach.vol.mL) ## not normal

bartlett.test(structure.study.final$cum.Leach.vol.mL ~ structure.study.final$Structure) ## meets HOV

fligner.test(structure.study.final$cum.Leach.vol.mL ~ structure.study.final$Structure) ## meets HOV

levene.test(structure.study.final$cum.Leach.vol.mL, structure.study.final$Structure) ## meets HOV



##########################just use KW test

kruskal.test(structure.study.final$cum.Leach.vol.mL ~ structure.study.final$Structure) ## difference!






#####Now small columns


small.columns.final

small.columns.final$log10.cum.Leach.vol.mL <- log10(small.columns.final$cum.Leach.vol.mL) ##log transform

##try rank transform

small.columns.final$rank.cum.Leach.vol.mL <- rank(small.columns.final$cum.Leach.vol.mL)

fligner.test(small.columns.final$rank.cum.Leach.vol.mL ~
               
               interaction(small.columns.final$Texture, small.columns.final$Plant.Influence))##HOV met

#ANOVA.small.Leach.volume <- aov(small.columns.final$small.columns.final$rank.cum.Leach.vol.mL ~

#                                 small.columns.final$Texture *

#                                small.columns.final$Plant.Influence)

#summary.small.columns.Leach.volume <- summary(ANOVA.small.Leach.volume)

#summary.small.columns.Leach.volume ## no sig





fligner.test(small.columns.final$log10.cum.Leach.vol.mL ~
               
               interaction(small.columns.final$Texture, small.columns.final$Plant.Influence))##HOV met

qqnorm(small.columns.final$log10.cum.Leach.vol.mL) ## normal

qqline(small.columns.final$log10.cum.Leach.vol.mL) ## normalish 

hist(small.columns.final$log10.cum.Leach.vol.mL)###stick with log10 trans


ANOVA.small.Leach.volume <- aov(small.columns.final$log10.cum.Leach.vol.mL ~
                                  
                                  small.columns.final$Texture *
                                  
                                  small.columns.final$Plant.Influence)

summary.small.columns.Leach.volume <- summary(ANOVA.small.Leach.volume)

summary.small.columns.Leach.volume ## no sig


Tukey.small.columns.Leach.volume <- TukeyHSD(ANOVA.small.Leach.volume)

Tukey.small.columns.Leach.volume


Letters.small.Leach.volume <- multcompLetters4(ANOVA.small.Leach.volume,Tukey.small.columns.Leach.volume)

Letters.small.Leach.volume ### all A


###########more prep for plots************














####Tall first#######

tall.columns.no.SS


tall.columns.no.SS.stats.mass <- summarySE(tall.columns.no.SS, measurevar="cum.TMX.micrg" ,
                                           
                                           groupvars=c("Texture","Plant.Influence",
                                                       
                                                       "Time.day"), na.rm =TRUE) ## 
tall.columns.no.SS.stats.mass


tall.columns.no.SS.stats.mass1 <- summarySE(tall.columns.no.SS, measurevar="cum.Leach.vol.mL",
                                            
                                            groupvars=c("Texture","Plant.Influence",
                                                        
                                                        "Time.day"), na.rm =TRUE) ##make df out of cum leachate

#tall.columns.no.SS.stats.mass1 <- subset(tall.columns.no.SS.stats.mass1,

#                                        select = c("cum.Leach.vol.mL", "sd", "se", "ci"))

colnames(tall.columns.no.SS.stats.mass1) <-  c("Texture", "Plant.Influence", "Time.day", "N1",
                                               
                                               "cum.Leach.vol.mL", "sd1",
                                               
                                               "se1", "ci1")  ##preserve imp common variables and change stats                                                        

tall.columns.no.SS.stats.mass1


##now merge dataframes


tall.columns.no.SS.stats.mass <- merge(tall.columns.no.SS.stats.mass,
                                       
                                       tall.columns.no.SS.stats.mass1,
                                       
                                       by = c("Texture", "Plant.Influence", "Time.day"))

tall.columns.no.SS.stats.mass ## good




### need to add column with cum.inf.cm

tall.columns.no.SS.stats.mass$Cum.inf.cm <- ifelse(tall.columns.no.SS.stats.mass$Time.day == 3, 0.9,
                                                   ifelse(tall.columns.no.SS.stats.mass$Time.day == 6, 1.8,
                                                          ifelse(tall.columns.no.SS.stats.mass$Time.day == 9, 2.7,
                                                                 ifelse(tall.columns.no.SS.stats.mass$Time.day == 12, 3.6,
                                                                        ifelse(tall.columns.no.SS.stats.mass$Time.day == 15, 4.5,
                                                                               ifelse(tall.columns.no.SS.stats.mass$Time.day == 18, 5.4,
                                                                                      ifelse(tall.columns.no.SS.stats.mass$Time.day == 21, 6.3,
                                                                                             ifelse(tall.columns.no.SS.stats.mass$Time.day == 24, 7.2,
                                                                                                    ifelse(tall.columns.no.SS.stats.mass$Time.day == 27, 8.1,
                                                                                                           ifelse(tall.columns.no.SS.stats.mass$Time.day == 30, 9,
                                                                                                                  ifelse(tall.columns.no.SS.stats.mass$Time.day == 31, 9.9,
                                                                                                                         ifelse(tall.columns.no.SS.stats.mass$Time.day == 33, 18.9,
                                                                                                                                NA
                                                                                                                         ))))))))))))




tall.columns.no.SS.stats.mass

Letters.tall.leach.mass

##Clay:No Plant    Clay:Plant    Sand:Plant Sand:No Plant 
######"a"           "a"           "b"           "b" 



######## add in sig letters and make all exclusions = NA

tall.columns.no.SS.stats.mass$Letters <-  ifelse(tall.columns.no.SS.stats.mass$Texture == "Clay" & 
                                                   
                                                   tall.columns.no.SS.stats.mass$Time.day == "33", "a",
                                                 
                                                 ifelse (tall.columns.no.SS.stats.mass$Texture == "Sand" & 
                                                           
                                                           tall.columns.no.SS.stats.mass$Time.day == "33", "b",
                                                         
                                                         NA)
)

tall.columns.no.SS.stats.mass

Letters.tall.Leach.vol

##Clay:No Plant    Clay:Plant Sand:No Plant    Sand:Plant 
##"a"          "ab"          "ab"           "b"


tall.columns.no.SS.stats.mass$Vol.diff <-  ifelse(tall.columns.no.SS.stats.mass$Texture == "Clay" & 
                                                    
                                                    tall.columns.no.SS.stats.mass$Time.day == "33" &
                                                    
                                                    tall.columns.no.SS.stats.mass$Plant.Influence == "No Plant", "*",
                                                  
                                                  ifelse (tall.columns.no.SS.stats.mass$Texture == "Sand" & 
                                                            
                                                            tall.columns.no.SS.stats.mass$Time.day == "33" &
                                                            
                                                            tall.columns.no.SS.stats.mass$Plant.Influence == "Plant", "*",
                                                          
                                                          NA)
)

tall.columns.no.SS.stats.mass  #####added significance star for leachate volume


##now make factorial combo (eg. clay No Plant)

tall.columns.no.SS.stats.mass$Combo <- paste(tall.columns.no.SS.stats.mass$Texture,
                                             
                                             tall.columns.no.SS.stats.mass$Plant.Influence)

tall.columns.no.SS.stats.mass ##good





#########now structure study

structure.study

structure.study.stats.mass <- summarySE(structure.study, measurevar="cum.TMX.micrg",
                                        
                                        groupvars=c("Structure",
                                                    
                                                    "Time.day"), na.rm =TRUE) ## stat summary

structure.study.stats.mass


structure.study.stats.mass1 <- summarySE(structure.study, measurevar="cum.Leach.vol.mL",
                                         
                                         groupvars=c("Structure",
                                                     
                                                     "Time.day"), na.rm =TRUE) ##make df out of cum leachate

structure.study.stats.mass1


colnames(structure.study.stats.mass1) <-  c("Structure", "Time.day", "N1",
                                            
                                            "cum.Leach.vol.mL", "sd1",
                                            
                                            "se1", "ci1")  ##preserve imp common variables and change stats                                                        

structure.study.stats.mass1


##now merge dataframes


structure.study.stats.mass <- merge(structure.study.stats.mass,
                                    
                                    structure.study.stats.mass1,
                                    
                                    by = c("Structure", "Time.day"))

structure.study.stats.mass ## good



structure.study.stats.mass$Letters <-  ifelse(structure.study.stats.mass$Structure == "NS" & 
                                                
                                                structure.study.stats.mass$Time.day == "33", "a",
                                              
                                              ifelse (structure.study.stats.mass$Structure == "S" & 
                                                        
                                                        structure.study.stats.mass$Time.day == "33", "b",
                                                      
                                                      NA)
)


structure.study.stats.mass$Vol.diff <- ifelse(structure.study.stats.mass$Structure == "NS" & 
                                                
                                                structure.study.stats.mass$Time.day == "33", "*",
                                              
                                              ifelse (structure.study.stats.mass$Structure == "S" & 
                                                        
                                                        structure.study.stats.mass$Time.day == "33", "*",
                                                      
                                                      NA)
)

structure.study.stats.mass

###added sig for Volume with *

structure.study.stats.mass ####good!

structure.study.stats.mass$Cum.inf.cm <- ifelse(structure.study.stats.mass$Time.day == 3, 0.9,
                                                ifelse(structure.study.stats.mass$Time.day == 6, 1.8,
                                                       ifelse(structure.study.stats.mass$Time.day == 9, 2.7,
                                                              ifelse(structure.study.stats.mass$Time.day == 12, 3.6,
                                                                     ifelse(structure.study.stats.mass$Time.day == 15, 4.5,
                                                                            ifelse(structure.study.stats.mass$Time.day == 18, 5.4,
                                                                                   ifelse(structure.study.stats.mass$Time.day == 21, 6.3,
                                                                                          ifelse(structure.study.stats.mass$Time.day == 24, 7.2,
                                                                                                 ifelse(structure.study.stats.mass$Time.day == 27, 8.1,
                                                                                                        ifelse(structure.study.stats.mass$Time.day == 30, 9,
                                                                                                               ifelse(structure.study.stats.mass$Time.day == 31, 9.9,
                                                                                                                      ifelse(structure.study.stats.mass$Time.day == 33, 18.9,
                                                                                                                             NA
                                                                                                                      ))))))))))))

structure.study.stats.mass



#######now small columns......mass

small.columns

small.columns.stats.mass <- summarySE(small.columns, measurevar="cum.TMX.micrg",
                                      
                                      groupvars=c("Texture", "Plant.Influence",
                                                  
                                                  "Time.day"), na.rm =TRUE)


small.columns.stats.mass

small.columns.stats.mass1 <- summarySE(small.columns,
                                       
                                       measurevar="cum.Leach.vol.mL",
                                       
                                       groupvars=c("Texture","Plant.Influence",
                                                   
                                                   "Time.day"), na.rm =TRUE) ##make df out of cum leachate

small.columns.stats.mass1

colnames(small.columns.stats.mass1) <-  c("Texture", "Plant.Influence", "Time.day", "N1",
                                          
                                          "cum.Leach.vol.mL", "sd1",
                                          
                                          "se1", "ci1")  ##preserve imp common variables and change stats                                                        

small.columns.stats.mass1


##now merge dataframes


small.columns.stats.mass <- merge(small.columns.stats.mass,
                                  
                                  small.columns.stats.mass1,
                                  
                                  by = c("Texture", "Plant.Influence", "Time.day"))

small.columns.stats.mass ## good

small.columns.stats.mass$Cum.inf.cm <- ifelse(small.columns.stats.mass$Time.day == 3, 0.9,
                                              ifelse(small.columns.stats.mass$Time.day == 6, 1.8,
                                                     ifelse(small.columns.stats.mass$Time.day == 9, 2.7,
                                                            ifelse(small.columns.stats.mass$Time.day == 12, 3.6,
                                                                   ifelse(small.columns.stats.mass$Time.day == 15, 4.5,
                                                                          ifelse(small.columns.stats.mass$Time.day == 18, 5.4,
                                                                                 ifelse(small.columns.stats.mass$Time.day == 21, 6.3,
                                                                                        ifelse(small.columns.stats.mass$Time.day == 24, 7.2,
                                                                                               ifelse(small.columns.stats.mass$Time.day == 27, 8.1,
                                                                                                      ifelse(small.columns.stats.mass$Time.day == 30, 9,
                                                                                                             ifelse(small.columns.stats.mass$Time.day == 31, 9.9,
                                                                                                                    ifelse(small.columns.stats.mass$Time.day == 33, 18.9,
                                                                                                                           NA
                                                                                                                    ))))))))))))

small.columns.stats.mass

Letters.small.mass 

##### all = "a"

small.columns.stats.mass$Letters <- ifelse(small.columns.stats.mass$Time.day == "33", "a", NA)

small.columns.stats.mass

##now make factorial combo (eg. clay No Plant)

small.columns.stats.mass$Combo <- paste(small.columns.stats.mass$Texture,
                                        
                                        small.columns.stats.mass$Plant.Influence)

small.columns.stats.mass ##good





#######now Stats on Cumualtive ET***************************************************************************

##Tall columns first (no Structure analysis SS)

tall.columns.no.SS.final ##working in df where final event was isolated (looking at cumualtive ET in mm)


hist(tall.columns.no.SS.final$ET.cum.mm) ## ~ normal 

qqnorm(structure.study.final$ET.cum.mm) ## ~ normal

qqline(structure.study.final$ET.cum.mm) ## ~ normal



##just use same analyses as leachate volume for consistency!  ##

tall.columns.no.SS.final$log10.ET.cum.mm <- log10(tall.columns.no.SS.final$ET.cum.mm)## log 10

hist(tall.columns.no.SS.final$log10.ET.cum.mm) ## normal 

qqnorm(tall.columns.no.SS.final$log10.ET.cum.mm) ##  normal

qqline(tall.columns.no.SS.final$log10.ET.cum.mm) ## normal



boxplot(tall.columns.no.SS.final$log10.ET.cum.mm ~
          
          tall.columns.no.SS.final$Texture * tall.columns.no.SS.final$Plant.Influence,
        
        ylab = "Cumulative ET (mm)", xlab = "Trt")## 

fligner.test(tall.columns.no.SS.final$log10.ET.cum.mm ~
               
               interaction(tall.columns.no.SS.final$Texture, tall.columns.no.SS.final$Plant.Influence)) 
###HOV met


ANOVA.tall.ET <- aov(tall.columns.no.SS.final$log10.ET.cum.mm ~
                       
                       tall.columns.no.SS.final$Texture * tall.columns.no.SS.final$Plant.Influence)

summary.tall.ET <- summary(ANOVA.tall.ET)

summary.tall.ET ## texture and plant influence. No interaction 

Tukey.tall.ET <- TukeyHSD(ANOVA.tall.ET)

Tukey.tall.ET

Letters.tall.ET <- multcompLetters4(ANOVA.tall.ET, Tukey.tall.ET)

Letters.tall.ET

#    Plant No Plant 
#      'a"      "b" 


#Sand:Plant Sand:No Plant    Clay:Plant Clay:No Plant           ##ET Sand  significantly >  ET Clay
#     "a"          "ab"          "ab"           "b" 







########now the structure study***********************************

structure.study.final

boxplot(structure.study.final$ET.cum.mm ~
          
          structure.study.final$Structure, ylab = "Cumulative ET (mm)",
        
        xlab = "TRT") ## 

hist(structure.study.final$ET.cum.mm) ## normal 

qqnorm(structure.study.final$ET.cum.mmL) ## normal

qqline(structure.study.final$ET.cum.mm) ## normal

bartlett.test(structure.study.final$ET.cum.mm ~ structure.study.final$Structure) ## meets HOV

fligner.test(structure.study.final$ET.cum.mm ~ structure.study.final$Structure) ## meets HOV

levene.test(structure.study.final$ET.cum.mm, structure.study.final$Structure) ## meets HOV



##########################just use KW test for consistency

kruskal.test(structure.study.final$ET.cum.mm ~ structure.study.final$Structure) ## difference!







#####Now small columns


small.columns.final

small.columns.final$log10.ET.cum.mm <- log10(small.columns.final$ET.cum.mm) ##log transform






fligner.test(small.columns.final$log10.ET.cum.mm ~
               
               interaction(small.columns.final$Texture, small.columns.final$Plant.Influence))##HOV met

qqnorm(small.columns.final$log10.ET.cum.mm) ## not normal

qqline(small.columns.final$log10.ET.cum.mm) ## not normal 

hist(small.columns.final$log10.ET.cum.mm)###s


ANOVA.small.ET <- aov(small.columns.final$log10.ET.cum.mm ~
                        
                        small.columns.final$Texture *
                        
                        small.columns.final$Plant.Influence)

summary.small.columns.ET <- summary(ANOVA.small.ET)

summary.small.columns.ET ## no sig


Tukey.small.columns.ET <- TukeyHSD(ANOVA.small.ET)

Tukey.small.columns.ET


Letters.small.ET <- multcompLetters4(ANOVA.small.ET,Tukey.small.columns.ET)

Letters.small.ET ### all A



##...maybe try Rank transform anyway

#small.columns.final$rank.ET.cum.mm <- rank(small.columns.final$ET.cum.mm)

#ANOVA.small.ET.rank <- aov(small.columns.final$ET.cum.mm ~

#                                 small.columns.final$Texture *

#                                small.columns.final$Plant.Influence)

#summary.small.columns.ET.rank <- summary(ANOVA.small.ET.rank)

#summary.small.columns.ET.rank ## no sig

####so really doesn't matter what test is used!!! no sig





###########more prep for plots************

tall.columns.no.SS


tall.columns.no.SS.stats.ET <- summarySE(tall.columns.no.SS, measurevar="ET.cum.mm" ,
                                         
                                         groupvars=c("Texture","Plant.Influence",
                                                     
                                                     "Time.day"), na.rm =TRUE) ## 
tall.columns.no.SS.stats.ET


#Leachate_full$Leach.Vol.mL[Leach.Vol.mL == 0] <- NA   ###replace 0 with NA in volumne leachate again

#to be compatible with previous plots in"Leachate_full" analysis

tall.columns.no.SS.stats.ET1 <- summarySE(tall.columns.no.SS, measurevar="cum.Leach.vol.mL",
                                          
                                          groupvars=c("Texture","Plant.Influence",
                                                      
                                                      "Time.day"), na.rm =TRUE) ##make df out of cum leachate


colnames(tall.columns.no.SS.stats.ET1) <-  c("Texture", "Plant.Influence", "Time.day", "N1",
                                             
                                             "cum.Leach.vol.mL", "sd1",
                                             
                                             "se1", "ci1")  ##preserve imp common variables and change stats                                                        

tall.columns.no.SS.stats.ET1


##now merge dataframes


tall.columns.no.SS.stats.ET <- merge(tall.columns.no.SS.stats.ET,
                                     
                                     tall.columns.no.SS.stats.ET1,
                                     
                                     by = c("Texture", "Plant.Influence", "Time.day"))

tall.columns.no.SS.stats.ET ## good



### need to add column with cum.inf.cm

tall.columns.no.SS.stats.ET$Cum.inf.cm <- ifelse(tall.columns.no.SS.stats.ET$Time.day == 3, 0.9,
                                                 ifelse(tall.columns.no.SS.stats.ET$Time.day == 6, 1.8,
                                                        ifelse(tall.columns.no.SS.stats.ET$Time.day == 9, 2.7,
                                                               ifelse(tall.columns.no.SS.stats.ET$Time.day == 12, 3.6,
                                                                      ifelse(tall.columns.no.SS.stats.ET$Time.day == 15, 4.5,
                                                                             ifelse(tall.columns.no.SS.stats.ET$Time.day == 18, 5.4,
                                                                                    ifelse(tall.columns.no.SS.stats.ET$Time.day == 21, 6.3,
                                                                                           ifelse(tall.columns.no.SS.stats.ET$Time.day == 24, 7.2,
                                                                                                  ifelse(tall.columns.no.SS.stats.ET$Time.day == 27, 8.1,
                                                                                                         ifelse(tall.columns.no.SS.stats.ET$Time.day == 30, 9,
                                                                                                                ifelse(tall.columns.no.SS.stats.ET$Time.day == 31, 9.9,
                                                                                                                       ifelse(tall.columns.no.SS.stats.ET$Time.day == 33, 18.9,
                                                                                                                              NA
                                                                                                                       ))))))))))))




tall.columns.no.SS.stats.ET

Letters.tall.ET

#Sand:Plant Sand:No Plant    Clay:Plant Clay:No Plant           ##ET Sand  significantly >  ET Clay
#     "a"          "ab"          "ab"           "b" 

##add sig letters for ETfrom above^^^^

tall.columns.no.SS.stats.ET$Letters <-  ifelse(tall.columns.no.SS.stats.ET$Texture == "Clay" & 
                                                 
                                                 tall.columns.no.SS.stats.ET$Time.day == "33" &
                                                 
                                                 tall.columns.no.SS.stats.ET$Plant.Influence == "No Plant", "b",
                                               
                                               ifelse(tall.columns.no.SS.stats.ET$Texture == "Sand" & 
                                                        
                                                        tall.columns.no.SS.stats.ET$Time.day == "33" &
                                                        
                                                        tall.columns.no.SS.stats.ET$Plant.Influence == "Plant", "a",
                                                      
                                                      ifelse(tall.columns.no.SS.stats.ET$Texture == "Clay" & 
                                                               
                                                               tall.columns.no.SS.stats.ET$Time.day == "33" &
                                                               
                                                               tall.columns.no.SS.stats.ET$Plant.Influence == "Plant", "ab",
                                                             
                                                             ifelse(tall.columns.no.SS.stats.ET$Texture == "Sand" & 
                                                                      
                                                                      tall.columns.no.SS.stats.ET$Time.day == "33" &
                                                                      
                                                                      tall.columns.no.SS.stats.ET$Plant.Influence == "No Plant", "ab",
                                                                    
                                                                    NA)
                                                      )))

tall.columns.no.SS.stats.ET                                               

###add sig to Leachate Volume                                               

Letters.tall.Leach.vol

##Clay:No Plant    Clay:Plant Sand:No Plant    Sand:Plant 
##"a"          "ab"          "ab"           "b"


tall.columns.no.SS.stats.ET$Vol.diff <-  ifelse(tall.columns.no.SS.stats.ET$Texture == "Clay" & 
                                                  
                                                  tall.columns.no.SS.stats.ET$Time.day == "33" &
                                                  
                                                  tall.columns.no.SS.stats.ET$Plant.Influence == "No Plant", "*",
                                                
                                                ifelse (tall.columns.no.SS.stats.ET$Texture == "Sand" & 
                                                          
                                                          tall.columns.no.SS.stats.ET$Time.day == "33" &
                                                          
                                                          tall.columns.no.SS.stats.ET$Plant.Influence == "Plant", "*",
                                                        
                                                        NA)
)


tall.columns.no.SS.stats.ET


##now make factorial combo (eg. clay No Plant)

tall.columns.no.SS.stats.ET$Combo <- paste(tall.columns.no.SS.stats.ET$Texture,
                                           
                                           tall.columns.no.SS.stats.ET$Plant.Influence)

tall.columns.no.SS.stats.ET ##good









#########now structure study

structure.study

structure.study.stats.ET <- summarySE(structure.study, measurevar="ET.cum.mm",
                                      
                                      groupvars=c("Structure",
                                                  
                                                  "Time.day"), na.rm =TRUE) ## stat summary

structure.study.stats.ET


structure.study.stats.ET1 <- summarySE(structure.study, measurevar="cum.Leach.vol.mL",
                                       
                                       groupvars=c("Structure",
                                                   
                                                   "Time.day"), na.rm =TRUE) ##make df out of cum leachate

structure.study.stats.ET1


colnames(structure.study.stats.ET1) <-  c("Structure", "Time.day", "N1",
                                          
                                          "cum.Leach.vol.mL", "sd1",
                                          
                                          "se1", "ci1")  ##preserve imp common variables and change stats                                                        

structure.study.stats.ET1


##now merge dataframes


structure.study.stats.ET <- merge(structure.study.stats.ET,
                                  
                                  structure.study.stats.ET1,
                                  
                                  by = c("Structure", "Time.day"))

structure.study.stats.ET ## good



structure.study.stats.ET$Letters <-  ifelse(structure.study.stats.ET$Structure == "NS" & 
                                              
                                              structure.study.stats.ET$Time.day == "33", "a",
                                            
                                            ifelse (structure.study.stats.ET$Structure == "S" & 
                                                      
                                                      structure.study.stats.ET$Time.day == "33", "b",
                                                    
                                                    NA)
)


structure.study.stats.ET$Vol.diff <- ifelse(structure.study.stats.ET$Structure == "NS" & 
                                              
                                              structure.study.stats.ET$Time.day == "33", "*",
                                            
                                            ifelse (structure.study.stats.ET$Structure == "S" & 
                                                      
                                                      structure.study.stats.ET$Time.day == "33", "*",
                                                    
                                                    NA)
)

structure.study.stats.ET

###added sig for Volume with *

structure.study.stats.ET ####good!

structure.study.stats.ET$Cum.inf.cm <- ifelse(structure.study.stats.ET$Time.day == 3, 0.9,
                                              ifelse(structure.study.stats.ET$Time.day == 6, 1.8,
                                                     ifelse(structure.study.stats.ET$Time.day == 9, 2.7,
                                                            ifelse(structure.study.stats.ET$Time.day == 12, 3.6,
                                                                   ifelse(structure.study.stats.ET$Time.day == 15, 4.5,
                                                                          ifelse(structure.study.stats.ET$Time.day == 18, 5.4,
                                                                                 ifelse(structure.study.stats.ET$Time.day == 21, 6.3,
                                                                                        ifelse(structure.study.stats.ET$Time.day == 24, 7.2,
                                                                                               ifelse(structure.study.stats.ET$Time.day == 27, 8.1,
                                                                                                      ifelse(structure.study.stats.ET$Time.day == 30, 9,
                                                                                                             ifelse(structure.study.stats.ET$Time.day == 31, 9.9,
                                                                                                                    ifelse(structure.study.stats.ET$Time.day == 33, 18.9,
                                                                                                                           NA
                                                                                                                    ))))))))))))

structure.study.stats.ET






#######now small columns

small.columns

small.columns.stats.ET <- summarySE(small.columns, measurevar="ET.cum.mm",
                                    
                                    groupvars=c("Texture", "Plant.Influence",
                                                
                                                "Time.day"), na.rm =TRUE)


small.columns.stats.ET

small.columns.stats.ET1 <- summarySE(small.columns,
                                     
                                     measurevar="cum.Leach.vol.mL",
                                     
                                     groupvars=c("Texture","Plant.Influence",
                                                 
                                                 "Time.day"), na.rm =TRUE) ##make df out of cum leachate

small.columns.stats.ET1

colnames(small.columns.stats.ET1) <-  c("Texture", "Plant.Influence", "Time.day", "N1",
                                        
                                        "cum.Leach.vol.mL", "sd1",
                                        
                                        "se1", "ci1")  ##preserve imp common variables and change stats                                                        

small.columns.stats.ET1


##now merge dataframes


small.columns.stats.ET <- merge(small.columns.stats.ET,
                                
                                small.columns.stats.ET1,
                                
                                by = c("Texture", "Plant.Influence", "Time.day"))

small.columns.stats.ET ## good

small.columns.stats.ET$Cum.inf.cm <- ifelse(small.columns.stats.ET$Time.day == 3, 0.9,
                                            ifelse(small.columns.stats.ET$Time.day == 6, 1.8,
                                                   ifelse(small.columns.stats.ET$Time.day == 9, 2.7,
                                                          ifelse(small.columns.stats.ET$Time.day == 12, 3.6,
                                                                 ifelse(small.columns.stats.ET$Time.day == 15, 4.5,
                                                                        ifelse(small.columns.stats.ET$Time.day == 18, 5.4,
                                                                               ifelse(small.columns.stats.ET$Time.day == 21, 6.3,
                                                                                      ifelse(small.columns.stats.ET$Time.day == 24, 7.2,
                                                                                             ifelse(small.columns.stats.ET$Time.day == 27, 8.1,
                                                                                                    ifelse(small.columns.stats.ET$Time.day == 30, 9,
                                                                                                           ifelse(small.columns.stats.ET$Time.day == 31, 9.9,
                                                                                                                  ifelse(small.columns.stats.ET$Time.day == 33, 18.9,
                                                                                                                         NA
                                                                                                                  ))))))))))))

small.columns.stats.ET

Letters.small.ET 

##### all = "a"

small.columns.stats.ET$Letters <- ifelse(small.columns.stats.ET$Time.day == "33", "a", NA)

small.columns.stats.ET

##now make factorial combo (eg. clay No Plant)

small.columns.stats.ET$Combo <- paste(small.columns.stats.ET$Texture,
                                      
                                      small.columns.stats.ET$Plant.Influence)

small.columns.stats.ET ##good





##make or insert dataframe that has "stairstep" rainfall


cum.rain <- read.csv("C:/Users/Jesse/Desktop/Neonicotinoids_15/Data/Column-Study/Raw-Data-Outputs/cum_rain.csv")
attach(cum.rain)





######done with preplot analysis



#####attach column pics


library(png)

column.60.pic <- readPNG("C:/Users/Jesse/Desktop/Neonicotinoids_15/Pics/Columns-60.png")

column.60.pic.grob <- rasterGrob(column.60.pic, interpolate = T)

column.20.pic <- readPNG("C:/Users/Jesse/Desktop/Neonicotinoids_15/Pics/columns-20.png")

column.20.pic.grob <- rasterGrob(column.20.pic, interpolate = T)

#######


#######attach picture grobs of panel letters (because the grid functions make it a pain in the ass to 
#adjust multipanel figures)

a <- readPNG("C:/Users/Jesse/Desktop/Neonicotinoids_15/Pics/leachate_pic_a.png")

b <- readPNG("C:/Users/Jesse/Desktop/Neonicotinoids_15/Pics/leachate_pic_b.png")

c <- readPNG("C:/Users/Jesse/Desktop/Neonicotinoids_15/Pics/leachate_pic_c.png")

d <- readPNG("C:/Users/Jesse/Desktop/Neonicotinoids_15/Pics/leachate_pic_d.png")

e <- readPNG("C:/Users/Jesse/Desktop/Neonicotinoids_15/Pics/leachate_pic_e.png")

f <- readPNG("C:/Users/Jesse/Desktop/Neonicotinoids_15/Pics/leachate_pic_f.png")


#######################################################################

#######plotting a CUmualtive ET (mm) vs t seems most appropriate
####keep the same the theme with cumulative rainfall on inverse y axis



##actually, shoulf probably convert ET to cm to be consistent with rainfall on iverse y

##tall columns
tall.columns.no.SS.stats.ET$ET.cum.cm <- tall.columns.no.SS.stats.ET$ET.cum.mm/10

tall.columns.no.SS.stats.ET$se.cm <- tall.columns.no.SS.stats.ET$se/10

tall.columns.no.SS.stats.ET$sd.cm <- tall.columns.no.SS.stats.ET$sd/10


##tall columns + structure study
structure.study.stats.ET$ET.cum.cm <- structure.study.stats.ET$ET.cum.mm/10

structure.study.stats.ET$se.cm <- structure.study.stats.ET$se/10

structure.study.stats.ET$sd.cm <- structure.study.stats.ET$sd/10


##small columns

small.columns.stats.ET$ET.cum.cm <- small.columns.stats.ET$ET.cum.mm/10

small.columns.stats.ET$se.cm <- small.columns.stats.ET$se/10

small.columns.stats.ET$sd.cm <- small.columns.stats.ET$sd/10




####lets add in leachate volume as cm also... volume (mL)/area(cm)

###remember that: se <- SE for ET (mm)
#                 se1 <- SE for cumulative leachate volume (mL)
#                 se.cm <- SE for ET (cm)
#                 se1.cm <- SE for cumulative leachate volume (cm)


##tall columns
tall.columns.no.SS.stats.ET$cum.Leach.cm <- tall.columns.no.SS.stats.ET$cum.Leach.vol.mL/area

tall.columns.no.SS.stats.ET$se1.cm <- tall.columns.no.SS.stats.ET$se1/area

##tall columns + structure study
structure.study.stats.ET$cum.Leach.cm <- structure.study.stats.ET$cum.Leach.vol.mL/area

structure.study.stats.ET$se1.cm <- structure.study.stats.ET$se1/area

##small columns

small.columns.stats.ET$cum.Leach.cm <- small.columns.stats.ET$cum.Leach.vol.mL/area

small.columns.stats.ET$se1.cm <- small.columns.stats.ET$se/10



####get cumulative ET and ET rate (daily) for table in paper

#isolate day 33

#tall columns w/oit structureless trt
ET.table.tall <- tall.columns.no.SS.stats.ET[tall.columns.no.SS.stats.ET$Time.day == 33,]

ET.table.tall$ET.rate.cm <- ET.table.tall$ET.cum.cm/33

ET.table.tall$ET.rate.cm.se <- ET.table.tall$se.cm/33

ET.table.tall$ET.rate.cm.sd <- ET.table.tall$sd.cm/33

#export to an excel file

write.xlsx(ET.table.tall , "C:/Users/Jesse/Desktop/Neonicotinoids_15/Data/Column-Study/Raw-Data-Outputs/ET.table.tall.xlsx")




#structure study

structure.study.stats.ET

ET.table.structure <- structure.study.stats.ET[structure.study.stats.ET$Time.day == 33,]  

ET.table.structure$ET.rate.cm <- ET.table.structure$ET.cum.cm/33

ET.table.structure$ET.rate.cm.se <- ET.table.structure$se.cm/33

ET.table.structure$ET.rate.cm.sd <- ET.table.structure$sd.cm/33

#export to an excel file

write.xlsx(ET.table.structure , "C:/Users/Jesse/Desktop/Neonicotinoids_15/Data/Column-Study/Raw-Data-Outputs/ET.table.structure.xlsx")






ET.table.small <- small.columns.stats.ET[small.columns.stats.ET$Time.day == 33,]

ET.table.small$ET.rate.cm <- ET.table.small$ET.cum.cm/33

ET.table.small$ET.rate.cm.se <- ET.table.small$se.cm/33

ET.table.small$ET.rate.cm.sd <- ET.table.small$sd.cm/33

#export to an excel file

write.xlsx(ET.table.small, "C:/Users/Jesse/Desktop/Neonicotinoids_15/Data/Column-Study/Raw-Data-Outputs/ET.table.small.xlsx")



###make a P-ET figure to show how much ET tapped into storage per treatment through time


P.minus.ET.rain <- read.csv("C:/Users/Jesse/Desktop/Neonicotinoids_15/Data/Column-Study/Raw-Data-Outputs/cum.rain.for.P-ET.csv")

attach(P.minus.ET.rain)

P.minus.ET.rain







###first isolate each treatment's et rate in cm day^-1

tall.clay.no.plants.et.rate <- subset(ET.table.tall, Texture == "Clay" &         ##tall clay no plant
                                        
                                        Plant.Influence == "No Plant", ET.rate.cm)
tall.clay.no.plants.et.rate 

tall.clay.plants.et.rate <- subset(ET.table.tall, Texture == "Clay" &         ##tall clay plant
                                        
                                        Plant.Influence == "Plant", ET.rate.cm)
tall.clay.plants.et.rate

tall.sand.no.plants.et.rate <- subset(ET.table.tall, Texture == "Sand" &      ##tall Sand no plant
                                     
                                     Plant.Influence == "No Plant", ET.rate.cm)
tall.sand.no.plants.et.rate

tall.sand.plants.et.rate <- subset(ET.table.tall, Texture == "Sand" &       ##tall Sand plant
                                        
                                        Plant.Influence == "Plant", ET.rate.cm)
tall.sand.plants.et.rate 

tall.structure.et.rate <- subset(ET.table.structure, Structure == "S", ET.rate.cm)  ##tall structure

tall.structure.et.rate

tall.no.structure.et.rate <- subset(ET.table.structure, Structure == "NS", ET.rate.cm)  ##tall  no structure

tall.no.structure.et.rate

###now small columns

small.clay.no.plants.et.rate <- subset(ET.table.small, Texture == "Clay" &         ##tall clay no plant
                                        
                                        Plant.Influence == "No Plant", ET.rate.cm)
small.clay.no.plants.et.rate 

small.clay.plants.et.rate <- subset(ET.table.small, Texture == "Clay" &         ##tall clay plant
                                     
                                     Plant.Influence == "Plant", ET.rate.cm)
small.clay.plants.et.rate

small.sand.no.plants.et.rate <- subset(ET.table.small, Texture == "Sand" &      ##tall Sand no plant
                                        
                                        Plant.Influence == "No Plant", ET.rate.cm)
small.sand.no.plants.et.rate

small.sand.plants.et.rate <- subset(ET.table.small, Texture == "Sand" &       ##tall Sand plant
                                     
                                     Plant.Influence == "Plant", ET.rate.cm)
small.sand.plants.et.rate


####calcualte ET, P-ET from calculated et rate above ^^^ in tables(cm/day)




##assign et rates to treatments

##again we are assuming that the rate stayed constant throughout exp



##nested ifelse()????



P.minus.ET.rain$et_cm.d <- ifelse(P.minus.ET.rain$Size_cm == "60" & P.minus.ET.rain$Texture == "Clay" &         ##tall clay no plant
         
                                  P.minus.ET.rain$Plant.Influence == "No Plant", tall.clay.no.plants.et.rate,
         
                           ifelse(P.minus.ET.rain$Size_cm == "60" & P.minus.ET.rain$Texture == "Clay" &         ##tall clay no plant
                                     
                                  P.minus.ET.rain$Plant.Influence == "Plant", tall.clay.plants.et.rate,
                                  
                           ifelse(P.minus.ET.rain$Size_cm == "60" & P.minus.ET.rain$Texture == "Sand" &         ##tall clay no plant
                                     
                                     P.minus.ET.rain$Plant.Influence == "No Plant", tall.sand.no.plants.et.rate,
                                  
                           ifelse(P.minus.ET.rain$Size_cm == "60" & P.minus.ET.rain$Texture == "Sand" &         ##tall clay no plant
                                    
                                    P.minus.ET.rain$Plant.Influence == "Plant", tall.sand.plants.et.rate,
                                  
                                  
                                  
                           
                                  
                            ifelse(P.minus.ET.rain$Structure == "NS", tall.no.structure.et.rate,
                                  
                                  
                                    
                                    
                                  
                                  
                                  
                           ifelse(P.minus.ET.rain$Size_cm == "20" & P.minus.ET.rain$Texture == "Clay" &         ##tall clay no plant
                                    
                                    P.minus.ET.rain$Plant.Influence == "No Plant", small.clay.no.plants.et.rate,
                                  
                           ifelse(P.minus.ET.rain$Size_cm == "20" & P.minus.ET.rain$Texture == "Clay" &         ##tall clay no plant
                                    
                                    P.minus.ET.rain$Plant.Influence == "Plant", small.clay.plants.et.rate,
                                  
                           ifelse(P.minus.ET.rain$Size_cm == "20" & P.minus.ET.rain$Texture == "Sand" &         ##tall clay no plant
                                    
                                    P.minus.ET.rain$Plant.Influence == "No Plant", small.sand.no.plants.et.rate,
                                  
                           ifelse(P.minus.ET.rain$Size_cm == "20" & P.minus.ET.rain$Texture == "Sand" &         ##tall clay no plant
                                    
                                    P.minus.ET.rain$Plant.Influence == "Plant", small.sand.plants.et.rate,
                                  
                                  NA)
        ))))))))


View(P.minus.ET.rain)


##wtf why is the NS trt not accepting the rate?????? 
###no idea why...

###append/fix manually?


###but this works ... only works if you make NA missing in excel file. Whatever
P.minus.ET.rain[P.minus.ET.rain$Structure == "NS",7] <- tall.no.structure.et.rate






####now calculate ET

P.minus.ET.rain$ET <- as.numeric(P.minus.ET.rain$et_cm.d) * P.minus.ET.rain$time.day 

##had to use as.numeric()....it treats that column as a factor??

P.minus.ET.rain



##P-ET

P.minus.ET.rain$P.ET <- P.minus.ET.rain$P.cm - P.minus.ET.rain$ET

View(P.minus.ET.rain)

##good



####subset for plots

###tall no structure study

tall.columns.P.minus.ET <- P.minus.ET.rain[P.minus.ET.rain$Size_cm == "60", ]

tall.columns.P.minus.ET <- subset(tall.columns.P.minus.ET, Structure == "S" | Texture == "Sand")

tall.columns.P.minus.ET

###probably want to make a texture/plant influence combo...makes graphing easier in ggplot

tall.columns.P.minus.ET$Combo <- paste(tall.columns.P.minus.ET$Texture, tall.columns.P.minus.ET$Plant.Influence)

tall.columns.P.minus.ET


##tall strcuture study included

structure.study.P.minus.ET <- subset(P.minus.ET.rain, Size_cm == "60" &
                            
                            Structure == "S" |
                            
                            Structure == "NS")


structure.study.P.minus.ET <- structure.study.P.minus.ET[structure.study.P.minus.ET$Plant.Influence != "No Plant",]
structure.study.P.minus.ET

##small columns

small.columns.P.minus.ET <- P.minus.ET.rain[P.minus.ET.rain$Size_cm == "20", ]

small.columns.P.minus.ET

###probably want to make a texture/plant influence combo...makes graphing easier in ggplot

small.columns.P.minus.ET$Combo <- paste(small.columns.P.minus.ET$Texture, small.columns.P.minus.ET$Plant.Influence)

small.columns.P.minus.ET








######now start plotting!!!!!!!!
cleanup <- theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background =  element_blank(),
                 axis.line = element_line(color = "black"))

# write a simple function to add footnote...found this online
# make.plot.label <- function(footnoteText =
#                            format(Sys.time(), "%d %b %Y"),
#                          size = 1.5, color = grey(.5))
# {
#   require(grid)
#   pushViewport(viewport())
#   grid.text(label = footnoteText ,
#             x = unit(1,"npc") - unit(2, "mm"),
#             y = unit(2, "mm"),
#             just = c("left", "top"),
#             gp = gpar(cex = size, col = color))
#   popViewport()
# }
# 


# make.plot.label("Footnote", color = "black") ##then use something like this to add figure letter label


###try getting P-ET figure set up############################################################################

tall.columns.P.minus.ET.graph  <- ggplot(tall.columns.P.minus.ET, aes(x = time.day,
                                                                      
                                                                      y = P.ET,
                                                                      
                                                                      ymin = -5,
                                                                      
                                                                      ymax = 12.5
                                                                      
                                                                      
                                                    )) +
  
  geom_line(size = 1, aes(linetype = Combo, color = Combo)) + ###linesize = 1, linetype = 2
  
  scale_linetype_manual(name = "",
                        
                        values = c(2, 1, 2, 1),
                        
                        labels = c("Clay, No Plant", "Clay, Plant",
                                   
                                   "Sand, No Plant", "Sand, Plant")) +
  
  scale_color_manual(name = "",
                     
                     values = c('red3', "red3", 'blue1','blue1'), 
                     
                     labels = c("Clay, No Plant", "Clay, Plant",
                                
                                "Sand, No Plant", "Sand, Plant")) +
  
  ylab("P-ET (cm)") +
  
  # ggtitle("a)") +
  
  #make.plot.label("a)", color = "black")
  
  scale_y_continuous(#breaks = c(-5, 0, 5, 10),
                     
                     expand = c(0,0), trans="reverse", position = "right") +
  
  
  theme_classic() +
  
  theme(plot.margin = unit(c(6,16.5,-33,54), units="points"),
        
        # plot.title = element_text(size = 16, face = "bold", hjust = -0.5, vjust = -0.25),
        
        axis.title.y = element_text(vjust = 0.3, size = 12, face = "bold", hjust = 0.5),
        
        axis.text.y = element_text(size = 15),
        
        panel.border = element_rect(color = "black",
                                    
                                    fill = NA, size = 3),
        
        axis.ticks.y = element_line(size = 1, color = "black"),
        
        legend.margin = margin(0, 7, 10, 1, "point"), ##remove margin on legend
        
        legend.text = element_text( size = 16),  #size = 16
        
        legend.position = "", #c(0.43, 0.15)
        
        legend.key=element_blank(),  ##removed gray border around legend symbols
        
        legend.key.size = unit(10, "mm"),
        
        #legend.background = element_rect(color = "black", size = 1), ##legend box
        
        legend.key.height = unit(7, "mm"),
        
        legend.key.width = unit(9, "mm"),
        
        #                                         legend.box.spacing = unit(10, "mm"),
        
        legend.text.align = 0,
        
        legend.direction = "horizontal")
  
        #axis.title.y = element_text(vjust = 1.5, hjust = 2.5))




tall.columns.P.minus.ET.graph

Leach.tall.TMX.ppb.graph.with.ET <- ggplot(tall.columns.no.SS.stats, aes(Time.day, TMX_PPB,
                                                                      
                                                                      shape = Combo,
                                                                      
                                                                      ymax= TMX_PPB + se,
                                                                      
                                                                      xmax = Time.day)) + ##use colour or shape to add in factor
  
  geom_point(aes( shape = Combo, color = Combo, size = Combo, stroke = 1.5)) +
  
  geom_errorbar(aes(ymin = TMX_PPB - se, ymax = TMX_PPB + se)) +
  
  xlab("Time (d)") +
  
  ylab(expression(TMX ~ (mu * g ~ L^{-1}))) +
  
  scale_shape_manual(name = "",
                     
                     values = c(1, 19, 1, 19),
                     
                     labels = c("Clay, No Plant", "Clay, Plant",
                                
                                "Sand, No Plant", "Sand, Plant")) +
  
  scale_color_manual(name = "",
                     
                     values = c('red3', "red3", 'blue1','blue1'), 
                     
                     labels = c("Clay, No Plant", "Clay, Plant",
                                
                                "Sand, No Plant", "Sand, Plant")) +
  
  scale_size_manual(name = "",
                    
                    values = c(3.5, 3.5, 3.5, 3.5),
                    
                    labels = c("Clay, No Plant", "Clay, Plant",
                               
                               "Sand, No Plant", "Sand, Plant")) +
  
  geom_text(aes(label = Letters, y = TMX_PPB, fontface = "italic"),
            
            position = position_dodge(width=2.5),
            
            hjust= 2.2, vjust = -1) +  #####add in signficance letter--adjusted for errror bars
  
  
  ##try to add a series of short line segments to make a legend for top graph
  
  geom_segment(aes(x= 0, y = 18, xend = 5, yend = 18), color = "red3", linetype = 2, size = 1) + #Clay, No Plant

  geom_segment(aes(x= 0, y = 15, xend = 5, yend = 15), color = "red3", linetype = 1, size = 1) + #Clay, Plant

  geom_segment(aes(x= 0, y = 12, xend = 5, yend = 12), color = "blue1", linetype = 2, size = 1) + #Sand, No Plant

  geom_segment(aes(x= 0, y = 9, xend = 5, yend = 9), color = "blue1", linetype = 1, size = 1) + #Sand, Plant


  
  
  cleanup +
  
  theme(plot.margin = unit(c(-1,46.5,1,2), units="points"), ##good legend size
        
        #                                       panel.border = element_rect(color = "black",
        
        #                                                                 fill = NA, size = 1.75),
        
        axis.text.y = element_text(size = 16),  ## change font size of y axis
        
        axis.title.y = element_text(size = 16),
        
        axis.text.x = element_text(size = 16),  ## change font size of y axis
        
        axis.title.x = element_text(size = 16, face = "bold"),
        
        legend.margin = margin(0, 7, 10, 1, "point"), ##remove margin on legend
        
        legend.text = element_text( size = 15),
        
         legend.position = c(0.4, 0.75), #c(0.43, 0.95)
        
        legend.key=element_blank(),  ##removed gray border around legend symbols
        
        legend.key.size = unit(7, "mm"),
        
        #legend.background = element_rect(color = "black", size = 1), #box around legend
        
        legend.key.height = unit(5, "mm"),
        
        legend.key.width = unit(5, "mm"),
        
        #                                         legend.box.spacing = unit(10, "mm"),
        
        legend.text.align = 0,
        
        legend.direction = "vertical",
        
        axis.line.x = element_line(color = "black", size = 1.75),
        
        
        axis.ticks.x= element_line(size = 1.75, color = "black"),
        
        axis.line.y = element_line(color = "black", size = 1.75),
        
        axis.ticks.y = element_line(size = 1.75, color = "black"),
        
        axis.ticks.length = unit(2, "mm")) +
  
  ## add secondary axis to bottom plot and remove ticks
  
  scale_y_continuous(sec.axis = sec_axis(~. * 1, labels = NULL,
                                         
                                         breaks = c(18)))

Leach.tall.TMX.ppb.graph.with.ET




Leach.tall.TMX.ppb.graph.comb.with.ET <- grid.arrange(tall.columns.P.minus.ET.graph, Leach.tall.TMX.ppb.graph.with.ET,
                                                   
                                                   heights = c(4/11, 7/11))

Leach.tall.TMX.ppb.graph.comb.with.ET


Leach.tall.TMX.ppb.graph.with.ET
#pdf("C:/Users/Jesse/Desktop/Neonicotinoids_15/Data/Column-Study/Raw-Data-Outputs/Leach.tall.TMX.ppb.graph.with.ET.pdf")
dev.off()

#ggsave(file = "C:/Users/Jesse/Desktop/Neonicotinoids_15/Data/Column-Study/Raw-Data-Outputs/Leach.tall.TMX.ppb.graph.comb.with.ET.pdf", Leach.tall.TMX.ppb.graph.comb.with.ET)


#########work on structure study


structure.study.P.minus.ET.graph <- ggplot(structure.study.P.minus.ET, aes(x = time.day,
                                                                                
                                                                                y = P.ET,
                                                                                
                                                                                ymin = -6,
                                                                                
                                                                                ymax = 12.5
                                                                                
                                                                                
)) +
  
  geom_line(size = 1, aes(linetype = Structure, color = Structure)) + ###linesize = 1, linetype = 2
  
  scale_linetype_manual(name = "",
                        
                        values = c(2, 1),
                        
                        labels = c("Unstructured Clay", "Structured Clay")) +
  
  scale_color_manual(name = "",
                     
                     values = c('black', "black"),
                     
                     labels = c("Unstructured Clay", "Structured Clay")) +
  
  ylab("P-ET (cm)") +
  
  scale_y_continuous(#breaks = c(-5, 0, 5, 10),
    
    expand = c(0,0), trans="reverse", position = "right") +
  
  
  theme_classic() +
  
  theme(plot.margin = unit(c(6,16.5,-33,54), units="points"),
        
        axis.title.y = element_text(vjust = 0.3, size = 12, face = "bold", hjust = 0.5),
        
        axis.text.y = element_text(size = 15),
        
        panel.border = element_rect(color = "black",
                                    
                                    fill = NA, size = 3),
        
        axis.ticks.y = element_line(size = 1, color = "black"),
        
        legend.margin = margin(0, 7, 10, 1, "point"), ##remove margin on legend
        
        legend.text = element_text( size = 16),  #size = 16
        
        legend.position = "", #c(0.43, 0.15)
        
        legend.key=element_blank(),  ##removed gray border around legend symbols
        
        legend.key.size = unit(10, "mm"),
        
        legend.background = element_rect(color = "black", size = 1),
        
        legend.key.height = unit(7, "mm"),
        
        legend.key.width = unit(7, "mm"),
        
        #                                         legend.box.spacing = unit(10, "mm"),
        
        legend.text.align = 0,
        
        legend.direction = "horizontal")

#axis.title.y = element_text(vjust = 1.5, hjust = 2.5))



structure.study.P.minus.ET.graph


structure.study.TMX.ppb.graph.with.ET <- ggplot(structure.study.stats, aes(Time.day, TMX_PPB,
                                                                        
                                                                        shape = Structure,
                                                                        
                                                                        ymax= TMX_PPB + se,
                                                                        
                                                                        xmax = Time.day)) + ##use colour or shape to add in factor
  
  geom_point(aes(shape = Structure, color = Structure, size = Structure, stroke = 1.5)) +
  
  geom_errorbar(aes(ymin = TMX_PPB - se, ymax = TMX_PPB + se)) +
  
  xlab("") +
  
  ylab(expression(TMX ~ (mu * g ~ L^{-1}))) +
  
  scale_shape_manual(name = "",
                     
                     values = c(1, 19),
                     
                     labels = c("Unstructured Clay", "Structured Clay")) +
  
  scale_color_manual(name = "",
                     
                     values = c('black', "black"), 
                     
                     labels = c("Unstructured Clay", "Structured Clay")) +
  
  scale_size_manual(name = "",
                    
                    values = c(3.5, 3.5),
                    
                    labels = c("Unstructured Clay", "Structured Clay")) +
  
  geom_text(aes(label = Letters, y = TMX_PPB, fontface = "italic"),
            
            position = position_dodge(width=0.9),
            
            hjust= 2.25, vjust = 0.1) +  #####add in signficance letter--adjusted for errror bars
  
  ##try to add a series of short line segments to make a legend for top graph
  
  geom_segment(aes(x= 0, y = 15.5, xend = 5, yend = 15.5), color = "black", linetype = 2, size = 1) + #Clay, No Plant
  
  geom_segment(aes(x= 0, y = 11.5, xend = 5, yend = 11.5), color = "black", linetype = 1, size = 1) + #Clay, Plant
  
  
  
  cleanup +
  
  theme(plot.margin = unit(c(-1,46.5,1,2), units="points"), ##good legend size
        
        #                                       panel.border = element_rect(color = "black",
        
        #                                                                 fill = NA, size = 1.75),
        
        axis.text.y = element_text(size = 16),  ## change font size of y axis
        
        axis.title.y = element_text(size = 16),
        
        axis.text.x = element_text(size = 16),  ## change font size of y axis
        
        axis.title.x = element_text(size = 16),
        
        legend.margin = margin(0, 7, 10, 0, "point"), ##remove margin on legend
        
        legend.text = element_text(size = 15),
        
        legend.position = c(0.45, 0.75),
        
        legend.key=element_blank(),  ##removed gray border around legend symbols
        
        legend.key.size = unit(7, "mm"),
        
        #legend.background = element_rect(color = "black", size = 1.25), #legend.box
        
        legend.key.height = unit(7, "mm"),
        
        legend.key.width = unit(7, "mm"),
        
        #                                         legend.box.spacing = unit(10, "mm"),
        
        legend.text.align = 0,
        
        axis.line.x = element_line(color = "black", size = 1.75),
        
        
        axis.ticks.x= element_line(size = 1.75, color = "black"),
        
        axis.line.y = element_line(color = "black", size = 1.75),
        
        axis.ticks.y = element_line(size = 1.75, color = "black"),
        
        axis.ticks.length = unit(2, "mm")) +
  
  ## add secondary axis to bottom plot and remove ticks
  
  scale_y_continuous(sec.axis = sec_axis(~. * 1, labels = NULL,
                                         
                                         breaks = c(18)))

structure.study.TMX.ppb.graph.with.ET


structure.study.TMX.ppb.graph.comb.with.ET <- grid.arrange(structure.study.P.minus.ET.graph, structure.study.TMX.ppb.graph.with.ET,
                                                        
                                                        heights = c(4/11, 7/11))

structure.study.TMX.ppb.graph.comb.with.ET


structure.study.TMX.ppb.graph.with.ET
#pdf("C:/Users/Jesse/Desktop/Neonicotinoids_15/Data/Column-Study/Raw-Data-Outputs/structure.study.TMX.ppb.graph.with.ET.pdf")
dev.off()

#ggsave(file = "C:/Users/Jesse/Desktop/Neonicotinoids_15/Data/Column-Study/Raw-Data-Outputs/structure.study.TMX.ppb.graph.comb.with.ET.pdf", structure.study.TMX.ppb.graph.comb.with.ET)



#################now small columns


Leach.small.P.minus.ET.graph <-  ggplot(small.columns.P.minus.ET, aes(x = time.day,
                                                                     
                                                                     y = P.ET,
                                                                     
                                                                     ymin = -5,
                                                                     
                                                                     ymax = 12.5
                                                                     
                                                                     
)) +
  
  geom_line(size = 1, aes(linetype = Combo, color = Combo)) + ###linesize = 1, linetype = 2
  
  scale_linetype_manual(name = "",
                        
                        values = c(2, 1, 2, 1),
                        
                        labels = c("Clay, No Plant", "Clay, Plant",
                                   
                                   "Sand, No Plant", "Sand, Plant")) +
  
  scale_color_manual(name = "",
                     
                     values = c('red3', "red3", 'blue1','blue1'), 
                     
                     labels = c("Clay, No Plant", "Clay, Plant",
                                
                                "Sand, No Plant", "Sand, Plant")) +
  
  ylab("P-ET (cm)") +
  
  scale_y_continuous(#breaks = c(-5, 0, 5, 10),
    
    expand = c(0,0), trans="reverse", position = "right") +
  
  
  theme_classic() +
  
  theme(plot.margin = unit(c(6,16.5,-33,54), units="points"),
        
        axis.title.y = element_text(vjust = 0.3, size = 12, face = "bold", hjust = 0.5),
        
        axis.text.y = element_text(size = 15),
        
        panel.border = element_rect(color = "black",
                                    
                                    fill = NA, size = 3),
        
        axis.ticks.y = element_line(size = 1, color = "black"),
        
        legend.margin = margin(0, 7, 10, 1, "point"), ##remove margin on legend
        
        legend.text = element_text( size = 16),  #size = 16
        
        legend.position = "", #c(0.43, 0.15)
        
        legend.key=element_blank(),  ##removed gray border around legend symbols
        
        legend.key.size = unit(7, "mm"),
        
        legend.background = element_rect(color = "black", size = 1),
        
        legend.key.height = unit(5, "mm"),
        
        legend.key.width = unit(5, "mm"),
        
        #                                         legend.box.spacing = unit(10, "mm"),
        
        legend.text.align = 0,
        
        legend.direction = "horizontal")

#axis.title.y = element_text(vjust = 1.5, hjust = 2.5))




Leach.small.P.minus.ET.graph


Leach.small.TMX.ppb.graph.with.ET <- ggplot(small.columns.stats, aes(Time.day, TMX_PPB,
                                                                  
                                                                  shape = Combo,
                                                                  
                                                                  ymax= TMX_PPB + se,
                                                                  
                                                                  xmax = Time.day)) + ##use colour or shape to add in factor
  
  geom_point(aes(shape = Combo, color = Combo, size = Combo, stroke = 1.5)) +
  
  geom_errorbar(aes(ymin = TMX_PPB - se, ymax = TMX_PPB + se)) +
  
  xlab("") + #"Time (d)"
  
  ylab(expression(TMX ~ (mu * g ~ L^{-1}))) +
  
  scale_shape_manual(name = "",
                     
                     values = c(1, 19, 1, 19),
                     
                     labels = c("Clay, No Plant", "Clay, Plant",
                                
                                "Sand, No Plant", "Sand, Plant")) +
  
  scale_color_manual(name = "",
                     
                     values = c('red3', "red3", 'blue1','blue1'), 
                     
                     labels = c("Clay, No Plant", "Clay, Plant",
                                
                                "Sand, No Plant", "Sand, Plant")) +
  
  scale_size_manual(name = "",
                    
                    values = c(3.5, 3.5, 3.5, 3.5),
                    
                    labels = c("Clay, No Plant", "Clay, Plant",
                               
                               "Sand, No Plant", "Sand, Plant")) +
  
  geom_text(aes(label = Letters, y = TMX_PPB, fontface = "italic"),
            
            position = position_dodge(width=0.9),
            
            hjust= 2.5) +  #####add in signficance letter--adjusted for errror bars
  
  geom_segment(aes(x= 5, y = 17.5, xend = 25, yend = 17.5), color = "black", linetype = 2, size = 1.5,
               
               data = small.columns.stats) +
  
  #geom_abline(intercept = 17.5, slope = 0, size = 1.5, color = "firebrick1", linetype = 5) +
  
  annotate("text", x = 15, y = 25, label = "Aquatic Life Benchmark") + ##add benchamrk
  
  ##try to make this x axis line up with others in this column...add whiteline with similar length as above
  
  geom_segment(aes(x= 0, y = 50, xend = 5, yend = 50), color = "white", linetype = 2, size = 1) + #Clay, No Plant
  
  
  cleanup +
  
  theme(plot.margin = unit(c(-1,46.5,1,2), units="points"), ##good legend size
        
        #                                       panel.border = element_rect(color = "black",
        
        #                                                                 fill = NA, size = 1.75),
        
        axis.text.y = element_text(size = 16),  ## change font size of y axis
        
        axis.title.y = element_text(size = 16),
        
        axis.text.x = element_text(size = 16),  ## change font size of y axis
        
        axis.title.x = element_text(size = 16, face = "bold"),
        
        legend.margin = margin(0, 5, 10, 1, "point"), ##remove margin on legend
        
        legend.text = element_text( size = 16),
        
        legend.position = "none",
        
        legend.key=element_blank(),  ##removed gray border around legend symbols
        
        legend.key.size = unit(10, "mm"),
        
        legend.background = element_rect(color = "black", size = 1.25),
        
        legend.key.height = unit(7, "mm"),
        
        legend.key.width = unit(7, "mm"),
        
        #                                         legend.box.spacing = unit(10, "mm"),
        
        legend.text.align = 0,
        
        axis.line.x = element_line(color = "black", size = 1.75),
        
        
        axis.ticks.x= element_line(size = 1.75, color = "black"),
        
        axis.line.y = element_line(color = "black", size = 1.75),
        
        axis.ticks.y = element_line(size = 1.75, color = "black"),
        
        axis.ticks.length = unit(2, "mm")) +
  
  ## add secondary axis to bottom plot and remove ticks
  
  scale_y_continuous(sec.axis = sec_axis(~. * 1, labels = NULL,
                                         
                                         breaks = c(90))) +
  
  scale_x_continuous(breaks = c(0, 10, 20, 30))


Leach.small.TMX.ppb.graph.with.ET


Leach.small.TMX.ppb.graph.comb.with.ET <- grid.arrange(Leach.small.P.minus.ET.graph, Leach.small.TMX.ppb.graph.with.ET,
                                                    
                                                    heights = c(4/11, 7/11))

Leach.small.TMX.ppb.graph.comb.with.ET

Leach.small.TMX.ppb.graph.with.ET
pdf("C:/Users/Jesse/Desktop/Neonicotinoids_15/Data/Column-Study/Raw-Data-Outputs/Leach.small.TMX.ppb.graph.with.ET.pdf")
dev.off()

#ggsave(file = "C:/Users/Jesse/Desktop/Neonicotinoids_15/Data/Column-Study/Raw-Data-Outputs/Leach.small.TMX.ppb.graph.comb.with.ET.pdf", Leach.small.TMX.ppb.graph.comb.with.ET)


#####################now TMX mass to be combined with P-ET figures

tall.columns.no.SS.stats.mass

Leach.tall.TMX.ppb.mass.graph.grid.with.ET <- ggplot(tall.columns.no.SS.stats.mass, aes(cum.Leach.vol.mL, cum.TMX.micrg,
                                                                                
                                                                                shape = Combo,
                                                                                
                                                                                ymax= cum.TMX.micrg + se,
                                                                                
                                                                                xmax = cum.Leach.vol.mL)) + ##use colour or shape to add in factor
  
  geom_errorbar(aes(ymin = cum.TMX.micrg - se, ymax = cum.TMX.micrg + se),
                
                colour="black", size = 0.75) + ### caps won't show
  
  geom_errorbarh(aes(xmin = cum.Leach.vol.mL - se1, xmax = cum.Leach.vol.mL + se1),
                 
                 colour="black", size = 0.75) +
  
  geom_point(aes(shape = Combo, color = Combo, size = Combo, stroke = 1.5)) +
  
  xlab("Volume Leached (mL)") +
  
  ylab(TMX ~ Leached ~ (mu*g)) +
  
  scale_shape_manual(name = "",
                     
                     values = c(1, 19, 1, 19),
                     
                     labels = c("Clay, No Plant", "Clay, Plant",
                                
                                "Sand, No Plant", "Sand, Plant")) +
  
  scale_color_manual(name = "",
                     
                     values = c('red3', "red3", 'blue1','blue1'), 
                     
                     labels = c("Clay, No Plant", "Clay, Plant",
                                
                                "Sand, No Plant", "Sand, Plant")) +
  
  scale_size_manual(name = "",
                    
                    values = c(3.5, 3.5, 3.5, 3.5),
                    
                    labels = c("Clay, No Plant", "Clay, Plant",
                               
                               "Sand, No Plant", "Sand, Plant")) +
  
  geom_text(aes(label = Letters, y = cum.TMX.micrg, fontface = "italic"),
            
            
            
            hjust= -1.25, vjust = -0.5) +  #####add in signficance letter--adjusted for errror bars
  
  geom_text(aes(label = Vol.diff, y = cum.TMX.micrg),
            
            position = position_dodge(width =0.9),
            
            hjust= -2.25, vjust = 0.5, size = 6, fontface = "bold") +
  
  annotation_custom(column.60.pic.grob, xmin = -200, xmax = 800, ymin=2.5, ymax=52.5) + ##add in column pic
  
  # annotation_raster(column.60.pic, ymin = 10, ymax= 50, xmin = 500,xmax = 1300) +
  
  cleanup +
  
  theme(plot.margin = unit(c(5,36.5,1,27), units="points"), ##good legend size
        
        panel.border = element_rect(color = "black",
                                    
                                    fill = NA, size = 1.75),
        
        axis.text.y = element_text(size = 16),  ## change font size of y axis
        
        axis.title.y = element_text(size = 16),
        
        axis.text.x = element_text(size = 16),  ## change font size of y axis
        
        axis.title.x = element_text(size = 16,  face = "bold"),
        
        legend.margin = margin(0, 5, 10, 1, "point"), ##remove margin on legend
        
        legend.text = element_text( size = 16),
        
        legend.position = "none",
        
        legend.key=element_blank(),  ##removed gray border around legend symbols
        
        legend.key.size = unit(10, "mm"),
        
        legend.background = element_rect(color = "black", size = 1.25),
        
        legend.key.height = unit(7, "mm"),
        
        legend.key.width = unit(7, "mm"),
        
        #                                         legend.box.spacing = unit(10, "mm"),
        
        legend.text.align = 0,
        
        axis.line.x = element_line(color = "black", size = 1.75),
        
        
        axis.ticks.x= element_line(size = 1.75, color = "black"),
        
        axis.line.y = element_line(color = "black", size = 1.75),
        
        axis.ticks.y = element_line(size = 1.75, color = "black"),
        
        axis.ticks.length = unit(2, "mm")) +
  
  scale_x_continuous(limits = c(0, 4000))






Leach.tall.TMX.ppb.mass.graph.grid.with.ET







######now structure study***************************

structure.study.TMX.mass.graph.grid.with.ET <- ggplot(structure.study.stats.mass, aes(cum.Leach.vol.mL, cum.TMX.micrg,
                                                                              
                                                                              shape = Structure,
                                                                              
                                                                              ymax= cum.TMX.micrg + se,
                                                                              
                                                                              xmax = cum.Leach.vol.mL + se1)) + ##use colour or shape to add in factor
  
  geom_errorbar(aes(ymin = cum.TMX.micrg - se, ymax = cum.TMX.micrg + se),
                
                size = 0.75) +
  
  geom_errorbarh(aes(xmin = cum.Leach.vol.mL - se1, xmax = cum.Leach.vol.mL + se1),
                 
                 size = 0.75) +
  
  geom_point(aes(shape = Structure, color = Structure, size = Structure, stroke = 1.5)) +
  
  
  xlab("") +
  
  ylab(TMX ~ Leached ~ (mu*g)) +
  
  scale_shape_manual(name = "",
                     
                     values = c(1, 19),
                     
                     labels = c("Untructured Clay", "Structured Clay")) +
  
  scale_color_manual(name = "",
                     
                     values = c('black', "black"), 
                     
                     labels = c("Unstructured Clay", "Structured Clay")) +
  
  scale_size_manual(name = "",
                    
                    values = c(3.5, 3.5),
                    
                    labels = c("Unstructured Clay", "Structured Clay")) +
  
  geom_text(aes(label = Letters, y = cum.TMX.micrg, fontface = "italic"),
            
            position = position_dodge(width=0.9),
            
            hjust= -1.8, vjust =-0.5) +  #####add in signficance letter--adjusted for errror bars
  
  geom_text(aes(label = Vol.diff, y = cum.TMX.micrg),
            
            position = position_dodge(width =0.9),
            
            hjust= -2.75, vjust = 1, size = 6, fontface = "bold") +
  
  annotation_custom(column.60.pic.grob, xmin = -200, xmax = 800, ymin=-1, ymax=49) + ##add in column pic
  
  cleanup +
  
  theme(plot.margin = unit(c(5,36.5,1,27), units="points"), ##good legend size
        
        panel.border = element_rect(color = "black",
                                    
                                    fill = NA, size = 1.75),
        
        axis.text.y = element_text(size = 16),  ## change font size of y axis
        
        axis.title.y = element_text(size = 16),
        
        axis.text.x = element_text(size = 16),  ## change font size of y axis
        
        axis.title.x = element_text(size = 16),
        
        legend.margin = margin(0, 5, 10, 1, "point"), ##remove margin on legend
        
        legend.text = element_text(size = 16),
        
        legend.position = "none",
        
        legend.key=element_blank(),  ##removed gray border around legend symbols
        
        legend.key.size = unit(10, "mm"),
        
        legend.background = element_rect(color = "black", size = 1.25),
        
        legend.key.height = unit(7, "mm"),
        
        legend.key.width = unit(7, "mm"),
        
        #                                         legend.box.spacing = unit(10, "mm"),
        
        legend.text.align = 0,
        
        axis.line.x = element_line(color = "black", size = 1.75),
        
        
        axis.ticks.x= element_line(size = 1.75, color = "black"),
        
        axis.line.y = element_line(color = "black", size = 1.75),
        
        axis.ticks.y = element_line(size = 1.75, color = "black"),
        
        axis.ticks.length = unit(2, "mm")) +
  
  scale_x_continuous( limits = c(0,4000))


structure.study.TMX.mass.graph.grid.with.ET




########Small columns    TMX mass***************************************

small.columns.stats.mass


Leach.small.TMX.mass.graph.grid.with.ET <- ggplot(small.columns.stats.mass, aes(cum.Leach.vol.mL, cum.TMX.micrg,
                                                                        
                                                                        shape = Combo,
                                                                        
                                                                        ymax= cum.TMX.micrg + se,
                                                                        
                                                                        xmax = cum.Leach.vol.mL + se1)) +
  
  geom_errorbar(aes(ymin = cum.TMX.micrg - se, ymax = cum.TMX.micrg + se),
                
                size = 0.75) +
  
  geom_errorbarh(aes(xmin = cum.Leach.vol.mL - se1, xmax = cum.Leach.vol.mL + se1),
                 
                 size = 0.75) +
  
  geom_point(aes(shape = Combo, color = Combo, size = Combo, stroke = 1.5)) +
  
  xlab("") + #"Volume Leached (mL)"
  
  ylab(TMX ~ Leached ~ (mu*g)) +
  
  scale_shape_manual(name = "",
                     
                     values = c(1, 19, 1, 19),
                     
                     labels = c("Clay, No Plant", "Clay, Plant",
                                
                                "Sand, No Plant", "Sand, Plant")) +
  
  scale_color_manual(name = "",
                     
                     values = c('red3', "red3", 'blue1','blue1'), 
                     
                     labels = c("Clay, No Plant", "Clay, Plant",
                                
                                "Sand, No Plant", "Sand, Plant")) +
  
  scale_size_manual(name = "",
                    
                    values = c(3.5, 3.5, 3.5, 3.5),
                    
                    labels = c("Clay, No Plant", "Clay, Plant",
                               
                               "Sand, No Plant", "Sand, Plant")) +
  
  geom_text(aes(label = Letters, y = cum.TMX.micrg, fontface = "italic"),
            
            position = position_dodge(width=0.9),
            
            hjust= -1.4, vjust = -0.5) +  #####add in signficance letter--adjusted for errror bars
  
  #geom_abline(intercept = 17.5, slope = 0, size = 1.5, color = "firebrick1", linetype = 5) +
  
  annotation_custom(column.20.pic.grob, xmin = -200, xmax = 800, ymin = 45, ymax = 250) + ##add in column pic
  
  cleanup +
  
  theme(plot.margin = unit(c(5,36.5,1,21), units="points"), ##good legend size
        
        panel.border = element_rect(color = "black",
                                    
                                    fill = NA, size = 1.75),
        
        axis.text.y = element_text(size = 16),  ## change font size of y axis
        
        axis.title.y = element_text(size = 16),
        
        axis.text.x = element_text(size = 16),  ## change font size of y axis
        
        axis.title.x = element_text(size = 16, face = "bold"),
        
        legend.margin = margin(0, 5, 10, 1, "point"), ##remove margin on legend
        
        legend.text = element_text( size = 16),
        
        legend.position = "none",
        
        legend.key=element_blank(),  ##removed gray border around legend symbols
        
        legend.key.size = unit(10, "mm"),
        
        legend.background = element_rect(color = "black", size = 1.25),
        
        legend.key.height = unit(7, "mm"),
        
        legend.key.width = unit(7, "mm"),
        
        #                                         legend.box.spacing = unit(10, "mm"),
        
        legend.text.align = 0,
        
        axis.line.x = element_line(color = "black", size = 1.75),
        
        
        axis.ticks.x= element_line(size = 1.75, color = "black"),
        
        axis.line.y = element_line(color = "black", size = 1.75),
        
        axis.ticks.y = element_line(size = 1.75, color = "black"),
        
        axis.ticks.length = unit(2, "mm")) +
  
  scale_x_continuous(limits = c(0, 4000))




Leach.small.TMX.mass.graph.grid.with.ET





######arrange

# Leach.tall.TMX.ppb.graph.comb.grid
# 
# structure.study.TMX.ppb.graph.comb.grid
# 
# Leach.small.TMX.ppb.graph.comb.grid
# 
# Leach.tall.TMX.ppb.mass.graph.grid
# 
# structure.study.TMX.mass.graph.grid
# 
# Leach.small.TMX.mass.graph.grid
# 
# Leach.small.TMX.mass.graph.grid


TMX.mass.and.conc.Leachate.grid.with.ET <- grid.arrange(Leach.small.TMX.ppb.graph.comb.with.ET,###conc 
                                                
                                                structure.study.TMX.ppb.graph.comb.with.ET,##conc

                                                Leach.tall.TMX.ppb.graph.comb.with.ET,##conc
                                                
                                                Leach.small.TMX.mass.graph.grid.with.ET, ##mass

                                                structure.study.TMX.mass.graph.grid.with.ET, ##mass
                                                
                                                Leach.tall.TMX.ppb.mass.graph.grid.with.ET, ##mass

                                                layout_matrix = rbind(c(1, 4),
                                                                      c(2, 5),
                                                                      c(3, 6)
                                                )) 

TMX.mass.and.conc.Leachate.grid.with.ET

#ggsave(file = "C:/Users/Jesse/Desktop/Neonicotinoids_15/Data/Column-Study/Raw-Data-Outputs/TMX.mass.and.conc.Leachate.grid.with.ET.pdf", TMX.mass.and.conc.Leachate.grid.with.ET)
TMX.mass.and.conc.Leachate.grid.with.ET
dev.off()









































####TMX_ppb- Tall columns#################################################################################



tall.columns.no.SS.stats

Leach.tall.TMX.ppb.inf <- ggplot(cum.rain, aes(x = Day, y = Cumulative.rainfall.cm,
                                                               
                                                               ymin = 0,
                                                               
                                                               ymax = Cumulative.rainfall.cm,
                                                               
                                                               xmax = Day)) +
                                 geom_line(size = 1, linetype = 2) +
                                 
                                 scale_y_continuous(limits=c(25,0),
                                                    
                                                    expand=c(0,0), trans="reverse", position = "right") +
                                 theme_classic() +
                                 
                                 theme(plot.margin = unit(c(6,24.75,-37,50), units="points"),
                                       
                                       axis.title.y = element_text(vjust = 0.3, size = 14, face = "bold"),
                                    
                                       axis.text.y = element_text(size = 16),
                                          
                                       panel.border = element_rect(color = "black",
                                                                   
                                                                   fill = NA, size = 3),
                                       
                                       axis.ticks.y = element_line(size = 1, color = "black")) +
                                
                                 ylab("T. Rainfall (cm)") +
  
                                 theme(axis.title.y = element_text(vjust = 0.25, hjust = -.7))

                                

Leach.tall.TMX.ppb.inf


Leach.tall.TMX.ppb.graph <- ggplot(tall.columns.no.SS.stats, aes(Time.day, TMX_PPB,
                                                                 
                                                                 shape = Combo,
                                                                 
                                                                 ymax= TMX_PPB + se,
                                                                 
                                                                 xmax = Time.day)) + ##use colour or shape to add in factor
                                   
                                   geom_point(aes( shape = Combo, color = Combo, size = Combo, stroke = 1.5)) +
                                    
                                   geom_errorbar(aes(ymin = TMX_PPB - se, ymax = TMX_PPB + se)) +
  
                                   xlab("Time (d)") +
                                     
                                   ylab(TMX ~ (ng~mL^{-1})) +
  
                                   #scale_y_log10() +
                                 
                                   scale_shape_manual(name = "",
                                     
                                                      values = c(1, 19, 1, 19),
                                                      
                                                      labels = c("Clay, No Plant", "Clay, Plant",
                                                                 
                                                                  "Sand, No Plant", "Sand, Plant")) +
  
                                   scale_color_manual(name = "",
                                     
                                                      values = c('red3', "red3", 'blue1','blue1'), 
     
                                                      labels = c("Clay, No Plant", "Clay, Plant",
             
                                                                 "Sand, No Plant", "Sand, Plant")) +
  
                                   scale_size_manual(name = "",
                    
                                                     values = c(3.5, 3.5, 3.5, 3.5),
                    
                                                     labels = c("Clay, No Plant", "Clay, Plant",
                    
                                                                "Sand, No Plant", "Sand, Plant")) +
                                   
                                   geom_text(aes(label = Letters, y = TMX_PPB, fontface = "italic"),
                                             
                                             position = position_dodge(width=0.9),
                                             
                                             hjust= -2.2) +  #####add in signficance letter--adjusted for errror bars
  
                                   cleanup +
                                     
                                   theme(plot.margin = unit(c(-1,56.5,1,1), units="points"), ##good legend size
                                                            
  #                                       panel.border = element_rect(color = "black",
                                                                     
    #                                                                 fill = NA, size = 1.75),
                                         
                                         axis.text.y = element_text(size = 16),  ## change font size of y axis
                                         
                                         axis.title.y = element_text(size = 16),
                                         
                                         axis.text.x = element_text(size = 16),  ## change font size of y axis
                                         
                                         axis.title.x = element_text(size = 16),
                                         
                                        legend.margin = margin(0, 5, 10, 1, "point"), ##remove margin on legend
                                         
                                         legend.text = element_text( size = 16),
                                         
                                         legend.position = c(0.13, 0.75),
                                         
                                         legend.key=element_blank(),  ##removed gray border around legend symbols
                                         
                                         legend.key.size = unit(10, "mm"),
                                        
                                         legend.background = element_rect(color = "black", size = 1.25),
                                         
                                         legend.key.height = unit(7, "mm"),
                                         
                                         legend.key.width = unit(7, "mm"),
                                       
#                                         legend.box.spacing = unit(10, "mm"),
                                      
                                         legend.text.align = 0,

                                         axis.line.x = element_line(color = "black", size = 1.75),


                                         axis.ticks.x= element_line(size = 1.75, color = "black"),

                                         axis.line.y = element_line(color = "black", size = 1.75),
    
                                         axis.ticks.y = element_line(size = 1.75, color = "black"),

                                         axis.ticks.length = unit(2, "mm")) 
         
                                      ## add secondary axis to bottom plot and remove ticks
                                     

                                    #scale_y_continuous(sec.axis = sec_axis(~. * 1, labels = NULL,
                                                                           
  #                                      breaks = c(18)))
                                                                           
                                                                           



                                                                         
                                    
####had to add in a break (18) at some high y value to maintain double y axis, also works to separate two different axes



Leach.tall.TMX.ppb.graph 

Leach.tall.TMX.ppb.graph.comb <- grid.arrange(Leach.tall.TMX.ppb.inf, Leach.tall.TMX.ppb.graph,
             
             heights = c(1/4, 3/4))

Leach.tall.TMX.ppb.graph.comb


#pdf("C:/Users/Jesse/Desktop/Neonicotinoids_15/Data/Column-Study/Raw-Data-Outputs/Leach.tall.TMX.ppb.graph.pdf")                                                          
#Leach.tall.TMX.ppb.graph  
#dev.off()


###have to use ggsave for grid obsjects
#ggsave(file = "C:/Users/Jesse/Desktop/Neonicotinoids_15/Data/Column-Study/Raw-Data-Outputs/Leach.tall.TMX.ppb.graph.comb.pdf", Leach.tall.TMX.ppb.graph.comb)





#####TMX ppb small columns*************************************************
small.columns.stats

Leach.small.TMX.ppb.inf <- ggplot(cum.rain, aes(x = Day, y = Cumulative.rainfall.cm,
                                                
                                                ymin = 0,
                                                
                                                ymax = Cumulative.rainfall.cm,
                                                
                                                xmax = Day)) +
  geom_line(size = 1, linetype = 2) +
  
  scale_y_continuous(limits=c(25,0),
                     
                     expand=c(0,0), trans="reverse", position = "right") +
  theme_classic() +
  
  theme(plot.margin = unit(c(6,24.75,-37,50), units="points"),
        
        axis.title.y = element_text(vjust = 0.3, size = 14, face = "bold"),
        
        axis.text.y = element_text(size = 16),
        
        panel.border = element_rect(color = "black",
                                    
                                    fill = NA, size = 3),
        
        axis.ticks.y = element_line(size = 1, color = "black")) +
  
  ylab("T. Rainfall (cm)") +
  
  theme(axis.title.y = element_text(vjust = 0.25, hjust = -.7))



Leach.small.TMX.ppb.inf


Leach.small.TMX.ppb.graph <- ggplot(small.columns.stats, aes(Time.day, TMX_PPB,
                                                                 
                                                                 shape = Combo,
                                                                 
                                                                 ymax= TMX_PPB + se,
                                                                 
                                                                 xmax = Time.day)) + ##use colour or shape to add in factor
  
  geom_point(aes(shape = Combo, color = Combo, size = Combo, stroke = 1.5)) +
  
  geom_errorbar(aes(ymin = TMX_PPB - se, ymax = TMX_PPB + se)) +
  
  xlab("Time (d)") +
  
  ylab(TMX ~ (ng~mL^{-1})) +
  
  scale_y_log10() +
  
  scale_shape_manual(name = "",
                     
                     values = c(1, 19, 1, 19),
                     
                     labels = c("Clay, No Plant", "Clay, Plant",
                                
                                "Sand, No Plant", "Sand, Plant")) +
  
  scale_color_manual(name = "",
                     
                     values = c('red3', "red3", 'blue1','blue1'), 
                     
                     labels = c("Clay, No Plant", "Clay, Plant",
                                
                                "Sand, No Plant", "Sand, Plant")) +
  
  scale_size_manual(name = "",
                    
                    values = c(3.5, 3.5, 3.5, 3.5),
                    
                    labels = c("Clay, No Plant", "Clay, Plant",
                               
                               "Sand, No Plant", "Sand, Plant")) +
  
  geom_text(aes(label = Letters, y = TMX_PPB, fontface = "italic"),
            
            position = position_dodge(width=0.9),
            
            hjust= -1.4) +  #####add in signficance letter--adjusted for errror bars
  
  geom_segment(aes(x= 15, y = 17.5, xend = 32, yend = 17.5), color = "firebrick1", linetype = 2, size = 1.5,
               
                   data = small.columns.stats) +
  
  #geom_abline(intercept = 17.5, slope = 0, size = 1.5, color = "firebrick1", linetype = 5) +
  
  annotate("text", x = 22, y = 23, label = "Aquatic Life Benchmark (Acute Toxicity)") + ##add benchamrk
  
  cleanup +
  
  theme(plot.margin = unit(c(-1,56.5,1,1), units="points"), ##good legend size
        
        #                                       panel.border = element_rect(color = "black",
        
        #                                                                 fill = NA, size = 1.75),
        
        axis.text.y = element_text(size = 16),  ## change font size of y axis
        
        axis.title.y = element_text(size = 16),
        
        axis.text.x = element_text(size = 16),  ## change font size of y axis
        
        axis.title.x = element_text(size = 16),
        
        legend.margin = margin(0, 5, 10, 1, "point"), ##remove margin on legend
        
        legend.text = element_text( size = 16),
        
        legend.position = c(0.13, 0.75),
        
        legend.key=element_blank(),  ##removed gray border around legend symbols
        
        legend.key.size = unit(10, "mm"),
        
        legend.background = element_rect(color = "black", size = 1.25),
        
        legend.key.height = unit(7, "mm"),
        
        legend.key.width = unit(7, "mm"),
        
        #                                         legend.box.spacing = unit(10, "mm"),
        
        legend.text.align = 0,
        
        axis.line.x = element_line(color = "black", size = 1.75),
        
        
        axis.ticks.x= element_line(size = 1.75, color = "black"),
        
        axis.line.y = element_line(color = "black", size = 1.75),
        
        axis.ticks.y = element_line(size = 1.75, color = "black"),
        
        axis.ticks.length = unit(2, "mm"),
        
        axis.title.y.right = element_text(size = 14)) 
  
  ## add secondary axis to bottom plot and remove ticks
  
  #scale_y_continuous(sec.axis = sec_axis(~. * 1, labels = NULL,
                                         
   #                                      breaks = c(85)))

Leach.small.TMX.ppb.graph


Leach.small.TMX.ppb.graph.comb <- grid.arrange(Leach.small.TMX.ppb.inf, Leach.small.TMX.ppb.graph,
                                              
                                              heights = c(1/4, 3/4))

Leach.small.TMX.ppb.graph.comb

?geom_segment()

?geom_abline()

?theme()

#pdf("C:/Users/Jesse/Desktop/Neonicotinoids_15/Data/Column-Study/Raw-Data-Outputs/Leach.small.TMX.ppb.graph.pdf")                                                          
Leach.small.TMX.ppb.graph  
dev.off()


###have to use ggsave for grid obsjects
#ggsave(file = "C:/Users/Jesse/Desktop/Neonicotinoids_15/Data/Column-Study/Raw-Data-Outputs/Leach.small.TMX.ppb.graph.comb.pdf", Leach.small.TMX.ppb.graph.comb)













#########TMX ppb structure

structure.study.stats



structure.study.TMX.ppb.inf <- ggplot(cum.rain, aes(x = Day, y = Cumulative.rainfall.cm,
                                                    
                                                    ymin = 0,
                                                    
                                                    ymax = Cumulative.rainfall.cm,
                                                    
                                                    xmax = Day)) +
  geom_line(size = 1, linetype = 2) +
  
  scale_y_continuous(limits=c(25,0),
                     
                     expand=c(0,0), trans="reverse", position = "right") +
  theme_classic() +
  
  theme(plot.margin = unit(c(6,24.75,-37,50), units="points"),
        
        axis.title.y = element_text(vjust = 0.3, size = 14, face = "bold"),
        
        axis.text.y = element_text(size = 16),
        
        panel.border = element_rect(color = "black",
                                    
                                    fill = NA, size = 3),
        
        axis.ticks.y = element_line(size = 1, color = "black")) +
  
  ylab("T. Rainfall (cm)") +
  
  theme(axis.title.y = element_text(vjust = 0.25, hjust = -.7))



structure.study.TMX.ppb.inf


structure.study.TMX.ppb.graph <- ggplot(structure.study.stats, aes(Time.day, TMX_PPB,
                                                             
                                                             shape = Structure,
                                                             
                                                             ymax= TMX_PPB + se,
                                                             
                                                             xmax = Time.day)) + ##use colour or shape to add in factor
  
  geom_point(aes(shape = Structure, color = Structure, size = Structure, stroke = 1.5)) +
  
  geom_errorbar(aes(ymin = TMX_PPB - se, ymax = TMX_PPB + se)) +
  
  xlab("Time (d)") +
  
  ylab(TMX ~ (ng~mL^{-1})) +
  
  scale_y_log10() +
  
  scale_shape_manual(name = "",
                     
                     values = c(1, 19),
                     
                     labels = c("Untructured Clay", "Structured Clay")) +
  
  scale_color_manual(name = "",
                     
                     values = c('black', "black"), 
                     
                     labels = c("Untructured Clay", "Structured Clay")) +
  
  scale_size_manual(name = "",
                    
                    values = c(3.5, 3.5),
                    
                    labels = c("Untructured Clay", "Structured Clay")) +
  
  geom_text(aes(label = Letters, y = TMX_PPB, fontface = "italic"),
            
            position = position_dodge(width=0.9),
            
            hjust= -1) +  #####add in signficance letter--adjusted for errror bars
  
  cleanup +
  
  theme(plot.margin = unit(c(-1.8,56.5,1,1), units="points"), ##good legend size
        
        #                                       panel.border = element_rect(color = "black",
        
        #                                                                 fill = NA, size = 1.75),
        
        axis.text.y = element_text(size = 16),  ## change font size of y axis
        
        axis.title.y = element_text(size = 16),
        
        axis.text.x = element_text(size = 16),  ## change font size of y axis
        
        axis.title.x = element_text(size = 16),
        
        legend.margin = margin(0, 5, 10, 0, "point"), ##remove margin on legend
        
        legend.text = element_text(size = 16),
        
        legend.position = c(0.15, 0.75),
        
        legend.key=element_blank(),  ##removed gray border around legend symbols
        
        legend.key.size = unit(10, "mm"),
        
        legend.background = element_rect(color = "black", size = 1.25),
        
        legend.key.height = unit(7, "mm"),
        
        legend.key.width = unit(7, "mm"),
        
        #                                         legend.box.spacing = unit(10, "mm"),
        
        legend.text.align = 0,
        
        axis.line.x = element_line(color = "black", size = 1.75),
        
        
        axis.ticks.x= element_line(size = 1.75, color = "black"),
        
        axis.line.y = element_line(color = "black", size = 1.75),
        
        axis.ticks.y = element_line(size = 1.75, color = "black"),
        
        axis.ticks.length = unit(2, "mm")) 
  
  ## add secondary axis to bottom plot and remove ticks
  
  #scale_y_continuous(sec.axis = sec_axis(~. * 1, labels = NULL,
                                         
                                         #breaks = c(18)))

structure.study.TMX.ppb.graph


structure.study.TMX.ppb.graph.comb <- grid.arrange(structure.study.TMX.ppb.inf, structure.study.TMX.ppb.graph,
                                               
                                               heights = c(1/4, 3/4))

structure.study.TMX.ppb.graph.comb



#pdf("C:/Users/Jesse/Desktop/Neonicotinoids_15/Data/Column-Study/Raw-Data-Outputs/structure.study.TMX.ppb.graph.pdf")                                                          
structure.study.TMX.ppb.graph 
dev.off()


###have to use ggsave for grid obsjects
#ggsave(file = "C:/Users/Jesse/Desktop/Neonicotinoids_15/Data/Column-Study/Raw-Data-Outputs/structure.study.TMX.ppb.graph.comb.pdf", structure.study.TMX.ppb.graph.comb)








###work on mass transport*********************************************************

###tall columns

#limitsx <- aes(xmin = tall.columns.no.SS.stats.mass$cum.Leach.vol.mL - tall.columns.no.SS.stats.mass$se1,
              
              #xmax = tall.columns.no.SS.stats.mass$cum.Leach.vol.mL + tall.columns.no.SS.stats.mass$se1)

tall.columns.no.SS.stats.mass

Leach.tall.TMX.ppb.mass.graph <- ggplot(tall.columns.no.SS.stats.mass, aes(cum.Leach.vol.mL, cum.TMX.micrg,
                                                                 
                                                                 shape = Combo,
                                                                 
                                                                 ymax= cum.TMX.micrg + se,
                                                                 
                                                                 xmax = cum.Leach.vol.mL)) + ##use colour or shape to add in factor
  
  geom_errorbar(aes(ymin = cum.TMX.micrg - se, ymax = cum.TMX.micrg + se),
                
                 colour="black", size = 0.75) + ### caps won't show
  
  geom_errorbarh(aes(xmin = cum.Leach.vol.mL - se1, xmax = cum.Leach.vol.mL + se1),
                 
                  colour="black", size = 0.75) +
  
  geom_point(aes(shape = Combo, color = Combo, size = Combo, stroke = 1.5)) +
  
  xlab("Cumulative Leachate Volume (mL)") +
  
  ylab(TMX ~ Leached ~ (mu*g)) +
  
  scale_shape_manual(name = "",
                     
                     values = c(1, 19, 1, 19),
                     
                     labels = c("Clay, No Plant", "Clay, Plant",
                                
                                "Sand, No Plant", "Sand, Plant")) +
  
  scale_color_manual(name = "",
                     
                     values = c('red3', "red3", 'blue1','blue1'), 
                     
                     labels = c("Clay, No Plant", "Clay, Plant",
                                
                                "Sand, No Plant", "Sand, Plant")) +
  
  scale_size_manual(name = "",
                    
                    values = c(3.5, 3.5, 3.5, 3.5),
                    
                    labels = c("Clay, No Plant", "Clay, Plant",
                               
                               "Sand, No Plant", "Sand, Plant")) +
  
  geom_text(aes(label = Letters, y = cum.TMX.micrg, fontface = "italic"),
            
            
            
            hjust= -2.2, vjust = -0.5) +  #####add in signficance letter--adjusted for errror bars
  
  geom_text(aes(label = Vol.diff, y = cum.TMX.micrg),
            
            position = position_dodge(width =0.9),
            
            hjust= -2.2, vjust = 1.25, size = 6, fontface = "bold") +
  cleanup +
  
  theme(plot.margin = unit(c(5,56.5,1,5), units="points"), ##good legend size
        
        panel.border = element_rect(color = "black",
        
                                    fill = NA, size = 1.75),
        
        axis.text.y = element_text(size = 16),  ## change font size of y axis
        
        axis.title.y = element_text(size = 16),
        
        axis.text.x = element_text(size = 16),  ## change font size of y axis
        
        axis.title.x = element_text(size = 16),
        
        legend.margin = margin(0, 5, 10, 1, "point"), ##remove margin on legend
        
        legend.text = element_text( size = 16),
        
        legend.position = c(0.2, 0.8),
        
        legend.key=element_blank(),  ##removed gray border around legend symbols
        
        legend.key.size = unit(10, "mm"),
        
        legend.background = element_rect(color = "black", size = 1.25),
        
        legend.key.height = unit(7, "mm"),
        
        legend.key.width = unit(7, "mm"),
        
        #                                         legend.box.spacing = unit(10, "mm"),
        
        legend.text.align = 0,
        
        axis.line.x = element_line(color = "black", size = 1.75),
        
        
        axis.ticks.x= element_line(size = 1.75, color = "black"),
        
        axis.line.y = element_line(color = "black", size = 1.75),
        
        axis.ticks.y = element_line(size = 1.75, color = "black"),
        
        axis.ticks.length = unit(2, "mm")) 
  
  


  

#Leach.tall.TMX.ppb.mass.graph

pdf("C:/Users/Jesse/Desktop/Neonicotinoids_15/Data/Column-Study/Raw-Data-Outputs/Leach.tall.TMX.ppb.mass.graph1.pdf")
#Leach.tall.TMX.ppb.mass.graph
dev.off()
























######now structure study***************************

structure.study.TMX.mass.graph <- ggplot(structure.study.stats.mass, aes(cum.Leach.vol.mL, cum.TMX.micrg,
                                                                   
                                                                   shape = Structure,
                                                                   
                                                                   ymax= cum.TMX.micrg + se,
                                                                   
                                                                   xmax = cum.Leach.vol.mL + se1)) + ##use colour or shape to add in factor
 
  geom_errorbar(aes(ymin = cum.TMX.micrg - se, ymax = cum.TMX.micrg + se),
                
                size = 0.75) +
  
  geom_errorbarh(aes(xmin = cum.Leach.vol.mL - se1, xmax = cum.Leach.vol.mL + se1),
                 
                 size = 0.75) +
   
  geom_point(aes(shape = Structure, color = Structure, size = Structure, stroke = 1.5)) +
  
  
  xlab("Cumulative Leachate Volume (mL)") +
  
  ylab(TMX ~ Leached ~ (mu*g)) +
  
  scale_shape_manual(name = "",
                     
                     values = c(1, 19),
                     
                     labels = c("Unstructured Clay", "Structured Clay")) +
  
  scale_color_manual(name = "",
                     
                     values = c('black', "black"), 
                     
                     labels = c("Unstructured Clay", "Structured Clay")) +
  
  scale_size_manual(name = "",
                    
                    values = c(3.5, 3.5),
                    
                    labels = c("Unstructured Clay", "Structured Clay")) +
  
  geom_text(aes(label = Letters, y = cum.TMX.micrg, fontface = "italic"),
            
            position = position_dodge(width=0.9),
            
            hjust= -1.8, vjust =-0.5) +  #####add in signficance letter--adjusted for errror bars
  
  geom_text(aes(label = Vol.diff, y = cum.TMX.micrg),
            
            position = position_dodge(width =0.9),
            
            hjust= -2.2, vjust = 1.25, size = 6, fontface = "bold") +
  
  cleanup +
  
  theme(plot.margin = unit(c(5,56.5,1,5), units="points"), ##good legend size
        
                                              panel.border = element_rect(color = "black",
        
                                                                         fill = NA, size = 1.75),
        
        axis.text.y = element_text(size = 16),  ## change font size of y axis
        
        axis.title.y = element_text(size = 16),
        
        axis.text.x = element_text(size = 16),  ## change font size of y axis
        
        axis.title.x = element_text(size = 16),
        
        legend.margin = margin(0, 5, 10, 0, "point"), ##remove margin on legend
        
        legend.text = element_text(size = 16),
        
        legend.position = c(0.25, 0.85),
        
        legend.key=element_blank(),  ##removed gray border around legend symbols
        
        legend.key.size = unit(10, "mm"),
        
        legend.background = element_rect(color = "black", size = 1.25),
        
        legend.key.height = unit(7, "mm"),
        
        legend.key.width = unit(7, "mm"),
        
        #                                         legend.box.spacing = unit(10, "mm"),
        
        legend.text.align = 0,
        
        axis.line.x = element_line(color = "black", size = 1.75),
        
        
        axis.ticks.x= element_line(size = 1.75, color = "black"),
        
        axis.line.y = element_line(color = "black", size = 1.75),
        
        axis.ticks.y = element_line(size = 1.75, color = "black"),
        
        axis.ticks.length = unit(2, "mm")) 
  

structure.study.TMX.mass.graph


#pdf("C:/Users/Jesse/Desktop/Neonicotinoids_15/Data/Column-Study/Raw-Data-Outputs/structure.study.TMX.mass.graph.pdf")
structure.study.TMX.mass.graph
dev.off()




########Small columns    TMX mass***************************************

small.columns.stats.mass


Leach.small.TMX.mass.graph <- ggplot(small.columns.stats.mass, aes(cum.Leach.vol.mL, cum.TMX.micrg,
                                                             
                                                             shape = Combo,
                                                             
                                                             ymax= cum.TMX.micrg + se,
                                                             
                                                             xmax = cum.Leach.vol.mL + se1)) +
  
  geom_errorbar(aes(ymin = cum.TMX.micrg - se, ymax = cum.TMX.micrg + se),
                
                size = 0.75) +
  
  geom_errorbarh(aes(xmin = cum.Leach.vol.mL - se1, xmax = cum.Leach.vol.mL + se1),
                 
                 size = 0.75) +
  
  geom_point(aes(shape = Combo, color = Combo, size = Combo, stroke = 1.5)) +
  
                              xlab("Cumulative Leachate Volume (mL)") +
  
                              ylab(TMX ~ Leached ~ (mu*g)) +
  
                              scale_shape_manual(name = "",
                     
                                                 values = c(1, 19, 1, 19),
                     
                                                 labels = c("Clay, No Plant", "Clay, Plant",
                                
                                                            "Sand, No Plant", "Sand, Plant")) +
  
                              scale_color_manual(name = "",
                     
                                                 values = c('red3', "red3", 'blue1','blue1'), 
                     
                                                 labels = c("Clay, No Plant", "Clay, Plant",
                                
                                                 "Sand, No Plant", "Sand, Plant")) +
  
                              scale_size_manual(name = "",
                    
                                                values = c(3.5, 3.5, 3.5, 3.5),
                    
                                                labels = c("Clay, No Plant", "Clay, Plant",
                               
                                                           "Sand, No Plant", "Sand, Plant")) +
  
                              geom_text(aes(label = Letters, y = cum.TMX.micrg, fontface = "italic"),
            
                                            position = position_dodge(width=0.9),
            
                                            hjust= -1.4, vjust = -0.5) +  #####add in signficance letter--adjusted for errror bars
  
  #geom_abline(intercept = 17.5, slope = 0, size = 1.5, color = "firebrick1", linetype = 5) +
  
                              cleanup +
  
                              theme(plot.margin = unit(c(5,56.5,1,5), units="points"), ##good legend size
        
                                    panel.border = element_rect(color = "black",
        
                                                                fill = NA, size = 1.75),
        
                                    axis.text.y = element_text(size = 16),  ## change font size of y axis
        
                                    axis.title.y = element_text(size = 16),
        
                                    axis.text.x = element_text(size = 16),  ## change font size of y axis
        
                                    axis.title.x = element_text(size = 16),
        
                                    legend.margin = margin(0, 5, 10, 1, "point"), ##remove margin on legend
        
                                    legend.text = element_text( size = 16),
        
                                    legend.position = c(0.2, 0.85),
        
                                    legend.key=element_blank(),  ##removed gray border around legend symbols
        
                                    legend.key.size = unit(10, "mm"),
        
                                    legend.background = element_rect(color = "black", size = 1.25),
        
                                    legend.key.height = unit(7, "mm"),
        
                                    legend.key.width = unit(7, "mm"),
        
        #                                         legend.box.spacing = unit(10, "mm"),
        
                                    legend.text.align = 0,
        
                                    axis.line.x = element_line(color = "black", size = 1.75),
        
        
                                    axis.ticks.x= element_line(size = 1.75, color = "black"),
        
                                    axis.line.y = element_line(color = "black", size = 1.75),
        
                                    axis.ticks.y = element_line(size = 1.75, color = "black"),
        
                                    axis.ticks.length = unit(2, "mm"))
  
  
  

Leach.small.TMX.mass.graph


#pdf("C:/Users/Jesse/Desktop/Neonicotinoids_15/Data/Column-Study/Raw-Data-Outputs/Leach.small.TMX.mass.graph.pdf")
Leach.small.TMX.mass.graph
dev.off()



###########now arrange all figures in the grid???***************************************


Leach.tall.TMX.ppb.inf.grid <- ggplot(cum.rain, aes(x = Day, y = Cumulative.rainfall.cm,
                                                    
                                                    ymin = 0,
                                                    
                                                    ymax = Cumulative.rainfall.cm,
                                                    
                                                    xmax = Day)) +
  geom_line(size = 1, linetype = 2) +
  
  scale_y_continuous(limits=c(25,0),
                     
                     expand=c(0,0), trans="reverse", position = "right", breaks = c(0,10,20)) +
  theme_classic() +
  
  theme(plot.margin = unit(c(6,17,-31.5,54), units="points"),
        
        axis.title.y = element_text(vjust = 0.3, size = 12, face = "bold"),
        
        axis.text.y = element_text(size = 15),
        
        panel.border = element_rect(color = "black",
                                    
                                    fill = NA, size = 3),
        
        axis.ticks.y = element_line(size = 1, color = "black")) +
  
  ylab("Rain (cm)") +
  
  theme(axis.title.y = element_text(vjust = 1.5, hjust = 2.5))



Leach.tall.TMX.ppb.inf.grid


Leach.tall.TMX.ppb.graph.grid <- ggplot(tall.columns.no.SS.stats, aes(Time.day, TMX_PPB,
                                                                 
                                                                 shape = Combo,
                                                                 
                                                                 ymax= TMX_PPB + se,
                                                                 
                                                                 xmax = Time.day)) + ##use colour or shape to add in factor
  
  geom_point(aes( shape = Combo, color = Combo, size = Combo, stroke = 1.5)) +
  
  geom_errorbar(aes(ymin = TMX_PPB - se, ymax = TMX_PPB + se)) +
  
  xlab("") +
  
  ylab(expression(TMX ~ (mu * g ~ L^{-1}))) +
  
  scale_shape_manual(name = "",
                     
                     values = c(1, 19, 1, 19),
                     
                     labels = c("Clay, No Plant", "Clay, Plant",
                                
                                "Sand, No Plant", "Sand, Plant")) +
  
  scale_color_manual(name = "",
                     
                     values = c('red3', "red3", 'blue1','blue1'), 
                     
                     labels = c("Clay, No Plant", "Clay, Plant",
                                
                                "Sand, No Plant", "Sand, Plant")) +
  
  scale_size_manual(name = "",
                    
                    values = c(3.5, 3.5, 3.5, 3.5),
                    
                    labels = c("Clay, No Plant", "Clay, Plant",
                               
                               "Sand, No Plant", "Sand, Plant")) +
  
  geom_text(aes(label = Letters, y = TMX_PPB, fontface = "italic"),
            
            position = position_dodge(width=2.5),
            
            hjust= 2.2, vjust = -1) +  #####add in signficance letter--adjusted for errror bars
  
  cleanup +
  
  theme(plot.margin = unit(c(-1,46.5,1,2), units="points"), ##good legend size
        
        #                                       panel.border = element_rect(color = "black",
        
        #                                                                 fill = NA, size = 1.75),
        
        axis.text.y = element_text(size = 16),  ## change font size of y axis
        
        axis.title.y = element_text(size = 16),
        
        axis.text.x = element_text(size = 16),  ## change font size of y axis
        
        axis.title.x = element_text(size = 16),
        
        legend.margin = margin(0, 7, 10, 1, "point"), ##remove margin on legend
        
        legend.text = element_text( size = 16),
        
        legend.position = c(0.325, 0.75),
        
        legend.key=element_blank(),  ##removed gray border around legend symbols
        
        legend.key.size = unit(10, "mm"),
        
        legend.background = element_rect(color = "black", size = 1.25),
        
        legend.key.height = unit(7, "mm"),
        
        legend.key.width = unit(7, "mm"),
        
        #                                         legend.box.spacing = unit(10, "mm"),
        
        legend.text.align = 0,
        
        axis.line.x = element_line(color = "black", size = 1.75),
        
        
        axis.ticks.x= element_line(size = 1.75, color = "black"),
        
        axis.line.y = element_line(color = "black", size = 1.75),
        
        axis.ticks.y = element_line(size = 1.75, color = "black"),
        
        axis.ticks.length = unit(2, "mm")) +
  
  ## add secondary axis to bottom plot and remove ticks
  
  scale_y_continuous(sec.axis = sec_axis(~. * 1, labels = NULL,
                                         
                                         breaks = c(18)))

####had to add in a break (18) at some high y value to maintain double y axis, also works to separate two different axes



Leach.tall.TMX.ppb.graph.grid 

Leach.tall.TMX.ppb.graph.comb.grid <- grid.arrange(Leach.tall.TMX.ppb.inf.grid, Leach.tall.TMX.ppb.graph.grid,
                                              
                                              heights = c(1/4, 3/4))

Leach.tall.TMX.ppb.graph.comb.grid







#####TMX ppb small columns*************************************************
small.columns.stats

Leach.small.TMX.ppb.inf.grid <- ggplot(cum.rain, aes(x = Day, y = Cumulative.rainfall.cm,
                                                     
                                                     ymin = 0,
                                                     
                                                     ymax = Cumulative.rainfall.cm,
                                                     
                                                     xmax = Day)) +
  geom_line(size = 1, linetype = 2) +
  
  scale_y_continuous(limits=c(25,0),
                     
                     expand=c(0,0), trans="reverse", position = "right", breaks = c(0,10,20)) +
  theme_classic() +
  
  theme(plot.margin = unit(c(6,17,-32.5,54), units="points"),
        
        axis.title.y = element_text(vjust = 0.3, size = 12, face = "bold"),
        
        axis.text.y = element_text(size = 15),
        
        panel.border = element_rect(color = "black",
                                    
                                    fill = NA, size = 3),
        
        axis.ticks.y = element_line(size = 1, color = "black")) +
  
  ylab("Rain (cm)") +
  
  theme(axis.title.y = element_text(vjust = 1.5, hjust = 4)) 
  
 



Leach.small.TMX.ppb.inf.grid


Leach.small.TMX.ppb.graph.grid <- ggplot(small.columns.stats, aes(Time.day, TMX_PPB,
                                                             
                                                             shape = Combo,
                                                             
                                                             ymax= TMX_PPB + se,
                                                             
                                                             xmax = Time.day)) + ##use colour or shape to add in factor
  
  geom_point(aes(shape = Combo, color = Combo, size = Combo, stroke = 1.5)) +
  
  geom_errorbar(aes(ymin = TMX_PPB - se, ymax = TMX_PPB + se)) +
  
  xlab("Time (d)") +
  
  ylab(expression(TMX ~ (mu * g ~ L^{-1}))) +
  
  scale_shape_manual(name = "",
                     
                     values = c(1, 19, 1, 19),
                     
                     labels = c("Clay, No Plant", "Clay, Plant",
                                
                                "Sand, No Plant", "Sand, Plant")) +
  
  scale_color_manual(name = "",
                     
                     values = c('red3', "red3", 'blue1','blue1'), 
                     
                     labels = c("Clay, No Plant", "Clay, Plant",
                                
                                "Sand, No Plant", "Sand, Plant")) +
  
  scale_size_manual(name = "",
                    
                    values = c(3.5, 3.5, 3.5, 3.5),
                    
                    labels = c("Clay, No Plant", "Clay, Plant",
                               
                               "Sand, No Plant", "Sand, Plant")) +
  
  geom_text(aes(label = Letters, y = TMX_PPB, fontface = "italic"),
            
            position = position_dodge(width=0.9),
            
            hjust= 2.5) +  #####add in signficance letter--adjusted for errror bars
  
  geom_segment(aes(x= 5, y = 17.5, xend = 25, yend = 17.5), color = "firebrick1", linetype = 2, size = 1.5,
               
               data = small.columns.stats) +
  
  #geom_abline(intercept = 17.5, slope = 0, size = 1.5, color = "firebrick1", linetype = 5) +
  
  annotate("text", x = 15, y = 25, label = "Aquatic Life Benchmark") + ##add benchamrk
  
  cleanup +
  
  theme(plot.margin = unit(c(-1,46.5,1,2), units="points"), ##good legend size
        
        #                                       panel.border = element_rect(color = "black",
        
        #                                                                 fill = NA, size = 1.75),
        
        axis.text.y = element_text(size = 16),  ## change font size of y axis
        
        axis.title.y = element_text(size = 16),
        
        axis.text.x = element_text(size = 16),  ## change font size of y axis
        
        axis.title.x = element_text(size = 16, face = "bold"),
        
        legend.margin = margin(0, 5, 10, 1, "point"), ##remove margin on legend
        
        legend.text = element_text( size = 16),
        
        legend.position = "none",
        
        legend.key=element_blank(),  ##removed gray border around legend symbols
        
        legend.key.size = unit(10, "mm"),
        
        legend.background = element_rect(color = "black", size = 1.25),
        
        legend.key.height = unit(7, "mm"),
        
        legend.key.width = unit(7, "mm"),
        
        #                                         legend.box.spacing = unit(10, "mm"),
        
        legend.text.align = 0,
        
        axis.line.x = element_line(color = "black", size = 1.75),
        
        
        axis.ticks.x= element_line(size = 1.75, color = "black"),
        
        axis.line.y = element_line(color = "black", size = 1.75),
        
        axis.ticks.y = element_line(size = 1.75, color = "black"),
        
        axis.ticks.length = unit(2, "mm")) +
  
  ## add secondary axis to bottom plot and remove ticks
  
  scale_y_continuous(sec.axis = sec_axis(~. * 1, labels = NULL,
                                         
                                         breaks = c(85)))

Leach.small.TMX.ppb.graph.grid


Leach.small.TMX.ppb.graph.comb.grid <- grid.arrange(Leach.small.TMX.ppb.inf.grid, Leach.small.TMX.ppb.graph.grid,
                                               
                                               heights = c(1/4, 3/4))

Leach.small.TMX.ppb.graph.comb.grid
















#########TMX ppb structure

structure.study.stats



structure.study.TMX.ppb.inf.grid <- ggplot(cum.rain, aes(x = Day, y = Cumulative.rainfall.cm,
                                                         
                                                         ymin = 0,
                                                         
                                                         ymax = Cumulative.rainfall.cm,
                                                         
                                                         xmax = Day)) +
  geom_line(size = 1, linetype = 2) +
  
  scale_y_continuous(limits=c(25,0),
                     
                     expand=c(0,0), trans="reverse", position = "right", breaks = c(0,10,20)) +
  theme_classic() +
  
  theme(plot.margin = unit(c(6,17,-32.5,54), units="points"),
        
        axis.title.y = element_text(vjust = 0.3, size = 12, face = "bold"),
        
        axis.text.y = element_text(size = 15),
        
        panel.border = element_rect(color = "black",
                                    
                                    fill = NA, size = 3),
        
        axis.ticks.y = element_line(size = 1, color = "black")) +
  
  ylab("Rain (cm)") +
  
  theme(axis.title.y = element_text(vjust = 1.5, hjust = 4))



structure.study.TMX.ppb.inf.grid


structure.study.TMX.ppb.graph.grid <- ggplot(structure.study.stats, aes(Time.day, TMX_PPB,
                                                                   
                                                                   shape = Structure,
                                                                   
                                                                   ymax= TMX_PPB + se,
                                                                   
                                                                   xmax = Time.day)) + ##use colour or shape to add in factor
  
  geom_point(aes(shape = Structure, color = Structure, size = Structure, stroke = 1.5)) +
  
  geom_errorbar(aes(ymin = TMX_PPB - se, ymax = TMX_PPB + se)) +
  
  xlab("") +
  
  ylab(expression(TMX ~ (mu * g ~ L^{-1}))) +
  
  scale_shape_manual(name = "",
                     
                     values = c(1, 19),
                     
                     labels = c("Unstructured Clay", "Structured Clay")) +
  
  scale_color_manual(name = "",
                     
                     values = c('black', "black"), 
                     
                     labels = c("Unstructured Clay", "Structured Clay")) +
  
  scale_size_manual(name = "",
                    
                    values = c(3.5, 3.5),
                    
                    labels = c("Unstructured Clay", "Structured Clay")) +
  
  geom_text(aes(label = Letters, y = TMX_PPB, fontface = "italic"),
            
            position = position_dodge(width=0.9),
            
            hjust= 2.25) +  #####add in signficance letter--adjusted for errror bars
  
  cleanup +
  
  theme(plot.margin = unit(c(-1,46.5,1,2), units="points"), ##good legend size
        
        #                                       panel.border = element_rect(color = "black",
        
        #                                                                 fill = NA, size = 1.75),
        
        axis.text.y = element_text(size = 16),  ## change font size of y axis
        
        axis.title.y = element_text(size = 16),
        
        axis.text.x = element_text(size = 16),  ## change font size of y axis
        
        axis.title.x = element_text(size = 16),
        
        legend.margin = margin(0, 7, 10, 0, "point"), ##remove margin on legend
        
        legend.text = element_text(size = 16),
        
        legend.position = c(0.35, 0.9),
        
        legend.key=element_blank(),  ##removed gray border around legend symbols
        
        legend.key.size = unit(8, "mm"),
        
        legend.background = element_rect(color = "black", size = 1.25),
        
        legend.key.height = unit(7, "mm"),
        
        legend.key.width = unit(7, "mm"),
        
        #                                         legend.box.spacing = unit(10, "mm"),
        
        legend.text.align = 0,
        
        axis.line.x = element_line(color = "black", size = 1.75),
        
        
        axis.ticks.x= element_line(size = 1.75, color = "black"),
        
        axis.line.y = element_line(color = "black", size = 1.75),
        
        axis.ticks.y = element_line(size = 1.75, color = "black"),
        
        axis.ticks.length = unit(2, "mm")) +
  
  ## add secondary axis to bottom plot and remove ticks
  
  scale_y_continuous(sec.axis = sec_axis(~. * 1, labels = NULL,
                                         
                                         breaks = c(18)))

structure.study.TMX.ppb.graph.grid


structure.study.TMX.ppb.graph.comb.grid <- grid.arrange(structure.study.TMX.ppb.inf.grid, structure.study.TMX.ppb.graph.grid,
                                                   
                                                   heights = c(1/4, 3/4))

structure.study.TMX.ppb.graph.comb.grid






###work on mass transport*********************************************************

###tall columns

#limitsx <- aes(xmin = tall.columns.no.SS.stats.mass$cum.Leach.vol.mL - tall.columns.no.SS.stats.mass$se1,

#xmax = tall.columns.no.SS.stats.mass$cum.Leach.vol.mL + tall.columns.no.SS.stats.mass$se1)

tall.columns.no.SS.stats.mass

Leach.tall.TMX.ppb.mass.graph.grid <- ggplot(tall.columns.no.SS.stats.mass, aes(cum.Leach.vol.mL, cum.TMX.micrg,
                                                                           
                                                                           shape = Combo,
                                                                           
                                                                           ymax= cum.TMX.micrg + se,
                                                                           
                                                                           xmax = cum.Leach.vol.mL)) + ##use colour or shape to add in factor
  
  geom_errorbar(aes(ymin = cum.TMX.micrg - se, ymax = cum.TMX.micrg + se),
                
                colour="black", size = 0.75) + ### caps won't show
  
  geom_errorbarh(aes(xmin = cum.Leach.vol.mL - se1, xmax = cum.Leach.vol.mL + se1),
                 
                 colour="black", size = 0.75) +
  
  geom_point(aes(shape = Combo, color = Combo, size = Combo, stroke = 1.5)) +
  
  xlab("") +
  
  ylab(TMX ~ Leached ~ (mu*g)) +
  
  scale_shape_manual(name = "",
                     
                     values = c(1, 19, 1, 19),
                     
                     labels = c("Clay, No Plant", "Clay, Plant",
                                
                                "Sand, No Plant", "Sand, Plant")) +
  
  scale_color_manual(name = "",
                     
                     values = c('red3', "red3", 'blue1','blue1'), 
                     
                     labels = c("Clay, No Plant", "Clay, Plant",
                                
                                "Sand, No Plant", "Sand, Plant")) +
  
  scale_size_manual(name = "",
                    
                    values = c(3.5, 3.5, 3.5, 3.5),
                    
                    labels = c("Clay, No Plant", "Clay, Plant",
                               
                               "Sand, No Plant", "Sand, Plant")) +
  
  geom_text(aes(label = Letters, y = cum.TMX.micrg, fontface = "italic"),
            
            
            
            hjust= -1.25, vjust = -0.5) +  #####add in signficance letter--adjusted for errror bars
  
  geom_text(aes(label = Vol.diff, y = cum.TMX.micrg),
            
            position = position_dodge(width =0.9),
            
            hjust= -2.25, vjust = 0.5, size = 6, fontface = "bold") +
  
  annotation_custom(column.60.pic.grob, xmin = -200, xmax = 800, ymin=5, ymax=55) + ##add in column pic
  
 # annotation_raster(column.60.pic, ymin = 10, ymax= 50, xmin = 500,xmax = 1300) +
  
  cleanup +
  
  theme(plot.margin = unit(c(5,36.5,1,27), units="points"), ##good legend size
        
        panel.border = element_rect(color = "black",
                                    
                                    fill = NA, size = 1.75),
        
        axis.text.y = element_text(size = 16),  ## change font size of y axis
        
        axis.title.y = element_text(size = 16),
        
        axis.text.x = element_text(size = 16),  ## change font size of y axis
        
        axis.title.x = element_text(size = 16),
        
        legend.margin = margin(0, 5, 10, 1, "point"), ##remove margin on legend
        
        legend.text = element_text( size = 16),
        
        legend.position = "none",
        
        legend.key=element_blank(),  ##removed gray border around legend symbols
        
        legend.key.size = unit(10, "mm"),
        
        legend.background = element_rect(color = "black", size = 1.25),
        
        legend.key.height = unit(7, "mm"),
        
        legend.key.width = unit(7, "mm"),
        
        #                                         legend.box.spacing = unit(10, "mm"),
        
        legend.text.align = 0,
        
        axis.line.x = element_line(color = "black", size = 1.75),
        
        
        axis.ticks.x= element_line(size = 1.75, color = "black"),
        
        axis.line.y = element_line(color = "black", size = 1.75),
        
        axis.ticks.y = element_line(size = 1.75, color = "black"),
        
        axis.ticks.length = unit(2, "mm")) +
  
  scale_x_continuous(limits = c(0, 4000))






Leach.tall.TMX.ppb.mass.graph.grid







######now structure study***************************

structure.study.TMX.mass.graph.grid <- ggplot(structure.study.stats.mass, aes(cum.Leach.vol.mL, cum.TMX.micrg,
                                                                         
                                                                         shape = Structure,
                                                                         
                                                                         ymax= cum.TMX.micrg + se,
                                                                         
                                                                         xmax = cum.Leach.vol.mL + se1)) + ##use colour or shape to add in factor
  
  geom_errorbar(aes(ymin = cum.TMX.micrg - se, ymax = cum.TMX.micrg + se),
                
                size = 0.75) +
  
  geom_errorbarh(aes(xmin = cum.Leach.vol.mL - se1, xmax = cum.Leach.vol.mL + se1),
                 
                 size = 0.75) +
  
  geom_point(aes(shape = Structure, color = Structure, size = Structure, stroke = 1.5)) +
  
  
  xlab("") +
  
  ylab(TMX ~ Leached ~ (mu*g)) +
  
  scale_shape_manual(name = "",
                     
                     values = c(1, 19),
                     
                     labels = c("Untructured Clay", "Structured Clay")) +
  
  scale_color_manual(name = "",
                     
                     values = c('black', "black"), 
                     
                     labels = c("Unstructured Clay", "Structured Clay")) +
  
  scale_size_manual(name = "",
                    
                    values = c(3.5, 3.5),
                    
                    labels = c("Unstructured Clay", "Structured Clay")) +
  
  geom_text(aes(label = Letters, y = cum.TMX.micrg, fontface = "italic"),
            
            position = position_dodge(width=0.9),
            
            hjust= -1.8, vjust =-0.5) +  #####add in signficance letter--adjusted for errror bars
  
  geom_text(aes(label = Vol.diff, y = cum.TMX.micrg),
            
            position = position_dodge(width =0.9),
            
            hjust= -2.75, vjust = 1, size = 6, fontface = "bold") +
  
  annotation_custom(column.60.pic.grob, xmin = -200, xmax = 800, ymin=-1, ymax=49) + ##add in column pic
  
  cleanup +
  
  theme(plot.margin = unit(c(5,36.5,1,27), units="points"), ##good legend size
        
        panel.border = element_rect(color = "black",
                                    
                                    fill = NA, size = 1.75),
        
        axis.text.y = element_text(size = 16),  ## change font size of y axis
        
        axis.title.y = element_text(size = 16),
        
        axis.text.x = element_text(size = 16),  ## change font size of y axis
        
        axis.title.x = element_text(size = 16),
        
        legend.margin = margin(0, 5, 10, 1, "point"), ##remove margin on legend
        
        legend.text = element_text(size = 16),
        
        legend.position = "none",
        
        legend.key=element_blank(),  ##removed gray border around legend symbols
        
        legend.key.size = unit(10, "mm"),
        
        legend.background = element_rect(color = "black", size = 1.25),
        
        legend.key.height = unit(7, "mm"),
        
        legend.key.width = unit(7, "mm"),
        
        #                                         legend.box.spacing = unit(10, "mm"),
        
        legend.text.align = 0,
        
        axis.line.x = element_line(color = "black", size = 1.75),
        
        
        axis.ticks.x= element_line(size = 1.75, color = "black"),
        
        axis.line.y = element_line(color = "black", size = 1.75),
        
        axis.ticks.y = element_line(size = 1.75, color = "black"),
        
        axis.ticks.length = unit(2, "mm")) +
  
  scale_x_continuous( limits = c(0,4000))


structure.study.TMX.mass.graph.grid




########Small columns    TMX mass***************************************

small.columns.stats.mass


Leach.small.TMX.mass.graph.grid <- ggplot(small.columns.stats.mass, aes(cum.Leach.vol.mL, cum.TMX.micrg,
                                                                   
                                                                   shape = Combo,
                                                                   
                                                                   ymax= cum.TMX.micrg + se,
                                                                   
                                                                   xmax = cum.Leach.vol.mL + se1)) +
  
  geom_errorbar(aes(ymin = cum.TMX.micrg - se, ymax = cum.TMX.micrg + se),
                
                size = 0.75) +
  
  geom_errorbarh(aes(xmin = cum.Leach.vol.mL - se1, xmax = cum.Leach.vol.mL + se1),
                 
                 size = 0.75) +
  
  geom_point(aes(shape = Combo, color = Combo, size = Combo, stroke = 1.5)) +
  
  xlab("Volume Leached (mL)") +
  
  ylab(TMX ~ Leached ~ (mu*g)) +
  
  scale_shape_manual(name = "",
                     
                     values = c(1, 19, 1, 19),
                     
                     labels = c("Clay, No Plant", "Clay, Plant",
                                
                                "Sand, No Plant", "Sand, Plant")) +
  
  scale_color_manual(name = "",
                     
                     values = c('red3', "red3", 'blue1','blue1'), 
                     
                     labels = c("Clay, No Plant", "Clay, Plant",
                                
                                "Sand, No Plant", "Sand, Plant")) +
  
  scale_size_manual(name = "",
                    
                    values = c(3.5, 3.5, 3.5, 3.5),
                    
                    labels = c("Clay, No Plant", "Clay, Plant",
                               
                               "Sand, No Plant", "Sand, Plant")) +
  
  geom_text(aes(label = Letters, y = cum.TMX.micrg, fontface = "italic"),
            
            position = position_dodge(width=0.9),
            
            hjust= -1.4, vjust = -0.5) +  #####add in signficance letter--adjusted for errror bars
  
  #geom_abline(intercept = 17.5, slope = 0, size = 1.5, color = "firebrick1", linetype = 5) +
  
  annotation_custom(column.20.pic.grob, xmin = -200, xmax = 800, ymin = 45, ymax = 250) + ##add in column pic
  
  cleanup +
  
  theme(plot.margin = unit(c(5,36.5,1,21), units="points"), ##good legend size
        
        panel.border = element_rect(color = "black",
                                    
                                    fill = NA, size = 1.75),
        
        axis.text.y = element_text(size = 16),  ## change font size of y axis
        
        axis.title.y = element_text(size = 16),
        
        axis.text.x = element_text(size = 16),  ## change font size of y axis
        
        axis.title.x = element_text(size = 16, face = "bold"),
        
        legend.margin = margin(0, 5, 10, 1, "point"), ##remove margin on legend
        
        legend.text = element_text( size = 16),
        
        legend.position = "none",
        
        legend.key=element_blank(),  ##removed gray border around legend symbols
        
        legend.key.size = unit(10, "mm"),
        
        legend.background = element_rect(color = "black", size = 1.25),
        
        legend.key.height = unit(7, "mm"),
        
        legend.key.width = unit(7, "mm"),
        
        #                                         legend.box.spacing = unit(10, "mm"),
        
        legend.text.align = 0,
        
        axis.line.x = element_line(color = "black", size = 1.75),
        
        
        axis.ticks.x= element_line(size = 1.75, color = "black"),
        
        axis.line.y = element_line(color = "black", size = 1.75),
        
        axis.ticks.y = element_line(size = 1.75, color = "black"),
        
        axis.ticks.length = unit(2, "mm")) +
  
  scale_x_continuous(limits = c(0, 4000))




Leach.small.TMX.mass.graph.grid





######arrange

Leach.tall.TMX.ppb.graph.comb.grid

structure.study.TMX.ppb.graph.comb.grid

Leach.small.TMX.ppb.graph.comb.grid

Leach.tall.TMX.ppb.mass.graph.grid

structure.study.TMX.mass.graph.grid

Leach.small.TMX.mass.graph.grid

Leach.small.TMX.mass.graph.grid


TMX.mass.and.conc.Leachate.grid <- grid.arrange(Leach.tall.TMX.ppb.graph.comb.grid,##conc
                                                
                                                structure.study.TMX.ppb.graph.comb.grid,##conc
             
                                                Leach.small.TMX.ppb.graph.comb.grid,###conc                                   
             
                                                Leach.tall.TMX.ppb.mass.graph.grid, ##mass
                                                
                                                structure.study.TMX.mass.graph.grid, ##mass
             
                                                Leach.small.TMX.mass.graph.grid, ##mass
             
                                                layout_matrix = rbind(c(1, 4),
                                                                      c(2, 5),
                                                                      c(3, 6)
)) 

TMX.mass.and.conc.Leachate.grid

#ggsave(file = "C:/Users/Jesse/Desktop/Neonicotinoids_15/Data/Column-Study/Raw-Data-Outputs/TMX.mass.and.conc.Leachate.grid1.pdf", TMX.mass.and.conc.Leachate.grid)
TMX.mass.and.conc.Leachate.grid
dev.off()


#TMX.mass.and.conc.Leachate.grid<- grid.arrange(Leach.tall.TMX.ppb.graph.comb.grid,##conc

  #            structure.study.TMX.ppb.graph.comb.grid,##conc

  #            structure.study.TMX.ppb.graph.comb.grid,###conc                                   

   #           Leach.small.TMX.ppb.graph.comb.grid, ###conc

   #           Leach.tall.TMX.ppb.mass.graph.grid, ##mass

    #          structure.study.TMX.mass.graph.grid, ##mass

     #         Leach.small.TMX.mass.graph.grid, ##mass
             
      #        ncol = 2)

#TMX.mass.and.conc.Leachate.grid

#TMX.mass.and.conc.Leachate.grid <- grid.arrange(Leach.tall.TMX.ppb.graph.comb.grid,##conc
             
  #           structure.study.TMX.ppb.graph.grid,##conc
             
   #          ncol = 1)


TMX.mass.and.conc.Leachate.grid





































