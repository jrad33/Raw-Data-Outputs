rm(list = ls())

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
###Estimating ET for Column Experiment



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

Leachate_full$Leach.Vol.mL[Leach.Vol.mL == 0] <- NA   ###replace 0 with NA in volumne leachate

Leachate_full$TMX_PPB[TMX_PPB == 0] <- NA   ###replace 0 with NA in leachate conc

Leachate_full




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

Leachate_full$ET.rate.mm.d <- Leachate_full$ET.cum.mm/Leachate_full$dt


Leachate_full[Leachate_full$Event == 10,]
Leachate_full[Leachate_full$Event == 11,]
Leachate_full[Leachate_full$Event == 12,]























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




######done with preplot analysis


#####attach column pics


library(png)

column.60.pic <- readPNG("C:/Users/Jesse/Desktop/Neonicotinoids_15/Pics/Columns-60.png")

column.60.pic.grob <- rasterGrob(column.60.pic, interpolate = T)

column.20.pic <- readPNG("C:/Users/Jesse/Desktop/Neonicotinoids_15/Pics/columns-20.png")

column.20.pic.grob <- rasterGrob(column.20.pic, interpolate = T)

#######















######now start plotting!!!!!!!!


#######plotting a CUmualtive ET (mm) vs t seems most appropriate
####keep the same the theme with cumulative rainfall on inverse y axis



##actually, shoulf probably convert ET to mm to be consistent with rainfall on iverse y

##tall columns
tall.columns.no.SS.stats.ET$ET.cum.cm <- tall.columns.no.SS.stats.ET$ET.cum.mm/10

tall.columns.no.SS.stats.ET$se.cm <- tall.columns.no.SS.stats.ET$se/10

##tall columns + structure study
structure.study.stats.ET$ET.cum.cm <- structure.study.stats.ET$ET.cum.mm/10

structure.study.stats.ET$se.cm <- structure.study.stats.ET$se/10

##small columns

small.columns.stats.ET$ET.cum.cm <- small.columns.stats.ET$ET.cum.mm/10

small.columns.stats.ET$se.cm <- small.columns.stats.ET$se/10




#####cleanup themes for ggplot

cleanup <- theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background =  element_blank(),
                axis.line = element_line(color = "black"))








Leach.tall.ET.inf <- ggplot(tall.columns.no.SS.stats.ET, aes(x = Time.day, y= Cum.inf.cm,
                                                                    
                                                                    ymin = 0,
                                                                    
                                                                    ymax = Cum.inf.cm,
                                                                    
                                                                    xmax = Time.day)) +
  geom_line(size = 1, linetype = 2) +
  
  scale_y_continuous(limits=c(25,0),
                     
                     expand=c(0,0), trans="reverse", position = "right", breaks = c(0,10,20)) +
  theme_classic() +
  
  theme(plot.margin = unit(c(6,17,-35.5,45), units="points"),
        
        axis.title.y = element_text(vjust = 0.3, size = 12, face = "bold"),
        
        axis.text.y = element_text(size = 15),
        
        panel.border = element_rect(color = "black",
                                    
                                    fill = NA, size = 3),
        
        axis.ticks.y = element_line(size = 1, color = "black")) +
  
  ylab("Rain (cm)") +
  
  theme(axis.title.y = element_text(vjust = 1.5, hjust = 0.5))



Leach.tall.ET.inf


Leach.tall.ET.graph <- ggplot(tall.columns.no.SS.stats.ET, aes(Time.day, ET.cum.cm,
                                                                      
                                                                      shape = Combo,
                                                                      
                                                                      ymax= ET.cum.cm + se.cm,
                                                                      
                                                                      xmax = Time.day)) + ##use colour or shape to add in factor
  
  geom_point(aes( shape = Combo, color = Combo, size = Combo, stroke = 1.5)) +
  
  geom_errorbar(aes(ymin = ET.cum.cm - se.cm, ymax = ET.cum.cm + se.cm)) +
  
  xlab("Time (d)") +
  
  ylab("Cumulative ET (cm)") +
  
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
  
  geom_text(aes(label = Letters, y = ET.cum.cm, fontface = "italic"),
            
            position = position_dodge(width=0.5),
            
            hjust= -1.5, vjust = 0.5) +  #####add in signficance letter--adjusted for errror bars
  
  annotation_custom(column.60.pic.grob, xmin = -1, xmax = 11, ymin=3, ymax=15) + ##add in column pic
  
  cleanup +
  
  theme(plot.margin = unit(c(-1,46.5,1,2), units="points"), ##good legend size
        
        #                                       panel.border = element_rect(color = "black",
        
        #                                                                 fill = NA, size = 1.75),
        
        axis.text.y = element_text(size = 16),  ## change font size of y axis
        
        axis.title.y = element_text(size = 16),
        
        axis.text.x = element_text(size = 16),  ## change font size of y axis
        
        axis.title.x = element_text(size = 16),
        
        legend.margin = margin(0, 7, 10, 1, "point"), ##remove margin on legend
        
        legend.text = element_text( size = 18),
        
        legend.position = c(0.8, 0.2),
        
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
                                         
                                         breaks = c(15)))

####had to add in a break (18) at some high y value to maintain double y axis, also works to separate two different axes



Leach.tall.ET.graph  

Leach.tall.ET.graph.comb <- grid.arrange(Leach.tall.ET.inf, Leach.tall.ET.graph ,
                                                   
                                                   heights = c(1/4, 3/4))

Leach.tall.ET.graph.comb


pdf("C:/Users/Jesse/Desktop/Neonicotinoids_15/Data/Column-Study/Raw-Data-Outputs/Leach.tall.ET.graph.pdf")                                                          
Leach.tall.ET.graph 
dev.off()


###have to use ggsave for grid obsjects
ggsave(file = "C:/Users/Jesse/Desktop/Neonicotinoids_15/Data/Column-Study/Raw-Data-Outputs/Leach.tall.ET.graph.comb.pdf", Leach.tall.ET.graph.comb)







######now structure study***************************


structure.study.stats.ET



structure.study.ET.inf <- ggplot(structure.study.stats.ET, aes(x = Time.day, y = Cum.inf.cm,
                                                                      
                                                                      ymin = 0,
                                                                      
                                                                      ymax = Cum.inf.cm,
                                                                      
                                                                      xmax = Time.day)) +
  geom_line(size = 1, linetype = 2) +
  
  scale_y_continuous(limits=c(25,0),
                     
                     expand=c(0,0), trans="reverse", position = "right", breaks = c(0,10,20)) +
  theme_classic() +
  
  theme(plot.margin = unit(c(6,17,-32.5,45), units="points"),
        
        axis.title.y = element_text(vjust = 0.3, size = 12, face = "bold"),
        
        axis.text.y = element_text(size = 15),
        
        panel.border = element_rect(color = "black",
                                    
                                    fill = NA, size = 3),
        
        axis.ticks.y = element_line(size = 1, color = "black")) +
  
  ylab("Rain (cm)") +
  
  theme(axis.title.y = element_text(vjust = 1.5, hjust = 0.5))



structure.study.ET.inf


structure.study.ET.graph <- ggplot(structure.study.stats.ET, aes(Time.day, ET.cum.cm,
                                                                        
                                                                        shape = Structure,
                                                                        
                                                                        ymax= ET.cum.cm + se.cm,
                                                                        
                                                                        xmax = Time.day)) + ##use colour or shape to add in factor
  
  geom_point(aes(shape = Structure, color = Structure, size = Structure, stroke = 1.5)) +
  
  geom_errorbar(aes(ymin = ET.cum.cm - se.cm, ymax = ET.cum.cm + se.cm)) +
  
  xlab("Time (d)") +
  
  ylab("Cumulative ET (cm)") +
  
  scale_shape_manual(name = "",
                     
                     values = c(1, 19),
                     
                     labels = c("Unstructured Clay", "Structured Clay")) +
  
  scale_color_manual(name = "",
                     
                     values = c('black', "black"), 
                     
                     labels = c("Unstructured Clay", "Structured Clay")) +
  
  scale_size_manual(name = "",
                    
                    values = c(3.5, 3.5),
                    
                    labels = c("Unstructured Clay", "Structured Clay")) +
  
  geom_text(aes(label = Letters, y = ET.cum.cm, fontface = "italic"),
            
            position = position_dodge(width=0.9),
            
            hjust= -2.25) +  #####add in signficance letter--adjusted for errror bars
  
  annotation_custom(column.60.pic.grob, xmin = -1, xmax = 11, ymin=3, ymax=15) + ##add in column pic
  
  cleanup +
  
  theme(plot.margin = unit(c(-1,46.5,1,2), units="points"), ##good legend size
        
        #                                       panel.border = element_rect(color = "black",
        
        #                                                                 fill = NA, size = 1.75),
        
        axis.text.y = element_text(size = 16),  ## change font size of y axis
        
        axis.title.y = element_text(size = 16),
        
        axis.text.x = element_text(size = 16),  ## change font size of y axis
        
        axis.title.x = element_text(size = 16),
        
        legend.margin = margin(0, 7, 10, 0, "point"), ##remove margin on legend
        
        legend.text = element_text(size = 18),
        
        legend.position = c(0.8, 0.2),
        
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
                                         
                                         breaks = c(15)))

structure.study.ET.graph





structure.study.ET.graph.comb <- grid.arrange(structure.study.ET.inf, structure.study.ET.graph,
                                                        
                                                        heights = c(1/4, 3/4))

structure.study.ET.graph.comb




pdf("C:/Users/Jesse/Desktop/Neonicotinoids_15/Data/Column-Study/Raw-Data-Outputs/structure.study.ET.graph.pdf")                                                          
structure.study.ET.graph 
dev.off()


###have to use ggsave for grid obsjects
ggsave(file = "C:/Users/Jesse/Desktop/Neonicotinoids_15/Data/Column-Study/Raw-Data-Outputs/structure.study.ET.graph.comb.pdf", structure.study.ET.graph.comb)









######now small columns########################################################################


small.columns.stats.ET

Leach.small.ET.inf <- ggplot(small.columns.stats.ET, aes(x = Time.day, y = Cum.inf.cm,
                                                         
                                                         ymin = 0,
                                                         
                                                         ymax = Cum.inf.cm,
                                                         
                                                         xmax = Time.day)) +
  geom_line(size = 1, linetype = 2) +
  
  scale_y_continuous(limits=c(25,0),
                     
                     expand=c(0,0), trans="reverse", position = "right", breaks = c(0,10,20)) +
  theme_classic() +
  
  theme(plot.margin = unit(c(6,17,-32.5,45), units="points"),
        
        axis.title.y = element_text(vjust = 0.3, size = 12, face = "bold"),
        
        axis.text.y = element_text(size = 15),
        
        panel.border = element_rect(color = "black",
                                    
                                    fill = NA, size = 3),
        
        axis.ticks.y = element_line(size = 1, color = "black")) +
  
  ylab("Rain (cm)") +
  
  theme(axis.title.y = element_text(vjust = 1.5, hjust = 0.5)) 


Leach.small.ET.inf


Leach.small.ET.graph <- ggplot(small.columns.stats.ET, aes(Time.day, ET.cum.cm,
                                                                  
                                                                  shape = Combo,
                                                                  
                                                                  ymax= ET.cum.cm + se.cm,
                                                                  
                                                                  xmax = Time.day)) + ##use colour or shape to add in factor
  
  geom_point(aes(shape = Combo, color = Combo, size = Combo, stroke = 1.5)) +
  
  geom_errorbar(aes(ymin = ET.cum.cm - se.cm, ymax = ET.cum.cm + se.cm)) +
  
  xlab("Time (d)") +
  
  ylab("Cumulative ET (cm)") +
  
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
  
  geom_text(aes(label = Letters, y = ET.cum.cm, fontface = "italic"),
            
            position = position_dodge(width=0.9),
            
            hjust= -2.5) +  #####add in signficance letter--adjusted for errror bars
  
  #geom_segment(aes(x= 5, y = 17.5, xend = 25, yend = 17.5), color = "firebrick1", linetype = 2, size = 1.5,
               
   #            data = small.columns.stats.ET) +
  
  #geom_abline(intercept = 17.5, slope = 0, size = 1.5, color = "firebrick1", linetype = 5) +
  
  #annotate("text", x = 15, y = 25, label = "Aquatic Life Benchmark") + ##add benchamrk
  
  annotation_custom(column.20.pic.grob, xmin = 1, xmax = 9, ymin=2, ymax=11) + ##add in column pic
  
  cleanup +
  
  theme(plot.margin = unit(c(-1,46.5,1,2), units="points"), ##good legend size
        
        #                                       panel.border = element_rect(color = "black",
        
        #                                                                 fill = NA, size = 1.75),
        
        axis.text.y = element_text(size = 16),  ## change font size of y axis
        
        axis.title.y = element_text(size = 16),
        
        axis.text.x = element_text(size = 16),  ## change font size of y axis
        
        axis.title.x = element_text(size = 16),
        
        legend.margin = margin(0, 5, 10, 1, "point"), ##remove margin on legend
        
        legend.text = element_text( size = 18),
        
        legend.position = c(0.8, 0.2),
        
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
                                         
                                         breaks = c(10))) 
  
  #scale_y_continuous(breaks = c(0, 5, 10, 15))

  

Leach.small.ET.graph


Leach.small.ET.graph.comb <- grid.arrange(Leach.small.ET.inf, Leach.small.ET.graph,
                                                    
                                                    heights = c(1/4, 3/4))

Leach.small.ET.graph.comb






pdf("C:/Users/Jesse/Desktop/Neonicotinoids_15/Data/Column-Study/Raw-Data-Outputs/Leach.small.ET.graph.pdf")                                                          
Leach.small.ET.graph
dev.off()


###have to use ggsave for grid obsjects
ggsave(file = "C:/Users/Jesse/Desktop/Neonicotinoids_15/Data/Column-Study/Raw-Data-Outputs/Leach.small.ET.graph.comb.pdf", Leach.small.ET.graph.comb)




#######################now work on the grid (multipanel figure)######################################




########grid##################################################################################################


####Tall columns-Grid##########################

Leach.tall.ET.inf.grid <- ggplot(tall.columns.no.SS.stats.ET, aes(x = Time.day, y= Cum.inf.cm,
                                                             
                                                             ymin = 0,
                                                             
                                                             ymax = Cum.inf.cm,
                                                             
                                                             xmax = Time.day)) +
  geom_line(size = 1, linetype = 2) +
  
  scale_y_continuous(limits=c(25,0),
                     
                     expand=c(0,0), trans="reverse", position = "right", breaks = c(0,10,20)) +
  theme_classic() +
  
  theme(plot.margin = unit(c(6,17,-35.5,45), units="points"),
        
        axis.title.y = element_text(vjust = 0.3, size = 12, face = "bold"),
        
        axis.text.y = element_text(size = 15),
        
        panel.border = element_rect(color = "black",
                                    
                                    fill = NA, size = 3),
        
        axis.ticks.y = element_line(size = 1, color = "black")) +
  
  ylab("Rain (cm)") +
  
  theme(axis.title.y = element_text(vjust = 1.5, hjust = 0.5))



Leach.tall.ET.inf.grid


Leach.tall.ET.graph.grid <- ggplot(tall.columns.no.SS.stats.ET, aes(Time.day, ET.cum.cm,
                                                               
                                                               shape = Combo,
                                                               
                                                               ymax= ET.cum.cm + se.cm,
                                                               
                                                               xmax = Time.day)) + ##use colour or shape to add in factor
  
  geom_point(aes( shape = Combo, color = Combo, size = Combo, stroke = 1.5)) +
  
  geom_errorbar(aes(ymin = ET.cum.cm - se.cm, ymax = ET.cum.cm + se.cm)) +
  
  xlab("") +
  
  ylab("Cumulative ET (cm)") +
  
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
  
  geom_text(aes(label = Letters, y = ET.cum.cm, fontface = "italic"),
            
            position = position_dodge(width=0.5),
            
            hjust= -1.5, vjust = 0.5) +  #####add in signficance letter--adjusted for errror bars
  
  annotation_custom(column.60.pic.grob, xmin = -1, xmax = 11, ymin=3, ymax=15) + ##add in column pic
  
  cleanup +
  
  theme(plot.margin = unit(c(-1,46.5,1,2), units="points"), ##good legend size
        
        #                                       panel.border = element_rect(color = "black",
        
        #                                                                 fill = NA, size = 1.75),
        
        axis.text.y = element_text(size = 16),  ## change font size of y axis
        
        axis.title.y = element_text(size = 16),
        
        axis.text.x = element_text(size = 16),  ## change font size of y axis
        
        axis.title.x = element_text(size = 16),
        
        legend.margin = margin(0, 7, 10, 1, "point"), ##remove margin on legend
        
        legend.text = element_text( size = 18),
        
        legend.position = c(0.8, 0.2),
        
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
                                         
                                         breaks = c(15)))

####had to add in a break (18) at some high y value to maintain double y axis, also works to separate two different axes



Leach.tall.ET.graph.grid  

Leach.tall.ET.graph.comb.grid <- grid.arrange(Leach.tall.ET.inf.grid, Leach.tall.ET.graph.grid,
                                         
                                         heights = c(1/4, 3/4))

Leach.tall.ET.graph.comb.grid








######now structure study_Grid***************************


structure.study.stats.ET



structure.study.ET.inf.grid <- ggplot(structure.study.stats.ET, aes(x = Time.day, y = Cum.inf.cm,
                                                               
                                                               ymin = 0,
                                                               
                                                               ymax = Cum.inf.cm,
                                                               
                                                               xmax = Time.day)) +
  geom_line(size = 1, linetype = 2) +
  
  scale_y_continuous(limits=c(25,0),
                     
                     expand=c(0,0), trans="reverse", position = "right", breaks = c(0,10,20)) +
  theme_classic() +
  
  theme(plot.margin = unit(c(6,17,-32.5,45), units="points"),
        
        axis.title.y = element_text(vjust = 0.3, size = 12, face = "bold"),
        
        axis.text.y = element_text(size = 15),
        
        panel.border = element_rect(color = "black",
                                    
                                    fill = NA, size = 3),
        
        axis.ticks.y = element_line(size = 1, color = "black")) +
  
  ylab("Rain (cm)") +
  
  theme(axis.title.y = element_text(vjust = 1.5, hjust = 0.5))



structure.study.ET.inf.grid


structure.study.ET.graph.grid <- ggplot(structure.study.stats.ET, aes(Time.day, ET.cum.cm,
                                                                 
                                                                 shape = Structure,
                                                                 
                                                                 ymax= ET.cum.cm + se.cm,
                                                                 
                                                                 xmax = Time.day)) + ##use colour or shape to add in factor
  
  geom_point(aes(shape = Structure, color = Structure, size = Structure, stroke = 1.5)) +
  
  geom_errorbar(aes(ymin = ET.cum.cm - se.cm, ymax = ET.cum.cm + se.cm)) +
  
  xlab("") +
  
  ylab("Cumulative ET (cm)") +
  
  scale_shape_manual(name = "",
                     
                     values = c(1, 19),
                     
                     labels = c("Unstructured Clay", "Structured Clay")) +
  
  scale_color_manual(name = "",
                     
                     values = c('black', "black"), 
                     
                     labels = c("Unstructured Clay", "Structured Clay")) +
  
  scale_size_manual(name = "",
                    
                    values = c(3.5, 3.5),
                    
                    labels = c("Unstructured Clay", "Structured Clay")) +
  
  geom_text(aes(label = Letters, y = ET.cum.cm, fontface = "italic"),
            
            position = position_dodge(width=0.9),
            
            hjust= -2.25) +  #####add in signficance letter--adjusted for errror bars
  
  annotation_custom(column.20.pic.grob, xmin = 1, xmax = 9, ymin=2, ymax=11) + ##add in column pic
  
  cleanup +
  
  theme(plot.margin = unit(c(-1,46.5,1,2), units="points"), ##good legend size
        
        #                                       panel.border = element_rect(color = "black",
        
        #                                                                 fill = NA, size = 1.75),
        
        axis.text.y = element_text(size = 16),  ## change font size of y axis
        
        axis.title.y = element_text(size = 16),
        
        axis.text.x = element_text(size = 16),  ## change font size of y axis
        
        axis.title.x = element_text(size = 16),
        
        legend.margin = margin(0, 7, 10, 0, "point"), ##remove margin on legend
        
        legend.text = element_text(size = 18),
        
        legend.position = c(0.8, 0.2),
        
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
                                         
                                         breaks = c(15)))

structure.study.ET.graph.grid





structure.study.ET.graph.comb.grid <- grid.arrange(structure.study.ET.inf.grid, structure.study.ET.graph.grid,
                                              
                                              heights = c(1/4, 3/4))

structure.study.ET.graph.comb.grid




######now small columns########################################################################


small.columns.stats.ET

Leach.small.ET.inf.grid <- ggplot(small.columns.stats.ET, aes(x = Time.day, y = Cum.inf.cm,
                                                         
                                                         ymin = 0,
                                                         
                                                         ymax = Cum.inf.cm,
                                                         
                                                         xmax = Time.day)) +
  geom_line(size = 1, linetype = 2) +
  
  scale_y_continuous(limits=c(25,0),
                     
                     expand=c(0,0), trans="reverse", position = "right", breaks = c(0,10,20)) +
  theme_classic() +
  
  theme(plot.margin = unit(c(6,17,-32.5,45), units="points"),
        
        axis.title.y = element_text(vjust = 0.3, size = 12, face = "bold"),
        
        axis.text.y = element_text(size = 15),
        
        panel.border = element_rect(color = "black",
                                    
                                    fill = NA, size = 3),
        
        axis.ticks.y = element_line(size = 1, color = "black")) +
  
  ylab("Rain (cm)") +
  
  theme(axis.title.y = element_text(vjust = 1.5, hjust = 0.5)) 


Leach.small.ET.inf.grid


Leach.small.ET.graph.grid <- ggplot(small.columns.stats.ET, aes(Time.day, ET.cum.cm,
                                                           
                                                           shape = Combo,
                                                           
                                                           ymax= ET.cum.cm + se.cm,
                                                           
                                                           xmax = Time.day)) + ##use colour or shape to add in factor
  
  geom_point(aes(shape = Combo, color = Combo, size = Combo, stroke = 1.5)) +
  
  geom_errorbar(aes(ymin = ET.cum.cm - se.cm, ymax = ET.cum.cm + se.cm)) +
  
  xlab("Time (d)") +
  
  ylab("Cumulative ET (cm)") +
  
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
  
  geom_text(aes(label = Letters, y = ET.cum.cm, fontface = "italic"),
            
            position = position_dodge(width=0.9),
            
            hjust= -2.5) +  #####add in signficance letter--adjusted for errror bars
  
  #geom_segment(aes(x= 5, y = 17.5, xend = 25, yend = 17.5), color = "firebrick1", linetype = 2, size = 1.5,
  
  #            data = small.columns.stats.ET) +
  
  #geom_abline(intercept = 17.5, slope = 0, size = 1.5, color = "firebrick1", linetype = 5) +
  
  #annotate("text", x = 15, y = 25, label = "Aquatic Life Benchmark") + ##add benchamrk
  
  annotation_custom(column.20.pic.grob, xmin = 1, xmax = 9, ymin=2, ymax=11) + ##add in column pic
  
  cleanup +
  
  theme(plot.margin = unit(c(-1,46.5,1,2), units="points"), ##good legend size
        
        #                                       panel.border = element_rect(color = "black",
        
        #                                                                 fill = NA, size = 1.75),
        
        axis.text.y = element_text(size = 16),  ## change font size of y axis
        
        axis.title.y = element_text(size = 16),
        
        axis.text.x = element_text(size = 16),  ## change font size of y axis
        
        axis.title.x = element_text(size = 16, face = "bold"),
        
        legend.margin = margin(0, 5, 10, 1, "point"), ##remove margin on legend
        
        legend.text = element_text( size = 18),
        
        legend.position = c(0.8, 0.2),
        
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
                                         
                                         breaks = c(10)))

Leach.small.ET.graph.grid


Leach.small.ET.graph.comb.grid <- grid.arrange(Leach.small.ET.inf.grid, Leach.small.ET.graph.grid,
                                          
                                          heights = c(1/4, 3/4))

Leach.small.ET.graph.comb.grid




###############arrange


##"grobs"

Leach.tall.ET.graph.comb.grid
structure.study.ET.graph.comb.grid
Leach.small.ET.graph.comb.grid


blank1 <- textGrob("")
blank2 <- textGrob("")
blank3 <- textGrob("")



ET.grid <- grid.arrange(Leach.tall.ET.graph.comb.grid,
                        
                        structure.study.ET.graph.comb.grid,
                        
                        Leach.small.ET.graph.comb.grid,
                        
                        blank1,
                        
                        blank2,
                        
                        blank3,
                        
                        layout_matrix = rbind(c(1, 1, 1, 4),
                                              c(2, 2, 2, 5),
                                              c(3, 3, 3, 6)
                                              
))

ET.grid 

ggsave(file = "C:/Users/Jesse/Desktop/Neonicotinoids_15/Data/Column-Study/Raw-Data-Outputs/ET.grid.pdf", ET.grid)
ET.grid
dev.off()




