Leachate_full <- read.csv("C:/Users/Jesse/Desktop/Neonicotinoids_15/Data/Column-Study/Raw-Data-Outputs/Leachate_full_R.csv")
attach(Leachate_full)

library(dplyr)
library(ggplot2)
library(Rmisc)
library(gridExtra)
library(grid)
library(multcompView)
library(lawstat)
#library(dplyr)

Leachate_full

str(Leachate_full)

head(Leachate_full)

#colnames(Leachate_Final) <- "Sample.ID"

Leachate_Final$Sample.ID

Leachate_full

###first clean data source

###leachate volumes <0.3 = 0

#if(Leachate_full$Leach.Vol.mL <= 0.3) {
#  Leachate_full$Leach.Vol.mL <- 0
#}
  
Leachate_full$Leach.Vol.mL[Leach.Vol.mL <= 0.3] <- 0 #### this works!

Leachate_full$Leach.Vol.mL[Leach.Vol.mL == 0] <- NA   ###replace 0 with NA in volumne leachate

Leachate_full$TMX_PPB[TMX_PPB == 0] <- NA   ###replace 0 with NA in leachate conc

Leachate_full




#Leachate_full$cum.Leach.vol.mL <- for(i in  Leachate_full$Leach.Vol.mL) {
  
 #                                 if(Leachate_full$Sample.ID == "identity" & Leachate_full$Time.day == "identity") {
                                    
  #                                cumsum(i)  
                                    
   #                               } else{NA
                                      
    #                              }
       #  
        #                          }




cum.na <- function(x) { 
  
  x[which(is.na(x))] <- 0 
  
  return(cumsum(x))  ####### use this to treat na as zero in cumulative calculation. IMPORTANT!!!!!
} 

Leachate_full$cum.Leach.vol.mL <- ave(Leachate_full$Leach.Vol.mL,
                                      
                                      Leachate_full$Sample.ID, FUN=cum.na) 

### cumulative function but NAs not removed in calculation, "KF-60 cm-LS/SS-5" leached less than 400 mL total!

### ave() is a shorthand way of treating a function by group (for staced data)----just adjust function within body^^



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


?summarySE()



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

View(tall.columns.no.SS)

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





################Stats on short columns ---tmx transport

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

##now make factorial combo (eg. clay No Plant)

tall.columns.no.SS.stats.mass$Combo <- paste(tall.columns.no.SS.stats.mass$Texture,
                                   
                                             tall.columns.no.SS.stats.mass$Plant.Influence)

tall.columns.no.SS.stats.mass ##good





#########now structure study

View(structure.study)

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

small.columns.stats.mass$Letters <- ifelse(small.columns.stats$Time.day == "33", "a", NA)

small.columns.stats.mass

##now make factorial combo (eg. clay No Plant)

small.columns.stats.mass$Combo <- paste(small.columns.stats.mass$Texture,
                                             
                                        small.columns.stats.mass$Plant.Influence)

small.columns.stats.mass ##good



######done with preplot analysis


######now start plotting!!!!!!!!

####TMX_ppb

cleanup <- theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background =  element_blank(),
                 axis.line = element_line(color = "black"))

tall.columns.no.SS.stats

Leach.tall.TMX.ppb.inf <- ggplot(tall.columns.no.SS.stats, aes(x = Time.day, y= Cum.inf.cm,
                                                               
                                                               ymin = 0,
                                                               
                                                               ymax = Cum.inf.cm,
                                                               
                                                               xmax = Time.day)) +
                                 geom_line(size = 1, linetype = 2) +
                                 
                                 scale_y_continuous(limits=c(25,0),
                                                    
                                                    expand=c(0,0), trans="reverse", position = "right") +
                                 theme_classic() +
                                 
                                 theme(plot.margin = unit(c(6,21,-37,52), units="points"),
                                       
                                       axis.title.y = element_text(vjust = 0.3, size = 14, face = "bold"),
                                    
                                       axis.text.y = element_text(size = 16),
                                          
                                       panel.border = element_rect(color = "black",
                                                                   
                                                                   fill = NA, size = 1.75),
                                       
                                       axis.ticks.y = element_line(size = 1.75, color = "black")) +
                                
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
                                     
                                   theme(plot.margin = unit(c(0,56.5,1,2), units="points"), ##good legend size
                                                            
                                         panel.border = element_rect(color = "black",
                                                                     
                                                                     fill = NA, size = 1.75),
                                         
                                         axis.text.y = element_text(size = 16),  ## change font size of y axis
                                         
                                         axis.title.y = element_text(face = "bold", size = 16),
                                         
                                         axis.text.x = element_text(size = 16),  ## change font size of y axis
                                         
                                         axis.title.x = element_text(size = 16, face = "bold"),
                                         
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
                                          
                                         axis.ticks.y = element_line(size = 1.75, color = "black"),
         
                                         axis.ticks.x= element_line(size = 1.75, color = "black"))


Leach.tall.TMX.ppb.graph 

Leach.tall.TMX.ppb.graph.comb <- grid.arrange(Leach.tall.TMX.ppb.inf, Leach.tall.TMX.ppb.graph,
             
             heights = c(1/4, 3/4))

Leach.tall.TMX.ppb.graph.comb








pdf("C:/Users/Jesse/Desktop/Neonicotinoids_15/Data/Column-Study/Raw-Data-Outputs/Leach.tall.TMX.ppb.graph.pdf")                                                          
Leach.tall.TMX.ppb.graph  
dev.off()


###have to use ggsave for grid obsjects
ggsave(file = "C:/Users/Jesse/Desktop/Neonicotinoids_15/Data/Column-Study/Raw-Data-Outputs/Leach.tall.TMX.ppb.graph.comb.pdf", Leach.tall.TMX.ppb.graph.comb)


?geom_point
?geom_smooth
?pch


balls <- ggplot(tall.columns.no.SS.stats, aes(Time.day, TMX_PPB)) +
  
                geom_point(size = 3)
balls

tall.columns.no.SS.stats # 

small.columns.stats.mass ###NA not the problem, when plotting in base r




