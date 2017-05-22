rm(list = ls())

Leachate_Final <- read.csv("C:/Users/Jesse/Desktop/Neonicotinoids_15/Data/Column-Study/Raw-Data-Outputs/Final-Leachate-R.csv")

attach(Leachate_Final)



library(lawstat)

View(Leachate_Final)

hist(Leachate_Final$TMX_PPB)

short_columns <- Leachate_Final[Leachate_Final$Size_cm == "20", c("Sample.ID","Structure", "Texture", "Plant.Influence", "TMX_PPB")] ## subsetting 20 cm columns

short_columns

tall_columns <- Leachate_Final[Leachate_Final$Size_cm == "60", c("Sample.ID","Structure", "Texture", "Plant.Influence", "TMX_PPB")] #subsetting 60 cm columns

tall_columns

structural_analysis <- Leachate_Final[ Leachate_Final$Size_cm == "60" & Leachate_Final$Structure == "S" | Leachate_Final$Structure == "NS" | Leachate_Final != "No Plant",  c("Size_cm", "Sample.ID","Structure", "Texture", "Plant.Influence", "TMX_PPB")]

# subseting the structure study

structural_analysis

structural_analysis <- na.omit(structural_analysis) # omiting NA

structural_analysis

structural_analysis1 <- subset(structural_analysis, structural_analysis$Plant.Influence == "Plant") ## exclude "No Plant"

structural_analysis1

hist(structural_analysis1$TMX_PPB) # distribution in structure study

qqnorm(structural_analysis1$TMX_PPB) ## not normal

qqline(structural_analysis1$TMX_PPB)

bartlett.test(structural_analysis1$TMX_PPB ~ structural_analysis1$Structure)  ## no HOV

fligner.test(structural_analysi1s$TMX_PPB ~ structural_analysis1$Structure)

levene.test(structural_analysis1$TMX_PPB, structural_analysis1$Structure)

boxplot(structural_analysis1$TMX_PPB ~ structural_analysis1$Structure, font.main = 2)

boxplot(structural_analysis1$TMX_PPB ~ structural_analysis1$Structure, ylab = "TMX ng mL^-1", xlab = "TRT")

View(structural_analysis1)

kruskal.test(structural_analysis1$TMX_PPB ~ structural_analysis1$Structure)## no difference


####Try to remove outliers?

structural_analysis1

hist(structural_analysis1$TMX_PPB)

subset(structural_analysis1, TMX_PPB < 0.86 | TMX_PPB > 0.88)

## have to make cut offs just below and just above value--- can't set the cutoffs at that value

subset(structural_analysis1, TMX_PPB < 0.87 | TMX_PPB > 0.87) ### doesn't work   


subset(structural_analysis2,  TMX_PPB != 0.87) ### doesn't work  

structural_analysis1[structural_analysis1$TMX_PPB != 0.87, c("Size_cm", "Sample.ID","Structure", "Texture", "Plant.Influence", "TMX_PPB")] ## doesn't work

###structural_analysis2 <- subset(structural_analysis1, TMX_PPB != 0.87)


structural_analysis2 <- subset(structural_analysis1, Sample.ID != "KF-60 cm-LS-Exp-4") ## or this--call the obs code and omit




structural_analysis2 ## good, now redo analysis (remove outlier in Structured clay "KF-60 cm-LS-Exp-4")

boxplot(structural_analysis2$TMX_PPB ~ structural_analysis2$Structure, ylab = "TMX ng mL^-1", xlab = "TRT")

kruskal.test(structural_analysis2$TMX_PPB ~ structural_analysis2$Structure) ## difference not noticeable? P =0.052

##try to run rank transformation and do one way anova

structural_analysis2$Rank.TMX_PPB <- rank(structural_analysis2$TMX_PPB)

aov.str <- aov(structural_analysis2$Rank.TMX_PPB ~ structural_analysis2$Structure)

summary.str <- summary(aov.str)

summary.str ## p=0.038 however the boxplot of raw data would challange this in a paper---more appropriate to remove both outliers--see below


hist(structural_analysis2$TMX_PPB) ## not normal 

qqnorm(structural_analysis2$TMX_PPB) ## not normal

qqline(structural_analysis2$TMX_PPB) ## not normal

bartlett.test(structural_analysis2$TMX_PPB ~ structural_analysis2$Structure) ## meets HOV

fligner.test(structural_analysis2$TMX_PPB ~ structural_analysis2$Structure) ## meets HOV

levene.test(structural_analysis2$TMX_PPB, structural_analysis2$Structure) ## meets HOV




#################################################
structural_analysis3 <- subset(structural_analysis2, Sample.ID != "KF-60 cm-LS/SS-4") ### remove outlier in Structured and unstructured clay

## outlier in unstructured soil "KF-60 cm-LS/SS-4" likely due to edge flow development. 

## remove outlier "KF-60 cm-LS-Exp-4" --- showed lowest cumulative leaching pattern in terms of water volume, [TMX], and showed a decreasing trend in leachate concentration when analyzed by container during the final event

library(markdown)

structural_analysis3

boxplot(structural_analysis3$TMX_PPB ~ structural_analysis3$Structure, ylab = "TMX ng mL^-1", xlab = "TRT") ## difference not noticeable

hist(structural_analysis3$TMX_PPB) ## not normal 

qqnorm(structural_analysis3$TMX_PPB) ## not normal

qqline(structural_analysis3$TMX_PPB) ## not normal

bartlett.test(structural_analysis3$TMX_PPB ~ structural_analysis3$Structure) ## meets HOV

fligner.test(structural_analysis3$TMX_PPB ~ structural_analysis3$Structure) ## meets HOV

levene.test(structural_analysis3$TMX_PPB, structural_analysis3$Structure) ## meets HOV


###Try transforming

structural_analysis3$log.TMX <- log10(structural_analysis3$TMX_PPB) # log base 10

hist(structural_analysis3$log.TMX) ##  not really normal

qqnorm(structural_analysis3$log.TMX) ## not really normal

qqline(structural_analysis3$log.TMX) ## not really normal

##########################just use KW test

kruskal.test(structural_analysis3$TMX_PPB ~ structural_analysis3$Structure) ## difference!







## what if we just remove the outler in unstructured soil?

structural_analysis4 <- subset(structural_analysis4, Sample.ID != "KF-60 cm-LS/SS-4") ## removed "KF-60 cm-LS/SS-4"

boxplot(structural_analysis4$TMX_PPB ~ structural_analysis4$Structure, ylab = "TMX ng mL^-1", xlab = "TRT") ## difference not apparent

kruskal.test(structural_analysis4$TMX_PPB ~ structural_analysis4$Structure) ## difference not apparent

hist(structural_analysis4$TMX_PPB) ## 

qqnorm(structural_analysis4$TMX_PPB) ## 

qqline(structural_analysis4$TMX_PPB) ## 

bartlett.test(structural_analysis$TMX_PPB ~ structural_analysis4$Structure) ## no HOV

fligner.test(structural_analysis4$TMX_PPB ~ structural_analysis4$Structure) ## no HOV

levene.test(structural_analysis4$TMX_PPB, structural_analysis4$Structure) ## no HOV

Leachate_Final




###Stats on 20 cm columns first

short_columns

## problem outliers >>> "NK-20 cm-LS-Exp-4," or column 4-15, only leached 2 containers (950 mL) for final>> hydraulic connection issue? 
## also took abnormally long time to yield that drainage

short_columns

short_columns <- subset(short_columns, Sample.ID != "NK-20 cm-LS-Exp-4")

short_columns


### Note these are two way interactions 


boxplot(short_columns$TMX_PPB) ## 


hist(short_columns$TMX_PPB) ## 

qqnorm(short_columns$TMX_PPB) ## 

qqline(short_columns$TMX_PPB) ## ~normal


bartlett.test(short_columns$TMX_PPB ~ interaction(short_columns$Texture, short_columns$Plant.Influence)) ## meets HOV

fligner.test(short_columns$TMX_PPB ~ interaction(short_columns$Texture, short_columns$Plant.Influence)) ## no HOV

levene.test(short_columns$TMX_PPB ~ short_columns$Texture * short_columns$Plant.Influence) ## 

aov <- aov(short_columns$TMX_PPB ~ short_columns$Texture * short_columns$Plant.Influence)

summary <- summary(aov)

summary ## Texture difference, Plant effect not important

boxplot(short_columns$TMX_PPB ~ short_columns$Texture * short_columns$Plant.Influence, ylab = "TMX (ng mL^-1)", xlab = "Trt")

Tukey <- TukeyHSD(aov)

Tukey

library(multcompView)

Letters.small <- multcompLetters4(aov,Tukey) ##Sand:Plant Sand:No Plant Clay:No Plant    Clay:Plant 
#####################                           "a"          "ab"          "ab"           "b"
Letters.small


####should really transform 20 cm columns since no HOV fligners

short_columns$log10.TMX <- log10(short_columns$TMX_PPB)

hist(short_columns$log10.TMX) ## normal

qqnorm(short_columns$log10.TMX) ## 

qqline(short_columns$log10.TMX) ##normal


fligner.test(short_columns$log10.TMX ~ interaction(short_columns$Texture, short_columns$Plant.Influence))## meets HOV

aov_small_columns <- aov(short_columns$log10.TMX ~ short_columns$Texture * short_columns$Plant.Influence)

summary_small_columns <- summary(aov_small_columns)

summary_small_columns ## Texture difference, Plant effect not important

#boxplot(short_columns$TMX_PPB ~ short_columns$Texture * short_columns$Plant.Influence, ylab = "TMX (ng mL^-1)", xlab = "Trt")

Tukey_small_columns <- TukeyHSD(aov_small_columns)

Tukey_small_columns

library(multcompView)

Letters.small1 <- multcompLetters4(aov_small_columns,Tukey_small_columns)

Letters.small1

###same result as with untransformed data

######just note that log trNSFORMED DATA was used on 20 cm columns



#### what if no outlier in Clay no plant

#short_columns1 <- subset(short_columns, Sample.ID != "KF-20 cm-LS-Ctr-4") 

#boxplot(short_columns1$TMX_PPB ~ short_columns1$Texture * short_columns1$Plant.Influence, ylab = "TMX (ng mL^-1)", xlab = "Trt")

#aov1 <- aov(short_columns1$TMX_PPB ~ short_columns1$Texture * short_columns$Plant.Influence)

#summary1 <- summary(aov2)

#summary1








### now the rest of the 60 cm columns---- Two way interaction again

tall_columns

tall_columns1 <- subset(tall_columns, Structure != "NS"| Texture == "Sand")

tall_columns1 

boxplot(tall_columns1$TMX_PPB ~ tall_columns1$Texture * tall_columns1$Plant.Influence, ylab = "TMX (ng mL^-1)", xlab = "Trt")

## remove outlier "KF-60 cm-LS-Exp-4" --- showed lowest cumulative leaching pattern in terms of water volume, [TMX], and showed a decreasing trend in leachate concentration when analyzed by container during the final event

tall_columns2 <- subset(tall_columns1, Sample.ID != "KF-60 cm-LS-Exp-4")

tall_columns2

boxplot(tall_columns2$TMX_PPB ~ tall_columns2$Texture * tall_columns2$Plant.Influence, ylab = "TMX (ng mL^-1)", xlab = "Trt")

hist(tall_columns2$TMX_PPB) ## 

qqnorm(tall_columns2$TMX_PPB) ## 

qqline(tall_columns2$TMX_PPB) ## definitely not normal

 tall_columns2$rank.TMX <- rank(tall_columns2$TMX_PPB)## rank transform data to normality
 

 boxplot(tall_columns2$rank.TMX ~ tall_columns2$Texture * tall_columns2$Plant.Influence, ylab = "TMX (ng mL^-1)", xlab = "Trt")## not intuitive
 
 bartlett.test(tall_columns2$rank.TMX ~ interaction(tall_columns2$Texture, tall_columns2$Plant.Influence)) ## HOV met
 
 fligner.test(tall_columns2$rank.TMX ~ interaction(tall_columns2$Texture, tall_columns2$Plant.Influence)) ## HOV met
 
 levene.test(tall_columns2$rank.TMX ~ tall_columns2$Texture * tall_columns2$Plant.Influence) ## 
 
 ANOVA.tall <- aov(tall_columns2$rank.TMX ~ tall_columns2$Texture * tall_columns2$Plant.Influence)
 
 summary.tall <- summary(ANOVA.tall)
 
 summary.tall ## Texture difference, Plant effect not important
 
 Tukey.tall <- TukeyHSD(ANOVA.tall)
 
 Tukey.tall
 
 Letters.tall.leach <- multcompLetters4(ANOVA.tall, Tukey.tall)
 
tall_columns2


Leachate_Final


