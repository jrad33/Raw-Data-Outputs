rm(list = ls())

#install.packages("data.table")
library(data.table)

#install.packages("xlsx")
library(xlsx)
library(dplyr)
library(ggplot2)
library(Rmisc)
library(gridExtra)
library(grid)
library(multcompView)
library(lawstat)

#JBR
##Ponded infiltration tests run on columns using 300 mL increments (0.9 cm).
##tests run in November 2015

Column.infiltration <- read.csv("C:/Users/Jesse/Desktop/Neonicotinoids_15/Data/Column-Study/Raw-Data-Outputs/Infiltration_R.csv")
attach <- Column.infiltration

Column.infiltration





###################################reference code for for loops within column...didn't really work though####

#hours <- c()
#for(i in 0:length(Column.infiltration$time.elapsed_hr)) {
  #print(Column.infiltration$time.elapsed_hr[i] - Column.infiltration$time.elapsed_hr[i-1])
  #print(Column.infiltration$time.elapsed_hr[i])
  
#if(Column.infiltration$time.elapsed_hr == 0){
# hours[i] <- 0
    #}
    #else{       
    #   hours[i] <- Column.infiltration$time.elapsed_hr[i] - Column.infiltration$time.elapsed_hr[i-1]
         #}
  
    #}



    #Column.infiltration$dt_hr <- diff(Column.infiltration$time.elapsed_hr[1],
                                  
    #Column.infiltration$time.elapsed_hr[2:length((Column.infiltration$time.elapsed_hr))])

    #diff(Column.infiltration$time.elapsed_hr, lag = 2)

    #Column.infiltration



##convert time from seconds to hr


#########################################################################end of ref code#######################





###initial calculations

Column.infiltration$time.elapsed_hr <- Column.infiltration$time.elapsed_seconds/(60^2)

Column.infiltration$time.elapsed_hr

Column.infiltration



##calculate dt for each time step (hr)


Column.infiltration$dt_hr <- ave(Column.infiltration$time.elapsed_hr, ####calculates dt per time step modifying ave()
                                 
                                 factor(Column.infiltration$Column.ID), 
                                 
                                 FUN = function(x) c(0, diff(x)))  ##weird syntax makes a function that first 

 ##makes column = 0, then runs diff()


Column.infiltration$dt_hr

Column.infiltration


###infiltration rate

h <- 0.9 ### cm of head added per increment

dt <- Column.infiltration$dt_hr

Column.infiltration$i_cm.hr <- h/dt ###i in cm per hr

Column.infiltration


?cumsum()
?ave()


###calculate I (cm)

Column.infiltration$I <- h * (Column.infiltration$Point-1)

Column.infiltration

###calculate dI (cm)



##sqrt of t (cm)

Column.infiltration$sqrt.t <- sqrt(Column.infiltration$time.elapsed_hr) 

Column.infiltration



Column.infiltration$dsqrt.t <- ave(Column.infiltration$sqrt.t, ####calculates dsqrt.t per time step modifying ave()
                                
                                factor(Column.infiltration$Column.ID), 
                                
                                FUN = function(x) c(0, diff(x)))

Column.infiltration





###calculate dI (cm)

Column.infiltration$dI <- ave(Column.infiltration$I, ####calculates dsqrt.t per time step modifying ave()
                                   
                                   factor(Column.infiltration$Column.ID), 
                                   
                                   FUN = function(x) c(0, diff(x)))
Column.infiltration


###calculate I/sqrt.t (some funky units... cm hr^-0.5)

Column.infiltration$I_sqrt.t_ratio <- Column.infiltration$I/Column.infiltration$sqrt.t

Column.infiltration


###calculate dI/dsqrt.t (some funky units... cm hr^-0.5)

Column.infiltration$dI_dsqrt.t_ratio <- Column.infiltration$dI/Column.infiltration$dsqrt.t

Column.infiltration$dI_dsqrt.t_ratio

Column.infiltration
















############# based on equations from Vandervaere et al (2000) and using Phillip's 1-D soln to get A :

#fit lm to each rep: either as I = S*sqrt(t) + At

#  (I/sqrt(t)) = S +A*sqrt(t) .... plotting (I/sqrt(t)) a function of sqrt(t)

# where A ~ Ks/3  or Ks = 3A and A ~ slope of line Philip (1990)




#or the differential form: (dI/dsqrt(t)) = S + 2A*sqrt(t) .... plotting (dI/dsqrt(t)) as a function of dsqrt(t)
 
## with the differential form A ~ 2/3*Ks or Ks = A/0.67 and A ~ slope of line Philip (1990)


plot(Column.infiltration$sqrt.t, Column.infiltration$dI_dsqrt.t_ratio)


###20 cm sand (2_5) observation 4 is an outlier.. too quick either an error or there was a large pref flow pathway
#### 20 

plot(Column.infiltration$sqrt.t, Column.infiltration$I_sqrt.t_ratio) 

#### should be positive relationship, problem??


Column.infiltration

########################################################################################################


####use for loop to fit linear model, extract slope, and calcualte Ks*****






###start with (I/t^0.5 v. t^0.5)*******


# R extracts the numeric value of "factor" levels

#for reach column assigned (e.g Sand == 2), so need to treat columns as characters in the loop

##then reassign as factors afterwards....

Column.infiltration$Texture <- as.character(Column.infiltration$Texture)
Column.infiltration$Size <- as.character(Column.infiltration$Size)
Column.infiltration$Structure <- as.character(Column.infiltration$Structure)

lm_df_1 <- data.frame("Column.ID" = levels(Column.infiltration$Column.ID),
                      "Texture" = NA,
                      "Size" = NA,
                      "Structure" = NA,
                      "Intercept" = NA,
                      "Slope" = NA)                ###make dummy df with levels()

keep_cols <- c("Texture", "Size", "Structure")     ##vector of column names that we want to keep from orginal



for (i in levels(Column.infiltration$Column.ID)) { ## levels() is sweet... lets you treat data as factor
  
  temp_df_1 <- Column.infiltration[Column.infiltration$Column.ID == i, ]
  
  ##temporary df allows for lm() and coef() to be fed to the loop sequentially, if we just used the OG df then
  ##loop would freeze at some value
  #lm_df_1[lm_df_1$Column.ID == i, 2] <- temp_df_1[levels(temp_df_1$Texture), 2]
  
  
  ###makes df with column written with column.id and two regression coefficients.. coef() reads the output vector
  lm_df_1[lm_df_1$Column.ID == i, c(5,6)] <- coef(lm(temp_df_1$I_sqrt.t_ratio ~ temp_df_1$sqrt.t)) 
  
  
  ## R extracts the numeric value of "factor" levels
  
  #for reach column assigned (e.g Sand == 2), so need to treat columns as characters in the loop
  
  lm_df_1[lm_df_1$Column.ID == i, keep_cols] <- temp_df_1[1, keep_cols] 
  
}


lm_df_1





####convert originial df columns back to factors (useful for summary statistic analyses)

Column.infiltration$Texture <- as.factor(Column.infiltration$Texture)
Column.infiltration$Size <- as.factor(Column.infiltration$Size)
Column.infiltration$Structure <- as.factor(Column.infiltration$Structure)

####convert  columns back to factors (useful for summary statistic analyses)

lm_df_1$Texture <- as.factor(lm_df_1$Texture)
lm_df_1$Size <- as.factor(lm_df_1$Size)
lm_df_1$Structure <- as.factor(lm_df_1$Structure)

str(lm_df_1)





#ifelse(temp_df_1$Texture == "Clay", lm_df_1$Texture[i] == "Clay",
       
#     ifelse(temp_df_1$Texture == "Sand", lm_df_1$Texture[i] == "Clay",
              
#      ifelse(temp_df_1$Size == "20 cm", lm_df_1$Size[i] == "20 cm",
                     
#     ifelse(temp_df_1$Size == "60 cm", lm_df_1$Size[i] == "60 cm",
                            
#     ifelse(temp_df_1$Structure == "Structure", lm_df_1$Structure[i] == "Structure",
                                   
#    ifelse(temp_df_1$Structure == "Structureless", lm_df_1$Structure[i] == "Structureless",
                                   
#     ifelse(temp_df_1$Structure == "NA", lm_df_1$Structure[i] == "NA", NA)
                                   
#))))))





lm_df_1 ### negative slopes in 2_8, 3-1, 4_2, and 4_8 is screwed

##toss 2_8 and 4_2... negative slope means ~ bad contact and exp error according to (Vandervaere, 2000)



### also according to this guy, you can take column 3_1 and fit to the linear zone (exclude first point)


plot(Column.infiltration$I_sqrt.t_ratio[Column.infiltration$Column.ID == "3_1"] ~ ###find linear zone for 3_1
       Column.infiltration$sqrt.t[Column.infiltration$Column.ID == "3_1"])

Column.infiltration[Column.infiltration$Column.ID == "3_1", ] ###rows 90-95

column.3_1 <- Column.infiltration[Column.infiltration$Column.ID == "3_1", ]

column.3_1.fix <- column.3_1[3:8,] ##rows 3:8 in here

column.3_1.fix


plot(column.3_1.fix$I_sqrt.t_ratio ~ ##good
       column.3_1.fix$sqrt.t)

ceof <- coef(lm(column.3_1.fix$I_sqrt.t_ratio ~
           column.3_1.fix$sqrt.t))

lm_df_1$Intercept[lm_df_1$Column.ID == "3_1"] <- ceof[1] ##change intercept

lm_df_1$Slope[lm_df_1$Column.ID == "3_1"] <- ceof[2] ##change slope

lm_df_1


###calculate Ks where A ~ Ks/3  or Ks = 3A and A ~ slope of line Philip (1990)

lm_df_1$Ks <- lm_df_1$Slope*3 ## cm/hr

lm_df_1$Ks.cm.day <-lm_df_1$Ks *24 

###why do 2_8 and 4_2 have a negative slope?

## postive slopes here seem to be matched with small number of points and high dI values





###toss 2_8 and 4_2 due to negative slope (possible exp. error or something wrong with column packing?)
##also toss 4_15 due to odd experimental conditions (was also dropped for leachate analysis)


lm_df_1 <- subset(lm_df_1, Column.ID != "2_8" & Column.ID != "4_2" & Column.ID != "4_15")

lm_df_1

### high variation in 20 cm sand (4-15) and repacked (structureless clay)







##summarize
sum1 <- summarySE(lm_df_1, measurevar = "Ks", groupvars = c("Texture", "Size", "Structure"))
sum1

sum1.5 <- summarySE(lm_df_1, measurevar = "Ks.cm.day", groupvars = c("Texture", "Size", "Structure"))


sum1.5 

###values are high! and variable ....but they show important differences (e.g. the sands show the highest Ks,

#                                                                         unstructured clay is lowest etc.)

#export to an excel file

write.xlsx(sum1.5 , "C:/Users/Jesse/Desktop/Neonicotinoids_15/Data/Column-Study/Raw-Data-Outputs/Infiltration_R_table_output.xlsx")






#######now the differential function... (dI/dsqrt(t)/dsqrt(t)) vs sqrt(t))#########################



Column.infiltration$Texture <- as.character(Column.infiltration$Texture)
Column.infiltration$Size <- as.character(Column.infiltration$Size)
Column.infiltration$Structure <- as.character(Column.infiltration$Structure)

lm_df_2 <- data.frame("Column.ID" = levels(Column.infiltration$Column.ID),
                      "Texture" = NA,
                      "Size" = NA,
                      "Structure" = NA,
                      "Intercept" = NA,
                      "Slope" = NA)                ###make dummy df with levels()

keep_cols <- c("Texture", "Size", "Structure")     ##vector of column names that we want to keep from orginal



for (i in levels(Column.infiltration$Column.ID)) { ## levels() is sweet... lets you treat data as factor
  
  temp_df_2 <- Column.infiltration[Column.infiltration$Column.ID == i, ]
  
  ##temporary df allows for lm() and coef() to be fed to the loop sequentially, if we just used the OG df then
  ##loop would freeze at some value
  #lm_df_1[lm_df_1$Column.ID == i, 2] <- temp_df_1[levels(temp_df_1$Texture), 2]
  
  
  ###makes df with column written with column.id and two regression coefficients.. coef() reads the output vector
  lm_df_2[lm_df_1$Column.ID == i, c(5,6)] <- coef(lm(temp_df_2$dI_dsqrt.t_ratio ~ temp_df_2$sqrt.t)) 
  
  
  ## R extracts the numeric value of "factor" levels
  
  #for reach column assigned (e.g Sand == 2), so need to treat columns as characters in the loop
  
  lm_df_2[lm_df_2$Column.ID == i, keep_cols] <- temp_df_2[1, keep_cols] 
  
}


lm_df_2

###### many negative slopes possible because:- dI remains constant,
##                                           - steady state was never reached
##                                           -few point measuremnts 
#                                             
##                                                 ???















