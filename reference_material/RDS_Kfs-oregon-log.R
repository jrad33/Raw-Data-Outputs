
#install.packages(c("stats", "zoo", "survival",
                   
#                "lattice", "psych", "latticeExtra", "zoo", "hydroGOF", "qualityTools"))

require(stats); require(graphics); require(zoo) #; require(splinefun)
library(survival); library(plyr); library(lattice); library(psych)
library(latticeExtra); library(zoo); library(hydroGOF); library(Hmisc);  library(qualityTools)
# library(minpack.lm); library(reshape2)


### Sorptivity Method Comparison                
### Ryan D. Stewart, January 2013
###

rm(list=ls())

setwd("D:/Research/OSU/Chile Project/Data/Infiltration-Data/Beerkan/Lewisburg")
#setwd("C:/Documents and Settings/stewarry/My Documents/Research/Chile Project/Data/Infiltration-Data/Beerkan/Lewisburg")
#setwd('/Users/drupp/Rupp/Chile/R-Scripts')
# setwd("C:/Users/Ryan/Documents/Soil_shrinkage")


########
# INPUT DATA
########


p = 2.96 # fitting parameter
q = 2.16 # fitting parameter

Ksaggrmax = 3.6 # cm/hr
Kscrackmax = 2900 # cm/hr
Ksaggrmax3D = 0.1 # cm/hr
Kscrackmax3D = 190 # cm/hr
Bcrack = 1. # Parameter for crack geometries
Baggr = 1 # Parameter for pore geometries

Ks1D = 34 #cm/hr - This is Kip+ia
Ks3D = 1
oneDcoef = 0.35 # Coefficient to convert C2 into Ks
alph = 0.4
phimax = 0.604
phimin = 0.19

rdisk = 4.8 #cm
hf = 10 #cm

### Output
#outputfilename = "fieldsorptivity_output.csv"

############
###########
# Code starts here
###


### New Ks model!

# 1D model

Ufunc = seq(0,1,.01)
phipedon = phimax - 1 + (1 - (phimax-phimin))^(1/3) # Phi Pedon assuming isotropic shrinkage
phi_crack = (phipedon-phimin)*((1-Ufunc^q)/(1+p*Ufunc^q))
phi_sub = (phimax-phipedon)*((1-Ufunc^q)/(1+p*Ufunc^q))

# 1D model

Ksfunc1 = Kscrackmax*(phi_crack/(1-phi_sub))*((1-Ufunc^q)/(1+p*Ufunc^q))^Bcrack
Ksfunc2 = Ksaggrmax*(1-(phi_crack/(1-phi_sub)))*((p+1)/(p+Ufunc^-q))^Baggr
Ksfunc = Ksfunc1+Ksfunc2

# 3D model

Ksfunc3D1 = Kscrackmax3D*(phi_crack/(1-phi_sub))*((1-Ufunc^q)/(1+p*Ufunc^q))^Bcrack
Ksfunc3D2 = Ksaggrmax3D*(1-(phi_crack/(1-phi_sub)))*((p+1)/(p+Ufunc^-q))^Baggr
Ksfunc3D = Ksfunc3D1+Ksfunc3D2


## Read data files #############################################################                                  

read.csv("C:/Users/Jesse/Desktop/Neonicotinoids_15/Data/Column-Study/Raw-Data-Outputs/reference_material/field_sorptivity_gravimetric.csv")
z <- read.table(filename, as.is=TRUE, sep = ",", header = FALSE, strip.white = TRUE, comment.char = "#")

column.names = c("Ring #","Date", "C1 (I = C1(t^0.5) + C2(t))","C1 (I*(t^-0.5) v. t^0.5)","C1 (dI*(dt^-0.5) v. t^0.5)","Short-time (S = I(t^0.5))", "C2 (I = C1(t^0.5) + C2(t))","C2 (I*(t^-0.5) v. t^0.5)","C2 (dI*(dt^-0.5) v. t^0.5)", "Time0", "Time1", "Time2", "Time3", "Time4", "Time5", "Time6", "Time7", "Time8", "Time9", "Time10")
colnames(z) = column.names

# Get sorptivity for each approach
Ring = z[,1]
Date = z[,2]
S1 = z[,3]
S2 = z[,4]
S3 = z[,5]
S4 = z[,6]
K1 = z[,7]
K2 = z[,8]
K3 = z[,9]
t0 = z[,10]
t1 = z[,11]
t2 = z[,12]
t3 = z[,13]
t4 = z[,14]
t5 = z[,15]
t6 = z[,16]
t7 = z[,17]
t8 = z[,18]
t9 = z[,19]
t10 = z[,20]
gwc_8cm = z[,21]
vwc_8cm = z[,22]
moistureratio_8cm = z[,23]
specvolume_8cm = z[,24]
bulkdensity_8cm = z[,25]
voidratio_8cm = z[,26]
degreesat_8cm = z[,27]
gwc_12cm = z[,28]
vwc_12cm = z[,29]
moistureratio_12cm = z[,30]
specvolume_12cm = z[,31]
bulkdensity_12cm = z[,32]
voidratio_12cm = z[,33]
degreesat_12cm = z[,34]
u_avg = z[,35]
vwc_avg = z[,36]
moistureratio_avg = z[,37]
specvolume_avg = z[,38]
bulkdensity_avg = z[,39]
voidratio_avg = z[,40]
degreesat_avg = z[,41]

Size <- length(Ring)

a = 1.025 # G&A correction factor


# umax = max(gwc_avg)+.013 # maximum gravimetric water content (g/g)

soliddens = 2.65 #Mg/m3

waterdens = 1.0 #Mg/m3

umax = phimax/(1-phimax)*(waterdens/soliddens)

U = u_avg/umax

###########

# Infiltration Models

rdisc = 4.8 #cm

Iactual = c(seq(0,10,1)) #cm3
Iactual = Iactual*100/(pi*(rdisc^2)) #cm
Iactual_mat = matrix(data = Iactual, nrow = Size, ncol = 11, byrow = TRUE, dimnames = NULL)
Iactual_mat = round(Iactual_mat, 2)

time = z[,10:20]  # sec
sqrt_time = time^0.5  #sec^0.5
x = sqrt_time
is.na(x)
sqrt_time = ifelse(is.na(x), 0, x)
time = sqrt_time^2


# Calculate coefficients using [R]

### Smiles (I/t^0.5 v. t^0.5)

Smiles_vec = rep(1,Size)

I = Iactual_mat #cm
t05 = x #sec^0.5

Sm = I/t05
is.na(Sm)
Sm = ifelse(is.na(Sm), 0, Sm)


for (i in 1:Size) {
  
  Smiles_vec[i] = lm(Sm[i,2:11] ~ t05[i,2:11]) 
    
}

Smiles <- data.frame(matrix(unlist(Smiles_vec), nrow=Size, byrow=T))

#Vander_mat = matrix(data = Vander_vec, nrow = Size, ncol = 1, byrow = TRUE, dimnames = NULL)

Sm_int = Smiles[,1]   #cm s^-0.5

Sm_C1 = Sm_int * 60  #cm h^-0.5

Sm_slope = Smiles[,2]  #cm s^-1

Sm_C2 = Sm_slope * 3600 #cm h^-1

Ks_total <- Sm_C2/oneDcoef   #cm h^-1


### Data Frame

RemoveOct <- c('9/23/2011','10/20/2011', '11/28/2011', '2/6/2012', '3/8/2012', '4/9/2012', '5/8/2012', '6/3/2012','7/5/2012', '8/7/2012', '9/6/2012', '10/4/2012', '11/6/2012', '12/10/2012', '1/6/2013', '2/4/2013', '3/4/2013') 


# Date[Date %in% '10/3/2012'] = '10/4/2012'  ### Note: This is used to recombine 10/3 and 10/4 data into a single set. Could be removed to separate them again. 


Sorp_theta <- data.frame(degreesat_avg, Date, Ring, Sm_C1, Sm_C2, Ks_total)

Sorp_theta1 <- Sorp_theta

# Sorp_theta1 <- subset(Sorp_theta, Date %in% RemoveOct)


Sorp_theta1$Date <- as.Date(Sorp_theta1$Date, "%m/%d/%Y")

Sorp_theta_pos <- Sorp_theta1[Sorp_theta1$Sm_C1 > 0,]

Sm_C1_stats = ddply(Sorp_theta_pos, c('degreesat_avg'), function(echs) c(count=nrow(echs), mean=mean(echs$Sm_C1), geomean=geometric.mean(echs$Sm_C1), median=median(echs$Sm_C1), sd=sd(echs$Sm_C1)))

Sm_C2_stats = ddply(Sorp_theta_pos, c('degreesat_avg'), function(echs) c(count=nrow(echs), mean=mean(echs$Sm_C2), geomean=geometric.mean(echs$Sm_C2),median=median(echs$Sm_C2), sd=sd(echs$Sm_C2)))

Sm_Ktot_stats = ddply(Sorp_theta_pos, c('degreesat_avg'), function(echs) c(count=nrow(echs), mean=mean(echs$Ks_total), geomean=geometric.mean(echs$Ks_total),median=median(echs$Ks_total), sd=sd(echs$Ks_total)))

#### Compare Sorptivity Calculations
### Against Date of sampling for all data!!!



Sm_Date = ddply(Sorp_theta_pos, c('Date'), function(echs) c(count=nrow(echs), mean=mean(echs$Sm_C1), geomean=geometric.mean(echs$Sm_C1),median=median(echs$Sm_C1), sd=sd(echs$Sm_C1)))

C2_Date = ddply(Sorp_theta_pos, c('Date'), function(echs) c(count=nrow(echs), mean=mean(echs$Sm_C2), geomean=geometric.mean(echs$Sm_C2),median=median(echs$Sm_C2), sd=sd(echs$Sm_C2)))

K_Date = ddply(Sorp_theta_pos, c('Date'), function(echs) c(count=nrow(echs), mean=mean(echs$Ks_total), geomean=geometric.mean(echs$Ks_total),median=median(echs$Ks_total), sd=sd(echs$Ks_total)))


Date1 <- as.Date("8/1/2011", "%m/%d/%Y")

Date2 <- as.Date("4/1/2013", "%m/%d/%Y")

Date3 <- as.Date("9/1/2011", "%m/%d/%Y")

theXDates <- c(Date1,Date2)

Polydate1 <- as.Date("5/20/2012", "%m/%d/%Y")

Polydate2 <- as.Date("10/20/2012", "%m/%d/%Y")


Polyx <- c(Polydate1,Polydate1,Polydate2,Polydate2)

Polyy <- c(0.1,100,100,0.1)

Polycolor <- rgb(190, 190, 190, alpha=150, maxColorValue=255)



PolyxK <- c(Polydate1,Polydate1,Polydate2,Polydate2)

PolyyK <- c(1,5000,5000,1)

Polycolor <- rgb(190, 190, 190, alpha=150, maxColorValue=255)

theYlab212 = expression(paste(bold(S~~bgroup("(",cm~~hr^{-0.5},")")))) 
theXlab212 = expression(paste(bold("Date of Sampling")))

theYlab222 = expression(paste(bold(C~~bgroup("(",cm~~hr^{-1},")")))) 
theXlab222 = expression(paste(bold("Date of Sampling")))

win.graph(width = 18, height = 18)
par(mfrow = c(2, 1))
par(oma=c(0,0,1,0)) 
par(mar=c(6,6,2,3))

yearDate <- as.Date(paste(rep(2011:2013, each = 12), rep(1:12, 3), 1, sep = "-"))
yearDate <- zoo(c(0,0), yearDate)
plot(Sorp_theta_pos$Sm_C1 ~ Sorp_theta_pos$Date, xlab = theXlab212,
     log="y", xaxt="n", ann=FALSE, pch=21, lwd=1, 
     ylab = theYlab212, 
     ylim=c(.1,100), xaxp=c(Date1,Date2,21), yaxs='i', xaxs='i', xlim=theXDates, 
     cex.lab=1.5, cex.axis=1.3,)
axis(1, at = time(yearDate), labels = FALSE)
tt <- time(yearDate)
ix <- seq(1, length(tt)) #every month
fmt <- "%b-%Y" # format for axis labels
labs <- format(tt[ix], fmt)
axis(side = 1, at = tt[ix], labels = labs,  cex.axis = 1.2, format= "%m/%Y", las = 2)  
mtext(theYlab212, side=2, line = 3, cex=1.5)
#   times <- time(yearDate)
#   ticks <- seq(times[1], times[length(times)], by = "months")
#   axis(1, at = ticks, labels = FALSE, tcl = -0.3)
polygon(Polyx, Polyy, col=Polycolor, border = "black", lty="dotted")
points(Sorp_theta_pos$Sm_C1 ~ Sorp_theta_pos$Date, pch=21, col="black" )
#points(Sorp_theta_pos$W_C1 ~ Sorp_theta_pos$Date, pch=21, col="black")
points(Sm_Date$mean ~ Sm_Date$Date, pch=21, col="black", bg="lightblue", cex=2,)
text(Date3,50,"a)",cex=1.5,font=2)



yearDate <- as.Date(paste(rep(2011:2013, each = 12), rep(1:12, 3), 1, sep = "-"))
yearDate <- zoo(c(0,0), yearDate)
plot(Sorp_theta_pos$Sm_C2 ~ Sorp_theta_pos$Date,  xlab = theXlab222, 
     log="y", xaxt="n", ann=FALSE, pch=21, lwd=1, 
     ylab = theYlab222,
     ylim=c(1,5000), xaxp=c(Date1,Date2,21), yaxs='i', xaxs='i', 
     xlim=theXDates, cex.lab=1.5, cex.axis=1.3, )
axis(1, at = time(yearDate), labels = FALSE)
#axix(2, labels = TRUE)
tt <- time(yearDate)
ix <- seq(1, length(tt)) #every month
fmt <- "%b-%Y" # format for axis labels
labs <- format(tt[ix], fmt)
axis(side = 1, at = tt[ix], labels = labs,  cex.axis = 1.2, format= "%m/%Y", las = 2)  
mtext(theYlab222, side=2, line = 3, cex=1.5)
#   times <- time(yearDate)
#   ticks <- seq(times[1], times[length(times)], by = "months")
#   axis(1, at = ticks, labels = FALSE, tcl = -0.3)
polygon(PolyxK, PolyyK, col=Polycolor, border = "black", lty="dotted")
points(Sorp_theta_pos$Sm_C2 ~ Sorp_theta_pos$Date, pch=21, col="black" )
points(C2_Date$mean ~ C2_Date$Date, pch=21, col="black", bg="lightsalmon", cex=2,)
text(Date3,2000,"b)",cex=1.5,font=2)


### Subset to exclude cracked summertime values

WinterDate <- c('9/23/2011','10/20/2011', '11/28/2011', '2/6/2012', '3/8/2012', '4/9/2012', '5/8/2012', '10/4/2012', '11/6/2012', '12/10/2012', '1/6/2013', '2/4/2013', '3/4/2013') 

SummerDate <- c('6/3/2012', '7/5/2012', '8/7/2012', '9/6/2012', '10/3/2012')

SameLocation <- c('9/23/2011','10/20/2011', '11/28/2011', '2/6/2012', '3/8/2012', '4/9/2012', '5/8/2012')


Sorp_sub <- subset(Sorp_theta, Date %in% WinterDate)

Sorp_sameloc <- subset(Sorp_theta, Date %in% SameLocation)

Sorp_sameloc_pos <- Sorp_sameloc[Sorp_sameloc$Sm_C1 > 0,]



Size_sub <- length(Sorp_sub$Ring)

Sorp_sub$Date <- as.Date(Sorp_sub$Date, "%m/%d/%Y")
                     
Sorp_positive <- Sorp_sub[Sorp_sub$Sm_C1 > 0,]

Lew_Sm_posistats = ddply(Sorp_positive, c('degreesat_avg'), function(echs) c(count=nrow(echs), mean=mean(echs$Sm_C1), geomean=geometric.mean(echs$Sm_C1),median=median(echs$Sm_C1), sd=sd(echs$Sm_C1)))

Lew_C2_posistats = ddply(Sorp_positive, c('degreesat_avg'), function(echs) c(count=nrow(echs), mean=mean(echs$Sm_C2), geomean=geometric.mean(echs$Sm_C2),median=median(echs$Sm_C2), sd=sd(echs$Sm_C2)))

S_prime_Sm <- Lew_Sm_posistats$mean/((1-a*Lew_Sm_posistats$degreesat_avg)^0.5)



Sfun <- lm(S_prime_Sm ~ Lew_Sm_posistats$degreesat_avg)


S0_Sm <- mean(S_prime_Sm)


# Subsets on randomly placed rings (after 10/3/2012)

Sorp_subset <- subset(Sorp_theta, Sorp_theta[[2]] %in% c("10/4/2012","11/6/2012","12/10/2012","1/6/2013","2/4/2013","3/4/2013"), drop = TRUE)

Sorp_sub_Sm <- Sorp_subset[Sorp_subset$Sm_C1 > 0,]

Lew_Sm_stats = ddply(Sorp_sub_Sm, c('degreesat_avg'), function(echs) c(count=nrow(echs), mean=mean(echs$Sm_C1), geomean=geometric.mean(echs$Sm_C1),median=median(echs$Sm_C1), sd=sd(echs$Sm_C1)))

Lew_C2_stats = ddply(Sorp_sub_Sm, c('degreesat_avg'), function(echs) c(count=nrow(echs), mean=mean(echs$Sm_C2), geomean=geometric.mean(echs$Sm_C2),median=median(echs$Sm_C2), sd=sd(echs$Sm_C2)))

S_prime_Sm2 <- Lew_Sm_stats$mean/((1-Lew_Sm_stats$degreesat_avg)^0.5)

A0_Sm <- Lew_C2_stats$mean/((1-Lew_C2_stats$degreesat_avg)^2)


##### Input data without summer dates

z_filter <- z

z_filter <- subset(z, Date %in% WinterDate)


time_filter = z_filter[,10:20]
sqrt_time_filter = time_filter^0.5
x_filter = sqrt_time_filter
is.na(x_filter)
sqrt_time_filter = ifelse(is.na(x_filter), 0, x_filter)
time_filter = sqrt_time_filter^2

Size_filter <- length(z_filter[,1])

Iactual_filter = c(seq(0,10,1))
Iactual_filter = Iactual_filter*100/(pi*(rdisc^2))
Iactual_mat_filter = matrix(data = Iactual_filter, nrow = Size_filter, ncol = 11, byrow = TRUE, dimnames = NULL)
Iactual_mat_filter = round(Iactual_mat_filter, 2)

degreesat_avg_filter <- z_filter[,41]

ring_filter <- z_filter[,1]

moistfilter <- z_filter[,37]

ufilter <- z_filter[,35]

ufilter = ufilter/umax

degsatfilter <- z_filter[,41]

datefilter <- z_filter[,2]

ringfilter <- z_filter[,1]

uprime = ufilter*(((1-(phimax-phimin)*(((p+1)*ufilter^q)/(1+p*ufilter^q))-phimin)/(1-phimax)))


# Calculate coefficients using [R]

### Smiles (I/t^0.5 v. t^0.5)

Smiles_vec_filter = rep(1,Size_filter)

I_filter = Iactual_mat_filter #cm
t05_filter = x_filter #sec^0.5

Sm_filter = I_filter/t05_filter
is.na(Sm_filter)
Sm_filter = ifelse(is.na(Sm_filter), 0, Sm_filter)


for (i in 1:Size_filter) {
  
  Smiles_vec_filter[i] = lm(Sm_filter[i,] ~ t05_filter[i,]) 
  
}

Smiles_filter <- data.frame(matrix(unlist(Smiles_vec_filter), nrow=Size_filter, byrow=T))

#Vander_mat = matrix(data = Vander_vec, nrow = Size, ncol = 1, byrow = TRUE, dimnames = NULL)

Sm_int_filter = Smiles_filter[,1]

Sm_C1_filter = Sm_int_filter * 60    #cm h^-0.5

Sm_S0_filter = Sm_C1_filter / ((1-ufilter)^0.5) #cm h^-0.5
 
Sm_slope_filter = Smiles_filter[,2]

Sm_C2_filter = Sm_slope_filter * 3600   #cm h^-1

Ks_filter <- Sm_C2_filter/oneDcoef   #cm h^-1

Ks_gravmod <- lm(Ks_filter ~ ufilter)

Ksgraveq <- substitute(#italic(y) == b %.% italic(x)*","
  ~~italic(r)^2~"="~r2, 
  list(#b = format(coef(Kmodel)[1], digits = 2), 
    r2 = format(summary(Ks_gravmod)$r.squared, digits = 2)))

as.character(as.expression(Ksgraveq));

### 3D Infiltration using Wu et al. 1999

a_wu = 0.91 # a from Wu et al., 1999

b_wu = 0.17 # Beta from Haverkamp et al., 1994

d_wu = 1 # cm (ring insertion depth)

Gstar = d_wu + rdisk/2

lambda_loam = 25 # cm^-1

f_loam = lambda_loam/Gstar

K3DWu <- Sm_C2_filter/(a_wu*f_loam) #cm/hr

K3DWu[K3DWu<0] <- NA

KWuframe <- data.frame(K3DWu,ufilter)

K3DWu_mean <- ddply(KWuframe, c('ufilter'), function(echs) c(count=nrow(echs),
         mean=mean(echs$K3DWu),  geomean=geometric.mean(echs$K3DWu), median=median(echs$K3DWu), sd=sd(echs$K3DWu)))

## To determine K values
Umod = K3DWu_mean$ufilter
phicrackmod = (phipedon-phimin)*((1-Umod^q)/(1+p*Umod^q))
phisubmod = (phimax-phipedon)*((1-Umod^q)/(1+p*Umod^q))
phiaggrmod = (phimax-phimin)*((Umod)^q*(p+1)/(1+p*(Umod)^q))+phimin

Kfsobs = log10(K3DWu_mean$geomean)
# # parameter fitting using levenberg marquart algorithm
# # initial guess for parameters
# ssq=function(parms){
#   
#   # parameters from the parameter estimation routine
#   k1=parms[1]
#   k2=parms[2]
#   b1=parms[3]
#   b2=parms[4]
#   
#   outdf = k1*(phicrackmod/(1-phisubmod))*((1-Umod^q)/(1+p*Umod^q))^b1+
#     (1-(phicrackmod/(1-phisubmod)))*((p+1)/(p+Umod^-q))^b2*k2
#   obsdf = K3DWu_mean$geomean
#   # Evaluate predicted vs experimental residual
#   preddf=melt(outdf,id.var="Umod",variable.name="K",value.name="Kfs")
#   expdf=melt(obsdf,id.var="uratio",variable.name="K",value.name="Kfs")
#   ssqres=preddf$Kfs-expdf$Kfs
#   
#   # return predicted vs experimental residual
#   return(ssqres)
#   
# }
# parms=c(k1=5000,k2=0.01,b1=1,b2=1)
# # fitting
# fitval=nls.lm(par=parms,fn=ssq,lower=c(0,0,1,1),upper=c(Inf,Inf,1,1))
# 
# parest=as.list(coef(fitval))

kfunc <- function(Umod,k1,k2){     
  b1 = 2
  b2 = 1
  phicrackmod = (phipedon-phimin)*((1-Umod^q)/(1+p*Umod^q))
  phisubmod = (phimax-phipedon)*((1-Umod^q)/(1+p*Umod^q))
  phiaggrmod = (phimax-phimin)*((Umod)^q*(p+1)/(1+p*(Umod)^q))+phimin
  Ksmod = log10(k1*(phicrackmod/(1-phisubmod))*((1-Umod^q)/(1+p*Umod^q))^b1+
    (1-(phicrackmod/(1-phisubmod)))*((p+1)/(p+Umod^-q))^b2*k2)
      # ((p+1)/((p+Umod^-q)*(1-phisubmod)))^b2*k2)
}

fitval=nls(Kfsobs~kfunc(Umod,k1,k2),start=list(k1=11,k2=0.3),algorithm="port",lower=c(0,0),upper=c(Inf,Inf))
parest=as.list(coef(fitval))

# Calculate 95% Confidence Intervals

interval <- confint(fitval, level = 0.95)

### Ks model with optimized parameters!
Kscrackpar = as.numeric(parest[1]) # cm/hr
Ksaggrpar = as.numeric(parest[2]) # cm/hr
Bcrack = 2
Baggr = 1

# Kscrackpar = Kscrackmax3D
# Ksaggrpar = Ksaggrmax3D
# Bcrack = 1
# Baggr = 1

Ksfunc3D1 = Kscrackpar*(phi_crack/(1-phi_sub))*((1-Ufunc^q)/(1+p*Ufunc^q))^Bcrack
Ksfunc3D2 = Ksaggrpar*(1-(phi_crack/(1-phi_sub)))*((p+1)/(p+Ufunc^-q))^Baggr
Ksfunc3D = Ksfunc3D1+Ksfunc3D2



# To calculate r^2 value:

Ksmod = Kscrackpar*(phicrackmod/(1-phisubmod))*((1-Umod^q)/(1+p*Umod^q))^Bcrack+
  (1-(phicrackmod/(1-phisubmod)))*((p+1)/(p+Umod^-q))^Baggr*Ksaggrpar

Kmodel <- lm(log10(K3DWu_mean$geomean) ~ log10(Ksmod))
Keq <- substitute(#italic(y) == b %.% italic(x)*","
  ~~italic(r)^2~"="~r2, 
  list(#b = format(coef(Kmodel)[1], digits = 2), 
    r2 = format(summary(Kmodel)$r.squared, digits = 2)))


### K from sorptivity

Ks_Sm_filter <- Sm_C1_filter^2 / (0.8*hf*(1-degreesat_avg_filter))   #cm h^-1

Kssm <- data.frame(Ks_Sm_filter, degreesat_avg_filter)

Kssm_mean <- ddply(Kssm, c('degreesat_avg_filter'), function(echs) c(count=nrow(echs), mean=mean(echs$Ks_Sm_filter), geomean=geometric.mean(echs$Ks_Sm_filter),median=median(echs$Ks_Sm_filter), sd=sd(echs$Ks_Sm_filter)))

Kssmfun <- lm(Kssm_mean$mean~Kssm_mean$degreesat_avg_filter)

Kssmeq <- substitute(#italic(y) == b %.% italic(x)*","
  ~~italic(r)^2~"="~r2, 
  list(#b = format(coef(Kmodel)[1], digits = 2), 
    r2 = format(summary(Kssmfun)$r.squared, digits = 2)))

as.character(as.expression(Kssmeq));

Kssmfunlog <- lm(log10(Kssm_mean$mean) ~ (Kssm_mean$degreesat_avg_filter))

S0Sm <- data.frame(Sm_S0_filter,degreesat_avg_filter)

SSm_mean <- ddply(S0Sm, c('degreesat_avg_filter'), function(echs) c(count=nrow(echs), mean=mean(echs$Sm_S0_filter), median=median(echs$Sm_S0_filter), sd=sd(echs$Sm_S0_filter)))

SSmfun <- lm(SSm_mean$mean~SSm_mean$degreesat_avg_filter)


# 
# #### Short-time Approx
# 
# White_vec_filter = rep(1,Size_filter)
# 
# 
# for (i in 1:Size_filter) {
#   
#   White_vec_filter[i] = lm(Iactual_mat_filter[i,1:5] ~ 0 + x_filter[i,1:5]) 
#   
# }
# 
# White_filter <- data.frame(matrix(unlist(White_vec_filter), nrow=Size_filter, byrow=T))
# 
# W_C1_filter = (White_filter[,1])
# 
# W_C1_filter = W_C1_filter * 60   #cm h^-0.5
# 
# W_S0_filter = W_C1_filter / (1-degreesat_avg_filter)   #cm h^-0.5
# 
# Ks_W_filter <- W_C1_filter^2 / (0.8*hf*(1-degreesat_avg_filter))   #cm h^-1
#  
# KsW <- data.frame(Ks_W_filter, degreesat_avg_filter)
# 
# Ksw_mean <- ddply(KsW, c('degreesat_avg_filter'), function(echs) c(count=nrow(echs), mean=mean(echs$Ks_W_filter), geomean=geometric.mean(echs$Ks_W_filter),median=median(echs$Ks_W_filter), sd=sd(echs$Ks_W_filter)))
# 
# KsWfun <- lm(Ksw_mean$mean~Ksw_mean$degreesat_avg_filter)
# 
# KsWfunlog <- lm(log10(Ksw_mean$mean) ~ (Ksw_mean$degreesat_avg_filter))
# 
# S0W <- data.frame(W_S0_filter,degreesat_avg_filter)
# 
# SW_mean <- ddply(S0W, c('degreesat_avg_filter'), function(echs) c(count=nrow(echs), mean=mean(echs$W_S0_filter),geomean=geometric.mean(echs$W_S0_filter), median=median(echs$W_S0_filter), sd=sd(echs$W_S0_filter)))
# 
# SWfun <- lm(SW_mean$mean~SW_mean$degreesat_avg_filter)
# 
# ## Long Time Approx
# 
# ### I = C + Kt
# 
# long_vec = rep(1,Size_filter)
# 
# longtime = as.matrix(z_filter[,10:20],)
# longtime = subset(longtime,!is.na(longtime[,6]))
# 
# ulongframe = data.frame(ufilter,z_filter[,15])
# ulongframe = subset(ulongframe,!is.na(ulongframe[,2]))
# ulong <- ulongframe[,1]
# 
# Sizelong = length(longtime[,1])
# 
# ILT = Iactual_mat_filter #cm
# ILT = ILT[1:Sizelong,]
# is.na(ILT)
# ILT = ifelse(is.na(ILT), 0, ILT)
# 
# 
# for (i in 1:Sizelong) {
#   
#   long_vec[i] = lm(ILT[i,5:11]~longtime[i,5:11])
#   
# }
# 
# long <- data.frame(matrix(unlist(long_vec), nrow=Sizelong, byrow=T))
# 
# LT_int = long[,1]   #cm
# 
# LT_slope = long[,2]  #cm s^-1
# 
# Ks_LT <- LT_slope*3600   #cm h^-1
# 
# Ks_LT = ifelse(is.na(Ks_LT), 0, Ks_LT)
# 
# Ks_LT <- data.frame(ulong,Ks_LT)
# 
# Ks_LT <- Ks_LT[Ks_LT$Ks_LT > 0,]
# 
# Ks_LT_stats = ddply(Ks_LT, c('ulong'), function(echs) c(count=nrow(echs), 
#                 mean=mean(echs$Ks_LT), geomean=geometric.mean(echs$Ks_LT), median=median(echs$Ks_LT), sd=sd(echs$Ks_LT)))

## Wintertime data

z_winter <- data.frame(ring_filter, degreesat_avg_filter, Sm_C1_filter, Sm_S0_filter, Sm_C2_filter, Ks_filter)

z_winter_pos <- z_winter[z_winter$Sm_C1_filter > 0,]

C2_mean <- ddply(z_winter_pos, c('degreesat_avg_filter'), function(echs) c(count=nrow(echs), mean=mean(echs$Sm_C2_filter), geomean=geometric.mean(echs$Sm_C2_filter),median=median(echs$Sm_C2_filter), sd=sd(echs$Sm_C2_filter)))

C2fun <- lm(C2_mean$mean ~ C2_mean$degreesat_avg_filter)

KS_mean <- ddply(z_winter_pos, c('degreesat_avg_filter'), function(echs) c(count=nrow(echs), mean=mean(echs$Ks_filter),geomean=geometric.mean(echs$Ks_filter), median=median(echs$Ks_filter), sd=sd(echs$Ks_filter)))

KSfun <- lm(KS_mean$mean ~ KS_mean$degreesat_avg_filter)

KSfunlog <- lm(log10(KS_mean$mean) ~ (KS_mean$degreesat_avg_filter))

### Soil "constants"

C2_Sm_filter = mean(Sm_C2_filter)  #cm/hr

S0_Sm_filter = mean(Sm_S0_filter) #cm/hr^0.5

A0_5 <- (1/alph)*(Sm_C2_filter-alph*Ks1D)/(((1-ufilter^q)/(1+p*(ufilter^q))^2)) #cm/hr

A5 <- data.frame(A0_5,ufilter,datefilter)  #cm/hr

A5$A0_5[A5$A0_5 < 0] <- 0

A0_5_mean <- mean(A5$A0_5)  #cm/hr

A5_group <- ddply(A5, ~ufilter,mean=mean(A0_5), sd=sd(A0_5)) #cm/hr

A5_mean <- ddply(A5, ~ufilter,summarise,mean=mean(A0_5), sd=sd(A0_5)) #cm/hr

A5_mean_calc <- mean(A5$A0_5)

A5_mean_calc2 <- mean(A5_mean$mean)

A5_sd <- mean(A5_mean$sd) # Note, this is the standard deviation of the means (possibly the better indicator, due to lognormal dist of individual points)

A5_sd_calc <- sd(A5$A0_5)#Note, this is the standard deviation of all measurements



K_gravity = data.frame(ufilter,Ks_filter) #cm/hr

K_grav_mean <- ddply(K_gravity, c('ufilter'), function(echs) c(count=nrow(echs),
                 mean=mean(echs$Ks_filter), geomean=geometric.mean(echs$Ks_filter),median=median(echs$Ks_filter), sd=sd(echs$Ks_filter)))



##### Input summer data

z_summer <- subset(z, Date %in% SummerDate)

# z_summer <- z

time_summer = z_summer[,10:20]
sqrt_time_summer = time_summer^0.5
x_summer = sqrt_time_summer
is.na(x_summer)
sqrt_time_summer = ifelse(is.na(x_summer), 0, x_summer)
time_summer = sqrt_time_summer^2

Size_summer <- length(z_summer[,1])

Iactual_summer = c(seq(0,10,1))
Iactual_summer = Iactual_summer*100/(pi*(rdisc^2))
Iactual_mat_summer = matrix(data = Iactual_summer, nrow = Size_summer, ncol = 11, byrow = TRUE, dimnames = NULL)
Iactual_mat_summer = round(Iactual_mat_summer, 2)

degreesat_avg_summer <- z_summer[,41]

ring_summer <- z_summer[,1]

moistsummer <- z_summer[,37]

usummer <- z_summer[,35]

usummer = usummer/umax

usummerprime = usummer*(((1-(phimax-phimin)*(((p+1)*usummer^q)/(1+p*usummer^q))-phimin)/(1-phimax)))


degsatsummer <- z_summer[,41]

datesummer <- z_summer[,2]

ringsummer <- z_summer[,1]



# Calculate coefficients using [R]

### Smiles (I/t^0.5 v. t^0.5)

Smiles_vec_summer = rep(1,Size_summer)

I_summer = Iactual_mat_summer #cm
t05_summer = x_summer #sec^0.5

Sm_summer = I_summer/t05_summer
is.na(Sm_summer)
Sm_summer = ifelse(is.na(Sm_summer), 0, Sm_summer)


for (i in 1:Size_summer) {
  
  Smiles_vec_summer[i] = lm(Sm_summer[i,] ~ t05_summer[i,]) 
  
}

Smiles_summer <- data.frame(matrix(unlist(Smiles_vec_summer), nrow=Size_summer, byrow=T))

#Vander_mat = matrix(data = Vander_vec, nrow = Size, ncol = 1, byrow = TRUE, dimnames = NULL)

Sm_int_summer = Smiles_summer[,1]

Sm_C1_summer = Sm_int_summer * 60    #cm h^-0.5

Sm_S0_summer = Sm_C1_summer / ((1-usummer)^0.5) #cm h^-0.5

Sm_slope_summer = Smiles_summer[,2]

Sm_C2_summer = Sm_slope_summer * 3600   #cm h^-1

Ks_summer <- Sm_C2_summer/alph   #cm h^-1

Ks_gravmod_summer <- lm(Ks_summer ~ usummer)

Ksgraveq_summer <- substitute(#italic(y) == b %.% italic(x)*","
  ~~italic(r)^2~"="~r2, 
  list(#b = format(coef(Kmodel)[1], digits = 2), 
    r2 = format(summary(Ks_gravmod_summer)$r.squared, digits = 2)))

as.character(as.expression(Ksgraveq_summer));

Ks_Sm_summer <- Sm_C1_summer^2 / (0.8*hf*(1-degreesat_avg_summer))   #cm h^-1

Kssm_summer <- data.frame(Ks_Sm_summer, degreesat_avg_summer)

Kssm_mean_summer <- ddply(Kssm_summer, c('degreesat_avg_summer'), function(echs) c(count=nrow(echs), mean=mean(echs$Ks_Sm_summer), geomean=geometric.mean(echs$Ks_Sm_summer),median=median(echs$Ks_Sm_summer), sd=sd(echs$Ks_Sm_summer)))

Kssmfun_summer <- lm(Kssm_mean_summer$mean~Kssm_mean_summer$degreesat_avg_summer)

Kssmeq_summer <- substitute(#italic(y) == b %.% italic(x)*","
  ~~italic(r)^2~"="~r2, 
  list(#b = format(coef(Kmodel)[1], digits = 2), 
    r2 = format(summary(Kssmfun_summer)$r.squared, digits = 2)))

as.character(as.expression(Kssmeq_summer));

Kssmfunlog_summer <- lm(log10(Kssm_mean_summer$mean) ~ (Kssm_mean_summer$degreesat_avg_summer))

S0Sm_summer <- data.frame(Sm_S0_summer,degreesat_avg_summer)

SSm_mean_summer <- ddply(S0Sm_summer, c('degreesat_avg_summer'), function(echs) c(count=nrow(echs), mean=mean(echs$Sm_S0_summer), geomean=geometric.mean(echs$Sm_S0_summer),median=median(echs$Sm_S0_summer), sd=sd(echs$Sm_S0_summer)))

SSmfun_summer <- lm(SSm_mean_summer$mean~SSm_mean_summer$degreesat_avg_summer)

z_summer2 <- data.frame(ring_summer, degreesat_avg_summer, Sm_C1_summer, Sm_S0_summer, Sm_C2_summer, Ks_summer)

z_summer_pos <- z_summer2[z_summer2$Sm_C1_summer > 0,]

C2_mean_summer <- ddply(z_summer_pos, c('degreesat_avg_summer'), function(echs) c(count=nrow(echs), mean=mean(echs$Sm_C2_summer), geomean=geometric.mean(echs$Sm_C2_summer),median=median(echs$Sm_C2_summer), sd=sd(echs$Sm_C2_summer)))

C2fun_summer <- lm(C2_mean_summer$mean ~ C2_mean_summer$degreesat_avg_summer)

KS_mean_summer <- ddply(z_summer_pos, c('degreesat_avg_summer'), function(echs) c(count=nrow(echs), mean=mean(echs$Ks_summer), geomean=geometric.mean(echs$Ks_summer),median=median(echs$Ks_summer), sd=sd(echs$Ks_summer)))

KSfun_summer <- lm(KS_mean_summer$mean ~ KS_mean_summer$degreesat_avg_summer)

KSfunlog_summer <- lm(log10(KS_mean_summer$mean) ~ (KS_mean_summer$degreesat_avg_summer))

### Soil "constants"

C2_Sm_summer = mean(Sm_C2_summer)  #cm/hr

S0_Sm_summer = mean(Sm_S0_summer) #cm/hr^0.5


A0_5_summer <- (1/alph)*(Sm_C2_summer-alph*Ks1D)/(((1-usummer^q)/(1+p*(usummer^q))^2)) #cm/hr

A5_summer <- data.frame(A0_5_summer,usummer)  #cm/hr

A5_summer$A0_5_summer[A5_summer$A0_5_summer < 0] <- 0

#A5 <- A5[A5$degsatsummer < 0.96 ,]

A0_5_mean_summer <- mean(A5_summer$A0_5_summer)  #cm/hr

A5_group_summer <- ddply(A5_summer, ~usummer,mean=mean(A0_5_summer), sd=sd(A0_5_summer)) #cm/hr

A5_mean_summer <- ddply(A5_summer, ~usummer,summarise,mean=mean(A0_5_summer), sd=sd(A0_5_summer)) #cm/hr

A5_mean_calc_summer <- mean(A5_summer$A0_5_summer)

A5_mean_calc2_summer <- mean(A5_mean_summer$mean)

A5_sd_summer <- mean(A5_mean_summer$sd) # Note, this is the standard deviation of the means (possibly the better indicator, due to lognormal dist of individual points)

A5_sd_calc_summer <- sd(A5_summer$A0_5_summer)#Note, this is the standard deviation of all measurements



K_gravity_summer = data.frame(usummer,Ks_summer) #cm/hr

K_grav_mean_summer <- ddply(K_gravity_summer, c('usummer'), function(echs) c(count=nrow(echs),
              mean=mean(echs$Ks_summer), geomean=geometric.mean(echs$Ks_summer),median=median(echs$Ks_summer), sd=sd(echs$Ks_summer)))




A5fun <- lm(A5_mean$mean~A5_mean$ufilter)


A5cor = cor(A0_5, z_filter[,41])
A5cor = signif(A5cor, digits = 4)
print(A5cor)
A5int = A5fun$coefficients[[1]]
A5int = signif(A5int, digits = 3)
A5slope = A5fun$coefficients[[2]]
A5slope = signif(A5slope, digits = 3)
A5R = format(summary(A5fun)$r.squared, digits = 3)


A5eq <- substitute(~~italic(r)^2~"="~r2, 
                   list(r2 = A5R))


### Plot Ks (this is 1D Philip)
# 
# 
win.graph(width = 26, height = 24)
par(mfrow = c(1, 1))
par(oma=c(0,0,1,0)) 
par(mar=c(4,6,0.5,1))

theYticks=c(0,100,10)
ylimsI=c(0,100)
theXticks=c(0,15.2,1)
xlimsI=c(0,16)

theXlab400 = expression(paste(bold("Water Content,"~u/u[max])))
theYlab400 = expression(paste(bold(K[s]~~bgroup("(",cm~~hr^{-1},")"))))

plot(K_gravity$Ks_filter ~ K_gravity$ufilter, xlab = theXlab400, ylab = theYlab400, type="p", log="y",
     pch=21, cex = 1.3, col="gray40", bg="lightgray",  ylim=c(1,4000), xaxp=c(0,1,10), yaxs='i', xaxs='i', xlim = c(0,1), cex.lab=1.8, cex.axis=1.5, family="Arial")
# abline(h=A0_5_mean,lty=2,lwd=2)
# abline(lm(log10(K_grav_mean$mean) ~ K_grav_mean$ufilter),lty=2,lwd=2)
lines(Ksfunc~Ufunc,lwd=2)
lines(Ksfunc1~Ufunc,lwd=2,lty=2)
lines(Ksfunc2~Ufunc,lwd=2,lty=3)
# points(K_gravity$Ks_filter ~ K_gravity$ufilter,pch=21, cex = 1.3, col="gray40", bg="lightgray",)
# points(K_gravity_summer$Ks_summer ~ K_gravity_summer$usummer,pch=21, cex = 1.3, col="black", bg="white",)
# points(K_grav_mean_summer$geomean ~ K_grav_mean_summer$usummer,pch=21, cex = 2, col="black", bg="darkslategray",)
points(K_grav_mean$geomean ~ K_grav_mean$ufilter,pch=21, cex = 2, col="black", bg="lightsalmon",)
# points(Ks_LT_stats$geomean ~ Ks_LT_stats$ulong,pch=21, cex = 2, col="black", bg="lightsalmon",)
legend(0.11, 30, c("Individual", "Geom. Mean"), pch=c(21,21),pt.bg=c("lightgray","lightsalmon"),col = c("gray40","black"),pt.cex=c(1.6,2.2),
       text.col = c("black","black"),
       #        title="Depth", title.col="black",
       cex=c(1.8), bg = "white")


### Plot A0


win.graph(width = 26, height = 24)
par(mfrow = c(1, 1))
par(oma=c(0,0,1,0)) 
par(mar=c(4,6,0,1))

theYticks=c(0,100,10)
ylimsI=c(0,100)
theXticks=c(0,15.2,1)
xlimsI=c(0,16)

theXlab4 = expression(paste(bold("Water Content,"~u/u[max])))
theYlab4 = expression(paste(bold(A[0]~"from 1D solution"~~bgroup("(",cm~~hr^{-1},")"))))

plot(A5$A0_5 ~ A5$ufilter, xlab = theXlab4, ylab = theYlab4, type="p", log="y",
     pch=21, cex = 1.3, col="black", bg="lightgray",  ylim=c(10,20000), xaxp=c(0,1,10), 
     yaxs='i', xaxs='i', xlim = c(0,1), cex.lab=1.8, cex.axis=1.5,)
# abline(h=A0_5_mean,lty=2,lwd=2)
# abline(lm(log10(A5_mean$mean) ~ A5_mean$ufilter),lty=2,lwd=2)
# abline(lm(log10(A5$A0_5) ~ A5$ufilter),lty=1.6)
points(A5_summer$A0_5_summer ~ A5_summer$usummer, pch=21, cex = 1.3, col = "black", bg = "white")
points(A5_mean_summer$mean ~ A5_mean_summer$usummer, cex=2, pch=21, col = "black", bg = "darkslategray")
points(A5$A0_5 ~ A5$ufilter,pch=21, cex = 1.3, col="black", bg="lightgray",)
points(A5_mean$mean ~ A5_mean$ufilter, cex=2, pch=21, col="black", bg="lightsalmon",)
legend(0.11, 100, c("Individual", "Arith. Mean"), pch=c(21,21),pt.bg=c("gray","lightsalmon"),col = c("black","black"),pt.cex=c(1.6,2.2),
       text.col = c("black","black"),
       #        title="Depth", title.col="black",
       cex=c(1.8), bg = "white")
text(0.035,15200,"a)",cex=1.6)


## Plot K3D as function of water content (Wu and Pan)


win.graph(width = 26, height = 24)
par(mfrow = c(1, 1))
par(oma=c(0,0,1,0)) 
par(mar=c(4,6,0.5,1))

theYticks=c(0,100,10)
ylimsI=c(0,100)
theXticks=c(0,15.2,1)
xlimsI=c(0,16)

theXlab400 = expression(paste(bold("Water Content,"~u/u[max])))
theYlab401 = expression(paste(bold(K[s]~~bgroup("(",mm~~hr^{-1},")"))))

plot(10*K3DWu ~ ufilter, xlab = theXlab400, ylab = theYlab401, type="p", log="y",
     pch=22, cex = 0.8, col="black", bg="lightgray",  ylim=c(1,2000), xaxp=c(0,1,10), 
     yaxs='i', xaxs='i', xlim = c(0,1), cex.lab=1.8, cex.axis=1.5, )
# abline(h=A0_5_mean,lty=2,lwd=2)
# abline(lm(log10(K_grav_mean$mean) ~ K_grav_mean$ufilter),lty=2,lwd=2)
lines(10*Ksfunc3D~Ufunc,lwd=2)
lines(10*Ksfunc3D1~Ufunc,lwd=2,lty=2)
abline(h=170,lwd=2,lty=6)
lines(10*Ksfunc3D2~Ufunc,lwd=2,lty=3)
points(10*K3DWu ~ ufilter, pch=22,cex=0.8,col="black",bg="lightgray")
points(10*K3DWu_mean$geomean ~ K3DWu_mean$ufilter,pch=22, cex = 1.8, col="black", bg="lightsalmon",)
text(0.17,1450,Keq,cex=1.7)
text(0.03,1500,c("b)"),cex=1.7)
legend(0.11, 10, c("Individual", "Geo. Mean"), pch=c(22,22),pt.bg=c("gray","lightsalmon"),col = c("black","black"),pt.cex=c(1.6,2),
       text.col = c("black","black"),
       #        title="Depth", title.col="black",
       cex=c(1.8), bg = "white")

### Compare 1D and 3D K values

# 
# win.graph(width = 26, height = 24)
# par(mfrow = c(1, 1))
# par(oma=c(0,0,1,0))
# par(mgp=c(4,1,0))
# par(mar=c(6,7,2,1.5))
# 
# theYlabKcomp = expression(paste(bold(K[s]~~"From 1D solution ("~cm~~hr^-1~")"))) 
# theXlabKcomp = expression(paste(bold(K[s]~~"From 3D solution ("~cm~~hr^-1~")"))) 
# 
# plot(Kframe$Ks_filter ~ Kframe$K3D, xlab = theXlabKcomp, ylab = theYlabKcomp, type="p", log="xy",
#      pch=21, cex = 1.3, col="black",  bg="red", yaxs='i', xaxs='i', cex.lab=1.8, cex.axis=1.5, 
#      family="Arial", xlim=c(1,5000), ylim=c(1,5000))
# abline(Klogmodel)
# lines(c(1,5000),c(1,5000), lty="dashed")
# text(8, 1500, Keq, cex = 1.8)
# 
# 
# 
# win.graph(width = 26, height = 24)
# par(mfrow = c(1, 1))
# par(oma=c(0,0,1,0))
# par(mgp=c(4,1,0))
# par(mar=c(6,7,2,1.5))
# 
# theYlabKcomp = expression(paste(bold(K[s]~~"From 1D solution ("~cm~~hr^-1~")"))) 
# theXlabKcomp = expression(paste(bold(K[s]~~"From 3D solution ("~cm~~hr^-1~")"))) 
# 
# plot(Kwetframe$Ks_filter ~ Kwetframe$K3D, xlab = theXlabKcomp, ylab = theYlabKcomp, type="p", log="xy",
#      pch=21, cex = 1.3, col="black",  bg="red", yaxs='i', xaxs='i', cex.lab=1.8, cex.axis=1.5, 
#      family="Arial", xlim=c(1,5000), ylim=c(1,5000))
# abline(Kwetlogmodel)
# lines(c(1,5000),c(1,5000), lty="dashed")
# text(8, 1500, Kweteq, cex = 1.8)
# 

## Add A0 term into K3D
# 
# A3DU <- (((Sm_C2_filter-gammy*Sm_C1_filter^2/(rdisc*phimax*(1-uprime)))/
#             ((((2-Bates)/3)*(1-(uprime)^BCeta))+(uprime)^BCeta))-Ks3D)*(((1+p*ufilter^q)/(1-ufilter^q))^2) #cm/hr
# 
# 
# A3Dframe <- data.frame(A3DU, ufilter)
# 
# A3Dframe0 <- A3Dframe
# 
# A3Dframe0$A3DU[A3Dframe0$A3DU < 0] <- 0
# 
# A3Dframe <- A3Dframe[A3Dframe$A3DU > 0,]
# 
# A3Dmean <- mean(A3Dframe0$A3DU) #cm/hr
# 
# A3D_mean <- ddply(A3Dframe, ~ufilter,summarise,mean=mean(A3DU), sd=sd(A3DU)) #cm/hr
# 
# A3D_mean_calc2 <- mean(A3D_mean$mean)
# 
# A3Dn <- length(A3D_mean$mean)
# 
# A3D_sd <- mean(A3D_mean$sd) # Note, this is the standard deviation of the means (possibly the better indicator, due to lognormal dist of individual points)
# 
# A3D_sd_calc <- sd(A3Dframe$A3DU)#Note, this is the standard deviation of all measurements
# 
# A_3D_mean <- ddply(A3Dframe, c('ufilter'), function(echs) c(count=nrow(echs),
#                                                              mean=mean(echs$A3DU), median=median(echs$A3DU), sd=sd(echs$A3DU)))
# 
# A3Derror <- qnorm(0.975)*A3D_sd_calc/sqrt(A3Dn)
# 
# A3Dupper <- A3Dmean+A3Derror #Note, these points are log-normally distributed.. confidence intervals may not work
# 
# A3Dlower <- A3Dmean-A3Derror
# 
# ## Summer
# 
# 
# A3DUsummer <- (((Sm_C2_summer-gammy*Sm_C1_summer^2/(rdisc*phimax*(1-usummerprime)))/((((2-Bates)/3)*(1-(usummerprime)^BCeta))+(usummerprime)^BCeta))-Ks3D)*(((1+p*usummer^q)/(1-usummer^q))^2) #cm/hr
# 
# A3Dframesummer <- data.frame(A3DUsummer, usummer)
# 
# A3Dframe0summer <- A3Dframesummer
# 
# A3Dframe0summer$A3DUsummer[A3Dframe0summer$A3DUsummer < 0] <- 0
# 
# A3Dframesummer <- A3Dframesummer[A3Dframesummer$A3DUsummer > 0,]
# 
# A3Dmeansummer <- mean(A3Dframe0summer$A3DUsummer) #cm/hr
# 
# A3D_meansummer <- ddply(A3Dframesummer, ~usummer,summarise,mean=mean(A3DUsummer), sd=sd(A3DUsummer)) #cm/hr
# 
# A3D_mean_calc2summer <- mean(A3D_meansummer$mean)
# 
# A3Dnsummer <- length(A3D_meansummer$mean)
# 
# A3D_sdsummer <- mean(A3D_meansummer$sd) # Note, this is the standard deviation of the means (possibly the better indicator, due to lognormal dist of individual points)
# 
# A3D_sd_calcsummer <- sd(A3Dframesummer$A3DUsummer)#Note, this is the standard deviation of all measurements
# 
# A_3D_meansummer <- ddply(A3Dframesummer, c('usummer'), function(echs) c(count=nrow(echs),
#                                                                         mean=mean(echs$A3DUsummer), median=median(echs$A3DUsummer), sd=sd(echs$A3DUsummer)))
# 
# 
# ## Plot A03D as function of water content
# 
# 
# win.graph(width = 26, height = 24)
# par(mfrow = c(1, 1))
# par(oma=c(0,0,1,0)) 
# par(mar=c(4,6,0.5,1))
# 
# theYticks=c(0,100,10)
# ylimsI=c(0,100)
# theXticks=c(0,15.2,1)
# xlimsI=c(0,16)
# 
# theXlab900 = expression(paste(bold("Water Content,"~u/u[max])))
# theYlab901 = expression(paste(bold(A[0]~"from 3D solution"~~bgroup("(",cm~~hr^{-1},")"))))
# 
# plot(A3Dframe$A3DU ~ A3Dframe$ufilter, xlab = theXlab900, ylab = theYlab901, type="p", log="y",
#      pch=21, cex = 1.3, col="black", bg="lightgray",  ylim=c(10,40000), xaxp=c(0,1,10), yaxs='i', xaxs='i', xlim = c(0,1), cex.lab=1.8, cex.axis=1.5,)
# # abline(h=A0_5_mean,lty=2,lwd=2)
# # abline(lm(log10(K_grav_mean$mean) ~ K_grav_mean$ufilter),lty=2,lwd=2)
# points(A3Dframesummer$A3DUsummer ~ A3Dframesummer$usummer, pch = 21, cex = 1.3, col ="black", bg="white")
# points(A_3D_mean$mean ~ A_3D_mean$ufilter,pch=21, cex = 2, col="black", bg="lightsalmon",)
# points(A_3D_meansummer$mean ~ A_3D_meansummer$usummer,pch=21, cex = 2, col="black", bg="darkslategray",)
# # abline(h=A3Dupper,lty=2,lwd=2)
# # abline(h=A3Dlower,lty=2,lwd=2)
# legend(0.11, 100, c("Individual", "Arith. Mean"), pch=c(21,21),pt.bg=c("lightgray","lightsalmon"),col = c("black","black"),pt.cex=c(1.6,2.2),
#        text.col = c("black","black"),
#        #        title="Depth", title.col="black",
#        cex=c(1.8), bg = "white")
# text(0.035,30000,"b)",cex=1.6)
# 
# 
# ### Plot A0 for 1D and 3D models
# 
# 
# win.graph(width = 26, height = 24)
# par(mfrow = c(1, 1))
# par(oma=c(0,0,1,0)) 
# par(mar=c(4,6,0,1))
# 
# theXlab4 = expression(paste(bold("Water Content,"~u/u[max])))
# theYlab4 = expression(paste(bold(A[0]~"from 1D solution"~~bgroup("(",cm~~hr^{-1},")"))))
# 
# plot(A5$A0_5 ~ A5$ufilter, xlab = theXlab4, ylab = theYlab4, type="p", log="y",
#      pch=21, cex = 1.3, col="black", bg="lightgray",  ylim=c(10,40000), xaxp=c(0,1,10), yaxs='i', xaxs='i', xlim = c(0,1), cex.lab=1.8, cex.axis=1.5, family="Arial")
# # abline(h=A0_5_mean,lty=2,lwd=2)
# # abline(lm(log10(A5_mean$mean) ~ A5_mean$ufilter),lty=2,lwd=2)
# # abline(lm(log10(A5$A0_5) ~ A5$ufilter),lty=1.6)
# points(A3Dframesummer$A3DUsummer ~ A3Dframesummer$usummer, pch = 23, cex = 1.3, col ="black", bg="white")
# points(A5_summer$A0_5_summer ~ A5_summer$usummer, pch=21, cex = 1.3, col = "black", bg = "white")
# points(A5$A0_5 ~ A5$ufilter,pch=21, cex = 1.3, col="black", bg="lightgray",)
# points(A3Dframe$A3DU ~ A3Dframe$ufilter,pch=23, cex = 1.2, col="black", bg="lightgray",)
# points(A_3D_meansummer$mean ~ A_3D_meansummer$usummer,pch=23, cex = 2, col="black", bg="darkslategray",)
# points(A5_mean_summer$mean ~ A5_mean_summer$usummer, cex=2, pch=21, col = "black", bg = "darkslategray")
# points(A_3D_mean$mean ~ A_3D_mean$ufilter, cex=1.7, pch=23, col="black", bg="lightsalmon",)
# points(A5_mean$mean ~ A5_mean$ufilter, cex=2, pch=21, col="black", bg="lightsalmon",)
# legend(0.11, 100, c("Individual", "Arith. Mean"), pch=c(21,21),pt.bg=c("lightgray","lightsalmon"),col = c("black","black"),pt.cex=c(1.6,2.2),
#        text.col = c("black","black"),
#        #        title="Depth", title.col="black",
#        cex=c(1.8), bg = "white")

### Histograms

## For 1D model

A0filt = A0_5[A0_5>0]
# 
# 
# win.graph(width = 24, height = 18)
# par(mfrow = c(1, 1))
# par(oma=c(0,0,1,0)) 
# par(mar=c(6,6,2,1.5))
# 
# theYticks=c(0,100,10)
# ylimsI=c(0,100)
# theXticks=c(0,15.2,1)
# xlimsI=c(0,16)
# 
theYlab41 = expression(paste(bold("Frequency"))) 
theXlab41 = expression(paste(bold("log "~A[0])))
# 
binwidth=c(1)
bins=seq(-8,4,by=binwidth)
# hist(log(A0filt),breaks=bins, main=" ",col="lightgray",ylim=c(0,70), xlab = theXlab41, ylab = theYlab41, cex.axis=1.5,cex.lab=1.8)

#Obtain the mean and standard deviation of the data
Mu = mean(log(A0filt));
Sigma = sd(log(A0filt));

population_x <- seq(
  qnorm(0.001, Mu, Sigma), 
  qnorm(0.999, Mu, Sigma), 
  length.out = 60
)

binwidth=c(1)
bins=seq(-5,20,by=binwidth)
# hist(log(A5$A0_5),breaks=bins, main=" ",col="lightgray",ylim=c(0,35), xlab = theXlab41, ylab = theYlab41, cex.axis=1.5,cex.lab=1.8)
# box()

# # 
# theYlab49 = expression(paste(bold("Quantiles for log-normal distribution"))) 
# theXlab49 = expression(paste(bold("Quantiles for "~A[0])))


# win.graph(width = 32, height = 18)
# par(mfrow = c(1, 2))
# par(oma=c(0,0,1,0)) 
# par(mar=c(5,5,2,1.5))
# 
# hist(log(A0filt),breaks=bins, main=" ",col="lightgray",xlim=c(2,12),ylim=c(0,50),
#      xlab = theXlab41, ylab = theYlab41, cex.axis=1.5,cex.lab=1.8)
# lines(population_x, 145 * dnorm(population_x, Mu, Sigma) * binwidth, 
#       col = "black", lwd=2)
# text(2.2,48,"a)",cex=2.)
# 
# qqPlot(A0filt, "log-normal", xlab=theXlab49, xlim=c(0,15000),ylim=c(0,15000),
#        ylab=theYlab49, cex.lab=1.5,cex.axis=1.3)
# text(400,14500,"b)",cex=2.)
z_plot <- z_filter[z_filter$Date %in% c("5/8/2012"),]
# z_C1 <- z_plot[,4]/60
# z_C2 <- z_plot[,7]/3600

timeplot = z_plot[,10:20]  # sec
sqrt_timeplot = timeplot^0.5  #sec^0.5
xplot = sqrt_timeplot
# is.na(xplot)
sqrt_time = ifelse(is.na(xplot), 0, xplot)
# time = sqrt_time^2


# Calculate coefficients using [R]

### Smiles (I/t^0.5 v. t^0.5)

Smilesplot_vec = rep(1,12)

Iplot = Iactual_mat[1:12,] #cm
t05plot = xplot #sec^0.5

Smplot = Iplot/t05plot
# is.na(Sm)
Smplot = ifelse(is.na(Smplot), 0, Smplot)


for (i in 1:12) {
  
  Smilesplot_vec[i] = lm(Smplot[i,2:11] ~ t05plot[i,2:11]) 
  
}

Smilesplot <- data.frame(matrix(unlist(Smilesplot_vec), nrow=12, byrow=T))

Smplot_int = Smilesplot[,1]   #cm s^-0.5

Smplot_slope = Smilesplot[,2]  #cm s^-1

Iplotmodel = Smplot_int+Smplot_slope*t05plot

win.graph(width = 26, height = 24)
par(mfrow = c(1, 1))
par(oma=c(0,0,1,0)) 
par(mar=c(5,5,2,1.5))


z_plot <- z_plot[,c(10:20)]
zp1 <- as.numeric(matrix(z_plot[1,], ncol=1))
zp2 <- as.numeric(matrix(z_plot[2,], ncol=1))
zp3 <- as.numeric(matrix(z_plot[3,], ncol=1))
zp4 <- as.numeric(matrix(z_plot[4,], ncol=1))
zp5 <- as.numeric(matrix(z_plot[5,], ncol=1))
zp6 <- as.numeric(matrix(z_plot[6,], ncol=1))
zp7 <- as.numeric(matrix(z_plot[7,], ncol=1))
zp8 <- as.numeric(matrix(z_plot[8,], ncol=1))
zp9 <- as.numeric(matrix(z_plot[9,], ncol=1))
zp10 <- as.numeric(matrix(z_plot[10,], ncol=1))
zp11 <- as.numeric(matrix(z_plot[11,], ncol=1))
zp12 <- as.numeric(matrix(z_plot[12,], ncol=1))
zplot <- data.frame(Iactual,zp1,zp2,zp3,zp4,zp5,zp6,zp7,zp8,zp9,zp10,zp11,zp12)
zplot2 <- zplot[,2:13]^0.5
# zplot2 <- data.frame(Iactual,zplot2)
zplot3 <- Iactual/zplot2
theXlabzplot = expression(paste(bold("Time (s)"))) 
theYlabzplot = expression(paste(bold("Infiltration (cm)")))

plot(zplot$Iactual ~ zplot$zp1, xlab = theXlabzplot, ylab = theYlabzplot, type="p", 
     pch=21, cex = 1.3, col="black",  bg="darkgreen", ylim=c(0,15),xaxp=c(0,5000,5), yaxs='i', xaxs='i', xlim = c(0,5000), cex.lab=1.8, cex.axis=1.5, family="Arial")
lines(zplot$Iactual ~ zplot$zp1,lwd=1.5,col="darkgreen")
lines(zplot$Iactual ~ zplot$zp2,lwd=1.5,col="orange")
lines(zplot$Iactual ~ zplot$zp3,lwd=1.5,col="burlywood3")
lines(zplot$Iactual ~ zplot$zp4,lwd=1.5,col="red")
lines(zplot$Iactual ~ zplot$zp5,lwd=1.5,col="blue")
lines(zplot$Iactual ~ zplot$zp6,lwd=1.5,col="lightblue")
lines(zplot$Iactual ~ zplot$zp7,lwd=1.5,col="purple")
lines(zplot$Iactual ~ zplot$zp8,lwd=1.5,col="brown")
lines(zplot$Iactual ~ zplot$zp9,lwd=1.5,col="green")
lines(zplot$Iactual ~ zplot$zp10,lwd=1.5,col="gray")
lines(zplot$Iactual ~ zplot$zp11,lwd=1.5,col="gold")
lines(zplot$Iactual ~ zplot$zp12,lwd=1.5,col="pink")
points(zplot$Iactual ~ zplot$zp1, pch=21, cex = 1.3, col="black",  bg="darkgreen")
points(zplot$Iactual ~ zplot$zp2, pch=21, cex = 1.3, col="black",  bg="orange")
points(zplot$Iactual ~ zplot$zp3, pch=21, cex = 1.3, col="black",  bg="burlywood3")
points(zplot$Iactual ~ zplot$zp4, pch=21, cex = 1.3, col="black",  bg="red")
points(zplot$Iactual ~ zplot$zp5, pch=21, cex = 1.3, col="black",  bg="blue")
points(zplot$Iactual ~ zplot$zp6, pch=21, cex = 1.3, col="black",  bg="lightblue")
points(zplot$Iactual ~ zplot$zp7, pch=21, cex = 1.3, col="black",  bg="purple")
points(zplot$Iactual ~ zplot$zp8, pch=21, cex = 1.3, col="black",  bg="brown")
points(zplot$Iactual ~ zplot$zp9, pch=21, cex = 1.3, col="black",  bg="green")
points(zplot$Iactual ~ zplot$zp11, pch=21, cex = 1.3, col="black",  bg="gold")
points(zplot$Iactual ~ zplot$zp12, pch=21, cex = 1.3, col="black",  bg="pink")
points(zplot$Iactual ~ zplot$zp10, pch=21, cex = 1.3, col="black",  bg="gray")
text(300,14,"a)",cex=1.7)

# legend(4400, 11.3, c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"), 
#        pt.bg = c("darkgreen", "orange", "burlywood3", "red", "blue", "lightblue", "purple",
#                "brown", "green", "gray","gold","pink"),
#        col="black",
#      text.col = "black", pch = 21, title="Ring", title.col="black",
#      cex=1.3, bg = "white")

win.graph(width = 26, height = 24)
par(mfrow = c(1, 1))
par(oma=c(0,0,1,0)) 
par(mar=c(5,6,2,1.5))

theXlabzplot2 = expression(paste(bold("Time"^0.5~"(s"^0.5~")"))) 
theYlabzplot2 = expression(paste(bold("Infiltration*Time"^-0.5~"(cm"~s^-0.5~")")))

plot(zplot3[,1] ~ zplot2[,1], xlab = c(""), ylab = theYlabzplot2, type="p", 
     pch=21, cex = 1.3, col="black",  bg="darkgreen", ylim=c(0,.4),, yaxs='i', xaxs='i', xlim = c(0,100), cex.lab=1.8, cex.axis=1.5)
lines(Iplotmodel[1,] ~ t05plot[1,],lwd=2,col="darkgreen")
lines(Iplotmodel[2,] ~ t05plot[2,],lwd=2,col="orange")
lines(Iplotmodel[3,] ~ t05plot[3,],lwd=2,col="burlywood3")
lines(Iplotmodel[4,] ~ t05plot[4,],lwd=2,col="red")
lines(Iplotmodel[5,] ~ t05plot[5,],lwd=2,col="blue")
lines(Iplotmodel[6,] ~ t05plot[6,],lwd=2,col="lightblue")
lines(Iplotmodel[7,] ~ t05plot[7,],lwd=2,col="purple")
lines(Iplotmodel[8,] ~ t05plot[8,],lwd=2,col="brown")
lines(Iplotmodel[9,] ~ t05plot[9,],lwd=2,col="green")
lines(Iplotmodel[10,] ~ t05plot[10,],lwd=2,col="gold")
lines(Iplotmodel[11,] ~ t05plot[11,],lwd=2,col="pink")
lines(Iplotmodel[12,] ~ t05plot[12,],lwd=2,col="gray")
# lines(zplot3[,1] ~ zplot2[,1],lwd=1.5,col="darkgreen")
# lines(zplot3[,2] ~ zplot2[,2],lwd=1.5,col="orange")
# lines(zplot3[,3] ~ zplot2[,3],lwd=1.5,col="burlywood3")
# lines(zplot3[,4] ~ zplot2[,4],lwd=1.5,col="red")
# lines(zplot3[,5] ~ zplot2[,5],lwd=1.5,col="blue")
# lines(zplot3[,6] ~ zplot2[,6],lwd=1.5,col="lightblue")
# lines(zplot3[,7] ~ zplot2[,7],lwd=1.5,col="purple")
# lines(zplot3[,8] ~ zplot2[,8],lwd=1.5,col="brown")
# lines(zplot3[,9] ~ zplot2[,9],lwd=1.5,col="green")
# lines(zplot3[,10] ~ zplot2[,10],lwd=1.5,col="gray")
# lines(zplot3[,11] ~ zplot2[,11],lwd=1.5,col="gold")
# lines(zplot3[,12] ~ zplot2[,12],lwd=1.5,col="pink")
points(zplot3[,1] ~ zplot2[,1], pch=21, cex = 1.5, col="black",  bg="darkgreen")
points(zplot3[,2] ~ zplot2[,2], pch=21, cex = 1.5, col="black",  bg="orange")
points(zplot3[,3] ~ zplot2[,3], pch=21, cex = 1.5, col="black",  bg="burlywood3")
points(zplot3[,4] ~ zplot2[,4], pch=21, cex = 1.5, col="black",  bg="red")
points(zplot3[,5] ~ zplot2[,5], pch=21, cex = 1.5, col="black",  bg="blue")
points(zplot3[,6] ~ zplot2[,6], pch=21, cex = 1.5, col="black",  bg="lightblue")
points(zplot3[,7] ~ zplot2[,7], pch=21, cex = 1.5, col="black",  bg="purple")
points(zplot3[,8] ~ zplot2[,8], pch=21, cex = 1.5, col="black",  bg="brown")
points(zplot3[,9] ~ zplot2[,9], pch=21, cex = 1.5, col="black",  bg="green")
points(zplot3[,10] ~ zplot2[,10], pch=21, cex = 1.5, col="black",  bg="gold")
points(zplot3[,11] ~ zplot2[,11], pch=21, cex = 1.5, col="black",  bg="pink")
points(zplot3[,12] ~ zplot2[,12], pch=21, cex = 1.5, col="black",  bg="gray")
mtext(side = 1, line = 3.6, theXlabzplot2,cex=1.9)
# text(4,0.38,"b)",cex=1.7)
