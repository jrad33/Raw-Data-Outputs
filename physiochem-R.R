install.packages("xtable")
library(xtable)
##try to make physiochemical tables in R

chemprop <- read.csv("C:/Users/Jesse/Desktop/Neonicotinoids_15/Data/Column-Study/Raw-Data-Outputs/Chem-prop-R.csv")

attach(chemprop)

chemprop

options(xtable.floating = FALSE)

options(xtable.timestamp = "")


print(xtable(chemprop), type = "html")



Tab

View(Tab)
