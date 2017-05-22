library(ggplot2)

str(mpg)

qplot(displ, hwy, data = mpg, color = drv, geom = c("point", "smooth"), facets = .~ drv) ### color by drv variable in mpg


##### in facets = columns ~ rows

