CEC <- read.csv("C:/Users/Jesse/Desktop/Neonicotinoids_15/Data/Column-Study/Raw-Data-Outputs/CEC-R.csv")

attach(CEC)

CEC

hist(CEC$CEC)

qqnorm(CEC$CEC) ## ~normal

qqline(CEC$CEC)

bartlett.test(CEC$CEC ~ CEC$Soil)  ## Appears to meet HOV

fligner.test(CEC$CEC ~ CEC$Soil)## Appears to be HOV assumption

boxplot(CEC$CEC ~ CEC$Soil, font.main = 2)

boxplot(CEC$CEC ~ CEC$Soil, ylab = "CEC (meq 1000 g^-1)", xlab = "Soil")


ANOVA1 = aov(CEC$CEC ~ CEC$Soil)

summary(ANOVA1) # significant p< 0.05

Tuk <- TukeyHSD(ANOVA1)

Tuk

par(mar = c(2,8,2,2))

plot(TukeyHSD(ANOVA1), las = 1)






library(multcompView)

multcompLetters4(ANOVA1, Tuk) ##### easy way to get sig letters report!


multcompl






exp_letters1
#Notice lowest mean treatments gets a "e"
#Ordered letters
multcompLetters2(y ~ treatments, exp_tukey$treatments[,"p adj"], experiment)
multcompLetters2(y ~ treatments, exp_tukey$treatments[,"p adj"], experiment, reversed = TRUE)


str



multcompLetters(x, compare = "<", threshold = 0.05, Letters = c(letters, LETTERS, "."), reversed = FALSE)

multcompLetters2(formula, x, data, ...)

multcompLetters3(z, y, x, data, ...)

multcompLetters4(object, comp, ...)

