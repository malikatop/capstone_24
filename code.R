# Malika Top
# Data Science Capstone Spring 2024

# Load the data
bangl <- read.csv("/Users/malikatop/Downloads/bangl_data.csv")

# Rename covariates
colnames(bangl)[2:22] = c("age_visit", "sex", "gest_weeks", 
                               "deliv_type", "birth_order", "drinkwater_cups", 
                               "hosp_child_yn", "pica_yn", "education", 
                               "educa_spouse", "smokenv", "home_emo", 
                               "home_avoid", "home_careg", "home_env", 
                               "home_play", "home_stim", "energy", "pb_ln_2", 
                               "mn_ln_2", "as_ln_2")


# Graphics to assess relationship between metal and CS @ various calcium percentiles
percentiles <- quantile(bangl$calc, probs = c(0.25, 0.5, 0.75))

breaks <- c(-Inf, percentiles, Inf)
labels <- c("25th Percentile", "50th Percentile", "75th Percentile", "< 25% and > 75%")

bangl$calc_cat <- cut(bangl$calc, breaks = breaks, labels=labels, include.lowest = TRUE)

calc_filtered <- bangl[bangl$calc_cat!="< 25% and > 75%",]
par(mfrow = c(2,2))
pb <- ggplot(calc_filtered, aes(x = pb_ln_2, y = y, color = calc_cat)) +
  geom_smooth(method = "lm", se = FALSE)+
  labs(x = "Cord Blood Pb Levels", y = "CS", color = "Calcium Percentile") +
  scale_color_manual(values = c("25th Percentile" = "blue", "50th Percentile" = "green", "75th Percentile" = "red" 
                               )) +
  theme_minimal()

mn <- ggplot(calc_filtered, aes(x = mn_ln_2, y = y, color = calc_cat)) +
  geom_smooth(method = "lm", se = FALSE)+
  labs(x = "Cord Blood Mn Levels", y = "CS", color = "Calcium Percentile") +
  scale_color_manual(values = c("25th Percentile" = "blue", "50th Percentile" = "green", "75th Percentile" = "red" 
  )) +
  theme_minimal()

as <- ggplot(calc_filtered, aes(x = as_ln_2, y = y, color = calc_cat)) +
  geom_smooth(method = "lm", se = FALSE)+
  labs(x = "Cord Blood As Levels", y = "CS", color = "Calcium Percentile") +
  scale_color_manual(values = c("25th Percentile" = "blue", "50th Percentile" = "green", "75th Percentile" = "red" 
  )) +
  theme_minimal()

# Correlation between calcium and various metals
# install.packages("ggplot2")
library(ggplot2)
corr_bangl <- data.frame(bangl$y, bangl$as_ln_2, bangl$pb_ln_2, bangl$mn_ln_2,
                         bangl$calc, bangl$age_visit)
colnames(corr_bangl)[1:6] = c("CS", "As", "Pb", 
                          "Mn", "Calcium", "Age Visit")

cor <- cor(corr_bangl)
library(corrplot)
par(mfrow=c(2,2))

corplot <- corrplot(cor, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)
pb
as
mn



# The independent variable is the metal exposures
# represented as a matrix to represent it as a mixture
A <- cbind(bangl$pb_ln_2,bangl$as_ln_2,bangl$mn_ln_2) 

# The moderator variable of interest is calcium
Z <- bangl$calc
Zc <- c(scale(Z,  center=TRUE, scale=FALSE))


# The dependent variable is cognitive score
Y <- bangl$y

bangl_sub <- data.frame(A, Z, Y)
colnames(bangl_sub)[1:3] = c("pb_ln", "as_ln", "mn_ln")


model <- lm(y~A+Z+quad_age+bangl$sex+bangl$education+bangl$educa_spouse
            +bangl$prot+quad_home+bangl$smokenv+bangl$drinkwater_cups+A*Z
            +Z*quad_age)
# Calculate HOME score sum
bangl$row_sum <- rowSums(bangl[ , c(13,18)], na.rm=TRUE)

# Rename matrix columns to be more informative
colnames(A)  <- c("Pb","As","Mn")

# Initial linear regression model with interaction between calcium and metals
model1 <- lm(Y~A+Z+bangl$age_visit+bangl$sex+bangl$education+bangl$educa_spouse
             +bangl$prot+bangl$row_sum+bangl$smokenv+bangl$drinkwater_cups+A*Z)
summary(model1)
coef(summary(model1))

# Secondary linear regression model with interaction between calcium and metals, 
# calcium and age visit
model2 <- lm(Y~A+Z+bangl$age_visit+bangl$sex+bangl$education+bangl$educa_spouse
             +bangl$prot+bangl$row_sum+bangl$smokenv+bangl$drinkwater_cups+A*Z
             + Z*bangl$age_visit)
summary(model2)
coef(summary(model2))
gvlma(model2)





###############Trying BKMR-CMA########################
source("source_BKMR_CMA.R")
# install.packages("bkmr")
library(bkmr)

# GQMS
# install.packages("gWQS")
library(gWQS)
# A is a matrix of the metal exposures, Arsenic, Manganese, Lead
A <- cbind(bangl$pb_ln,bangl$as_ln,bangl$mn_ln)
m <- bangl$calc
y <- bangl$y
bangl$row_sum <- rowSums(bangl[ , c(13,18)], na.rm=TRUE)
X <- bangl[c("age_visit", "sex", "education", "educa_spouse", 
             "prot", "row_sum", "smokenv", "drinkwater_cups")]


E.M <- NULL
E.Y <- bangl$age_visit

Z.M <- cbind(A,E.M) 
Z.Y <- cbind(A,E.Y) 
Zm.Y <- cbind(Z.Y,m)

colnames(Z.M)  <- colnames(A)  <- c("Pb","As","Mn")
colnames(Z.Y)  <- c("Pb","As","Mn","age")
colnames(Zm.Y) <- c("Pb","As","Mn","age","calc")

set.seed(1)
fit.y <- kmbayes(y=y, Z=Zm.Y, X=X, iter=5000, verbose=TRUE, varsel=FALSE) 
save(fit.y,file="bkmr_y.RData")

set.seed(2)
fit.y.TE <- kmbayes(y=y, Z=Z.Y, X=X, iter=5000, verbose=TRUE, varsel=FALSE) 
save(fit.y.TE,file="bkmr_y_TE.RData")

set.seed(3)
fit.m <- kmbayes(y=m, Z=Z.M, X=X, iter=5000, verbose=TRUE, varsel=FALSE) 
save(fit.m,file="bkmr_m.RData")


##### load models 
load("bkmr_y.RData")
load("bkmr_y_TE.RData")
load("bkmr_m.RData")

## mean level of confounders
X.predict <- matrix(colMeans(X),nrow=1)


## the change in exposure for which you want to estimate the mediation effects

## We will consider a change in all exposures from their 25th to 75th percentiles 
## fixing age at testing to its 10th and 90th percentiles 
## However, this contrast can be anything. 

## *** if modifiers are considered, you should fix the levels of the modifiers ***

astar.age10 <- c(apply(A, 2, quantile, probs=0.25), quantile(E.Y, probs=0.1))
astar.age90 <- c(apply(A, 2, quantile, probs=0.25), quantile(E.Y, probs=0.9))

a.age10 <- c(apply(A, 2, quantile, probs=0.75), quantile(E.Y, probs=0.1))
a.age90 <- c(apply(A, 2, quantile, probs=0.75), quantile(E.Y, probs=0.9))


## the index of the MCMC iterations to be used for inference 
sel<-seq(100,5001,by=15)

TE.age10 <- TE.bkmr(a=a.age10, astar=astar.age10, fit.y.TE=fit.y.TE, X.predict=X.predict, alpha=0.05, sel=sel, seed=122)

## look at the posterior mean, median, and 95% CI for TE
TE.age10$est



## repeat, now fixing age at testing to its 90th percentile 
TE.age90 <- TE.bkmr(a=a.age90, astar=astar.age90, fit.y.TE=fit.y.TE, X.predict=X.predict, sel=sel, seed=122)

## look at the posterior mean, median, and 95% CI for TE
TE.age90$est

CDE.age10 <- CDE.bkmr(a=a.age10, astar=astar.age10, m.quant=c(0.1,0.5,0.75), fit.y=fit.y, alpha=0.05, sel=sel, seed=777)

## look at the posterior mean, median, and 95% CI for the CDEs 
CDE.age10$est

mediationeffects.age90 <- mediation.bkmr(a.Y=a.age90, astar.Y=astar.age90, astar.M=astar, fit.m=fit.m, fit.y=fit.y, fit.y.TE=fit.y.TE,
                                         X.predict.M=X.predict, X.predict.Y=X.predict, alpha=0.05, sel=sel, seed=22, K=1000)
## save this object
save(mediationeffects.age90, file="mediationeffects_age90.RData")

## look at the posterior mean, median, and 95% CI for the TE, NDE, and NIE
mediationeffects.age90$est





