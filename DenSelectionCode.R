#load packages
library(lme4)
library(survival)
library(Hmisc)
library(ggplot2)
library(AICcmodavg)
library(rafalib)
library(MuMIn)
library(MASS)
library(pROC)
library(car)
library(caret)
library(psych)
library(dplyr)
library(tidyr)
library(broom)
library(tidyverse)
library(magrittr)
library(irr)
library(splitstackshape)
library(jtools)
library(ggstance)
library(ggh4x)
library(grid)
library(gridExtra)
library(jtools)
library(interactions)
library(BAMMtools)
library(questionr)

####################################################################
##Habitat Selection at Third Order
####################################################################

#import data
Used <- read.csv("Used.csv")
Available <- read.csv("Available.csv")
Recency <- read.csv("Recency.csv")

#Join recency of use to Used Dens
Recency <- Recency %>%
  left_join(Used, by = "SiteID") %>%
  dplyr::select(SiteID, Used_this_year)

#replace NAs with zeroes:
Used[is.na(Used)] <- 0
Available[is.na(Available)] <- 0

#Make one file of both used and available sites
AllThirdDens <- rbind(Used, Available)

##First, assess whether there's a difference in covariates between recently used 
#and not-recently-used dens
R <- AllThirdDens %>%
  right_join(Recency, by = "SiteID") %>%
  dplyr::filter(Use == 1)

R2 <- 
  R %>%
  dplyr::select(Used_this_year, SLP, EAST, SOUTH, DIST_WATER, DIST_BLDG, DIST_RD, DIST_ED, densrd50, densrd25, densrd100, densrd200, densed25, densed50, densed100, densed200, densbldg25, densbldg50, densbldg100, densbldg200, AG_25, NAT_25, ANTH_25, GRASS_25, AG_50, NAT_50, ANTH_50, GRASS_50, AG_100, NAT_100, ANTH_100, GRASS_100, AG_200, NAT_200, ANTH_200, GRASS_200, WAT_200, WAT_25, WAT_50, WAT_100) %>% #select quanttaive variables and grouping variable (Use)
  gather(key = variable, value = value, -Used_this_year) %>% #plonks em in one list
  group_by(Used_this_year, variable) %>% 
  summarise(value = list(value)) %>%
  spread(Used_this_year, value) %>% 
  group_by(variable) %>%
  mutate(p_value = t.test(unlist(Y), unlist(N))$p.value,
         t_value = t.test(unlist(Y), unlist(N))$statistic)
#Write these results to a table
Recencythird_ttesttable <- R2[,-c(2:3)]
Recencythird_ttesttable
Recencythirdtabledf <- as.data.frame(Recencythird_ttesttable)
Recencythirdtabledf <- Recencythirdtabledf %>%
  dplyr::select(variable, p_value, t_value)
Recencythirdtabledf[is.na(Recencythirdtabledf)] <- 0
Recencythirdtabledf
write.csv(Recencythirdtabledf, file = "third_RECENCY.csv")


##Generate summary statistics for covariate values
summarythird <- AllThirdDens %>% 
  dplyr::group_by(Use) %>% 
  dplyr::summarise_if(is.numeric, .funs=list(mean, sd, min, max)) %>%
  tidyr::pivot_longer(cols = -Use, 
                      names_to = c('.value', 'variable'), 
                      names_sep = '_fn')
summarythird$variable[summarythird$variable==1] <- "mean"
summarythird$variable[summarythird$variable==2] <- "sd"
summarythird$variable[summarythird$variable==3] <- "min"
summarythird$variable[summarythird$variable==4] <- "max"
summarythird <- summarythird%>%
  dplyr::select("Use",
                "variable",
                "SLP", 
                "SOUTH", 
                "EAST", 
                "DIST_WATER", 
                "DIST_BLDG", 
                "DIST_RD", 
                "DIST_ED", 
                "densrd50", 
                "densed25", 
                "AG_25", 
                "ANTH_25", 
                "WAT_200", 
                "NAT_25", 
                "GRASS_25")

#Summarise t-test results for covariate values
AllThirdDens$Use2 <- ifelse(AllThirdDens$Use == 1, "Used", "Avail") 
Y <- 
  AllThirdDens %>%  #specify dataframe
  dplyr::select(Use2, SLP, EAST, SOUTH, DIST_WATER, DIST_BLDG, DIST_RD, DIST_ED, densrd50, densed25, AG_25, NAT_25, ANTH_25, GRASS_25, WAT_200) %>% #select quanttaive variables and grouping variable (Use)
  gather(key = variable, value = value, -Use2) %>% 
  group_by(Use2, variable) %>% 
  summarise(value = list(value)) %>%
  spread(Use2, value) %>% 
  group_by(variable) %>%
  mutate(p_value = t.test(unlist(Avail), unlist(Used))$p.value,
         t_value = t.test(unlist(Avail), unlist(Used))$statistic)
#write this to a table
third_ttesttable <- Y[,-c(2:3)]
head(third_ttesttable)

##Now combine t-test results table and summary table
summarythirdUse <- data.frame(t(dplyr::filter(summarythird, Use == "1")))
colnames(summarythirdUse) <- as.character(summarythirdUse[2,])
summarythirdUse <- summarythirdUse[-c(1:2),]
summarythirdUse <- tibble::rownames_to_column(summarythirdUse, var = "variable") ##OK, we're happy with out use table

summarythirdAvail <- data.frame(t(dplyr::filter(summarythird, Use == "0")))
colnames(summarythirdAvail) <- as.character(summarythirdAvail[2,])
summarythirdAvail <- summarythirdAvail[-c(1:2),]
summarythirdAvail <- tibble::rownames_to_column(summarythirdAvail, var = "variable") ##OK, we're happy with out AVAIL table

#Join three tables (i.e., Use, Avail, T-test results)
quant_summary_third <- third_ttesttable %>%
  left_join(summarythirdAvail, by = "variable") %>%
  left_join(summarythirdUse, by = "variable") #in table, x represents available, and y represents use

#reorder columns
quant_summary_third <- 
  quant_summary_third %>%
  dplyr::select(mean.y, sd.y, min.y, max.y, mean.x, sd.x, min.x, max.x, t_value, p_value)

#rename columns
quant_summary_third <- 
  quant_summary_third %>%
  rename("Variable" = variable, "Mean" = mean.y, "SD" = sd.y, "Min" = min.y, "Max" = max.y, "MeanA" = mean.x, "SDA" = sd.x, "MinA" = min.x, "MaxA" = max.x, "t-stat" = t_value, "p-value" = p_value)

#reorder from smallest to largest p values
quant_summary_third <- 
  quant_summary_third %>%
  arrange(`p-value`)

#replace p values < 0.001 with "<0.001"
quant_summary_third$`p-value`[(quant_summary_third$`p-value` <= 0.001)] <- "<0.001"

#Write to CSV
write.csv(quant_summary_third, file = "third_quantitative.csv")

#For one binary variable (proximity to river), perform Chi Square test
tblRVYN <- table(AllThirdDens$Use2, AllThirdDens$RV_Y_N) 
chiRVYN <- chisq.test(tblRVYN)
chiRVYN #X-squared = 10.601, df = 1, p-value = 0.00113


##Begin modelling process
#Create River Valley binary variable
AllThirdDens$RV_Y_N <- ifelse(AllThirdDens$DIST_RIV < 75, "Y", "N")

#Develop univariate models
model1 <- glm(Use~SLP, data = AllThirdDens, family = binomial)
model2 <- glm(Use~SOUTH, data = AllThirdDens, family = binomial)
model3 <- glm(Use~EAST, data = AllThirdDens, family = binomial)
model4 <- glm(Use~densbldg200, data = AllThirdDens, family = binomial)
model5 <- glm(Use~densbldg100, data = AllThirdDens, family = binomial)
model6 <- glm(Use~densbldg50, data = AllThirdDens, family = binomial)
model7 <- glm(Use~densbldg25, data = AllThirdDens, family = binomial)
model8 <- glm(Use~densrd200, data = AllThirdDens, family = binomial)
model9 <- glm(Use~densrd100, data = AllThirdDens, family = binomial)
model10 <- glm(Use~densrd50, data = AllThirdDens, family = binomial)
model11 <- glm(Use~densrd25, data = AllThirdDens, family = binomial)
model12 <- glm(Use~densed200, data = AllThirdDens, family = binomial)
model13 <- glm(Use~densed100, data = AllThirdDens, family = binomial)
model14 <- glm(Use~densed50, data = AllThirdDens, family = binomial)
model15 <- glm(Use~densed25, data = AllThirdDens, family = binomial)
model15.1 <- glm(Use~RV_Y_N, data = AllThirdDens, family = binomial)


summary(model1); confint(model1) #CI = 0.007055534 0.08068467; retain
summary(model2); confint(model2) #CI = -0.1831946 0.3900530; don't retain
summary(model3); confint(model3) #CI = -0.1888189 0.4127517; don't retain
summary(model4); confint(model4) #CI = -511.868005 516.0699526; don't retain
summary(model5); confint(model5) #CI = -699.583944 554.3103714; don't retain
summary(model6); confint(model6) #CI = -591.698651 518.9987947; don't retain
summary(model7); confint(model7) #CI = -659.300270 616.9294213; don't retain
summary(model8); confint(model8) #CI = -53.412661  2.7709504; retain
summary(model9); confint(model9) #CI = -73.812967 -17.7556030
summary(model10); confint(model10) #CI =  -121.0150 -40.2633340; road density performs best at 50-m scale
summary(model11); confint(model11) #CI = -86.821561 -20.2784764
summary(model12); confint(model12) #CI = -36.440913 17.2058075; don't retain
summary(model13); confint(model13) #CI = -30.305853  7.3408659
summary(model14); confint(model14) #CI = -19.088904  6.1216345
summary(model15); confint(model15) #CI = -17.155146  0.7214427; edge density performs best at 25-m scale 
summary(model15.1); confint(model15.1) #CI = -2.576853 -0.6616790; retain

#check slope, road, and edge density for performance as quadratics
qodel1 <- glm(Use~I(SLP^2) + SLP, data = AllThirdDens, family = binomial)
qodel10 <- glm(Use~I(densrd50^2) + densrd50, data = AllThirdDens, family = binomial)
qodel15 <- glm(Use~I(densed25^2) + densed25, data = AllThirdDens, family = binomial)

model1$aic; qodel1$aic #slp better as linear
model10$aic; qodel10$aic #equal performance; retain linear
model15$aic; qodel15$aic #ed density better as linear

##Determine the best decay terms
#no decay term
model16 <- glm(Use~DIST_WATER, data = AllThirdDens, family = binomial)
model17 <- glm(Use~DIST_ED, data = AllThirdDens, family = binomial)
model18 <- glm(Use~DIST_RD, data = AllThirdDens, family = binomial)
model19 <- glm(Use~DIST_BLDG, data = AllThirdDens, family = binomial)

#decay term = 50
model16a <- glm(Use~exp(-50/DIST_WATER), data = AllThirdDens, family = binomial)
model17a <- glm(Use~exp(-50/DIST_ED), data = AllThirdDens, family = binomial)
model18a <- glm(Use~exp(-50/DIST_RD), data = AllThirdDens, family = binomial)
model19a <- glm(Use~exp(-50/DIST_BLDG), data = AllThirdDens, family = binomial)

#decay term = 100
model16b <- glm(Use~exp(-100/DIST_WATER), data = AllThirdDens, family = binomial)
model17b <- glm(Use~exp(-100/DIST_ED), data = AllThirdDens, family = binomial)
model18b <- glm(Use~exp(-100/DIST_RD), data = AllThirdDens, family = binomial)
model19b <- glm(Use~exp(-100/DIST_BLDG), data = AllThirdDens, family = binomial)

#decay term = 250
model16c <- glm(Use~exp(-250/DIST_WATER), data = AllThirdDens, family = binomial)
model17c <- glm(Use~exp(-250/DIST_ED), data = AllThirdDens, family = binomial)
model18c <- glm(Use~exp(-250/DIST_RD), data = AllThirdDens, family = binomial)
model19c <- glm(Use~exp(-250/DIST_BLDG), data = AllThirdDens, family = binomial)

#decay term = 500
model16d <- glm(Use~exp(-500/DIST_WATER), data = AllThirdDens, family = binomial)
model17d <- glm(Use~exp(-500/DIST_ED), data = AllThirdDens, family = binomial)
model18d <- glm(Use~exp(-500/DIST_RD), data = AllThirdDens, family = binomial)
model19d <- glm(Use~exp(-500/DIST_BLDG), data = AllThirdDens, family = binomial)

#decay term = 1000
model16e <- glm(Use~exp(-1000/DIST_WATER), data = AllThirdDens, family = binomial)
model17e <- glm(Use~exp(-1000/DIST_ED), data = AllThirdDens, family = binomial)
model18e <- glm(Use~exp(-1000/DIST_RD), data = AllThirdDens, family = binomial)
model19e <- glm(Use~exp(-1000/DIST_BLDG), data = AllThirdDens, family = binomial)

#decay term = 1500
model16f <- glm(Use~exp(-1500/DIST_WATER), data = AllThirdDens, family = binomial)
model17f <- glm(Use~exp(-1500/DIST_ED), data = AllThirdDens, family = binomial)
model18f <- glm(Use~exp(-1500/DIST_RD), data = AllThirdDens, family = binomial)
model19f <- glm(Use~exp(-1500/DIST_BLDG), data = AllThirdDens, family = binomial)

#decay term = 25
model16g <- glm(Use~exp(-25/DIST_WATER), data = AllThirdDens, family = binomial)
model17g <- glm(Use~exp(-25/DIST_ED), data = AllThirdDens, family = binomial)
model18g <- glm(Use~exp(-25/DIST_RD), data = AllThirdDens, family = binomial)
model19g <- glm(Use~exp(-25/DIST_BLDG), data = AllThirdDens, family = binomial)

#Use AIC to determine best decay term for each variable
AIC(model16); AIC(model16a); AIC(model16b); AIC(model16c); AIC(model16d);
AIC(model16e); AIC(model16f); AIC(model16g)#1000 for water

AIC(model17); AIC(model17a); AIC(model17b); AIC(model17c); AIC(model17d);
AIC(model17e); AIC(model17f); AIC(model17g) #25 for ed

AIC(model18); AIC(model18a); AIC(model18b); AIC(model18c); AIC(model18d);
AIC(model18e); AIC(model18f); AIC(model18g) #25 for road

AIC(model19); AIC(model19a); AIC(model19b); AIC(model19c); AIC(model19d)
AIC(model19e); AIC(model19f); AIC(model19g) # Don't retain this variable


#Add decay distance as a covariate to data table
AllThirdDens$WaterDecay <- (exp(-1000/AllThirdDens$DIST_WATER))
AllThirdDens$RdDecay <- (exp(-25/AllThirdDens$DIST_RD))
AllThirdDens$EdDecay <- (exp(-25/AllThirdDens$DIST_ED))

#Determine whether distance covariates perform better as linear or quadratic
AIC(glm(Use~WaterDecay, data = AllThirdDens, family = binomial)) #531.3213
AIC(glm(Use~I(WaterDecay^2), data = AllThirdDens, family = binomial)) #533.4113
#linear form is better

AIC(glm(Use~EdDecay, data = AllThirdDens, family = binomial)) #540.4111 
AIC(glm(Use~I(EdDecay^2), data = AllThirdDens, family = binomial)) #542.5852
#linear form is better

AIC(glm(Use~RdDecay, data = AllThirdDens, family = binomial)) #521.4851
AIC(glm(Use~I(RdDecay^2), data = AllThirdDens, family = binomial)) #529.2519
#linear form is better


##Determine at which scale land cover is strongest
#25-m diameter
model20 <- glm(Use~AG_25, data = AllThirdDens, family = binomial)
model21 <- glm(Use~ANTH_25, data = AllThirdDens, family = binomial)
model23 <- glm(Use~GRASS_25, data = AllThirdDens, family = binomial)
model24 <- glm(Use~NAT_25, data = AllThirdDens, family = binomial)

#50-m diameter
model20a <- glm(Use~AG_50, data = AllThirdDens, family = binomial)
model21a <- glm(Use~ANTH_50, data = AllThirdDens, family = binomial)
model23a <- glm(Use~GRASS_50, data = AllThirdDens, family = binomial)
model24a <- glm(Use~NAT_50, data = AllThirdDens, family = binomial)

#100-m diameter
model20b <- glm(Use~AG_100, data = AllThirdDens, family = binomial)
model21b <- glm(Use~ANTH_100, data = AllThirdDens, family = binomial)
model23b <- glm(Use~GRASS_100, data = AllThirdDens, family = binomial)
model24b <- glm(Use~NAT_100, data = AllThirdDens, family = binomial)

#200-m diameter
model20c <- glm(Use~AG_200, data = AllThirdDens, family = binomial)
model21c <- glm(Use~ANTH_200, data = AllThirdDens, family = binomial)
model23c <- glm(Use~GRASS_200, data = AllThirdDens, family = binomial)
model24c <- glm(Use~NAT_200, data = AllThirdDens, family = binomial)

AIC(model20); AIC(model20a); AIC(model20b); AIC(model20c) #AG
#all within 2 AIC; 50-diameter is best (AIC = 542.453)

AIC(model21); AIC(model21a); AIC(model21b); AIC(model21c) #ANTH
#25, 50, and 100-m all within 2 AIC; 50-diameter is best (AIC = 536.9162)

AIC(model23); AIC(model23a); AIC(model23b); AIC(model23c) #GRASS
#25-m is best and > 2 AIC from others (AIC = 529.9031)

AIC(model24); AIC(model24a); AIC(model24b); AIC(model24c) #ANTH
#25-m, 50-m, and 100-m all within 2 AIC; 50-m is best (AIC = 527.0419)

##Determine whether land cover variables perform better as quadratics
AIC(glm(Use~AG_50, data = AllThirdDens, family = binomial)) # 542.453
AIC(glm(Use~I(AG_50^2), data = AllThirdDens, family = binomial)) #541.4609
#equally predictive

AIC(glm(Use~NAT_25, data = AllThirdDens, family = binomial)) #527.7535
AIC(glm(Use~I(NAT_25^2), data = AllThirdDens, family = binomial)) #526.7138
#equally predictive

AIC(glm(Use~ANTH_25, data = AllThirdDens, family = binomial)) #537.0509
AIC(glm(Use~I(ANTH_25^2), data = AllThirdDens, family = binomial)) #538.004
#equally predictive

AIC(glm(Use~GRASS_25, data = AllThirdDens, family = binomial)) #529.9031
AIC(glm(Use~I(GRASS_25^2), data = AllThirdDens, family = binomial)) #531.7616
#equally predictive


##Begin multivariate modelling 
#make a vector of only covariates with P > 0.25 un univariate GLMs and selected
#densities and decay terms
ThirdDensQuant <- AllThirdDens%>%
  dplyr::select("SLP", 
                "WaterDecay", 
                "RdDecay", 
                "EdDecay",
                "densrd50", 
                "densed25",
                "AG_25", 
                "ANTH_25", 
                "NAT_25", 
                "GRASS_25")

#Assess colinearity (Pearson's)
cor(ThirdDensQuant[sapply(ThirdDensQuant, is.numeric)], method = c("pearson"))
rcorr(as.matrix(ThirdDensQuant))
#correlation between (Road density + road decay), (edge density + edge decay), 
#and (NAT_25 + ANTH_25)

##Use Mumin to determine all possible subsets
options(na.action = "na.fail")

# fit model with all parameters
all.parmsthird<-glm(Use~SLP + WaterDecay + EdDecay + RdDecay + densrd50 + densed25 + ANTH_25 + NAT_25 + GRASS_25 + AG_25 +RV_Y_N, data = AllThirdDens, family = binomial)

densubset <- expression(!(densed25 && EdDecay) + (NAT_25 && ANTH_25) + (densrd50 && RdDecay))
resultsthird<-dredge(all.parmsthird, subset= densubset)
resultsthird

##OK, build models are all within 2 AIC These are the models using 25 m riv decay
thirdmodel1 <- get.models(resultsthird, 1)[[1]]
thirdmodel2 <- get.models(resultsthird, 2)[[1]]
thirdmodel3 <- get.models(resultsthird, 3)[[1]]
thirdmodel4 <- get.models(resultsthird, 4)[[1]]
thirdmodel5 <- get.models(resultsthird, 5)[[1]]
thirdmodel6 <- get.models(resultsthird, 6)[[1]]
thirdmodel7 <- get.models(resultsthird, 7)[[1]]
thirdmodel8 <- get.models(resultsthird, 8)[[1]]
thirdmodel9 <- get.models(resultsthird, 9)[[1]]
thirdmodel10 <- get.models(resultsthird, 10)[[1]]
thirdmodel11 <- get.models(resultsthird, 11)[[1]]

summary(thirdmodel1) #all covariates associated with P < 0.10 = TOP main effects model
summary(thirdmodel2) #all covariates associated with P < 0.10 = TOP main effects model
summary(thirdmodel3) #doesn't meet all P < 0.10
summary(thirdmodel4) #doesn't meet all P < 0.10
summary(thirdmodel5) #doesn't meet all P < 0.10
summary(thirdmodel6) #doesn't meet all P < 0.10
summary(thirdmodel7) #doesn't meet all P < 0.10
summary(thirdmodel8) #doesn't meet all P < 0.10
summary(thirdmodel9) #doesn't meet all P < 0.10
summary(thirdmodel10) #doesn't meet all P < 0.10
summary(thirdmodel11) #doesn't meet all P < 0.10


##Add biologically plausible interactions (to top two main effects models) 
#Biologically plausible interactions were: distrd*densrd, disted*densed, distwater*RV_Y_N and all land cover pairs
#Top main effects model + interactions:
intxthird1 <- glm(Use~WaterDecay:RV_Y_N + WaterDecay + RV_Y_N + GRASS_25 + RdDecay, data = AllThirdDens, family = binomial)
AIC(thirdmodel1); AIC(intxthird1) #no interactions model has lower AIC, delta AIC < 2

intxthird2 <- glm(Use~WaterDecay + RV_Y_N + GRASS_25 + GRASS_25:RV_Y_N + RdDecay, data = AllThirdDens, family = binomial)
AIC(thirdmodel1); AIC(intxthird2) #no interactions model has lower AIC, delta AIC < 2

intxthird3 <- glm(Use~WaterDecay + RV_Y_N + GRASS_25 + RdDecay + RdDecay:RV_Y_N, data = AllThirdDens, family = binomial)
AIC(thirdmodel1); AIC(intxthird3) #interactions model has lower AIC, delta AIC < 2

#Second main effects model + interactions
intxthird1.1 <- glm(Use~WaterDecay:RV_Y_N + WaterDecay + RV_Y_N + GRASS_25 + densrd50, data = AllThirdDens, family = binomial)
AIC(thirdmodel2); AIC(intxthird1.1) #no interactions model has lower AIC, delta AIC < 2

intxthird2.1 <- glm(Use~WaterDecay + RV_Y_N + GRASS_25 + GRASS_25:RV_Y_N + densrd50, data = AllThirdDens, family = binomial)
AIC(thirdmodel2); AIC(intxthird2.1) #no interactions model has lower AIC, delta AIC < 2

intxthird3.1 <- glm(Use~WaterDecay + RV_Y_N + GRASS_25 + densrd50 + densrd50:RV_Y_N, data = AllThirdDens, family = binomial)
AIC(thirdmodel2); AIC(intxthird3.1) #no interactions model has lower AIC, delta AIC < 2

#Average model results
Thirdavg <- model.avg(thirdmodel1, thirdmodel2)
summary(Thirdavg)
confint(Thirdavg)
#calculate odds ratio using covariates
exp(Thirdavg$coefficients)
exp(confint(Thirdavg))

#Determine scaled/ centered coefficients
Thirdavg.1 <- model.avg(thirdmodel1, thirdmodel2, beta = "sd")
summary(Thirdavg.1)
confint(Thirdavg.1)


##Determine variance inflation factors
vif(thirdmodel1) #max VIF = 1.1
vif(thirdmodel2) #max VIF = 1.0

##Determine ROC AUC
roc.third <- roc(AllThirdDens$Use~fitted(thirdmodel1))
plot(roc.third)
auc(roc.third) # AUC = 0.74

roc.third.1 <- roc(AllThirdDens$Use~fitted(thirdmodel2))
plot(roc.third.1)
auc(roc.third.1) # AUC = 0.73


##Cross validation
df <- AllThirdDens%>%
  dplyr::select("Use2", "WaterDecay", "GRASS_25", "RdDecay", "RV_Y_N", "densrd50")

#create index matrix
index <- createDataPartition(df$Use2, p = 0.8, list = FALSE, times = 1)

#create train and test
train_df <- df[index,]
test_df <- df[-index,]

#convert outcome variable to type factor
train_df$Use2 <- as.factor(train_df$Use2)
test_df$Use2 <- as.factor(test_df$Use2)

#specify number of folds and train control
ctrlspecs <- trainControl(method = "cv", number = 5,
                          savePredictions = "all",
                          classProbs = TRUE)
#specify logistic regression model
modelo1 <- train(Use2~WaterDecay + RdDecay + GRASS_25 + RV_Y_N,
                 data = train_df,
                 method = "glm", family = binomial,
                 trControl = ctrlspecs)
print(modelo1)
summary(modelo1)
varImp(modelo1)

#apply model to test_df
#predict outcome using model from train_df applied to test_df 
predictions1 <- predict(modelo1, newdata = test_df)
predictions1

#create confusion matrix
confusionMatrix(data = predictions1, test_df$Use2) #kappa = 0.11, Accuracy = 0.74
#kappa and accuracy change slightly depending on seed set randomly by R

#repeat CV for second main effects model
modelo2 <- train(Use2~WaterDecay + densrd50 + GRASS_25 + RV_Y_N,
                 data = train_df,
                 method = "glm", family = binomial,
                 trControl = ctrlspecs)
print(modelo2) #Kappa = 0, accuracy = 0.75
summary(modelo2)
varImp(modelo2)
predictions2 <- predict(modelo1, newdata = test_df)
predictions2
confusionMatrix(data = predictions2, test_df$Use2)

##Get scaled + ceneterd covariate estimates for final models
export_summs(thirdmodel1, scale = TRUE,
             error_format = "[{conf.low}, {conf.high}]")
export_summs(thirdmodel2, scale = TRUE,
             error_format = "[{conf.low}, {conf.high}]")

####################################################################
##Repeat Analysis with den areas as units of replication
####################################################################

####################################################################
##Habitat Selection at Fourth Order
####################################################################
#import data
AllFourthDens <- read.csv("Den_Survey_Info.csv")

#Remove single NA value in Water column
AllFourthDens$water[is.na(AllFourthDens$water)] <- 0

#split into use sites and available sites
UseSites <- dplyr::filter(AllFourthDens, AllFourthDens$Use == "1")
AvailSites <- dplyr::filter(AllFourthDens, AllFourthDens$Use == "0")

#calculate mean East and South for available and use sites
AvgEastUse <- mean(na.omit(UseSites$EAST))
AvgEastAvail <- mean(na.omit(AvailSites$EAST))
AvgSouthUse <- mean(na.omit(UseSites$SOUTH))
AvgSouthAvail <- mean(na.omit(AvailSites$SOUTH))

#replace instances of 'NA' in East and South Index with average for Use or Avail
#sites, respectively
UseSites$EAST[is.na(UseSites$EAST)] <- AvgEastUse  # Replace with mean
AvailSites$EAST[is.na(AvailSites$EAST)] <- AvgEastAvail  # Replace with mean
UseSites$SOUTH[is.na(UseSites$SOUTH)] <- AvgSouthUse  # Replace with mean
AvailSites$SOUTH[is.na(AvailSites$SOUTH)] <- AvgSouthAvail  # Replace with mean

#Combine Use and Available sites into one dataframe
AllFourthDens <- rbind(UseSites, AvailSites)

##Generate Summary statistics for quantitative covariates
summary <- AllFourthDens %>% 
  dplyr::group_by(Use) %>% 
  dplyr::summarise_if(is.numeric, .funs=list(mean, sd, min, max)) %>%
  tidyr::pivot_longer(cols = -Use, 
                      names_to = c('.value', 'variable'), 
                      names_sep = '_fn')
summary$variable[summary$variable==1] <- "mean"
summary$variable[summary$variable==2] <- "sd"
summary$variable[summary$variable==3] <- "min"
summary$variable[summary$variable==4] <- "max"
summary <- summary[,-c(3:10, 12:13, 16:18,20,27:29)]

#Summarize t-test results for quantitative covariates
AllFourthDens$Use2 <- ifelse(AllFourthDens$Use == 1, "Used", "Avail")
AllFourthDens$Use2 <- factor(AllFourthDens$Use2, levels = c("Used", "Avail"))
X <- AllFourthDens %>% 
  dplyr::select(Use2, Use, avg_slope, EAST, SOUTH, avgCC, avgcc, tree, shrub, herb, inorg, water)
X <- AllFourthDens %>%
  gather(key = variable, value = value, -Use2) %>%
  group_by(Use2, variable) %>% 
  summarise(value = list(value)) %>%
  spread(Use2, value) %>% 
  group_by(variable) %>%
  mutate(p_value = t.test(unlist(Avail), unlist(Used))$p.value,
         t_value = t.test(unlist(Avail), unlist(Used))$statistic)
#Write to table
fourth_ttesttable <- X[,-c(2:3)]
head(fourth_ttesttable)

#Combine t-test results and summary table
summaryUse <- data.frame(t(dplyr::filter(summary, Use == "1")))
colnames(summaryUse) <- as.character(summaryUse[2,])
summaryUse <- summaryUse[-c(1:2),]
summaryUse <- tibble::rownames_to_column(summaryUse, var = "variable")

summaryAvail <- data.frame(t(dplyr::filter(summary, Use == "0")))
colnames(summaryAvail) <- as.character(summaryAvail[2,])
summaryAvail <- summaryAvail[-c(1:2),]
summaryAvail <- tibble::rownames_to_column(summaryAvail, var = "variable")

#Join all three tables
quant_summary_fourth <- fourth_ttesttable %>%
  left_join(summaryAvail, by = "variable") %>%
  left_join(summaryUse, by = "variable") #in table, x represents available, and y represents use

#reorder columns
quant_summary_fourth <- 
  quant_summary_fourth %>%
  select(mean.y, sd.y, min.y, max.y, mean.x, sd.x, min.x, max.x, t_value, p_value)

#rename columns
quant_summary_fourth <- 
  quant_summary_fourth %>%
  rename("Variable" = variable, "Mean" = mean.y, "SD" = sd.y, "Min" = min.y, 
         "Max" = max.y, "MeanA" = mean.x, "SDA" = sd.x, "MinA" = min.x, 
         "MaxA" = max.x, "t-stat" = t_value, "p-value" = p_value)

#Reorder from smallest to largest p values
quant_summary_fourth <- 
  quant_summary_fourth %>%
  arrange(`p-value`)

#Replace p values < 0.001 with "<0.001"
quant_summary_fourth$`p-value`[(quant_summary_fourth$`p-value` <= 0.001)] <- "<0.001"

#Replace variable names with final Variable names
quant_summary_fourth$Variable <- c("Hiding Cover", "Average Slope", "Shrub Cover", 
                                   "East Index", "Herb Cover", "Tree Cover", 
                                   "South Index", "Bare Ground Cover", 
                                   "Canopy Cover", "Water Cover")

#Write to CSV
write.csv(quant_summary_fourth, file = "fourth_quantitative.csv")

#Get mean and sd for aspect at fourth order
UseAngles <- UseSites%>%
  dplyr::select(asp)
UseAngles <- na.omit(UseAngles)

AvailAngles <- AvailSites%>%
  dplyr::select(asp)
AvailAngles <- na.omit(AvailAngles)

anglecir1 =  circular(UseAngles, type="angles", units="degrees",modulo="2pi", template='geographics')
anglecir2 =  circular(AvailAngles, type="angles", units="degrees",modulo="2pi", template='geographics')

summary(anglecir1)
summary(anglecir2)

circular::sd(anglecir1)
circular::sd(anglecir2)

##Determine differences in quantitative covariate values between recently used
#dens and those that weren't recently used
W <- 
  AllFourthDens %>%
  select(Used_this_year, avg_slope, EAST, SOUTH, avgCC, avgcc, tree, shrub, herb, inorg, water) %>%
  gather(key = variable, value = value, -Used_this_year) %>%
  group_by(Used_this_year, variable) %>% 
  summarise(value = list(value)) %>%
  spread(Used_this_year, value) %>% 
  group_by(variable) %>%
  mutate(p_value = t.test(unlist(Y), unlist(N))$p.value,
         t_value = t.test(unlist(Y), unlist(N))$statistic)

#Develop a table
Recencyfourth_ttesttable <- W[,-c(2:3)]
Recencyfourth_ttesttable
Recencyfourth_ttesttable <- as.data.frame(Recencyfourth_ttesttable)
Recencyfourth_ttesttable <- Recencyfourth_ttesttable %>%
  dplyr::select(p_value, t_value, variable)
write.csv(Recencyfourth_ttesttable, file = "fourth_recency_ttest.csv")
#No significant differences

##Summarise qualitative data
#Summarise mode
calculate_mode <- function(x) {
  uniqx <- unique (x)
  uniqx[which.max(tabulate(match(x,uniqx)))]
} #makes a function that should calculate mode

qualmodes1 <- AllFourthDens%>%
  group_by(Use) %>%
  summarise(garbage = calculate_mode(garbage))

qualmodes2 <- AllFourthDens%>%
  group_by(Use2) %>%
  summarise(soil = calculate_mode(soil))

qualmodes3 <- AllFourthDens%>%
  group_by(Use2) %>%
  summarise(manipulation = calculate_mode(manipulation))

qualmodes4 <- AllFourthDens%>%
  group_by(Use2) %>%
  summarise(mesopos = calculate_mode(mesopos))

qualmodes <- qualmodes1 %>%
  left_join(qualmodes2, by = "Use2") %>%
  left_join(qualmodes3, by = "Use2") %>%
  left_join(qualmodes4, by = "Use2")


##Use chi square tests to determine differences in qualitative variables between
#used and available sites 
tbl2mesopos <- table(AllFourthDens$Use2, AllFourthDens$mesopos) 
tbl2mesopos
chimesopos <- chisq.test(tbl2mesopos)

tbl2soil <- table(AllFourthDens$Use2, AllFourthDens$soil) 
tbl2soil
chisoil <- chisq.test(tbl2soil)

tbl2garbage <- table(AllFourthDens$Use2, AllFourthDens$garbage) 
tbl2garbage
chigarbage <- chisq.test(tbl2garbage)

tbl2manipulation <- table(AllFourthDens$Use2, AllFourthDens$manipulation) 
tbl2manipulation
chimanipulation <- chisq.test(tbl2manipulation)

chisquare<- c(chimesopos$statistic, chisoil$statistic, chigarbage$statistic, chimanipulation$statistic)
pvalue <- c(chimesopos$p.value, chisoil$p.value, chigarbage$p.value, chimanipulation$p.value)
variable <- c("mesopos", "soil", "garbage", "manipulation")
chiresults <- data.frame(variable,chisquare,pvalue, row.names = NULL)

#Merge two tables
qualmodessumm <- data.frame(t(qualmodes)) #first get the first table in shape
colnames(qualmodessumm) <- as.character(qualmodessumm[1,])
qualmodessumm <- qualmodessumm[-1,]
qualmodessumm <- tibble::rownames_to_column(qualmodessumm, var = "variable") ##OK, we're happy with out use table

qualmodessummUse <- qualmodessumm[,c(1,3)]
qualmodessummAvail <- qualmodessumm[,c(1:2)]

fourth_chisummary <- qualmodessummUse %>%
  left_join(qualmodessummAvail, by = "variable") %>%
  left_join(chiresults, by = "variable")

#Rename columns
fourth_chisummary <- 
  fourth_chisummary %>%
  rename("Variable" = variable, "Mode" = Used, "Mode2" = Avail, "|test stat|" = chisquare, "p-value" = pvalue)

#Reorder from smallest to largest p values
fourth_chisummary <- 
  fourth_chisummary %>%
  arrange(`p-value`)

#Replace any p values < 0.001 with "<0.001"
fourth_chisummary$`p-value`[(fourth_chisummary$`p-value` <= 0.001)] <- "<0.001"

#Replace variable names with final headings
fourth_chisummary$Variable <- c("Garbage Manipulation", "Mesoslope Position", "Soil Type", "Garbage Quantity")

#Write to CSV
write.csv(fourth_chisummary, file = "fourth_qualitative.csv")

##Assess differences between dens that were and were not used recently (qualitative covariates)
Use4 <- dplyr::filter(AllFourthDens, Use == "1")
Use4_split <- split(Use4, Use4$Used_this_year)
str(Use4_split)

#make new use vs non use dataframes
Yes4 <- Use4_split$`Y`
No4 <- Use4_split$`N`
Unk4 <- Use4_split$`UNK`

#Compare using chi square tests or fisher tests
tblmesopos <- table(Use4$Used_this_year, Use4$mesopos) 
fisher.test(tblmesopos)

tblsoil <- table(Use4$Used_this_year, Use4$soil) 
fisher.test(tblsoil,simulate.p.value=TRUE,B=1e7)

tblgarbage <- table(Use4$Used_this_year, Use4$garbage) 
fisher.test(tblgarbage)

tblmanipulation <- table(Use4$Used_this_year, Use4$manipulation) 
fisher.test(tblmanipulation)

##Begin modelling approach

#Build univariate models
fourth_mod1 <- glm(Use~avg_slope, data = AllFourthDens, family = binomial)
fourth_mod2 <- glm(Use~SOUTH, data = AllFourthDens, family = binomial)
fourth_mod3 <- glm(Use~EAST, data = AllFourthDens, family = binomial)
fourth_mod4 <- glm(Use~avgCC, data = AllFourthDens, family = binomial)
fourth_mod5 <- glm(Use~avgcc, data = AllFourthDens, family = binomial)
fourth_mod6 <- glm(Use~tree, data = AllFourthDens, family = binomial)
fourth_mod7 <- glm(Use~shrub, data = AllFourthDens, family = binomial)
fourth_mod8 <- glm(Use~herb, data = AllFourthDens, family = binomial)
fourth_mod9 <- glm(Use~inorg, data = AllFourthDens, family = binomial)
fourth_mod10 <- glm(Use~water, data = AllFourthDens, family = binomial)
fourth_mod11 <- glm(Use~soil, data = AllFourthDens, family = binomial)
fourth_mod12 <- glm(Use~mesopos, data = AllFourthDens, family = binomial)

summary(fourth_mod1); confint(fourth_mod1) #CI = 0.08214548  0.2086493, retain avg_slope
summary(fourth_mod2); confint(fourth_mod2) #CI = -0.6100876 0.2228572, don't retain south
summary(fourth_mod3); confint(fourth_mod3) #CI = 0.1424270 1.0444991, retain EAST
summary(fourth_mod4); confint(fourth_mod4) #CI = -0.006919326 0.01764401; don't retain Avg CC
summary(fourth_mod5); confint(fourth_mod5) #CI = 0.0443632  0.08192315; retain avg cc
summary(fourth_mod6); confint(fourth_mod6) #CI = -0.02206325 0.00198683; retain tree
summary(fourth_mod7); confint(fourth_mod7) #CI = 0.01376217  0.03532629; retain shrub
summary(fourth_mod8); confint(fourth_mod8) #CI = -0.01707128 -0.0002531265; retain herb
summary(fourth_mod9); confint(fourth_mod9) #CI = -0.004966525 0.01250596; don't retain bare ground
summary(fourth_mod10); confint(fourth_mod10) #CI = -0.06539192 0.1651048; don't retain water
summary(fourth_mod11); confint(fourth_mod11) #don't retain soil
summary(fourth_mod12); confint(fourth_mod12) #CI = variable; retain slope position

#Assess colinearity among variables
FourthDensQuant <- AllFourthDens[,-c(1:23,25:26, 29:35, 37:40, 47:53)]
cor(FourthDensQuant[sapply(FourthDensQuant, is.numeric)], method = c("pearson"))
rcorr(as.matrix(FourthDensQuant))
#herb and avgCC are correlated, but avgCC is not retained for main effects model

#Test all possible subsets
# fit model with all parameters
all.parmsfourth<-glm(Use~avg_slope + avgcc + EAST + tree + shrub + herb + mesopos, data = AllFourthDens, family = binomial)
resultsfourth<-dredge(all.parmsfourth)
resultsfourth

#save final three models within 2 AIC
finalfourthmod1 <- get.models(resultsfourth, 1)[[1]]
finalfourthmod2 <- get.models(resultsfourth, 2)[[1]]
finalfourthmod3 <- get.models(resultsfourth, 3)[[1]]

summary(finalfourthmod1) #all covariates associated with P < 0.10
summary(finalfourthmod2) #does not meet significance requirement
summary(finalfourthmod3) #all covariates associated with P < 0.10

#Test first and third main effects models for improvement using a series of 
#biologically plausible interaction terms as follows: (1) tree+shrub, (2) tree+herb,
#(3) shrub+herb, (4) avgslope+east, (5) avgslope+avgcc, (6) avgslp+mesopos, (7) avgcc+mespos
Intxfourthmod1 <- glm(Use~avgcc + avg_slope + EAST + herb + shrub + herb*shrub, data = AllFourthDens, family = binomial)
Intxfourthmod2 <- glm(Use~avgcc + avg_slope + EAST + herb + tree + herb*tree, data = AllFourthDens, family = binomial)
Intxfourthmod3 <- glm(Use~avgcc + avg_slope + EAST + herb + shrub + tree + shrub*tree, data = AllFourthDens, family = binomial)
Intxfourthmod4 <- glm(Use~avgcc + avg_slope*EAST + avg_slope + EAST + herb, data = AllFourthDens, family = binomial)
Intxfourthmod5 <- glm(Use~avg_slope*avgcc + avg_slope + avgcc + EAST + herb, data = AllFourthDens, family = binomial)
Intxfourthmod6 <- glm(Use~avgcc + avg_slope*mesopos + avg_slope + mesopos + EAST + herb, data = AllFourthDens, family = binomial)
Intxfourthmod7 <- glm(Use~avgcc*mesopos + avgcc + mesopos + avg_slope + EAST + herb, data = AllFourthDens, family = binomial)

#Develop interaction models based on third main effects model
Intxfourthmod8 <- glm(Use~avgcc + avg_slope + EAST + tree + shrub + tree*shrub, data = AllFourthDens, family = binomial)
Intxfourthmod9 <- glm(Use~avgcc + avg_slope + EAST + shrub + tree + herb + tree*herb, data = AllFourthDens, family = binomial)
Intxfourthmod10 <- glm(Use~avgcc + avg_slope + EAST + shrub + herb + shrub*herb, data = AllFourthDens, family = binomial)
Intxfourthmod11 <- glm(Use~avgcc + avg_slope*EAST + avg_slope + EAST + shrub, data = AllFourthDens, family = binomial)
Intxfourthmod12 <- glm(Use~avgcc*avg_slope + avgcc + avg_slope + EAST + shrub, data = AllFourthDens, family = binomial)
Intxfourthmod13 <- glm(Use~avgcc + avg_slope*mesopos + avg_slope + mesopos + EAST + shrub, data = AllFourthDens, family = binomial)
Intxfourthmod14 <- glm(Use~avgcc*mesopos + avgcc + mesopos + avg_slope + EAST + shrub, data = AllFourthDens, family = binomial)

#Compare interaction models to original models
AIC(finalfourthmod1); AIC(Intxfourthmod1) #no interactions better, > 2 AIC
AIC(finalfourthmod1); AIC(Intxfourthmod2) #no interactions better, > 2 AIC
AIC(finalfourthmod1); AIC(Intxfourthmod3) #no interactions better, > 2 AIC
AIC(finalfourthmod1); AIC(Intxfourthmod4) #no interactions better, < 2 AIC
AIC(finalfourthmod1); AIC(Intxfourthmod5) #no interactions better, < 2 AIC
AIC(finalfourthmod1); AIC(Intxfourthmod6) #no interactions better, > 2 AIC
AIC(finalfourthmod1); AIC(Intxfourthmod7) #no interactions better, > 2 AIC

AIC(finalfourthmod3); AIC(Intxfourthmod8) #no interactions better, < 2 AIC
AIC(finalfourthmod3); AIC(Intxfourthmod9) #no interactions better, < 2 AIC
AIC(finalfourthmod3); AIC(Intxfourthmod10) #no interactions better, < 2 AIC
AIC(finalfourthmod3); AIC(Intxfourthmod11) #interactions better, < 2 AIC
AIC(finalfourthmod3); AIC(Intxfourthmod12) #no interactions better, < 2 AIC
AIC(finalfourthmod3); AIC(Intxfourthmod13) #no interactions better, > 2 AIC
AIC(finalfourthmod3); AIC(Intxfourthmod14) #no interactions better, > 2 AIC
#no interaction terms improve model performance

#Average model results
Fourthavg <- model.avg(finalfourthmod1, finalfourthmod3)
summary(Fourthavg)
confint(Fourthavg)
#calculate odds ratio using covariates
exp(Fourthavg$coefficients)
exp(confint(Fourthavg))

#Determine scaled/ centered coefficients
Fourthavg.1 <- model.avg(finalfourthmod1, finalfourthmod3, beta = "sd")
summary(Fourthavg.1)
confint(Fourthavg.1)


#Determine VIFs
vif(finalfourthmod1) #max vif = 1.2
vif(finalfourthmod3) #max vif = 1.1

#Determine ROC AUC
roc.fourth <- roc(AllFourthDens$Use~fitted(finalfourthmod1))
plot(roc.fourth)
auc(roc.fourth) #AUC = 0.868

roc.fourth2 <- roc(AllFourthDens$Use~fitted(finalfourthmod3))
plot(roc.fourth2)
auc(roc.fourth2) #AUC = 0.8644

#Complete Cross Validation
df2 <- AllFourthDens%>%
  dplyr::select("Use2", "avg_slope", "EAST", "avgcc", "herb", "shrub")
#create index matrix
index2 <- createDataPartition(df2$Use2, p = 0.8, list = FALSE, times = 1)
#create train and test
train_df2 <- df2[index2,]
test_df2 <- df2[-index2,]

#convert outcome variable to type factor
train_df2$Use2 <- as.factor(train_df2$Use2)
test_df2$Use2 <- as.factor(test_df2$Use2)

#specify number of folds and train control and all that
ctrlspecs2 <- trainControl(method = "cv", number = 5,
                           savePredictions = "all",
                           classProbs = TRUE)
#specify logistic regression model
modelp1 <- train(Use2~avg_slope + EAST + avgcc + herb,
                 data = train_df2,
                 method = "glm", family = binomial,
                 trControl = ctrlspecs2)
print(modelp1)
summary(modelp1)
varImp(modelp1)

#apply model to test_df
#predict outcome using model from train_df applied to test_df 
predictions4.1 <- predict(modelp1, newdata = test_df2)
predictions4.1

#create confusion matrix
confusionMatrix(data = predictions4.1, test_df2$Use2) #kappa = 0.55, accuracy = 0.775

#Repeat CV for second main effects model
modelp2 <- train(Use2~avg_slope + EAST + avgcc + shrub,
                 data = train_df2,
                 method = "glm", family = binomial,
                 trControl = ctrlspecs2)
print(modelp2); summary(modelp2); varImp(modelp2)
predictions4.2 <- predict(modelp2, newdata = test_df2)
predictions4.2
confusionMatrix(data = predictions4.2, test_df2$Use2) #kappa = 0.45, accuracy = 0.725

##Get scaled and centered covariate estimates
export_summs(finalfourthmod1, scale = TRUE,
             error_format = "[{conf.low}, {conf.high}]")
export_summs(finalfourthmod3, scale = TRUE,
             error_format = "[{conf.low}, {conf.high}]")

####################################################################
##Analysis of Conflict, Den Proximity, Season, and Land Cover
####################################################################
#import dataset
jonreports <- read.csv("JonReports.csv")

#remove reports where coyote resp is unknown
conflict <- dplyr::filter(jonreports, jonreports$CoyoteResp !="Unknown" )

#remove reports that occur > 500m from path
conflict <- conflict %>%
  dplyr::filter(PathDist < 500)

#manipulate land cover variable
conflict$LandCover <- as.factor(conflict$SiteTypeRa)
conflict$LandCover <- ifelse(conflict$LandCover == "1", "Naturalized urban area", "Modified area")
conflict$LandCover <- factor(conflict$LandCover, levels = c("Naturalized urban area", "Modified area"))

#define threshold for conflict
conflict$benign <- ifelse(conflict$CoyoteResp >= "6",1,0)

#remove irrelevant columns
conflict <- conflict %>%
  dplyr::select(DenDist, PathDist, LandCover, Year, CoyoteResp, benign, Season, Month)

#Create decay term (alpha = 500)
conflict$DenDecay <- (exp(-500/conflict$DenDist))

conflict <- conflict %>%
  mutate(Season3 = case_when(Month == "4" ~ "Pup rearing",
                             Month == "5" ~ "Pup rearing",
                             Month == "6" ~ "Pup rearing",
                             Month == "7" ~"Pup rearing"))
conflict$Season3[is.na(conflict$Season3)] <- "Outside of pup rearing"

#Explore data using t-test
t.test(conflict$DenDist~conflict$benign) #t = 2.4535, df = 768.96, p-value = 0.01437
#conflict events occur closer to dens than conflict events

tblconflict <- table(conflict$benign, conflict$Season3) 
chisq.test(tblconflict) #X-squared = 65.315, df = 1, p-value = 6.384e-16
tblconflict
#conflict events occur more than expected during the pup rearing period

tblLC <- table(conflict$LandCover, conflict$benign)
chisq.test(tblLC) #X-squared = 36.117, df = 1, p-value = 1.859e-09
tblLC
#conflict events occur more than expected in natural areas

#Develop subsets based on Land Cover
ModData <- conflict%>%
  dplyr::filter(LandCover == "Modified area")
NatData <- conflict%>%
  dplyr::filter(LandCover == "Naturalized urban area")

#Develop a priori models
conf.model1 <- glm(benign~DenDecay*Season3, data = conflict, family = binomial)
conf.model2 <- glm(benign~DenDecay*Season3, data = NatData, family = binomial)
conf.model3 <- glm(benign~DenDecay*Season3, data = ModData, family = binomial)
conf.model4 <- glm(benign~Season3*LandCover, data = conflict, family = binomial)

summary(conf.model1)
summary(conf.model2)
summary(conf.model3)
summary(conf.model4)

#Summarise model results
export_summs(conf.model1, conf.model2, conf.model3, conf.model4,
             scale = TRUE,
             ci_level = 0.95,
             stars = NULL,
             error_format = '({statistic}, p = {p.value}), CI [{conf.low}, {conf.high}]',
             model.names = c("Den proximity and season (All)", 
                             "Den proximity and season (Naturalized urban area)",
                             "Den proximity and season (Modified area)",
                             "Land cover and season (All)"),
             coefs = c("Den distance decay term" = "DenDecay", 
                       "Pup rearing" = "Season3Pup rearing",
                       "Den distance * Pup rearing" = "DenDecay:Season3Pup rearing",
                       "Land cover" = "LandCover"),
             to.file = "xlsx",
             file.name = "ConflictBetas.xlsx")

#Determine Odds ratio and CIs
summ(conf.model1, exp = TRUE)
summ(conf.model2, exp = TRUE)
summ(conf.model3, exp = TRUE)
summ(conf.model4, exp = TRUE)


####################################################################
##Repeat for subset of conflict reports occurring in 2021
####################################################################


####################################################################
##Plot showing Third Order covariate selection (Figure 2 in text)
####################################################################
fakemodel <- glm(Use~SLP + NAT_25 + EdDecay, data = AllThirdDens, family = binomial)
Thirdplotsumms <- plot_summs(fakemodel, thirdmodel3,thirdmodel4, thirdmodel5,
                             thirdmodel6,thirdmodel7, thirdmodel8, thirdmodel9,
                             thirdmodel10, thirdmodel11, thirdmodel1, thirdmodel2,
                             coefs = c("Within 75 m of riverbank" = "RV_Y_NY",
                                       "Decay distance\n to water" = "WaterDecay",
                                       "Mowed grass (25-m)" = "GRASS_25",
                                       "Decay distance\n to road" = "RdDecay",
                                       "Road density (50-m)" = "densrd50",
                                       "Slope (%)" = "SLP",
                                       "Agricultural (25-m)" = "AG_25",
                                       "Anthropogenic (25-m)" = "ANTH_25",
                                       "Natural (25 m)" = "NAT_25",
                                       "Edge density (25-m)" = "densed25",
                                       "Decay distance\n to edge" = "EdDecay"),
                             scale = TRUE,
                             point.shape =  FALSE,
                             colors = c("white", "dodgerblue1", "dodgerblue1", "dodgerblue1", "dodgerblue1", "dodgerblue1", "dodgerblue1", "dodgerblue1","dodgerblue1","dodgerblue1", "dodgerblue4", "dodgerblue4"))
Thirdplotsumms1 <- Thirdplotsumms + xlim(-3.5,1.0) + ggtitle("Den sites") + 
  theme(legend.position = "none", 
        axis.title.x=element_text(family = "serif", colour = "white"),
        axis.title.y=element_text(family = "serif", angle = 90, colour = "black", face = "plain"),
        axis.text.y = element_text(family = "serif", colour = "black", face = "plain"),
        axis.text.x = element_text(family = "serif", colour = "black", face = "plain"),
        plot.title = element_text(hjust = 0.5, family="serif", colour = "black", face = "plain"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  labs(x = "Estimate", y = "Parameter") +
  force_panelsizes(rows = unit(4, "in"),
                   cols = unit(1.75, "in"))
Thirdplotsumms1

Thirdplotsumms2 <- Thirdplotsumms + xlim(-3.5,1.0) + ggtitle("Den sites") + 
  theme(legend.position = "none", 
        axis.title.x=element_text(family = "serif", colour = "black", face = "plain"),
        axis.title.y=element_text(family = "serif", angle = 90, colour = "black", face = "plain"),
        axis.text.y = element_text(family = "serif", colour = "black", face = "plain"),
        axis.text.x = element_text(family = "serif", colour = "black", face = "plain"),
        plot.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  labs(x = "Estimate", y = "Parameter") +
  force_panelsizes(rows = unit(4, "in"),
                   cols = unit(1.75, "in"))
Thirdplotsumms2
ggsave(Thirdplotsumms2, file = "thirdfigurejustsites.png", dpi = 300, height = 5, width = 4)

####################################################################
##Plot showing Fourth Order covariate selection (Figure 3 in text)
####################################################################
fourthsumms <- plot_summs(finalfourthmod2, finalfourthmod3, finalfourthmod1,
                          scale = TRUE,
                          coefs = c("Hiding\ncover (%)" = "avgcc",
                                    "Slope (\u00B0)" = "avg_slope",
                                    "East\nindex" = "EAST",
                                    "Herb\ncover (%)" = "herb",
                                    "Shrub\ncover (%)" = "shrub"),
                          point.shape =  FALSE,
                          colors = c("dodgerblue1", "dodgerblue4", "dodgerblue4"))
fourthsumms
Fourthsumms1 <- fourthsumms + 
  xlim(-1,2.25) + 
  theme(axis.ticks.y = element_blank(), 
        axis.text.y = element_text(family = "serif", colour = "black", face = "plain"),
        axis.text.x = element_text(family = "serif", colour = "black", face = "plain"),
        axis.title.x=element_text(family = "serif", colour = "black", face = "plain"),
        axis.title.y=element_text(family = "serif", colour = "black", angle = 90, face = "plain"),
        legend.position = "none", 
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  labs(x = "Estimate", y = "Parameter")
Fourthsumms1
ggsave(Fourthsumms1, file = "fourthorderfigure.png", dpi = 300, height = 3.0, width = 3.5)


####################################################################
##Plot showing interactions between den proximity, season, and conflict (Figure 4 in text)
####################################################################
plot1 <- interact_plot(conf.model1, pred = DenDecay, modx = Season3, interval = TRUE,
                       int.type = "confidence", int.width = .95)
plot1
Fplot1 <- plot1 +
  ggtitle("(A)") + 
  theme(axis.text.y = element_text(family = "serif", colour = "black", face = "plain"),
        axis.text.x = element_text(family = "serif", colour = "black", face = "plain"),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        element_blank(),
        plot.margin = unit(c(-0.25,0,0,0), "in"),
        legend.position = "none", 
        plot.title = element_text(hjust = 0.075, vjust = -12, family="serif", colour = "black", face = "plain"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  #force_panelsizes(rows = unit(1.5, "in"), cols = unit(1.5, "in")) +
  labs(x = "D", y = "L") +
  xlim(0,0.75)

Fplot1 #Final plot for first model

plot2 <- interact_plot(conf.model2, pred = DenDecay, modx = Season3, interval = TRUE,
                       int.type = "confidence", int.width = .95)
Fplot2 <- plot2 +
  ggtitle("(B)") + 
  theme(axis.text.y = element_text(family = "serif", colour = "black"),
        axis.text.x = element_text(family = "serif", colour = "black"),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        plot.margin = unit(c(-0.25,0,0,0), "in"),
        legend.position = "none", 
        plot.title = element_text(hjust = 0.075, vjust = -12, family="serif", colour = "black", face = "plain"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  #force_panelsizes(rows = unit(1.5, "in"), cols = unit(1.5, "in")) +
  labs(x = "D", y = "Likelihood of conflict encounter") +
  xlim(0,0.75)
Fplot2 #final plot for second model

plot3 <- interact_plot(conf.model3, pred = DenDecay, modx = Season3, interval = TRUE,
                       int.type = "confidence", int.width = .95)
Fplot3 <- plot3 +
  ggtitle("(C)") + 
  theme(axis.text.y = element_text(family = "serif", colour = "black"),
        axis.text.x = element_text(family = "serif", colour = "black"),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position = "none",
        plot.margin = unit(c(-0.25,0,0,0), "in"),
        legend.title = element_text(family = "serif", colour = "black"),
        legend.text = element_text(family = "serif", colour = "black"),
        plot.title = element_text(hjust = 0.075, vjust = -12, family="serif", colour = "black", face = "plain"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  #force_panelsizes(rows = unit(1.5, "in"), cols = unit(1.5, "in")) +
  labs(x = "Den decay distance", y = "Likelihood of conflict encounter") +
  xlim(0,0.75)

Fplot3 #Final plot for third model

plot4 <- interact_plot(conf.model4, pred = DenDecay, modx = LandCover, interval = TRUE,
                       int.type = "confidence", width = 0.95, colors = c("#00b159","#ffc425"))
Fplot4 <- plot4 +
  ggtitle("(D)") + 
  theme(axis.text.y = element_text(family = "serif", colour = "black"),
        axis.text.x = element_text(family = "serif", colour = "black"),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position = "none", 
        plot.margin = unit(c(-0.25,0,0,0), "in"),
        plot.title = element_text(hjust = 0.075, vjust = -12, family="serif", colour = "black", face = "plain"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  #force_panelsizes(rows = unit(1.5, "in"), cols = unit(1.5, "in")) +
  labs(x = "Den decay distance", y = "Likelihood of conflict encounter") +
  xlim(0,0.75)
Fplot4 #Final plot for fourth model

#Arrange plots
q = list(Fplot1, Fplot2, Fplot3, Fplot4) %>% map(~.x + labs(x=NULL, y=NULL))

#Add common axis labels
yleft <- textGrob(expression(paste("Likelihood of conflict encounter")), 
                  rot = 90, gp = gpar(fontfamily = "serif"))
bottom <- textGrob("Den distance decay", gp = gpar(fontfamily = "serif"))

conflictplot <- grid.arrange(grobs=q, ncol = 2, nrow = 2, 
                             left = yleft, bottom = bottom)
conflictplot

ggsave(conflictplot, file = "conflictfigure.png", dpi = 300, height = 5, width = 5)
