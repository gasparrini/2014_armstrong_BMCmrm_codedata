###############################################################################
# Updated version of the R code for the analysis in:
#
#   "Conditional Poisson models: a flexible alternative to conditional logistic
#     case cross-over analysis"
#   Ben Armstrong, Antonio Gasparrini, Aurelio Tobias
#   BMC Medical Research Methodology - 2014
#   http://www.ag-myresearch.com/bmcmrm2014b.html
#
# NB: the analysis is an exercise only. In particular, there is apoor control 
#   for temperature
#
# Original code by Ben Armstrong
# Corrections by Antonio Gasparrini
# Updated: 24 November 2014 
# For any problem with this code, please contact the authors
###############################################################################

library(foreign) # ENABLES READING THE DATA FILE, WHICH IS A STATA FORMAT

data <- read.dta("londondataset2002_2006.dta")
summary(data)
# SET THE DEFAULT ACTION FOR MISSING DATA TO na.exclude
# (MISSING EXCLUDED IN ESTIMATION BUT RE-INSERTED IN PREDICTION/RESIDUALS)
options(na.action="na.exclude")

# SCALE EXPOSURE
data$ozone10 <- data$ozone/10

# GENERATE MONTH AND YEAR
data$month  <- as.factor(months(data$date))
data$year   <- as.factor(format(data$date, format="%Y") )
data$dow    <- as.factor(weekdays(data$date))
data$stratum <- as.factor(data$year:data$month:data$dow)

data <- data[order(data$date),]

# FIT A CONDITIONAL POISSON MODEL WITH A YEAR X MONTH X DOW STRATA
library(gnm)
modelcpr1 <- gnm(numdeaths ~ ozone10 + temperature, data=data, family=poisson,
  eliminate=factor(stratum))
summary(modelcpr1) 

# ALLOW FOR OVERDISPERSION
modelcpr2 <- gnm(numdeaths ~ ozone10 + temperature, data=data, 
  family=quasipoisson, eliminate=factor(stratum))
summary(modelcpr2)

# ADD BRUMBACK AUTOCORRELATION ADJUSTMENT 
library(tsModel)   # FACILITATES GETTING LAGGED VALUES'
reslag1 <- Lag(residuals(modelcpr1,type="deviance"),1)
modelcpr3 <- gnm(numdeaths ~ ozone10 +  temperature + reslag1, data=data,
  family=quasipoisson, eliminate=factor(stratum))
summary(modelcpr3)

# ALLOW FOR AUTOCORRELATION AND OVERDISPERSION
library(tsModel)   # FACILITATES GETTING LAGGED VALUES'
reslag1 <- Lag(residuals(modelcpr1,type="deviance"),1)
modelcpr4 <- gnm(numdeaths ~ ozone10 +  temperature + reslag1, data=data,
  family=quasipoisson, eliminate=factor(stratum))
summary(modelcpr4)

# ILLUSTRATION OF ALLOWING FOR VARYING RATE DENOMINATORS
#  FOR THIS WE HAVE IMAGINED AVAILABILITY OF A RELEVANT POPULATION MEASURE CHANGING
#  AT SHORT TIME SCALES
data$population <- 3000000
logpop <- log(data$population)
modelcpr5 <- gnm(numdeaths ~ ozone10 + temperature, data=data, family=poisson,
  offset=logpop,eliminate=factor(stratum))
summary(modelcpr5) 

# FURTHER CODE FOR THE UNCONDITIONAL POISSON AND CONDITIONAL LOGISTIC (CASE CROSSOVER)
#  ANALYSES REPORTED IN THE TEXT

#  FIT UNCONDITIONAL POISSON MODEL
model_upr <- glm(numdeaths ~ ozone10  + temperature + factor(stratum),
  data=data, family=poisson)
summary(model_upr)

# FIT CONDITIONAL LOGISTIC MODEL

# EXPAND THE DATA IN A CASE-CROSSOVER FORMAT (AND EXCLUDE STRATA WITH 0)
# CODE IN THE FUNCTION funccmake()
source("funccmake.R")
dataexp <- funccmake(data$stratum,data$numdeaths,
  vars=cbind(ozone10=data$ozone10,temperature=data$temperature))
dataexp <- dataexp[dataexp$weights>0,]
Xexp <- as.matrix(dataexp)[,-seq(4)]

# RUN CLR
library(survival)
timeout <- as.numeric(factor(dataexp$stratum))
timein <- timeout-0.1
model_clr <- coxph(Surv(timein,timeout,status) ~ ozone10  + temperature,
  weights=weights, dataexp)
summary(model_clr)

#
