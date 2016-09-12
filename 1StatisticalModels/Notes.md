# Session 1: Statistaical Models for Complex Ecological Data

## Introduction

This session will introduce some of the statistical approaches that are available for dealing with the more complex data typical of many ecological studies. I assume experience of running basic statistical methods in R.

## Exercise 1: Metabolic Theory

For this first section of this session, we will be using the dataset from Hudson et al. (2014) on the field metabolic rates of birds and mammals. The data are estimates of field metabolic rate for individual birds and mammals, with associated estimates of body mass. Often there estimates for several individuals of a species, but sometimes only for one.

First, load the necessary packages and data:

```R
library(MResModelling)
data(HudsonFMR)
```

To make it easier to specify the models later, we will first create duplicates of the data columns with simpler names:

```R
HudsonFMR$mass <- HudsonFMR$M_kg
HudsonFMR$fmr<- HudsonFMR$FMR_kJ_per_day
```

To begin with, fit a simple linear model to the data. We fit both variables log transformed because we expect metabolic rates to follow a power-law relationship. Theory, which we will cover in the lecture tomorrow, suggests exponents of the relationship between mass and metabolic rate of either 2/3 or 3/4.

```R
model1 <- lm(log10(fmr)~log10(mass),data=HudsonFMR)
model1
# Call:
# lm(formula = log10(fmr) ~ log10(mass), data = HudsonFMR)
# 
# Coefficients:
# (Intercept)  log10(mass)  
#      2.9596       0.6528

summary(model1)
# Call:
#   lm(formula = log10(fmr) ~ log10(mass), data = HudsonFMR)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.03553 -0.14913  0.03124  0.16596  0.63925 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 2.959603   0.007681   385.3   <2e-16 ***
#   log10(mass) 0.652813   0.005660   115.3   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2456 on 1496 degrees of freedom
# Multiple R-squared:  0.8989,	Adjusted R-squared:  0.8989 
# F-statistic: 1.33e+04 on 1 and 1496 DF,  p-value: < 2.2e-16
```

From this you can see that the model fits very well, and we get a slope of 0.65, fairly close to the 2/3 expected by one branch of theory.

However, there are features of the data that mean that this simple model might not be the most appropriate.

Because datasets are now made available alongside many published papers and because there has been a rapid increase in computing power, there has been a growth in the use of synthetic analyses that make use of collations of existing data. Such analyses have the advantage of greater statistical power brought about by the increased size of datasets, and often greater ability to generalize patterns across different contexts. However, no matter how much we attempt to find datasets that have been collected in the same way, there are likely to be some differences in protocol.

Mixed-effects models are a powerful tool for dealing with heterogeneity in response variables that is not caused by the variables that you are asking questions about but that could otherwise bias the results and conclusions, including heterogeneity caused by assembling data from lots of different studies.

There are several R packages for running mixed-effects models. We will be using lme4, which is probably the most widely used.

The specification of models is very similar to the specification of simple linear models. We just need an extra component that describes the random effects. Remember from the lecture that we can fit random intercepts, which describe variation in the overall value of the response variable among different subsets of the data, and random slopes, which describe variation in the relationship between an explanatory variable and the response variable among the same subsets.

First, we will run a simple mixed-effects model with a random intercept describing variation in recorded metabolic rates among the studies from which data were taken:

```R
model2 <- lmer(log10(fmr)~log10(mass)+(1|Study),data=HudsonFMR)
# The '(1|Study)' describes a random intercept of study
```

It is possible that the slope of the relationship between mass and metabolic rate also varies among studies. 
