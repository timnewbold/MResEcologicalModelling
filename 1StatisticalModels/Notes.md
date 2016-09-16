# Session 1: Statistaical Models for Complex Ecological Data

## Introduction

This session will introduce some of the statistical approaches that are available for dealing with the more complex data typical of many ecological studies. I assume experience of running basic statistical methods in R.

## Exercise 1: Metabolic Theory - Mixed-effects models

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

model2
# Linear mixed model fit by REML ['lmerMod']
# Formula: log10(fmr) ~ log10(mass) + (1 | Study)
#    Data: HudsonFMR
# REML criterion at convergence: -1781.331
# Random effects:
#  Groups   Name        Std.Dev.
#  Study    (Intercept) 0.2288  
#  Residual             0.1141  
# Number of obs: 1498, groups:  Study, 126
# Fixed Effects:
# (Intercept)  log10(mass)  
#      2.9603       0.6806 

summary(model2)
# Fixed effects:
#             Estimate Std. Error t value
# (Intercept)  2.96033    0.02213  133.75
# log10(mass)  0.68061    0.01420   47.92
```

We can calculate pseudo-R<sup>2</sup>values, using a method proposed by Nakagawa & Schielzeth (2013)

```R
R2GLMER(model2)
# $conditional
# [1] 0.9799064
# 
# $marginal
# [1] 0.8991389
```

Accounting for the possibly different protocols used to collect the data, this model suggests a slightly steeper relationship between body mass and metabolic rate. In the information about the model that R returns, the section on random effects shows the variation in the response variable that is accounted for by variation among studies (a standard deviation of 0.2288) as well as the unexplained residual variation (standard deviation of 0.1141).

It is possible that the slope of the relationship between mass and metabolic rate also varies among studies. To test this possibility, we will now fit a model with a random slope of mass nested within study identity:

```R
model3 <- lmer(log10(fmr)~log10(mass)+(1+log10(mass)|Study),data=HudsonFMR)
# You can include additional random slopes to the left of the '|' if you have other fixed effects in your model

model3
# Linear mixed model fit by REML ['lmerMod']
# Formula: log10(fmr) ~ log10(mass) + (1 + log10(mass) | Study)
#    Data: HudsonFMR
# REML criterion at convergence: -1840.272
# Random effects:
#  Groups   Name        Std.Dev. Corr
#  Study    (Intercept) 0.2282       
#           log10(mass) 0.1366   0.43
#  Residual             0.1101       
# Number of obs: 1498, groups:  Study, 126
# Fixed Effects:
# (Intercept)  log10(mass)  
#      2.9176       0.6415  

summary(model3)
# Fixed effects:
#             Estimate Std. Error t value
# (Intercept)  2.91757    0.02646  110.26
# log10(mass)  0.64146    0.02031   31.59
```

Comparing the AIC values suggests that a model accounting for variation among studies in the slopes of the relationship gives a better fit to the data, and a slope estimate more like that obtained with a simple linear regression:

```R
AIC(model2,model3)
#        df       AIC
# model2  4 -1773.331
# model3  6 -1828.272
```

The data contain estimates for multiple species, often with multiple estimates for a single species. Metabolic rates might differ among species, which could influence our ability to detect the relationship between body mass and metabolic rate. So let's fit a model with a random intercept for species (note we have dropped the random slope here because there are insufficient data to create a reliable model that includes both random intercepts and the random slope term):

```R
# First create a column holding the species binomial
HudsonFMR$binomial <- paste(HudsonFMR$Genus,HudsonFMR$Species)

model4 <- lmer(log10(fmr)~log10(mass)+(1|Study)+(1|binomial),data=HudsonFMR)

model4
# Linear mixed model fit by REML ['lmerMod']
# Formula: log10(fmr) ~ log10(mass) + (1 | Study) + (1 | binomial)
#    Data: HudsonFMR
# REML criterion at convergence: -1903.492
# Random effects:
#  Groups   Name        Std.Dev.
#  binomial (Intercept) 0.19490 
#  Study    (Intercept) 0.09435 
#  Residual             0.10799 
# Number of obs: 1498, groups:  binomial, 133; Study, 126
# Fixed Effects:
# (Intercept)  log10(mass)  
#      2.9422       0.6587

# Note we are comparing with model here, which was the simpler (and nested) model without random slopes
AIC(model2,model4)
#        df       AIC
# model2  4 -1773.331
# model4  5 -1893.492
```

Clearly, including a random effect of species identity substantially improved the model fit. A lot of variation in the response variable (metabolic rate) is explained by variation among species. In fact a lot of the variation that we previously attributed to variation among studies is shown in this model to be attributable to variation among species (a potential difficulty for the model here is that study identity and species identity are obviously correlated; in other words, certain studies tended to focus on certain species).

Have a think about other ways that you could analyse these variables (for example, other explanatory variables that you could consider). For hints, you could have a look at the paper by Hudson et. al.

One thing to note from this exercise is that whichever model we used (even the simple linear regression), the answer that we obtained was relatively similar. This is probably because it is relatively easy to follow a standard protocol to estimate species metabolic rates and body masses, and because the variation in metabolic rates among species did not bias our estimates of the relationship between body mass and metabolic rate. This does not mean necessarily that it is appropriate to use the simple linear regression when we have reason to suspect some underlying variation or non-independence in the data. We will see an example later where failure to account for the hierarchical structure of the data could lead to totally the wrong conclusions.

## References

* Hudson, L.N., Isaac, N.B.J. & Reuman, D.C. (2013). The relationship between body mass and field metabolic rate among individual birds and mammals. <i> Journal of Animal Ecology</i> <b>82</b>: 1009-1020.
* Nakagawa, S. & Schielzeth, H. (2013). A general and simple method for obtaining R<sup>2</sup> from generalized linear mixed-effects models. <i>Methods in Ecology & Evolution</i> <b>4</b>: 133-142.

