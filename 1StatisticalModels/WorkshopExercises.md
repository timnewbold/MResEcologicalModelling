# Session 1: Statistical Models for Complex Ecological Data

## Introduction

This session will introduce some of the statistical approaches that are available for dealing with the more complex data typical of many ecological studies, as presented in the <a href="https://github.com/timnewbold/MResEcologicalModelling/blob/master/1StatisticalModels/Lecture1ApproachesStatisticalModelling.pdf">lecture</a> this morning. I assume experience of running basic statistical methods in R, but do let me know if you get stuck with anything.

Before you start you need to have R installed on your computer. You also need to install the version of RTools that is appropriate for your version of R; you can find the installers <a href="https://cran.r-project.org/bin/windows/Rtools">here</a>.

Then you will need to install a few R packages:

```R
install.packages("lme4")
install.packages("arm")
install.packages("emdbook")
install.packages("bbmle")
install.packages("R2WinBUGS")
# Alternatively, if you are using a Mac, install this package instead
install.packages("rjags")
# If the next line doesn't work, come to see me for a version
# of the StatisticalModels and MResModelling packages
install.packages("devtools")
library(devtools)
install.packages("survival")
install.packages("Formula")
install.packages("ggplot2")
install.packages("Hmisc")
install_github("timnewbold/StatisticalModels")
install_github("timnewbold/MResEcologicalModelling",subdir="MResModelling")
```

## Exercise 1: Functional responses - Linear Models and Maximum Likelihood Estimation

In this exercise, we will use data on predator functional responses (remember from the <a href="https://github.com/timnewbold/MResEcologicalModelling/blob/master/1StatisticalModels/Lecture1ApproachesStatisticalModelling.pdf">lecture</a> that these describe the relationship between prey density and number of prey eaten) of East African reed frogs, from Vonesh & Bolker (2005).

```R
library(emdbook)
data(ReedfrogFuncresp)
```

First, let's plot the functional response to inspect the data:

```R
plot(ReedfrogFuncresp$Initial,ReedfrogFuncresp$Killed)
```

It looks as though there might be a linear relationship between these variables. So let's try fitting a simple linear model:

```R
m1 <- lm(Killed~Initial,data=ReedfrogFuncresp)

m1
# Coefficients:
# (Intercept)      Initial  
#       2.727        0.276

summary(m1)
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  2.72651    1.96292   1.389    0.187    
# Initial      0.27603    0.03948   6.991 6.33e-06 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 5.04 on 14 degrees of freedom
# Multiple R-squared:  0.7773,    Adjusted R-squared:  0.7614 
# F-statistic: 48.88 on 1 and 14 DF,  p-value: 6.334e-06
```

This linear model seems to fit the data very well. However, as you heard in the <a href="https://github.com/timnewbold/MResEcologicalModelling/blob/master/1StatisticalModels/Lecture1ApproachesStatisticalModelling.pdf">lecture today</a>, and will hear much more about in the <a href="https://github.com/timnewbold/MResEcologicalModelling/blob/master/2SimpleEcologicalModels/Lecture2SimpleTheoreticalModels.pdf">lecture tomorrow</a>, there are different theoretical models describing functional responses. One of these (the Type II) functional response is a saturating function:

N<sub>killed</sub> = (aN)/(1+aHN),

where N is initial density, and a and H are the parameters 'attack rate' and 'handling time'.

It is not easy to transform the variables to make a linear relationship from this model, and even if we did do this the parameters would be difficult to interpret. So instead we will fit the model using maximum likelihood estimation. First though, as an experiment, let's refit the linear model using maximum likelihood.

```R
# We need to load the bblme package
library(bbmle)

# First we define the likelihood function
linearNLL <- function(killed,init,m,c,sd){
  
  # The following line calculates the expected y values
  # for a given slope and intercept
  y.pred <- m * init + c
  
  # The next line calculates the likelihood for this model
  # assuming that the residuals are a normal distribution
  # around 0 (as in a linear regression)
  suppressWarnings(-sum(dnorm(x = killed,mean = y.pred,sd = sd,log = TRUE)))
  
}

# Now we run the maximum likelihood estimation
m2 <- mle2(minuslogl = linearNLL,start = list(m=0.2,c=2,sd=1),
           data = list(init=ReedfrogFuncresp$Initial,
                       killed=ReedfrogFuncresp$Killed))

m2
# Coefficients:
#         m         c        sd 
# 0.2760248 2.7265615 4.7141473 
```

The parameters are of course very similar to those estimated by the linear regression. Before moving on to the more complex functional response model, let's plot the data and the fitted linear relationship:

```R
plot(ReedfrogFuncresp$Initial,ReedfrogFuncresp$Killed)

# Create a data frame with values of initial density to predict
preds <- data.frame(Initial=1:100)
# Predict number of prey killed
preds$Killed <- m2@coef['m']*preds$Initial+m2@coef['c']
# Plot (type="l" plots a line, and lwd=2 makes the line thick)
points(preds$Initial,preds$Killed,type="l",lwd=2,col="red")
```

Now we will fit the more complex 'Type II' functional response model to the data.

```R
# As before, we first need to define the likelihood function
binomLL <- function(killed,init,p){
  
  # The next two lines just specify that the first parameter we pass
  # will be the attack rate, and the second the handling time
  a = p[1]
  h = p[2]
  
  # This line defines the probability of a prey individual being killed
  # Note this is the equation for the number killed above, divided by
  # the intial number of prey
  pkilled <- a/(1+a*h*init)
  
  # The likelihood in this case is a binomial probability distribution
  # with the observed number of kills from a total potential number of kills
  # equal to the initial number of prey, and probability of being killed equal
  # from above
  -sum(dbinom(x = killed,size = init,prob = pkilled,log = TRUE))
  
}
# This line just defines the names for the parameters
parnames(binomLL) <- c("a","h")

# Now run the likelihood estimation
m3 <- mle2(minuslogl = binomLL,start = c(a=0.5,h=0.0125),
           data = list(init=ReedfrogFuncresp$Initial,
                       killed=ReedfrogFuncresp$Killed))

m3
# Coefficients:
#          a          h 
# 0.52630319 0.01664362 
# 
# Log-likelihood: -46.72
```

If we compare the AIC values of these models, we can see that the Type II functional response is a slightly better fit to the data (i.e. has a lower AIC value - see the <a href="https://github.com/timnewbold/MResEcologicalModelling/blob/master/1StatisticalModels/Lecture1ApproachesStatisticalModelling.pdf">lecture</a> for a reminder about AIC values):
```R
AIC(m2,m3)
#         AIC df
# 1 101.02429  3
# 2  97.44279  2
```

Finally, we will plot the fitted Type II functional response. Note that the model as defined above predicts the probability that a single prey individual is eaten, so we need to multiply by initial prey density to get the estimate of the total number of prey eaten.

```R
plot(ReedfrogFuncresp$Initial,ReedfrogFuncresp$Killed)

preds <- data.frame(Initial=1:100)
preds$Killed <- preds$Initial *
	m3@coef['a']/(1 + m3@coef['a'] * 
		m3@coef['h'] * preds$Initial)
points(preds$Initial,preds$Killed,type="l",lwd=2,col="red")
```

## Exercise 2: Land use impacts on Hymenoptera in Hawaii - Generalized linear models

For this section, you will be working with the data on hymenopterans in different land uses in Hawaii, which I described in the <a href="https://github.com/timnewbold/MResEcologicalModelling/blob/master/1StatisticalModels/Lecture1ApproachesStatisticalModelling.pdf">lecture</a>.

These data form one of the studies in the PREDICTS database, which I also discussed earlier, and which you will be working with later. This particular study sampled hymenopteran species in 754 different locations in three different land uses (primary vegetation, secondary vegetation and pastures).

First, load my R package associated with these sessions:

```R
library(MResModelling)
```

Now, load in the Hawaii data describing the abundance of all species:

```R
data(HawaiiHymenoptera)
```

We will work with the data for just one species:

```R
hh.sp <- droplevels(hh[(hh$Taxon_name_entered=="Pimpla punicipes"),])
```

The abundance data (given in the column called 'Measurement') have their own particular problems. Perhaps try plotting a histogram of these data, and thinking about what the problems might be. For now we will analyze species presence or absence. To do so we will create a new column with values of 1 where abundance is greater than zero, and values of 0 otherwise (i.e. where the species is absent):

```R
hh.sp$PresAbs <- ifelse(hh.sp$Measurement>0,1,0)
```

Before we run any models, we will plot the data as I did in the lecture:

```R
# The first line of code just tabulates the presences and absences by land use
presabs.counts <- table(hh.sp$Occur,hh.sp$LandUse)

# The next line divides the counts through by the column sums to calculate proportions
presabs.props <- sweep(presabs.counts,2,colSums(presabs.counts),'/')

# This line reorders the table, so that presences come before absences, and the land uses are in a logical order (primary vegetation followed by secondary vegetation followed by pasture)
presabs.props <- presabs.props[c(2,1),c(1,3,2)]

# Finally, we create a barplot of the proportions
barplot(presabs.props,col=c("#1b9e77","#d95f02"),las=1,ylab="Proportion of sites")
```

It looks from this plot that this species, <i>Pimpla punicipes</i>, is less likely to occur in secondary vegetation than in primary vegetation, and very unlikely to occur in pasture. But we will run a generalized linear model to check this:

```R
# Run a model of species presence or absence as a function of land use
m1 <- glm(PresAbs~LandUse,data=hh.sp,family=binomial)

print(summary(m1))

# Also, run a null model (with just an intercept) with which to compare our land-use model
m0 <- glm(PresAbs~1,data=hh.sp,family=binomial)
```

Looking at the output of model 1 seems to confirm the pattern in presence and absence observed in the data. We can plot an error bar to show the modelled result:

```R
# First, create a set of data with the land uses we want to show
nd <- data.frame(LandUse=c("Primary","Secondary","Pasture"))
# Now predict probability of presence for these land uses, including standard error
preds <- predict(object=m1,newdata=nd,se.fit=TRUE)
# Calculate the mean predicted value for each land use, back transforming the values using the inverse of the link function (for a binomial GLM, the default link function is logit)
y <- 1/(1+exp(-(preds$fit)))
# Now calculate the upper and lower confidence limits (multiplying by 1.96 gives 95% confidence limits)
yplus <- 1/(1+exp(-(preds$fit+1.96*preds$se.fit)))
yminus <- 1/(1+exp(-(preds$fit-1.96*preds$se.fit)))
# Finally, plot the error bar with the mean and confidence limits of the model predictions for each land use
errbar(x=nd$LandUse,y=y,yplus=yplus,yminus=yminus)
```

Now do an analysis of variance to compare the two models:

```R
anova(m0,m1,test="Chi")
```

This shows that the model including land use is significantly better than the null model. So land use does have a significant effect on the presence or absence of this species at this particular study location in Hawaii. 

Finally, though, let's find out how much of the variation in species presence or absence is explained by land use as a measure of fit of the model. Remember from the <a href="https://github.com/timnewbold/MResEcologicalModelling/blob/master/1StatisticalModels/Lecture1ApproachesStatisticalModelling.pdf">lecture</a>, that explained variation in a GLM is (null deviance - residual deviance)/null deviance.

```R
with(m1,(null.deviance-deviance)/null.deviance)
# Note that residual deviance is labelled simply 'deviance' in the model object.
```

This shows that our model explains only around 5% of the variation in species presence or absence. So land use has a significant effect on the species, but clearly there are other (unmeasured) factors that determine to a large extent whether or not the species is present and is detected at each location in the dataset.

Now we will work with data on the counts of all hymenopteran species at the 754 sampled sites in Hawaii:

```R
data(HawaiiHymenopteraSites)
```

We will work with the species richness data. First, produce a boxplot of species richness as a function of land use:

```R
boxplot(Species_richness~LandUse,data=hhs)
```

This suggests that there might be differences in species richness among the different land uses, but we produce some models to test this possibility.

But first, let's inspect the distribution of the species richness data:

```R
hist(hhs$Species_richness)
```

Clearly the species richness data are not normally distributed. Because species richness values are counts, it would be reasonable to start with a GLM with a Poisson error distribution.

```R
# Model species richness as a function of land use
m1 <- glm(Species_richness~LandUse,data=hhs,family=poisson)

print(summary(m1))

# Also run a null, intercept-only model
m0 <- glm(Species_richness~1,data=hhs,family=poisson)
```

As before, we will plot an error bar to show the modelled result:

```R
# First, create a set of data with the land uses we want to show
nd <- data.frame(LandUse=c("Primary","Secondary","Pasture"))
# Now predict probability of presence for these land uses, including standard error
preds <- predict(object=m1,newdata=nd,se.fit=TRUE)
# Calculate the mean predicted value for each land use, back transforming the values using the inverse of the link function (this time, because we ran a Poisson GLM, the model used a log link function)
y <- exp(preds$fit)
# Now calculate the upper and lower confidence limits (multiplying by 1.96 gives 95% confidence limits)
yplus <- exp(preds$fit+1.96*preds$se.fit)
yminus <- exp(preds$fit-1.96*preds$se.fit)
# Finally, plot the error bar with the mean and confidence limits of the model predictions for each land use
errbar(x=nd$LandUse,y=y,yplus=yplus,yminus=yminus)
```

So in this case, the model shows that species richness is, on average, slightly higher in secondary vegetation than in primary vegetation, but much lower in pasture than in either of the natural land-use tpyes.

Again, do an analysis of variance to compare the land-use model with the null model:

```R
anova(m0,m1,test="Chi")
```

As with the model of the presence or absence of the individual species, the model including land use is highly significantly better than the null model. But again, let's see how much of the variation in species richness land use explains:

```R
with(m1,(null.deviance-deviance)/null.deviance)
```

So land use explains nearly 17% of the variation in species richness. That's a lot for a single explanatory variable. 

## Exercise 3: Land use impacts on biodiversity - Mixed-effects models

In this exercise, we will use the full PREDICTS database to construct models relating the biodiversity of local ecological communities to land use. These data were described in Hudson et al. (2014), and released with Newbold et al. (2016). The data describe the total abundance and species richness sampled at over 18,000 locations ('sites'), in different land uses, around the world. The data were drawn from 573 published studies of land use impacts on biodiversity. 

Each of the underlying studies (of which the study in the last exercise was one) was conducted in a different geographical location, with different sampling protocols, and with different levels of sampling effort. It would be inappropriate to construct a simple linear or generalized linear model that ignores this hierarchical structure. Therefore, we will be using mixed-effects models.

First, load my StatisticalModels package, and read in the site-level data from Newbold et al. (2016):

```R
library(StatisticalModels)
data(PREDICTSSiteData)
```

Now, we will reorder the land uses into a somewhat intuitive order (by default they are ordered alphabetically, but it makes sense at least to have natural habitats - primary and secondary vegetation - first, and we will put urban last because it reasonable to expect this land use to have the lowest biodiversity).

```R
PREDICTSSites$LandUse <- factor(PREDICTSSites$LandUse,
                                levels=c("Primary Vegetation","Secondary Vegetation",
                                         "Plantation forest","Cropland","Pasture","Urban"))
```

The other variables that we will consider are human population density and distance to nearest road. The imported data table already contains log-transformed versions of these variables, rescaled to have values between 0 and 1: logHPD.rs and logDistRd.rs, respectively.

In this exercise, we will switch to using some functions I have written for easily creating mixed-effects models.

When we fitted mixed-effects models in Exercise 2 above, we used a very simple random-effects structure including only the study from which the data were taken. We certainly want to include study in these models too, because we expect sampled biodiversity to vary considerably among studies owing to differences in sampling methods and sampling effort (in this dataset, study identity designated as 'SS'):

```R
random1 <- GLMER(modelData = PREDICTSSites,responseVar = "LogAbund",fitFamily = "gaussian",
                 fixedStruct = "LandUse+poly(logHPD.rs,2)+poly(logDistRd.rs,2)",
                 randomStruct = "(1|SS)")
```

Note that we have included in this model land use, and also quadratic polynomial terms for both human population density and distance to nearest road. It is customary to compare different random-effects structures before comparing models with different fixed effects, and while doing so to use the most complex combination of fixed effects that will be considered (here land use plus human population density plus distance to nearest road).
 
In the PREDICTS data there is also sometimes spatial structuring of sites within studies (for example, if the authors of the original papers sampled sites arranged in spatial blocks within a landscape). Therefore, we might also consider a term to account for this structuring (there is one in the PREDICTS data: 'SSB').

If two factors that you want to include are overlapping then they are referred to as 'crossed'. For example, in the analysis of metabolic rates in Exercise 2, a species could be sampled in multiple studies, and therefore study and species identity were 'crossed' random effects. 

An alternative structure to the data is to have nested random effects. The case of spatial blocks within studies in the PREDICTS data is an example - a particular spatial block can only belong to one study. In this case, there are two ways to fit the random-effects structure. Sometimes the factors are specified in a way that doesn't account for their nestedness. For example, if there were 4 spatial blocks numbered 1 to 4 in one study, and 5 spatial blocks numbered 1 to 5 in a second study. In this case, the model has no way to know that block 1 in study 1 shouldn't be treated as being the same as block 1 in study 2. The hierarchical nature of the random effect must then be specified in the model as follows: (1|SS/SSB). Alternatively the hierarchical structure can be accounted for in the specification of the variables, for example by naming the spatial blocks 1.1, 1.2, 1.3, 1.4, 2.1, 2.2, 2.3, 2.4, and 2.5, where the first number indicates the identity of the study. In this case, the random effects can be specified either in the nested fashion (1|SS/SSB) or in the same way as crossed random effects (1|SS)+(1|SSB). This latter approach is easier to work with and is the way we will use with the PREDICTS data (where the nested random effects were specified in the data).

```R
random2 <- GLMER(modelData = PREDICTSSites,responseVar = "LogAbund",fitFamily = "gaussian",
                 fixedStruct = "LandUse+poly(logHPD.rs,2)+poly(logDistRd.rs,2)",
                 randomStruct = "(1|SS)+(1|SSB)")
```

If we compare the AIC values of these two models, we can see that the one including the effect of spatial structure of sites within studies is strongly supported over the one that includes only variation among studies:

```R
AIC(random1$model,random2$model)
#               df      AIC
# random1$model 12 34327.07
# random2$model 13 33758.41
```

If you look at the model output, in the column 'Std.Dev.' under 'Random effects:', you will see that study explained the greatest portion of the variation in abundance, but that the spatial structure of sites within studies also explained a substantial portion:


```R
random2$model
# Linear mixed model fit by REML ['lmerMod']
# Formula: LogAbund ~ LandUse + poly(logHPD.rs, 2) + poly(logDistRd.rs,  
#     2) + (1 | SS) + (1 | SSB)
#    Data: modelData
# REML criterion at convergence: 33732.41
# Random effects:
#  Groups   Name        Std.Dev.
#  SSB      (Intercept) 0.4410  
#  SS       (Intercept) 2.1926  
#  Residual             0.7766  
# Number of obs: 13197, groups:  SSB, 1531; SS, 428
```

Now we have identified our random-effects structure, we can select a combination of fixed effects that adequately describes the variation in our response variable (log-transformed total abundance in this case). One way to do this, as with simpler statistical models such as linear models, is to employ backward stepwise model selection. This entails dropping each term in turn and testing whether there is a significant reduction in the explained variation. My GLMERSelect function in the StatisticalModels package does this for you. The call for this function separates categorical fixed effects (fixedFactors) from continuous effects (fixedTerms). The fixedTerms parameter specifies that you want to start with quadratic polynomials (i.e. polynomials of order 2) of human population density and distance to nearest road. During model selection, simpler polynomial terms will be tested. The verbose=TRUE just means that the full details of the steps in the model selection will be printed on the screen.

```R
abundModelSelect <- GLMERSelect(modelData = PREDICTSSites,responseVar = "LogAbund",
                                fitFamily = "gaussian",fixedFactors = "LandUse",
                                fixedTerms = list(logHPD.rs=2,logDistRd.rs=2),
                                randomStruct = "(1|SS)+(1|SSB)",verbose = TRUE)
```

If you look at the output that was printed to screen, you will see that the effect of distance to nearest road was first simplified to a linear effect, and then dropped altogether, whereas the effect of human population density was retained as a quadratic polynomial. The effect of land use was retained. Viewing the table of statistics that is output by the function shows the same things (the row for the linear effect of human population density is blank because the quadratic effect was significantly better):

```R
abundModelSelect$stats
#                  terms        ChiSq     Df            P       dAIC
# 1              LandUse 1.113244e+02 5 , 11 2.150654e-22 101.324389
# 2    poly(logHPD.rs,2) 2.586414e+01 1 , 11 3.663110e-07  23.864142
# 3 poly(logDistRd.rs,2) 9.262984e-01 1 , 13 3.358266e-01  -1.073702
# 4    poly(logHPD.rs,1)           NA   <NA>           NA         NA
# 5 poly(logDistRd.rs,1) 5.733462e-03 1 , 12 9.396422e-01  -1.994267
```

We can plot an error bar showing the land-use effects from the final model using the PlotGLMERFactor in my StatisticalModels package:

```R
PlotGLMERFactor(model = abundModelSelect$model,data = abundModelSelect$data,
                responseVar = "Abundance",logLink = "e",catEffects = "LandUse",
                xtext.srt = 45)
```

The catEffects parameter specifies which factor(s) in the model to plot. The logLink="e" causes the response variable to be back-transformed from log<sub>e</sub> abundance to raw total abundance. The xtext.srt just rotates the x-axis labels by 45 degrees.

And we can plot the effect of human population density using PlotGLMERContinuous:

```R
PlotGLMERContinuous(model = abundModelSelect$model,data = abundModelSelect$data,
                    effects = "logHPD.rs",otherFactors = list(LandUse="Primary Vegetation"),
                    xlab = "(Log) Human population density (rescaled)",ylab = "Total abundance",
                    logLink = "e")
```

We have to specify all of the terms (factors and continuous effects) that were in the model to allow the function to plot. The 'otherFactors = list(LandUse="Primary Vegetation")' term specifies that we will plot predicted values for primary vegetation.

Now we will construct similar models for species richness. There is one extra complication introduced in modelling species richness: over-dispersion. It is common to model species richness assuming a Poisson distribution of errors. However, in a Poisson distribution, the variance of the values is equal to the mean of values. Observed species richness values commonly have a variance greater than the mean, a situation known as over-dispersion.

There are a number of solutions to over-dispersion. One is to use a more appropriate error distribution, such as the negative binomial distribution. There are some packages that can fit generalized linear mixed-effects models with a negative binomial distribution of errors. However, these packages are missing some of the functionality of the lme4 package, which we have been using so far. Another solution is to fit a random-intercept term with one level for each observation in the dataset (Rigby et al., 2008). In the case of the PREDICTS data we have been using, this would be a random intercept corresponding to site identity.

But before we get onto this, let's compare the study-only and study-plus-spatial-block random-effects structures that we considered for the models of total abundance. You will probably receive warnings about model convergence, but these aren't too serious (the convergence value is close to the accepted threshold, and the functions in the lme4 package are very cautious about convergence).

```R
randomR1 <- GLMER(modelData = PREDICTSSites,responseVar = "Species_richness",fitFamily = "poisson",
                 fixedStruct = "LandUse+poly(logHPD.rs,2)+poly(logDistRd.rs,2)",
                 randomStruct = "(1|SS)",REML = FALSE)
randomR2 <- GLMER(modelData = PREDICTSSites,responseVar = "Species_richness",fitFamily = "poisson",
                 fixedStruct = "LandUse+poly(logHPD.rs,2)+poly(logDistRd.rs,2)",
                 randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

AIC(randomR1$model,randomR2$model)
#                df      AIC
# randomR1$model 11 102591.5
# randomR2$model 12 100695.0
```

As with the models of total abundance, the random-effects structure that includes the effect of the spatial structure of sites is strongly favoured.

One way to test for over-dispersion is to compare the residual deviance to the residual degrees of freedom of a model. If the deviance is much larger than the degrees of freedom, this is an indication of over-dispersion (there other possible reasons though). There is a function in the StatisticalModels package that does this:

```R
GLMEROverdispersion(model = randomR2$model)
# $residDev
# [1] 31852.27
# 
# $residDF
# [1] 15621
# 
# $ratio
# [1] 2.039068
# 
# $P.ChiSq
# [1] 0
```

The large ratio of residual deviance to residual degrees of freedom indicates the presence of over-dispersion (confirmed by the significant chi-square test). Therefore, we will try a random-effects structure with a nested effect of site, within spatial block, within study (warning, this will take a little time to run!):

```R
randomR3 <- GLMER(modelData = PREDICTSSites,responseVar = "Species_richness",fitFamily = "poisson",
                  fixedStruct = "LandUse+poly(logHPD.rs,2)+poly(logDistRd.rs,2)",
                  randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

AIC(randomR2$model,randomR3$model)
#                df       AIC
# randomR2$model 12 100695.04
# randomR3$model 13  92211.73
```

The comparison of AIC values suggests that the random-effects structure including site is strongly favoured, and re-running the over-dispersion test shows that including an observation-level random effect has removed the over-dispersion (in fact, there is now under-dispersion):

```R
GLMEROverdispersion(model = randomR3$model)
# $residDev
# [1] 7448.711
# 
# $residDF
# [1] 15620
# 
# $ratio
# [1] 0.4768701
# 
# $P.ChiSq
# [1] 1
```

Now we have selected a random-effects structure for our species richness model, we can perform backward stepwise model selection (again, this might take a while to run!):

```R
richModelSelect <- GLMERSelect(modelData = PREDICTSSites,responseVar = "Species_richness",
                               fitFamily = "poisson",fixedFactors = "LandUse",
                               fixedTerms = list(logHPD.rs=2,logDistRd.rs=2),
                               randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",verbose = TRUE)
```

Looking at the statistics table from this model-selection routine, we can see that the quadratic polynomial for the effect of distance to nearest road was simplified to a linear effect, but that otherwise all terms were retained:

```R
richModelSelect$stats
#                  terms      ChiSq     Df            P        dAIC
# 1              LandUse 421.821312 5 , 12 5.864926e-89 411.8213121
# 2    poly(logHPD.rs,2)  24.853597 1 , 12 6.185353e-07  22.8535968
# 3 poly(logDistRd.rs,2)   1.190132 1 , 13 2.753029e-01  -0.8098678
# 4    poly(logHPD.rs,1)         NA   <NA>           NA          NA
# 5 poly(logDistRd.rs,1)   5.042986 1 , 12 2.472584e-02   3.0429855
```

Plotting the model shows the effects of land use, human population density and distance to nearest road:

```R
PlotGLMERFactor(model = richModelSelect$model,data = richModelSelect$data,
                responseVar = "Species richness",logLink = "e",catEffects = "LandUse",
                xtext.srt = 45)

PlotGLMERContinuous(model = richModelSelect$model,data = richModelSelect$data,
                    effects = c("logHPD.rs"),otherContEffects = c("logDistRd.rs"),
                    otherFactors = list(LandUse="Primary Vegetation"),
                    xlab = "(Log) Human population density (rescaled)",
                    ylab = "Species richness",logLink = "e")

PlotGLMERContinuous(model = richModelSelect$model,data = richModelSelect$data,
                    effects = c("logDistRd.rs"),otherContEffects = c("logHPD.rs"),
                    otherFactors = list(LandUse="Primary Vegetation"),
                    xlab = "(Log) Distance to road (rescaled)",
                    ylab = "Species richness",logLink = "e")
```

If you have time, you could consider experimenting with your own models. Perhaps you could include the effects of land-use intensity, or consider interactions among the explanatory variables. Let me know if you need a hand trying these things.

## Exercise 4: Metabolic Rates - Bayesian Models (If you have time, and are feeling adventurous)

For this exercise, we will be using the dataset from Hudson et al. (2014) on the field metabolic rates of birds and mammals. The data are estimates of field metabolic rate for individual birds and mammals, with associated estimates of body mass. Remember from the lecture that metabolic rates scale as a power function of body mass - much more of this tomorrow. In Hudson et al.'s data, there are often estimates for several individuals of a species, but sometimes only for one.

First, load the necessary packages and data:

```R
library(MResModelling)
data(HudsonFMR)
```

To make it easier to specify the models later, we will first create duplicates of the data columns with simpler names:

```R
HudsonFMR$mass <- HudsonFMR$M_kg
HudsonFMR$fmr <- HudsonFMR$FMR_kJ_per_day
```

The models will be run from within R, but call the WinBUGS program externally. You will need to download and install this program (instructions and the installers can be found <a href="http://www.mrc-bsu.cam.ac.uk/software/bugs/the-bugs-project-winbugs/">here</a>. The method used to install WinBUGS depends whether you have a 32-bit or 64-bit computer. If you need help, just give me a shout.

If you are running a computer with an operating system other than Windows, you will need to install JAGS instead. This software can be downloaded <a href="https://sourceforge.net/projects/mcmc-jags/">here</a>.

Before we start with the Bayesian models, let's run a linear model of (log) metabolic rate as a function of (log) body mass:

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

In order to run Bayesian models (using the Gibbs sampler), we will use the package 'R2WinBUGS':

```R
library(R2WinBUGS)
```

First, set up the x and y variables, and a variable for the number of observations:

```R
x <- log10(HudsonFMR$mass)
y <- log10(HudsonFMR$fmr)

n <- length(x)
```

We need to create a text file (.bug extension) that describes the model (likelihood function and prior probabilities for the parameters). We can do this using the 'sink' function in R. We will start with a model that has very weak prior probabilities (i.e. the data will dictate the posterior probabilities) - see the <a href="">lecture slides</a> for a reminder of prior and posterior probabilities. We will assume a normal prior distribution for the slope and intercept parameters with means of zero and very low precision (i.e. very high standard deviation). Note that the 'dnorm' function in WinBUGS, unlike its counterpart in R, uses a precision parameter (1/variance or 1/(standard deviation)<sup>2</sup>) instead of standard deviation itself. We will convert back to standard deviation in the model, for comparability with the linear model above. For the prior probability for the standard deviation parameter we will use a uniform distribution, because we need to prevent the sampler from trying values less than zero (this would lead to errors of course).

```R
sink("LinearRegressionWeakPriors.bug")
cat("
    model {
    	# Prior probabilities (note very low precision)
    	intercept ~ dnorm(0,0.001)
    	slope ~ dnorm(0,0.001)
    	sd ~ dunif(0,100)
    	# Likelihood
    	for (i in 1:n){
    		y[i] ~ dnorm(mu[i],tau)
    		mu[i] <- intercept + slope*x[i]
    	}    
    	# Calculate precision parameter from standard deviation
    	tau <- 1/(sd*sd)
    }
    ", fill=TRUE)
sink() 
```

Now we need to define the data and parameters, and set initial values for the parameters:

```R
# Define the data
mr.data <- list("x","y","n")
# Or if you are running JAGS on Mac
mr.data <- list('x'=x, 'y'=y,'n'=n)


# Define the parameters
params <- c("intercept","slope","sd")

# Set initial values for the parameters
inits <- list(list(intercept=0,slope=0,sd=1))
```

Now we can run the model. We will run just one parameter chain for now. Running multiple parameter chains (with different starting parameter values) can be useful to make sure that the optimizer doesn't get stuck in local minima in the likelihood surface. This is more likely with more complex, multi-parameter models.

```R
model2 <- bugs(data=mr.data, inits = inits, parameters = params, n.chains = 1,
               model="LinearRegressionWeakPriors.bug",
            bugs.directory = "C:/Users/ucbttne/Documents/WinBUGS14/")

# Or if you need to use JAGS on Mac
model2 <- jags.model(file="LinearRegressionWeakPriors.bug", data=mr.data, inits=inits,
                     n.chains = 1)
```

If you are running on WinBUGS, the model returns the posterior means of the sampling distributions for each parameter, the 'credible intervals' for each parameter (these are conceptually different to the 'confidence intervals' used in classical frequentist statistics), and some information about overall model fit (DIC is often used to compare Bayesian models; it is an analogue to the AIC used in standard statistical models).

```R
model2
# Inference for Bugs model at "LinearRegressionUninformativePriors.bug", fit using WinBUGS,
#  1 chains, each with 1200 iterations (first 200 discarded)
#  n.sims = 1000 iterations saved
#           mean  sd 2.5%  25%  50%  75% 97.5%
# intercept  3.0 0.0  2.9  3.0  3.0  3.0   3.0
# slope      0.7 0.0  0.6  0.6  0.7  0.7   0.7
# sd         0.2 0.0  0.2  0.2  0.2  0.2   0.3
# deviance  45.8 2.3 43.1 44.1 45.2 46.9  51.8
# 
# DIC info (using the rule, pD = Dbar-Dhat)
# pD = 2.9 and DIC = 48.7
# DIC is an estimate of expected predictive error (lower deviance is better).
```

The precision of the values printed in the overall table is not very high. If we display just the summary of the parameter values (which gives higher precision), you will see that the parameter estimates are very similar to those obtained by the linear model, above:

```R
model2$summary
#                 mean          sd      2.5%       25%     50%     75%      97.5%
# intercept  2.9599600 0.007610540  2.945000  2.955000  2.9600  2.9650  2.9740000
# slope      0.6529784 0.005656686  0.642395  0.649175  0.6531  0.6566  0.6645025
# sd         0.2458753 0.004657877  0.237000  0.242800  0.2457  0.2492  0.2553025
# deviance  45.9723500 2.537637046 43.149750 44.160000 45.2200 47.0650 52.4512494
```

If you are running on JAGS, you need to use a separate function to sample from the posterior parameter values. We will take 1000 samples from the posterior values in the MCMC chain:

```R
coef.samples <- jags.samples(model = model2,variable.names = c("intercept","slope","sd"),1000)
```

Then we need to manually calculate the posterior means, medians and credible intervals for each of the parameters:

```R
coef.estimates <- with(coef.samples,rbind(
  c(mean=mean(intercept),quantile(intercept,0.025),median=median(intercept),quantile(intercept,0.975)),
  c(mean=mean(slope),quantile(slope,0.025),median=median(slope),quantile(slope,0.975)),
  c(mean=mean(sd),quantile(sd,0.025),median=median(sd),quantile(sd,0.975))
))
row.names(coef.estimates) <- c("intercept","slope","sd")

# Now display these estimates
coef.estimates
#                mean      2.5%    median    97.5%
# intercept 2.9598977 2.9455276 2.9599261 2.974921
# slope     0.6531842 0.6415785 0.6531705 0.664649
# sd        0.2460361 0.2369341 0.2461339 0.255231
```

You will see that the parameter estimates are very similar, but not identical, to the estimates produced by BUGs.

Now let us suppose that we have information strongly suggesting a particular value of the slope parameter. For the sake of argument, let's assume that this is the value of 0.75 predicted by metabolic theory. We can refit the model, but this time inserting the value of 0.75, with a high precision, into the prior distribution for the slope parameter.

```R
sink("LinearRegressionStrongSlopePrior.bug")
cat("
    model {
        # Prior probabilities (note very low precision)
        intercept ~ dnorm(0,0.001)
		# We will use a prior probability distribution with a mean of 0.75
		# and a precision corresponding to a standard deviation of 0.005
        slope ~ dnorm(0.75,40000)
        sd ~ dunif(0,100)
        # Likelihood
        for (i in 1:n){
            y[i] ~ dnorm(mu[i],tau)
            mu[i] <- intercept + slope*x[i]
        }    
        # Calculate precision parameter from standard deviation
        tau <- 1/(sd*sd)
    }
    ", fill=TRUE)
sink() 

model3 <- bugs(data=mr.data, inits = inits, parameters = params, n.chains = 1,
               model="LinearRegressionStrongSlopePrior.bug",
            bugs.directory = "C:/Users/ucbttne/Documents/WinBUGS14/")

# Or if you need to use JAGS on Mac
model3 <- jags.model(file="LinearRegressionStrongSlopePrior.bug", data=mr.data, inits=inits,
                     n.chains = 1)
```

If you are running WinBUGS on Windows, print out the model:

```R
model3$summary
#                  mean           sd        2.5%      25%      50%      75%       97.5%
# intercept   3.0024470  0.007090661   2.9890000   2.9980   3.0030   3.0070   3.0160000
# slope       0.7089608  0.003823863   0.7015975   0.7063   0.7090   0.7115   0.7166025
# sd          0.2535648  0.004488158   0.2454975   0.2504   0.2535   0.2566   0.2624000
# deviance  140.5949000 12.702438976 117.3974990 131.5000 140.4000 148.8500 166.4199533
```

Alternatively if you are running JAGS, you will need to sample from the posterior parameter estimates manually:

```R
coef.samples <- jags.samples(model = model3,variable.names = c("intercept","slope","sd"),1000)

coef.estimates <- with(coef.samples,rbind(
  c(mean=mean(intercept),quantile(intercept,0.025),median=median(intercept),quantile(intercept,0.975)),
  c(mean=mean(slope),quantile(slope,0.025),median=median(slope),quantile(slope,0.975)),
  c(mean=mean(sd),quantile(sd,0.025),median=median(sd),quantile(sd,0.975))
))
row.names(coef.estimates) <- c("intercept","slope","sd")

# Now display these estimates
coef.estimates
#                mean      2.5%    median     97.5%
# intercept 3.0027276 2.9891472 3.0028965 3.0162356
# slope     0.7090535 0.7017819 0.7089133 0.7168327
# sd        0.2536759 0.2443110 0.2534583 0.2638492
```

You will see that the estimate of the slope parameter from this model falls somewhere between the estimate we obtained before and the mean of the prior probability distribution we used.

Our precision estimate was ridiculously high. Because the likelihood surface is quite steep, the data overwhelm the priors unless we use such a high precision on the priors. The influence of the priors would be greater if the likelihood profile were shallower, for example if the dataset being used were smaller. Nevertheless, you can hopefully see the dangers of using overly influential priors. Specifying informative priors can be useful if the study is building on previous work and thus if good quantitative prior estimates of a parameter can be used. But in extreme cases it could end up being pointless using the new dataset, and it could be very difficult to do anything but confirm our prior belief. Using prior probabilities to reflect a hunch about what a parameter's value should be is a very, very bad idea.

## References

* Hudson, L.N., Isaac, N.B.J. & Reuman, D.C. (2013). The relationship between body mass and field metabolic rate among individual birds and mammals. <i> Journal of Animal Ecology</i> <b>82</b>: 1009-1020.
* Hudson, L.N., Newbold, T., Contu, S., Hill, S.L.L., Lysenko, I., De Palma, A., Phillips, H.R.P., Senior, R.A., Bedford, F., Bennett, D., Booth, H., Choimes, S., Laginha Correia Pinto, D., Day, J., Echeverría-Londoño, S., Garon, M., Harrison, M.L.K., Ingram, D.I., Jung, M., Kemp, V., Kirkpatrick, L., Martin, C., Pan, Y., Robinson, A., White, H., [hundreds of data contributors], Collen, B., Ewers, R.M., Mace, G.M., Purves, D.W., Scharlemann, J.P.W. & Purvis, A. (2014). The PREDICTS database: a global database of how local terrestrial biodiversity responds to human impacts. <i>Ecology & Evolution</i> <b>4</b>: 4701–4735.
* Nakagawa, S. & Schielzeth, H. (2013). A general and simple method for obtaining R<sup>2</sup> from generalized linear mixed-effects models. <i>Methods in Ecology & Evolution</i> <b>4</b>: 133-142.
* Newbold, T., Hudson, L.N., Arnell, A.P., Contu, S., De Palma, A., Ferrier, S., Hill, S.L.L., Hoskins, A.J., Lysenko, I., Phillips, H.R.P., Burton, V.J., Chng, C.W.T., Emerson, S., Gao, D., Pask-Hale, G., Hutton, J., Jung, M., Sanchez-Ortiz, K., Simmons, B.I., Whitmee, S., Zhang, H., Scharlemann, J.P.W. & Purvis, A. (2016). Has land use pushed terrestrial biodiversity beyond the planetary boundary? A global assessment. <i>Science</i> <b>353</b>: 288-291.
* Rigby, R.A., Stasinopoulos, D.M. & Akantziliotou, C. (2008). A framework for modelling overdispersed count data, including the Poisson-shifted generalized inverse Gaussian distribution. <i>Computational Statistics & Data Analysis</i> <b>53</b>: 381-393.
* Vonesh, J.R. & Bolker, B.M. (2005). Compensatory larval responses shift trade-offs associated with predator-induced hatching plasticity. <i>Ecology</i> <b>86</b>: 1580-1591.

