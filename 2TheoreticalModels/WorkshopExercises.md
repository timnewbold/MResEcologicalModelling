# Session 2: Theoretical Ecological Models

## Introduction

In the first part of this session, we will be developing some simple competition and predation models, and simulating them in R. Then, in the second part, you will load and run the Madingley ecosystem model.

## Exercise 1: Lotka-Volterra Competition Model

In this exercise, we will be working with the simple, deterministic Lotka-Volterra competition model, simulating it using ordinary differential equations. We will explore the effect of different parameter combinations on model outcomes.

### Simple Logistic Growth

First, let's simulate simple logistic growth. Remember from the <a href="https://github.com/timnewbold/MResEcologicalModelling/blob/master/2SimpleEcologicalModels/Lecture2SimpleTheoreticalModels.pdf">lecture</a> that:

dN/dt = rN((K-N)/K)

To simulate this in R, we will use the deSolve package to solve an ordinary differential equation:

```R
install.packages("deSolve")
library(deSolve)
```

```R
# First, define the model. The parameters time, state and pars are included for use later.
# All we need to do here is define change in some state of the system, here change in the 
# abundance of the species (dN) as a function of the present value of that state (N) 
# and some other parameters, here the intrinsic growth rate (r) and the carrying capacity (K)
lgModel <- function(time, state, pars){
	with(as.list(c(state, pars)), {
		dN <- r * N * ((K - N) / K)
		list(dN)
	})
}



# Now we will set the parameters that will be fed into the function: 
# times, state values and parameters
# First, the current state of the system: the starting abundance of the species (N)
state <- c(N = 1)
# Second, the parameters that determine logistic growth: r and K
pars <- c(r = 0.02, K = 1000)
# Third, the time points for which we want to make a prediction
times <- seq(from = 0, to = 1000, by = 1)

# Now we can run and plot the model
# The 'ode' function takes the model we defined above, plus the state, time points
# and parameters, and fits the model as an ordinary differential equation
outLG <- ode(state, times, lgModel, pars)
plot(outLG)
```

Now try adjusting the parameters to explore the model behaviour.

### Adding competition

Remember from the <a href="https://github.com/timnewbold/MResEcologicalModelling/blob/master/2SimpleEcologicalModels/Lecture2SimpleTheoreticalModels.pdf">lecture</a> that the Lotka-Volterra model is simply a logistic growth model for two species, with extra terms to represent the competitive effect that each species has on the other:

dN<sub>1</sub>/dt = r<sub>1</sub>N<sub>1</sub>((K<sub>1</sub>-N<sub>1</sub>-&#945;<sub>12</sub>N<sub>2</sub>)/K<sub>1</sub>)

and:

dN<sub>2</sub>/dt = r<sub>2</sub>N<sub>2</sub>((K<sub>2</sub>-N<sub>2</sub>-&#945;<sub>21</sub>N<sub>1</sub>)/K<sub>2</sub>)

First, let's simulate a case where species 1 competitively excludes species 2, i.e. where K<sub>1</sub>/&#945;<sub>12</sub> > K<sub>2</sub> and K<sub>1</sub> > K<sub>2</sub>/&#945;<sub>21</sub> (see the <a href="https://github.com/timnewbold/MResEcologicalModelling/blob/master/2SimpleEcologicalModels/Lecture2SimpleTheoreticalModels.pdf">lecture slides</a>):

```R
# Set the carrying capacities, and the coefficients describing the strength of the competitive effect of each species on the other
K1 <- 200
a12 <- 0.25
K2 <- 80
a21 <- 0.5
```

We can plot the zero isoclines for the species under these parameters (see the <a href="https://github.com/timnewbold/MResEcologicalModelling/blob/master/2SimpleEcologicalModels/Lecture2SimpleTheoreticalModels.pdf">lecture slides</a> if you can't remember what these isoclines are):

```R
# First create a empty plot 
# (the 'expression(paste())' bit allows the sucscripts)
plot(
-9999, -9999, xlim = c(0, 200), ylim = c(0, 800), 
xlab = expression(paste("N"[1])), # (the 'expression(paste())' bit allows the sucscripts)
ylab = expression(paste("N"[2])), # (the 'expression(paste())' bit allows the sucscripts)
xaxs="i",yaxs="i" # These parameters remove white space at the ends of the axes
)

# Now plot the zero iscoline for species 1
segments(0, K1/a12, K1, 0, col = "red", lwd = 2)
# And for species 2
segments(0, K2, K2/a21, 0, col = "blue", lwd = 2)
```

Since the isoclines don't cross and the line for species 1 is always above and to the right of the line for species 2, we know that species 1 should always competitively exclude species 2 (see the <a href="https://github.com/timnewbold/MResEcologicalModelling/blob/master/2SimpleEcologicalModels/Lecture2SimpleTheoreticalModels.pdf">lecture slides</a> if you need a reminder). But we will simulate this just to check:

```R
# We will set the per-capita growth rates equal:
r1 = 0.02
r2 = 0.02

# Define the model
# Now there are two state values for each of the two species: dN1 and dN2
lvcModel <- function(time,state,pars){
	with(as.list(c(state,pars)), {
		dN1 <- r1 * N1 * ((K1 - N1 - a12 * N2) / K1)
		dN2 <- r2 * N2 * ((K2 - N2 - a21 * N1) / K2)
		list(c(dN1, dN2))
	})
}

# Set the starting conditions, parameters and time steps to be simulated
state <- c(N1 = 1, N2 = 1)
pars <- c(r1 = r1, r2 = r2, K1 = K1, K2 = K2, a12 = a12, a21 = a21)
times <- seq(from = 0, to = 1000, by = 1)

# Run and plot the model
outLVC <- ode(state, times, lvcModel, pars)
plot(outLVC)

# By default the two species are plotted side-by-side
# We can also plot them on the same plot manually
dev.off() # This resets the graphics device
dfLVC <- data.frame(outLVC)
plot(dfLVC$time, dfLVC$N1, type = "l", col = "red",
xlab = "Time", ylab = "N", ylim = c(0,max(c(dfLVC$N1,dfLVC$N2))))
points(dfLVC$time, dfLVC$N2, type = "l", col = "blue")
```

Try varying the starting conditions (but not the parameters). Species 1 should always outcompete species 2.

Now let's explore the case where the species co-exist stably, i.e. where K<sub>1</sub>/&#945;<sub>12</sub> > K<sub>2</sub> and K<sub>2</sub>/&#945;<sub>21</sub> > K<sub>1</sub> (again, see the <a href="https://github.com/timnewbold/MResEcologicalModelling/blob/master/2SimpleEcologicalModels/Lecture2SimpleTheoreticalModels.pdf">lecture slides</a> for an explanation of this):

```R
K1 <- 200
a12 <- 0.1
K2 <- 100
a21 <- 0.1
```
Plot the isoclines for these parameters (you will need to change the axis limits a bit). The lines should cross this time, and because population growth of species 2 is favoured when species 1 is more abundant, and vice versa, we expect stable coexistence of the two species.

Now run the model again. You should find that the two species stably coexist whatever initial abundances you set.

Finally, let's simulate the case where we expect unstable coexistence, or different outcomes depending on the starting conditions, i.e. where K<sub>2</sub> > K<sub>1</sub>/&#945;<sub>12</sub> and K<sub>1</sub> > K<sub>2</sub>/&#945;<sub>21</sub>:

```R
K1 <- 200
a12 <- 1.5
K2 <- 150
a21 <- 0.9
```

When you are simulating this model, you will probably need to increase the number of time steps in order to reach an equilibrium state. By varying the starting abundances of the species, you should be able to achieve cases where each species prevails.

## Exercise 2: Lotka-Volterra Predation Model

In this exercise, we will simulate predator-prey interactions using the deterministic Lotka-Volterra equations. As in Exercise 1, we will solve the model using ordinary differential equations.

First define the model describing changes in the abundance of both predators and prey, where N and P are the numbers of prey and predators, respectively, r is the intrinsic growth rate of the prey population, a is the attack rate of the predator on the prey, e is the conversion efficiency of energy during predation, and m is predator mortality (see the <a href="https://github.com/timnewbold/MResEcologicalModelling/blob/master/2SimpleEcologicalModels/Lecture2SimpleTheoreticalModels.pdf">lecture slides</a> for a remdinder of the equations for this model):

```R
lvpModel <- function(time,state,pars){
	with(as.list(c(state,pars)), {
		dN <- r * N - a * N * P
		dP <- e * a * N * P - m * P
		list(c(dN, dP))
	})
}
```
Now run the model:

```R
# Set the starting conditions and parameters
state <- c(N = 1000, P = 800)
pars <- c(r = 0.01,a=0.00005,e=0.1,m=0.01)
times <- seq(from = 0, to = 10000, by = 1)

# Run and plot the model
outLVP <- ode(state, times, lvpModel, pars)
plot(outLVP)
```

As before, we can plot both species on the same graph.

```R
dev.off()
with(data.frame(outLVP),{
	plot(time,N,type="l",log="y",ylim=c(10,10000),col="blue")
	points(time,P,type="l",col="red")
})
```

You can see that prey are always more abundant than predators (the conversion efficiency less than 1 means that energy is lost between trophic levels). Predators and prey undergo cycles of abundance, out of sync with one another (the classic predator-prey cycles).

Now try varying some of the parameters to explore the behaviour of the model.

## Exercise 3: Lotka-Volterra Competition Model With Stochasticity and Dispersal

Stochastic events and dispersal can be important in shaping the outcomes of interactions among species. We can include these processes into our simple theoretical models.

Let's start with the generalized form of the Lotka-Volterra competition model (based on exponential growth): dN<sub>i</sub>/dt = rN<sub>i</sub> - N&#931;<sub>j&#8712;J</sub>&#945;<sub>ij</sub>N.

We will switch to using a simulation-based approach, rather than ordinary differential equations so that we can incorporate the stochastic effects, i.e.: N<sub>t</sub> = N<sub>t-1</sub> + rN<sub>t-1</sub> - N<sub>t-1</sub>&#931;&#945;N<sub>t-1</sub>.

```R
# We will run the model for 50,000 time steps
times <- 0:50000

# We will simulate the dynamics of populations in two grid cells (represented by the matrices state1 and state2)
# Within each cell, we simulate 4 species
state1 <- matrix(nrow=length(times),ncol=4)
state1[1,] <- 1

state2 <- matrix(nrow=length(times),ncol=4)
state2[1,] <- 1

# Set the parameter for intrinsic growth rate (assumed to be the same for all species)
pars <- list(r=0.001)

# Now create parameters for the competition coefficients
pars$alpha <- matrix(0,nrow=4,ncol=4)
# Intraspecific competition (on the diagonal of the competition coefficients matrix)
# will be equal for all species
diag(pars$alpha) <- 1e-5
# Start with equal interspecific competition for all species, weaker than intraspecific competition
pars$alpha[lower.tri(pars$alpha)] <- 1e-7
pars$alpha[upper.tri(pars$alpha)] <- 1e-7
# Finally, make species 4 a superior competitor (i.e. having a stronger interspecific effect on 
# the other species
pars$alpha[1:3,4] <- 1e-5

# Now, run the model
for (t in 2:length(times)){
  for (s in 1:dim(state1)[2]){
    # The state in cell 1 in the current time step is the state in the previous time step plus...
    state1[t,s] <- state1[(t-1),s] + 
      # ... a component describing exponential growth in the absence of interactions, less... 
      pars$r * state1[(t-1),s] - 
      # ... the effects of competitive interactions
      state1[(t-1),s] * 
      # The next line describes the sum of the interaction coefficients of each other species on 
      # the focal species multiplied by the abundance of the other species
      sum(pars$alpha[s,1:dim(state1)[2]] * state1[(t-1),])
    # Repeat for cell 2
    state2[t,s] <- state2[(t-1),s] + pars$r * state2[(t-1),s] - state2[(t-1),s] * 
      sum(pars$alpha[s,1:dim(state2)[2]] * state2[(t-1),])
  }
  
}
```

Plotting this model shows that species 4 outcompetes the other 3 species:

```R
# First specify a 2-by-2 grid of plots
par(mfrow=c(2,2))

# Now plot the population dynamics for each species (we will plot the results for just one cell)
for (sp in 1:4){
	plot(state1[,sp],type="l",main=paste("Species ",sp,sep=""),xlab="Time",ylab="N",ylim=c(0,100))
}
```

Now we will introduce some stochasticity into the model, by assuming that a disturbance process occurs in each cell separately, in each time step, with a probability of 0.001. In the event of disturbance, the abundance of all species is reduced by 50%.

```R
# Re-run the model, but this time adding a random probability of disturbance for each cell, 
# and reducing abundance by 50% in the event of disturbance
for (t in 2:length(times)){
  for (s in 1:dim(state1)[2]){
    # This is the calculation of population growth and competition in cells 1 and 2, as before
    state1[t,s] <- state1[(t-1),s] + pars$r * state1[(t-1),s] - state1[(t-1),s] * 
      sum(pars$alpha[s,1:dim(state1)[2]] * state1[(t-1),])
    state2[t,s] <- state2[(t-1),s] + pars$r * state2[(t-1),s] - state2[(t-1),s] * 
      sum(pars$alpha[s,1:dim(state2)[2]] * state2[(t-1),])
  }
  
  # runif(1) draws a single number at random from a uniform distribution between 0 and 1
  # It is useful for simulating processes that have a certain probability of occurring
  # Here, we are simulating a random disturbance event that has a 0.1% chance of occurring,
  # and reducing all abundances by 50% if it does
  if (runif(1)<0.001){
    state1[t,] <- state1[t,]*0.5
  }
  if (runif(1)<0.001){
    state2[t,] <- state2[t,]*0.5
  }
  
}
```

Plotting this model shows less dominance by species 4 (the other species persist this time), although it still achieves greater abundance than the other species. And of course there is now random variation in abundances that wasn't seen before:

```R
par(mfrow=c(2,2))

for (sp in 1:4){
    plot(state1[,sp],type="l",main=paste("Species ",sp,sep=""),xlab="Time",ylab="N",ylim=c(0,100))
}
```

Now we will introduce dispersal.

```R
# First, we need a new set of parameters to describe dispersal rates
# We will assume a competition-colonization trade-off (i.e. that species 4 has a lower dispersal rate)
pars$d <- c(0.05,0.05,0.05,0.00001)

# Now re-run the model, but with dispersal of individuals between the two cells
for (t in 2:length(times)){
  for (s in 1:dim(state1)[2]){
    # This is the calculation of population growth and competition in cells 1 and 2, as before
    state1[t,s] <- state1[(t-1),s] + pars$r * state1[(t-1),s] - state1[(t-1),s] * 
      sum(pars$alpha[s,1:dim(state1)[2]] * state1[(t-1),])
    state2[t,s] <- state2[(t-1),s] + pars$r * state2[(t-1),s] - state2[(t-1),s] * 
      sum(pars$alpha[s,1:dim(state2)[2]] * state2[(t-1),])
    
    # First calculate the number of individuals that will emigrate from each cell
    cell1.emig <- pars$d[s] * state1[t,s]
    cell2.emig <- pars$d[s] * state2[t,s]
    
    # Now move the emigrating individuals between cells
    state1[t,s] <- state1[t,s] - cell1.emig
    state1[t,s] <- state1[t,s] + cell2.emig
    
    state2[t,s] <- state2[t,s] - cell2.emig
    state2[t,s] <- state2[t,s] + cell1.emig
    
  }
  
  # Apply disturbance, as before
  if (runif(1)<0.001){
    state1[t,] <- state1[t,]*0.5
  }
  if (runif(1)<0.001){
    state2[t,] <- state2[t,]*0.5
  }
}
```

Plotting the results from this model shows that species 4 is much less dominant, and species 1 to 3 maintain abundances that are nearly has high as species 4's abundance.

```R
par(mfrow=c(2,2))

for (sp in 1:4){
    plot(state1[,sp],type="l",main=paste("Species ",sp,sep=""),xlab="Time",ylab="N",ylim=c(0,100))
}
```

There are much more sophisticated and better models than this for simulating the effects of stochasticity and dispersal, but hopefully this has illustrated the effects that these processes can have and the approach that one might take to incorporate them into simple ecological models.

## Exercise 4: The Madingley Model

This session will be based entirely around the Madinlgey general ecosystem model (Harfoot et al., 2014). Of course, as you saw in the <a href="https://github.com/timnewbold/MResEcologicalModelling/blob/master/3ComplexEcologicalModels/Lecture3_ComplexEcologicalModels.pdf">lecture</a>, there are lots of other complex systems-level ecological models that could be used. But I am most familiar with the Madingley Model! 

The Madingley Model will only work on Windows machines at the moment, so get together in groups around a laptop with Windows (don't worry, there will be nothing on this in the exam)

Using the MadingleyR R package, you will download the model, run and plot basic simulations of the model, try altering the parameters of the model, and run scenarios of human impacts in the model.

### Setting up the model

First, you need to install a package for viewing the output files from the model. You can download the installer (Windows only unfortunately!) <a href="https://dl.dropboxusercontent.com/u/23999497/ScientificDataSet.msi">here</a>. To install, double-click on the downloaded file and follow the instructions.

You need to install an R package that works with NetCDF files (the format of many of the output files of the Madingley Model). Because the installer for this package is available only from my Dropbox, the installation process is a bit complicated.

First, you need to create a directory for downloading the installer, for example:

```R
installationDir <- "C:/Users/ucbttne/Documents/"
```

Now download the installer to this directory:

```R
download.file(url = "https://www.dropbox.com/s/nimp6yz1iifofz5/sds_0.1-5.zip?dl=0",
              destfile = paste(installationDir,"sds_0.1-5.zip",sep=""),method="wget")
```

Now install the package:

```R
install.packages(paste(installationDir,"sds_0.1-5.zip",sep=""), 
                 repos = NULL, type = "win.binary")
```

Now you need to install the R package that interfaces with the Madingley Model, the R package for plotting Madingley simulations, and some other required packages.

```R
install.packages("sp")
install.packages("rgdal")
install.packages("survival")
install.packages("Formula")
install.packages("ggplot2")
install.packages("acepack")
install.packages("latticeExtra")
install.packages("gridExtra")
install.packages("htmlTable")
install.packages("htmltools")
install.packages("data.table")
library(devtools)
install_github("timnewbold/MadingleyR",subdir = "MadingleyR")
install_github("timnewbold/MadingleyPlots",subdir = "MadingleyPlots")

library(MadingleyR)
library(MadingleyPlots)
```

Now, define a directory on your computer where all of the code, data and eventual simulations will be stored, for example:

```R
codeDir <- "F:/Temp"
```

Now you can download the Madingley code and necessary data to this location (this can take a while, depending on the speed of your internet connection):

```R
DownloadRelease(codeDir = codeDir)
GetMadingleyData(codeDir = codeDir)
```

### Exercise 4a: Running some basic simulations

In this exercise, we will run and plot some basic simulations from the Madingley Model.

To begin with, let's run a simulation for a single location.

The first step is to initialize the model with the parameters we want:

```R
init <- MadingleyInitialisation(duration = 50,
                                specific_locations = list(latitudes=2,longitudes=33))
params <- MadingleyParameters()
```

The duration=50 parameter means that the model will run for 50 years, and the specific_locations parameter defines the latitude and longitude of the location we wish to simulate. The empty call to MadingleyParameters means that default values for all of the parameters will be used (see <a href="https://github.com/timnewbold/MadingleyR/blob/master/EcologicalParameters.md">here</a> for the definitions and default values of the parameters, as described in Harfoot et al., 2014).

Having initialized the model, we can now run the specified simulation:

```R
madingleySim <- RunMadingley(codeDir = codeDir,init = init,params = params)
```

Finally, we can plot the time series of simulated biomass densities across the whole duration of this simulation:

```R
PlotTemporal(resultsDir = madingleySim$outputDir,plotName = "BiomassDensity")
```

You can see that plants (shown by the dark green line) have a higher biomass density than herbivores (light green), omnivores (purple) or carnivores (red), which is expected. But the estimated biomasses from a single simulation of a model as complex as the Madingley Model can be quite variable. Therefore, it is better to run multiple simulations, and then to calculate average estimated biomasses and confidence limits.

To do this we need to re-initialize the Madingley Model, specifying that the simulations be run in parallel (to make the most of all of the cores in your computer).

```R
init <- MadingleyInitialisation(duration = 50,parallel_simulations = "yes",
                                specific_locations = list(latitudes=2,longitudes=33))
```

Now when we run the model, we specify numSims=5 to tell the model to run 5 replicate simulations.

```R
madingleySim <- RunMadingley(codeDir = codeDir,init = init,params = params,numSims = 5)
```
When we plot the model this time, we will specify plotConfidence=TRUE to show the variation in estimated biomasses around the mean:

```R
PlotTemporal(resultsDir = madingleySim$outputDir,plotName = "BiomassDensity",
             plotConfidence = TRUE)
```

As well as plotting the simulated biomass densities across the whole length of the simulation, we can also plot properties of the ecosystem at the end of the simulation. Useful properties of the ecosystem include the relative biomasses of different trophic levels, which can be plotted as a biomass pyramid, or the relationship between the body mass of organisms and their abundance density (empirical evidence leads us to expect a negative power-law relationship between mass and density). There are functions for making both of these figures in MadingleyPlots (PlotPyramids and PlotMassDensity). We can plot a biomass pyramid from the simulation we just ran:

```R
PlotPyramids(resultsDir = madingleySim$outputDir,plotName = "BiomassPyramids")
``` 

Plotting the relationship between body mass and abundance density requires running the Madingley Model with a high level of output detail. This makes the simulations rather slow, so we won't do this right now, but if you want to try running this later, you just need to add output_detail="high" to the call to the MadingleyInitialisation function.

Now, let's try running a simulation for two different locations (the original cell in Uganda and a cell with much drier climatic conditions in northern Algeria), and compare the results. We will switch back to running just a single simulation in the interests of time, but you could try running an ensemble of simulations later if you are interested.

```R
init <- MadingleyInitialisation(duration = 50,parallel_cells = "yes",
                                specific_locations = list(latitudes=c(2,35),longitudes=c(33,1)))
madingleySim <- RunMadingley(codeDir = codeDir,init = init,params = params)
```

You can plot a map of the cells simulated using the following command:

```R
PlotCells(resultsDir = madingleySim$outputDir,map = "Africa")
```

Finally, plot the time series of estimated biomasses as before.

```R
PlotTemporal(resultsDir = madingleySim$outputDir,plotName = "BiomassDensity")
```

You will see that this ecosystem is estimated to have lower biomasses overall, especially of carnivores (the relatively low abundance of carnivores in ecosystems with low plant biomass matches theoretical expectations - see e.g. Post, 2004).

### Exercise 4b: Running simulations with altered ecological parameters

So far, we have been running simulations using the default ecological parameters. But it can be interesting to assess the effect on simulated dynamics of changing certain parameters.

Remember from yesterday's lecture that functional responses for a single predator or herbivore organism are composed of two main parameters: attack rates and handling times (the Madingley Model also contains parameters describing variation in attack rates and handling times across organisms, as a function of their body mass). Let's test the effect of increasing the handling time for all predator organisms in the model (the <a href="https://github.com/timnewbold/MadingleyR/blob/master/EcologicalParameters.md">default</a> handling-time scalar is 0.5; we will increase this by a factor of 4 to 2.0):

```R
init <- MadingleyInitialisation(duration = 50,
                                specific_locations = list(latitudes=2,longitudes=33))
params <- MadingleyParameters(pred_terr_handlingscalar = 2.0)
madingleySim <- RunMadingley(codeDir = codeDir,init = init,params = params)
PlotTemporal(resultsDir = madingleySim$outputDir,plotName = "BiomassDensity")
```

Any parameters not referred to in your call to MadingleyParameters will assume their default values.

You will see that there is now a much smaller biomass of carnivores and omnivores (because they now predate much more slowly), and consequently a higher biomass of herbivores and a lower biomass of plants.

Now try instead increasing the handling times for herbivory. The default value for the terrestrial handling time scalar is 0.7, so try increasing it to 2.8.

### Exercise 4c: Running simulations under scenarios of human impact

In this final exercise, we will look at simulating scenarios of human impact using the Madingley Model.

First, reset the ecological parameters to use all default values:

```R
params <- MadingleyParameters()
```

Now, re-initialize the model. This time, we need to set a burn-in period (argument 'burnin'), after which the scenario of human impact will be applied. We will stick with a total simulation length of 50 years, and adopt a 25-year burn-in.

```R
init <- MadingleyInitialisation(duration = 50,burnin = 25,parallel_simulations = "yes",
                                specific_locations = list(latitudes=2,longitudes=33))
```

Now when we run the model, we need to specify the scenarios we want, and labels to be applied to the output files (in case you run multiple scenarios, but we will just run one for now). In this scenario, we will remove a constant 80% of net primary production (which corresponds with a reduction in plant biomass). This is specified within the scenarios argument as npp="constant 0.8 0" (don't worry about the final number &#8210; this is used for types of scenario other than the constant-rate one we will use today).

```R
madingleySim <- RunMadingley(codeDir = codeDir,init = init,params = params,
                             numSims = 5,
                             scenarios = list(label="NPPConstant",
                                              npp="constant 0.8 0",
                                              temperature="no 0.0",
                                              harvesting="no 0.0"))
```

Plotting this model you will see that after the burn-in is finished and the removal of plant biomass begins, there is (not surprisingly!) a big decline in plant biomass. This leads to a reduction in the biomasses of the higher trophic levels, because there is a smaller total energy entering the ecosystem via the plants. Interestingly, there is a bigger decline in omnivores and carnivores than in herbivores, probably because it is harder for predators to find food, making them more sensitive to loss of plant biomass.

```R
PlotTemporal(resultsDir = madingleySim$outputDir,plotName = "BiomassDensity")
```

Now let's simulate a different type of human impact &#8210; direct hunting of animals. This is specified using the 'harvesting' component of the 'scenarios' argument. In this case the value refers not to a fraction but to the total biomass harvested (in kg per km<sup>2</sup>). We will assume that humans remove 65 tonnes of animal biomass per km<sup>2</sup> (this is probably a rather extreme scenario).

```R
madingleySim <- RunMadingley(codeDir = codeDir,init = init,params = params,
                             numSims = 5,
                             scenarios = list(label="NPPConstant",
                                              npp="no 0.0 0",
                                              temperature="no 0.0",
                                              harvesting="constant 65000"))
PlotTemporal(resultsDir = madingleySim$outputDir,plotName = "BiomassDensity")
```

You will see from the plot that carnivores are disproportionately impacted. Although herbivores, omnivores and carnivores are all harvested, carnivores emerge as being the most sensitive, probably again because of an interaction with food availability.

If you have finished all of the exercises, have a go at running your own scenarios of human impact. You could try exploring different levels of npp removal or harvesting rates, you could try applying both impacts simultaneously. If you need help running these scenarios, or if you want some different ideas, come and chat with me.

## References

* Harfoot, M.B.J., Newbold, T., Tittensor, D.P., Emmott, S., Hutton, J., Lyutsarev, V., Smith, M.J., Scharlemann, J.P.W. & Purves, D.W. (2014). Emergent global patterns of ecosystem structure and function from a mechanistic general ecosystem model. <i>PLoS Biology</i> <b>12</b>: e1001841.
* Post, D.M. (2004). The long and short of food-chain length. <i>Trends in Ecology & Evolution</i> <b>17</b>: 269-277.