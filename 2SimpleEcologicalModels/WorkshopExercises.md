# Session 2: Simple Theoretical Ecological Models

## Introduction

In this session, we will be developing some simple competition and predation models, and simulating them in R.

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
# Define the model
lgModel <- function(time, state, pars){
	with(as.list(c(state, pars)), {
		dN <- r * N * ((K - N) / K)
		list(dN)
	})
}

# Set the starting conditions, parameters and time steps for which we will make predictions
state <- c(N = 1)
pars <- c(r = 0.02, K = 1000)
times <- seq(from = 0, to = 1000, by = 1)

# Run and plot the model
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
K1 <- 200
a12 <- 0.25
K2 <- 80
a21 <- 0.5
```

We can plot the zero isoclines for the species under these parameters:

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

First define the model describing changes in the abundance of both predators and prey, where N and P are the numbers of prey and predators, respectively, r is the intrinsic growth rate of the prey population, a is the attack rate of the predator on the prey, e is the conversion efficiency of energy during predation, and m is predator mortality:

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
# We will simulate the dynamics of populations in two grid cells (represented by the matrices state1 and state2)
# Within each cell, we simulate 4 species
state1 <- matrix(nrow=length(times),ncol=4)
state1[1,] <- 1

state2 <- matrix(nrow=length(times),ncol=4)
state2[1,] <- 1

# We will run the model for 50,000 time steps
times <- 0:50000

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

## Additional Work

This is the end of the main exercises for the session. If you have reached this point with time to spare, you could try writing your own model. Or go back over the models from the exercises and explore the behaviour of the models under different conditions. You could try incorporating the effects of environmental change in your models, for example. If you need ideas, come and talk to me!
