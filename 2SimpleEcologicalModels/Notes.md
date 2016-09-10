# Session 2: Simple Theoretical Ecological Models

## Lotka-Volterra Competition Model

### Simple Logistic Growth

First, let's simulate simple logistic growth. Remember from the <a href="https://github.com/timnewbold/MResEcologicalModelling/blob/master/2SimpleEcologicalModels/Lecture2_SimpleEcologicalModels.pdf">lecture</a> that:

dN/dt = rN((K-N)/K)

To solve this in R, we will use the deSolve package:

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

# Set the starting conditions and parameters
state <- c(N = 1)
pars <- c(r = 0.02, K = 1000)
times <- seq(from = 0, to = 1000, by = 1)

# Run and plot the model
outLG <- ode(state, times, lgModel, pars)
plot(outLG)
```

Try adjusting the parameters to explore the model behaviour.

### Adding competition

Remember from the <a href="https://github.com/timnewbold/MResEcologicalModelling/blob/master/2SimpleEcologicalModels/Lecture2_SimpleEcologicalModels.pdf">lecture</a> that the Lotka-Volterra model is simply a logistic growth model for two species, with extra terms to represent the competitive effect that each species has on the other:

dN<sub>1</sub>/dt = r<sub>1</sub>N<sub>1</sub>((K<sub>1</sub>-N<sub>1</sub>-&#945;<sub>12</sub>N<sub>2</sub>)/K<sub>1</sub>)

and:

dN<sub>2</sub>/dt = r<sub>2</sub>N<sub>2</sub>((K<sub>2</sub>-N<sub>2</sub>-&#945;<sub>21</sub>N<sub>1</sub>)/K<sub>2</sub>)

First, let's simulate a case where species 1 competitively excludes species 2, i.e. where K<sub>1</sub>/&#945;<sub>12</sub> > K<sub>2</sub> and K<sub>1</sub> > K<sub>2</sub>/&#945;<sub>21</sub>:

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

Since the isoclines don't cross and the line for species 1 is always above and to the right of the line for species 2, we know that species 1 should always competitively exclude species 2. But we will simulate this just to check:

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

# Set the starting conditions and parameters
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

Now let's explore the case where the species co-exist stably, i.e. where K<sub>1</sub>/&#945;<sub>12</sub> > K<sub>2</sub> and K<sub>2</sub>/&#945;<sub>21</sub> > K<sub>1</sub>:

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

## Lotka-Volterra Predation Model

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

As before, we can plot both species on the same graph:

```R
with(data.frame(outLVP),{
	plot(time,N,type="l",log="y",ylim=c(10,10000),col="blue")
	points(time,P,type="l",col="red")
})

```

Now try varying some of the parameters to explore the behaviour of the model.

## Rosenzweig-MacArthur predation and competition model

The Lotka-Volterra predator-prey model assumes that, in the absence of predators, prey populations grow exponentially. However, competition among species mean that this is unrealistic. The Rosenzweig-MacArthur model includes both competitive and consumption effects of species on each other.
