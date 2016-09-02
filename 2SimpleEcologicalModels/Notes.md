# Session 2: Simple Theoretical Ecological Models

## Lotka-Volterra Competition Model

First, let's simulate simple logistic growth. Remember from the lecture that:

dN/dt = rN((K-N)/K)

To solve this in R, we will use the deSolve package:

```R
install.packages("deSolve")
library(deSolve)
```

```R
lgModel <- function(time, state, pars){
	with(as.list(c(state, pars)), {
		dN <- r * N * ((K - N) / K)
		list(dN)
	})
}

state <- c(N = 1)
pars <- c(r = 0.02, K = 1000)
times <- seq(from = 0, to = 1000, by = 1)

out <- ode(state, times, lgModel, pars)

```






```R
numTimeSteps <- 1000
N <- numeric(numTimeSteps)
r <- 100
K <- 1000
N[1] <- 1
for (i in 2:numTimeSteps){
	N[i] <- (K * N[i-1] * exp(r)) / (K + N[i-1] * (exp(r) - 1))
}
```

Simulate a case where species 1 competitively excludes species 2, i.e. where K<sub>1</sub>/&#945;<sub>12</sub> > K<sub>2</sub> and 

```R
K1 <- 100
```