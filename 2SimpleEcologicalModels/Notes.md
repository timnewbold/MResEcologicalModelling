# Session 2: Simple Theoretical Ecological Models

## Lotka-Volterra Competition Model

First, let's simulate simple logistic growth. Remember from the lecture that:

dN/dt = rN((K-N)/K)

Therefore:

N<sub>t+1</sub> = N<sub>t</sub> + rN<sub>t</sub>((K-N<sub>t</sub>)/K)

Simulate this in R, with a starting abundance (N) of 1:

```R
numTimeSteps <- 1000
N <- numeric(numTimeSteps)
r <- 10
K <- 1000
N[1] <- 1
for (i in 2:numTimeSteps){
	N[i] <- N[i-1] + r * N[i-1] * ((K - N[i-1])/K)
}
```

Simulate a case where species 1 competitively excludes species 2, i.e. where K<sub>1</sub>/&#945;<sub>12</sub> > K<sub>2</sub> and 

```R
K1 <- 100
```