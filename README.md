# `SpatialBinReg`: R package for spatial binary regression analysis under sparsity using flexible link functions.

## Overview

`RegBinEspacial` implements a spatial binary regression model for area data, using a class of flexible link functions that include an extra shape parameter that can conveniently adapt to the degree of asymmetry present in the data. The codes provided were written in R and also in C++ to optimize computational efficiency through the R package `Rcpp`. In this package, we illustrate our methodology by applying it to a motivating dataset on periodontal disease and also include a simulation study to investigate the robustness of our methods. We follow a Bayesian approach to perform parameter estimation and model comparison based on efficient Hamiltonian Monte Carlo samplers.

## A Hierarchical Spatial Model

Consider a spatial situation where we observe a binary response $y_{is}$ for subject $i$, at site $s$ within subject $i$. We assume that $Y_{is}\sim\mbox{Bernoulli}(p_{is})$ and for each individual $i$ the probability that $Y_{is}=1$ depends on a set of subject level covariates $x_i$ and on a neighbouring structure. The binary regression model with a spatial component is then given by

$$p_{is} = F_\lambda(x_{i}'\beta + \phi_{is}), i=1,\dots,n, ~s=1,\dots,m$$

where $F$ denotes a continuous cumulative distribution function (cdf), which can be any monotonically increasing function that maps an input in $\mathbb{R}$ onto the (0,1) interval and $F^{-1}$ is typically called a link function. $\lambda$ is the shape parameter, $x_{i}$ is the vector of covariates for subject $i$ (that do not vary across space), $\beta\in\mathbb{R}^k$ is the $k\times 1$ vector of covariate coefficients (fixed effects) and $\phi_{is}$ are spatially correlated random effects. 

We will be using the class of Power and Power Reverse link functions presented in the work of Bazán et al. (2017) to deal with unbalanced data scenarios.

The vignete is available  [here](https://github.com/alan-assuncao) (I need to change the link)

## Instalation
You can install the development version of `RegBinEspacial` from [GitHub](https://github.com/alan-assuncao/RegBinEspacial) with:

```R
devtools::install_github(repo = "https://github.com/alan-assuncao/RegBinEspacial")
```

## Example

The code provided here is provided for research purposes only. In the following, we simulate data from our model under the Cauchy Power link function and illustrate the use of our sampling method for the model parameters.

```R
rm(list=ls())# clear PC memory

###################################################################################
# Required packages
###################################################################################

library('BDgraph')	# to sample from G-wishart
library ('Rcpp')
library ('microbenchmark')
library ('rbenchmark')
library ('RcppArmadillo')

# R functions 

source("/home/alan/Documentos/TESE-ALAN/Artigo-regressao-binaria-espacial/ligacao-Cauchy-Potencia/funcoes-cauchy-potencia-R.R")
source("/home/alan/Documentos/TESE-ALAN/Artigo-regressao-binaria-espacial/ligacao-Cauchy-Potencia/funcoes-auxiliares.R")

#  Rcpp/C++ functions

sourceCpp("/home/alan/Documentos/TESE-ALAN/Artigo-regressao-binaria-espacial/ligacao-Cauchy-Potencia/hmcCpp-cauchy-potencia.cpp")

#################################################### Example of use #######################################

y=read.table("Y.txt",header=TRUE) # daddos sobre ocorrencia de perio

Y=as.matrix(y)

x=read.table("X.csv",header=TRUE,sep=',') # covariaveis

X=as.matrix(x)

pcov = length(X[1,])

m = dim(Y)[2]

n = dim(Y)[1]

W = matrix(W1,m,m) # adjacency matrix for random graphical structure

W_esparsa =W_sparsa(W) # calculating quantities for the sparse adjacency matrix W

D = diag(as.vector(W_esparsa$D_sparse)) # diagonal matrix with neighbors | random grid

rho = 0.9 # spatial correlation coefficient

S = D-rho*W #  CAR priori - 
#################################### HMC SETTINGS AND PREPARING THE SIMULATION #############################

SS = 1900       # chain size at the end of simulations
burn = 100      # burn-in
lagg =1         # lagg

SS. = (SS + burn)*lagg         # chain size
idx = seq(burn * lagg + 1, SS., by = lagg)

betainit = c(-0.3,0.3)
deltainit = -0.1

rbetadelta = 0

x1 = rnorm(n)
x2 = rnorm(n)
      
X=cbind(x1,x2)
      
pcov =length(X[1,]) # amount of covariates
      
betax=beta1*x1+beta2*x2
      
kap = m # kappa fixed for simulations
      
Omega_sim=rgwish( n = 1, adj =W, b = kap, D = S, threshold = 1e-8 ) # generating the precision matrix
      
Y=array(NA,c(n,m))
phi=array(NA,c(n,m))
      
for(k in 1:n)
{
  phi[k,] = mvrnorm(1, mu = rep(0, times = m), Sigma = solve(Omega_sim))

 p = F((betax[k]+phi[k,]),delta)
        
 Y[k,] = rbinom(m, size=1, prob=p)
}
     
para = c(deltainit,betainit)
     
map = optim(para, lpostbetadelta, gradbetadelta,phi, y = Y, X=X, control = list(fnscale = -1), method = 'BFGS', hessian = TRUE); map
      
G. = G(theta=c(round(map$par,2),phi),y=Y,X=X) # Fisher information matrix of the model
      
theta.current =c(map$par,phi)
      
D. <- length(theta.current)
      
theta      <- matrix( , SS., D.)
theta[1, ] <- theta.current  
      
######################################### INITIAL VALUES FOR PARAMETER CHAINS ######################
      
kgw = (m+n)
Sgw = S
      
Omegacpp = Omega_sim

rbetadelta = 0
      
################################################## SAMPLING PARAMETERS ##############################################################
tempo=system.time(
        for(r in 2:SS.) {
          
          # amostrando beta, xi e phi i=1,...,n
          
          B=hmcCpp(theta[r-1,],SS=1,burn=1,lag=1, epsilon=0.1, LF=22,epsphi = 0.001, Lphi = 22, M=G., y=Y,X=X, Omega=Omegacpp)
          
          theta[r,] = B$theta
          
          #  rbetaphi = rbetaphi + 2*B$"taxa-aceitacao phi"
          
          rbetadelta = rbetadelta + 2*B$"taxa-aceitacao beta-delta"
          
          phimc = matrix(theta[r,(2+pcov):(1+pcov+n*m)],nrow=n,ncol=m)  # spatial effects
          
          for (q in 1:n)
          {            
            Sgw = Sgw + phimc[q,]%*%t(phimc[q,]) # log-priori de phi  (implementing it this way is equivalent to calculating the trace of
# matrix multiplication)
          }
          
          Sgw = S+Sgw
          
          # sampling Omega
          
          Omegacpp=rgwish( n = 1, adj =W, b = kgw, D = Sgw, threshold = 1e-8 ) # generating the precision matrix
          
          print(r)
          
          Sgw = 0
        }         
)

 theta = theta[idx, ]
      
post = theta
      
deltapost = post[,1]
      
lambdapost = exp(deltapost)
      
betapost = post[,(2:(pcov+1))]
      
phipost = post[,((2+pcov):(1+pcov+n*m))]
      
#>   mean(deltapost);mean(betapost[,1]);mean(betapost[,2])
# [1] 0.672052
# [1] -0.6872023
# [1] 0.7101682
#>   exp(mean(deltapost))
# [1] 1.958251
```

## Documentation

## Reference
* Alves, J. S., Bazán, J. L. and Arellano-Valle, R. B. (2023) Flexible cloglog links for binomial regression models as
an alternative for imbalanced medical data. Biometrical Journal, 65, 2100325.
* Bazán, J., Torres-Avilés, F., Suzuki, A. K. and Louzada, F. (2017) Power and reversal power links for binary
regressions: An application for motor insurance policyholders. Applied Stochastic Models in Business and
Industry, 33, 22–34
* Mohammadi, A. and Wit, E. C. (2015) Bayesian structure learning in sparse gaussian graphical models.
* Girolami, M. and Calderhead, B. (2011) Riemann manifold Langevin and Hamiltonian Monte Carlo methods.
Journal of the Royal Statistical Society B, 73, 123–214.
