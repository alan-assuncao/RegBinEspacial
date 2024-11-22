# Bayesian inference under sparsity for spatial binary regression models.

We implement a spatial binary regression model for area data, using a class of flexible link functions that include an extra shape parameter that can be adapted according to the degree of skewness present in the data. The provided codes were written in R and also in C++ to optimize computational efficiency through the R package `Rcpp`.

## Overview
The implementations performed here derive from the results obtained in the paper, Alan S. Assunção, Ricardo S. Ehlers, Dipankar Bandyopadhyay "Bayesian inference under sparsity for spatial binary regression models". We illustrate our application to a motivating dataset on periodontal disease and also include a simulation study to investigate the robustness of our methods. We follow a Bayesian approach to perform parameter estimation and model comparison.

## A Hierarchical Spatial Model

Consider a spatial situation where we observe a binary response $y_{is}$ for subject $i$, at site $s$ within subject $i$. We assume that $Y_{is}\sim\mbox{Bernoulli}(p_{is})$ and for each individual $i$ the probability that $Y_{is}=1$ depends on a set of subject level covariates $x_i$ and on a neighbouring structure. The binary regression model with a spatial component is then given by

$$p_{is} = F_\lambda(x_{i}'\beta + \phi_{is}), i=1,\dots,n, ~s=1,\dots,m$$
$$\delta=log(\lambda) \sim U(-2,2) \mbox{ and } \beta \sim N_p(0,100\times I_p)$$
$$\phi_i \sim N_m(0,\Omega)\mbox{ } i=1,\ldots,n,$$
$$\Omega \sim \mbox{G-Wishart}_W(\kappa=m,S)$$

where $F_\lambda$ denotes a continuous cumulative distribution function (cdf), which can be any monotonically increasing function that maps an input in $\mathbb{R}$ onto the (0,1) interval and $F^{-1}$ is typically called a link function. $\lambda$ is the shape parameter, $x_{i}$ is the vector of covariates for subject $i$ (that do not vary across space), $\beta\in\mathbb{R}^k$ is the $k\times 1$ vector of covariate coefficients (fixed effects) and $\phi_{is}$ are spatially correlated random effects.

$\delta=log(\lambda) \in \mathbb{R}$, is a convenient parameterization for the shape parameter $\lambda$ to simplify Bayesian computation, as per Bazán et al. (2017), $S$ is the scale matrix and $\kappa$ the degrees of freedom. We will be using the class of Power and Power Reverse link functions presented in the work of [Bazán et al. (2017)](https://onlinelibrary.wiley.com/doi/pdf/10.1002/asmb.2215) to deal with unbalanced data scenarios.

## Main Functions

### Implementação do amostrador de Monte Carlo Hamiltoniano para $\beta$, $\delta$ e $\phi_i$ para $i=1,\ldots,n$

Para cada uma das funções de ligação adotadas nesta aplicação, há dois arquivos contendo funções para poder viabilizar a implementação dos modelos sob a respectiva função de ligação. O primeiro arquivo contém funções `R` que calcularão a matriz de informação do modelo sob a função de ligação desejada para o vetor de parâmetros $(\beta,\delta)$. A matriz de informação é necessária para a implementação do Monte Carlo Hamiloniano Riemann-Manifold, que fará a amostragem do vetor de parâmetros $(\beta,\delta)$. O segundo arquivo, possui extensão `.cpp` e possui funções que implementam, propriamente falando, os métodos de Monte Carlo Hamiltoniano para amostrar o vetor de parâmetros $(\beta,\delta)$ e os efeitos aleatórios espaciais $\phi_i$ $i=1,\ldots,n$. A descrição desses arquivos para cada função de ligação consta no quadro abaixo

Link | R file | Rcpp file
---  |---     |---
Power Cauchy          | funcoes-cauchy-potencia-R.R         | hmcCpp-cauchy-potencia.cpp
Reverse Power Cauchy  | funcoes-cauchy-reversa-potencia-R.R | hmcCpp-cauchy-reversa-potencia.cpp
Power Logistic        | funcoes-logistica-potencia-R.R      | hmcCpp-logistica-potencia.cpp
Reverse Power Logistic| funcoes-logistica-reversa-potencia-R.R |hmcCpp-logistica-reversa-potencia.cpp
Power Reverse Gumbel  | funcoes-gumbel-reversa-potencia-R.R | hmcCpp-gumbel-reversa-potencia.cpp
Reverse Power Reverse Gumbel | funcoes-gumbel-reversa-reversa-de-potencia-R.R | hmcCpp-gumbel-reversa-reversa-de-potencia.cpp
Power Normal          | funcoes-normal-potencia-R.R         | hmcCpp-normal-potencia.cpp
Reverse Power Normal  | funcoes-normal-reversa-potencia-R.R | hmcCpp-normal-reversa-potencia.cpp
Logit                 | funcoes-logito-R.R                  | hmcCpp-logito.cpp
Probit                | funcoes-probito-potencia-R.R        | hmcCpp-probito-potencia.cpp
Cloglog               | funcoes-cloglog-potencia-R.R        | hmcCpp-cloglog-potencia.cpp

Cada arquivo `R` acima é composto de quatro funções:
* `F` - implementa a função de distribuição acumulada que dá origem à respectiva função de ligação
* `lpostbetadelta` - implementa a log-posteriori do vetor de parâmetros $(\beta,\delta)$
* `gradbetadelta` - implementa o gradiente do vetor de parâmetros $(\beta,\delta)$ sob a respectiva função de ligação
* `G` - Calcula a matriz de informação do modelo sob a função de ligação especificada

Cada arquivo com funções `Rcpp` descrito no quadro é composto das três primeiras funções descritas acima, só que "traduzidas" para o `Rcpp`, junto com as demais funções a seguir:
* `lpostphi` - implementa a log-posteriori do vetor de efeitos aleatórios espaciais $\phi_i$ para $i=1,\ldots,n$
* `gradphi` - implementa o gradiente do vetor de efeitos aleatórios espaciais $\phi_i$ para $i=1,\ldots,n$
* `hmcCpp` - implementa métodos de Monte Carlo Hamiltoniano para amostrar o vetor de parâmetros $(\beta,\delta)$ e os efeitos aleatórios espaciais $\phi_i$ $i=1,\ldots,n$

### Sampling from G-Wishart

Para amostrar valores da distribuição G-Wishart, utilizaremos uma função disponível no pacote `R` `BDgraph` [Mohammadi, R., Massam, H., & Letac, G. (2021).](https://www.tandfonline.com/doi/full/10.1080/01621459.2021.1996377).

## Auxiliary functions
No arquivo **funcoes-auxiliares.R** contém duas funções:

* `W_sparsa` - Calculation of various quantities related to the adjacency matrix W. Return a list of results, namely: number of area units; number of adjacency pairs; adjacency pairs; and number of neighbors for each area unit
* `adjacency2` - A function that imports the adjacency matrix, formatting it and making it ready to be used. This function is designed to import the adjacency matrix when it has been saved in .csv format.

## Example

The code provided here is provided for research purposes only. In the following, we simulate data from our model under the Cauchy Power link function and illustrate the use of our sampling method for the model parameters. The simulated data can be found here in the simulated-data folder

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

y=read.table("Y.txt",header=TRUE) # binary responses

Y=as.matrix(y)

x=read.table("X.csv",header=TRUE,sep=',') # covariates

X=as.matrix(x)

pcov = length(X[1,])

m = dim(Y)[2]

n = dim(Y)[1]

W=adjacency2('/home/alan/Documentos/TESE-ALAN/Artigo-regressao-binaria-espacial/ligacao-Cauchy-Potencia/W_.csv',m)

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
      
kap = m # kappa fixed for simulations
      
Omega_sim=rgwish( n = 1, adj =W, b = kap, D = S, threshold = 1e-8 ) # generating the initial precision matrix

phi=mvrnorm(1, mu = rep(0, times =n), Sigma = solve(Omega_sim)) # generating the initial spatial random effects

betainit = rep(1,3) # rep(1,5)

deltainit = -0.1

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
          
          # sampling beta, xi e phi i=1,...,n
          
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
