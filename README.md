# `SpatialBinReg`: Pacote R para análise de regressão binária espacial sob esparsidade utilizando funções de ligação flexíveis.

## Overview

`RegBinEspacial` implements a spatial binary regression model for area data, using a class of flexible link functions that include an extra shape parameter that can conveniently adapt to the degree of asymmetry present in the data. The codes provided were written in R and also in C++ to optimize computational efficiency through the R package `Rcpp`. In this package, we illustrate our methodology by applying it to a motivating dataset on periodontal disease and also include a simulation study to investigate the robustness of our methods. We follow a Bayesian approach to perform parameter estimation and model comparison based on efficient Hamiltonian Monte Carlo samplers.

## A Hierarchical Spatial Model

Consider a spatial situation where we observe a binary response $y_{is}$ for subject $i$, at site $s$ within subject $i$. We assume that $Y_{is}\sim\mbox{Bernoulli}(p_{is})$ and for each individual $i$ the probability that $Y_{is}=1$ depends on a set of subject level covariates $x_i$ and on a neighbouring structure. The binary regression model with a spatial component is then given by

$$p_{is} = F_\lambda(x_{i}'\beta + \phi_{is}), i=1,\dots,n, ~s=1,\dots,m$$

where $F$ denotes a continuous cumulative distribution function (cdf), which can be any monotonically increasing function that maps an input in $\mathbb{R}$ onto the (0,1) interval and $F^{-1}$ is typically called a link function. $\lambda$ is the shape parameter, $x_{i}$ is the vector of covariates for subject $i$ (that do not vary across space), $\beta\in\mathbb{R}^k$ is the $k\times 1$ vector of covariate coefficients (fixed effects) and $\phi_{is}$ are spatially correlated random effects. 

We will be using the class of Power and Power Reverse link functions presented in the work of Bazán et al. (2017) to deal with unbalanced data scenarios.

The vignete is available  [here](https://github.com/alan-assuncao) (preciso mudar o link)

## Instalation
You can install the development version of `RegBinEspacial` from [GitHub](https://github.com/alan-assuncao/RegBinEspacial) with:

```R
devtools::install_github(repo = " ")
```

## Example

O código providenciado aqui está sendo fornecido apenas para fins de pesqusia

## Documentation

## Usage

## Reference
* Alves, J. S., Bazán, J. L. and Arellano-Valle, R. B. (2023) Flexible cloglog links for binomial regression models as
an alternative for imbalanced medical data. Biometrical Journal, 65, 2100325.
* Bazán, J., Torres-Avilés, F., Suzuki, A. K. and Louzada, F. (2017) Power and reversal power links for binary
regressions: An application for motor insurance policyholders. Applied Stochastic Models in Business and
Industry, 33, 22–34
* Mohammadi, A. and Wit, E. C. (2015) Bayesian structure learning in sparse gaussian graphical models.
* Girolami, M. and Calderhead, B. (2011) Riemann manifold Langevin and Hamiltonian Monte Carlo methods.
Journal of the Royal Statistical Society B, 73, 123–214.
