# COMMA ![my workflow](https://github.com/kimberlywebb/COMMA/actions/workflows/r.yml/badge.svg)

![ ](https://github.com/kimberlywebb/COMBO/blob/main/COMMA_hex_sticker.png?raw=true)

**COMMA:** **CO**rrecting **M**isclassified **M**ediation **A**nalysis

Overview
--------------------------------------------------

**COMMA** provides a set of functions for the analysis of mediation models with binary mediator misclassification. 

The two main parts are:

- Classification probability calculations
- Parameter estimation 


Classification probability calculations
--------------------------------------------------
The package allows users to compute the probability of the latent true mediators and the conditional probability of observing a mediator value given the latent true mediator, based on parameters estimated from the `COMMA_EM` and `COMMA_PVW` functions.


Parameter estimation 
--------------------------------------------------
Jointly estimate parameters from the true outcome and observation mechanisms, respectively, in a binary outcome misclassification model using the EM algorithm or MCMC. Parameters from the true outcome, first-stage observation, and second-stage observation mechanisms in a two-stage binary outcome misclassification model can also be estimated using the EM algorithm and MCMC.

Installation
--------------------------------------------------

``` r
# Install the development version from GitHub:
# install.packages("devtools")
devtools::install_github("kimberlywebb/COMMA")
```
