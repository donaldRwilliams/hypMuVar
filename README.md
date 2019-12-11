
hypMuVar: Bayesian Hypothesis Testing of Mean-Variance Relations
================================================================

``` r
library(hypMuVar)
```

### Overview

### Installation

``` r
# install.packages("devtools")
devtools::install_github("donaldRwilliams/hypMuVar")
```

### Random slopes for the location and scale:

``` r
fit <- melsm(fixed_location = rt ~  congruency,
             random_location = ~  congruency | ID,
             fixed_scale = sigma ~  congruency,
             random_scale = ~ congruency | ID,  
             mixture = "SSVS", k = 3,
             adapt = 5000, 
             iter = 10000, 
             data = stroop)
```

``` r
summary(fit)
#> hypMuVar: Bayesian Hypothesis Testing 
#>           of Mean--Variance Relations
#> Model: MELSM
#> Chains: 4 
#> Samples: 40000 
#> Mixture: SSVS 
#> Rho Test: muvar 
#> Credible Interval: 0.95 
#> ----
#> Call:
#> melsm.default(fixed_location = rt ~ congruency, random_location = ~congruency | 
#>     ID, fixed_scale = sigma ~ congruency, random_scale = ~congruency | 
#>     ID, mixture = "SSVS", k = 3, adapt = 5000, iter = 10000, 
#>     data = stroop)
#> ----
#> 
#> Random Effects Correlations
#> Posterior Distirbutions
#>                                                                     Post.mean Post.sd Cred.lb Cred.ub
#>  rho_01: location_(Intercept)_location_congruencyincongruent        -0.270    0.147   -0.557   0.014 
#>  rho_02: location_(Intercept)_scale_(Intercept)                      0.655    0.060    0.522   0.757 
#>  rho_03: location_(Intercept)_scale_congruencyincongruent           -0.365    0.105   -0.534  -0.131 
#>  rho_12: location_congruencyincongruent_scale_(Intercept)           -0.050    0.109   -0.317   0.131 
#>  rho_13: location_congruencyincongruent_scale_congruencyincongruent  0.735    0.082    0.536   0.858 
#>  rho_23: scale_(Intercept)_scale_congruencyincongruent              -0.564    0.076   -0.687  -0.400 
#> 
#> Inclusion Probabilities
#>                                                                     Pr.k1 Pr.k2 Pr.k3
#>  rho_01: location_(Intercept)_location_congruencyincongruent                         
#>  rho_02: location_(Intercept)_scale_(Intercept)                     0     1     0    
#>  rho_03: location_(Intercept)_scale_congruencyincongruent           0.018 0     0.982
#>  rho_12: location_congruencyincongruent_scale_(Intercept)           0.49  0.125 0.385
#>  rho_13: location_congruencyincongruent_scale_congruencyincongruent 0     1     0    
#>  rho_23: scale_(Intercept)_scale_congruencyincongruent                               
#> 
#> ----
#> 
#> Random Effects Standard Deviations
#>                                        Post.mean Post.sd Cred.lb Cred.ub
#>  tau_0: location_(Intercept)           0.100     0.007   0.087   0.113  
#>  tau_1: location_congruencyincongruent 0.027     0.004   0.018   0.036  
#>  tau_2: scale_(Intercept)              0.342     0.024   0.297   0.391  
#>  tau_3: scale_congruencyincongruent    0.234     0.020   0.196   0.276  
#> 
#> ----
#> 
#> Fixed Effects
#>                                         Post.mean Post.sd Cred.lb Cred.ub
#>  beta_0: location_(Intercept)            0.707    0.009    0.687   0.725 
#>  beta_1: location_congruencyincongruent  0.067    0.004    0.058   0.075 
#>  eta_0:  scale_(Intercept)              -1.784    0.033   -1.853  -1.724 
#>  eta_1:  scale_congruencyincongruent     0.155    0.026    0.103   0.205 
#> 
#> ----
#> note 
#> k1: 'spike' 
#> k2: positive slab 
#> k3: negative slab
```

#### Marginal Bayes Factors

k2 (positive component) vs. the compliment

``` r
bf <- marginal_bf(fit, H1 = "k2")

bf
#> hypMuVar: Bayesian Hypothesis Testing 
#>           of Mean--Variance Relations
#> Model: MELSM
#> Mixture: SSVS 
#> Rho Test: muvar 
#> Hypotheses
#> H1: rho = k2 
#> H2: rho = compliment 
#> -----
#>  Correlation                                                        BF_12
#>  rho_01: location_(Intercept)_location_congruencyincongruent           NA
#>  rho_02: location_(Intercept)_scale_(Intercept)                       Inf
#>  rho_03: location_(Intercept)_scale_congruencyincongruent           0.000
#>  rho_12: location_congruencyincongruent_scale_(Intercept)           0.286
#>  rho_13: location_congruencyincongruent_scale_congruencyincongruent   Inf
#>  rho_23: scale_(Intercept)_scale_congruencyincongruent                 NA
```

k2 (positive component) vs. k1(null component)

``` r
# note testing against the compliment is preferable
bf <- marginal_bf(fit, H1 = "k2", H2 = "k1")

bf
#> hypMuVar: Bayesian Hypothesis Testing 
#>           of Mean--Variance Relations
#> Model: MELSM
#> Mixture: SSVS 
#> Rho Test: muvar 
#> Hypotheses
#> H1: rho = k2 
#> H2: rho = k1 
#> -----
#>  Correlation                                                        BF_12
#>  rho_01: location_(Intercept)_location_congruencyincongruent           NA
#>  rho_02: location_(Intercept)_scale_(Intercept)                       Inf
#>  rho_03: location_(Intercept)_scale_congruencyincongruent           0.000
#>  rho_12: location_congruencyincongruent_scale_(Intercept)           0.255
#>  rho_13: location_congruencyincongruent_scale_congruencyincongruent   Inf
#>  rho_23: scale_(Intercept)_scale_congruencyincongruent                 NA
```

#### Posterior Model Probabilities

``` r
prob <- model_prob(fit)

prob
#> hypMuVar: Bayesian Hypothesis Testing 
#>           of Mean--Variance Relations
#> Model: MELSM
#> Mixture: SSVS 
#> Rho Test: muvar 
#> Components: 3 
#> ----
#> 
#>         Post.prob    BF.1i rho_02 rho_03 rho_12 rho_13
#> Model 1     0.482    1.000      2      3      1      2
#> Model 2     0.381    1.265      2      3      3      2
#> Model 3     0.119    4.051      2      3      2      2
#> Model 4     0.008   59.481      2      1      1      2
#> Model 5     0.006   85.274      2      1      2      2
#> Model 6     0.004  118.233      2      1      3      2
#> Model 7     0.000 2753.143      2      2      3      2
#> Model 8     0.000 3212.000      2      2      1      2
#> Model 9     0.000 3854.400      2      2      2      2
#> ----
#> 
#> rho_02: location_(Intercept)_scale_(Intercept)
#> rho_03: location_(Intercept)_scale_congruencyincongruent
#> rho_12: location_congruencyincongruent_scale_(Intercept)
#> rho_13: location_congruencyincongruent_scale_congruencyincongruent
```
