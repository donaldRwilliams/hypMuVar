
<style>
pre {
  overflow-x: auto;
}
pre code {
  word-wrap: normal;
  white-space: pre;
}
</style>
hypMuVar: Bayesian Hypothesis Testing of Mean-Variance Relations
----------------------------------------------------------------

Random slopes for the location and scale:

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
#>                                                                    
#>  rho_01: location_(Intercept)_location_congruencyincongruent       
#>  rho_02: location_(Intercept)_scale_(Intercept)                    
#>  rho_03: location_(Intercept)_scale_congruencyincongruent          
#>  rho_12: location_congruencyincongruent_scale_(Intercept)          
#>  rho_13: location_congruencyincongruent_scale_congruencyincongruent
#>  rho_23: scale_(Intercept)_scale_congruencyincongruent             
#>  Post.mean Post.sd Cred.lb Cred.ub
#>  -0.270    0.147   -0.557   0.014 
#>   0.655    0.060    0.522   0.757 
#>  -0.365    0.105   -0.534  -0.131 
#>  -0.050    0.109   -0.317   0.131 
#>   0.735    0.082    0.536   0.858 
#>  -0.564    0.076   -0.687  -0.400 
#> 
#> Inclusion Probabilities
#>                                                                     Pr.k1
#>  rho_01: location_(Intercept)_location_congruencyincongruent             
#>  rho_02: location_(Intercept)_scale_(Intercept)                     0    
#>  rho_03: location_(Intercept)_scale_congruencyincongruent           0.018
#>  rho_12: location_congruencyincongruent_scale_(Intercept)           0.49 
#>  rho_13: location_congruencyincongruent_scale_congruencyincongruent 0    
#>  rho_23: scale_(Intercept)_scale_congruencyincongruent                   
#>  Pr.k2 Pr.k3
#>             
#>  1     0    
#>  0     0.982
#>  0.125 0.385
#>  1     0    
#>             
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
