# mlesky : Maximum likelihood inference of effective population size through time using GMRF-skygrid model

This package is related to previous Bayesian implementations of the GMRF-skygrid model with the following notable differences: 

* The GMRF process takes place on the 2nd-order difference of log(Ne), which is more similar to the `skygrowth` model of Volz & Didelot 
* When computing the GMRF likelihood, the smoothing parameter (ie the precision of random walk) is fixed
* We use a novel cross-validation approach for selecting the smoothing parameter 

## Installation

In R, install the `devtools` package and run 
```
devtools::install_github('emvolz-phylodynamics/mleksy')
```

## Roadmap 

* Add non-parametric bootstrap as alternative for CI
* Add regression models for non-genetic time-series 
* Parameter to select order of differencing in GMRF 

## MRSA example 

Here we reanalyze dated phylogenies from __Volz & Didelot, Systematic Biology 2018__. 


```r
require(mlesky)

# load the tree
tree <- ape::read.tree(system.file('mrsa.nwk', package = 'mlesky'))

# run mleksy with defaults
(fit <- mlskygrid( tree ))
```

```
## mlskygrid fit
## 	Smoothing parameter tau = 1 
## 
## Estimated Ne(t): 
##    Time before most recent sample        2.5%        MLE      97.5%
## 1                      -44.068800   0.2822168   2.050882   14.90385
## 2                      -42.306048   0.6136106   3.682573   22.10089
## 3                      -40.543296   0.7438961   4.273402   24.54908
## 4                      -38.780544   1.6219125   9.159138   51.72277
## 5                      -37.017792   3.8761971  31.145086  250.24950
## 6                      -35.255040   5.0408578  55.874618  619.33367
## 7                      -33.492288   4.3411625  50.577077  589.25246
## 8                      -31.729536   3.4978649  39.848799  453.97031
## 9                      -29.966784   6.0190939  90.282889 1354.19054
## 10                     -28.204032   8.6188396 203.984287 4827.74843
## 11                     -26.441280  10.1180625 292.159643 8436.12672
...
```

```r
plot ( fit )
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-1.png)

Note there's a lot of uncertainty >20 years in the past when there are very few lineages remaining in the tree. We can focus the analysis on the portion of the tree that is more informative using the `NeStartTimeBeforePresent` parameter. 
We also us cross-validation to optimize the smoothing parameter: 

```r
(
fit <- mlskygrid( tree
  , tau = NULL, tau_lower = 1, tau_upper = 20
  , ncpu = 6
  , res = 50
  , NeStartTimeBeforePresent = 20)
)
```


```
## mlskygrid fit
## 	Smoothing parameter tau = 16.169103648728 
## 
## Estimated Ne(t): 
##    Time before most recent sample       2.5%        MLE      97.5%
## 1                           -20.0   2.145668   4.868885   11.04833
## 2                           -19.6   2.776083   5.733437   11.84125
## 3                           -19.2   3.556557   6.756595   12.83589
## 4                           -18.8   4.518919   7.979456   14.09003
## 5                           -18.4   5.708317   9.462214   15.68474
## 6                           -18.0   7.188612  11.288379   17.72630
## 7                           -17.6   9.055953  13.582237   20.37082
## 8                           -17.2  11.439543  16.504876   23.81310
## 9                           -16.8  14.512483  20.267514   28.30474
## 10                          -16.4  18.517597  25.166042   34.20150
## 11                          -16.0  23.786166  31.607349   42.00023
...
```

```r
plot( fit )
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2-1.png)

## HIV example 

We analyzed 399 HIV-1 sequences from Senegal between 1990 and 2014. 
All sequences are subtype CRF02_AG. `treedater` analysis shows a common ancestor around 1970 with LTT having rapid change in the early 1980s when the HIV epidemic was expanding. 

Here is the `mlesky` analysis 


```r
require(mlesky)

# load the tree 
tree <- ape::read.tree( system.file('sn02ag2.0.nwk', package='mlesky') )

# mlesky with defaults (tau = 1, res = 25)
(fit <- mlskygrid( tree ))
```

```
## mlskygrid fit
## 	Smoothing parameter tau = 1 
## 
## Estimated Ne(t): 
##    Time before most recent sample         2.5%          MLE       97.5%
## 1                      -42.668213    0.2645066     2.239153    18.95531
## 2                      -40.961485    0.7530669     8.095022    87.01669
## 3                      -39.254756    0.8960033    14.657942   239.79296
## 4                      -37.548028    0.8643560    12.256399   173.79335
## 5                      -35.841299    0.7560536     6.119641    49.53353
## 6                      -34.134570    0.9710828     3.536995    12.88287
## 7                      -32.427842    7.5995030    10.317788    14.00838
## 8                      -30.721113   61.1948002    77.950490    99.29404
## 9                      -29.014385  215.3944139   278.510600   360.12148
## 10                     -27.307656  442.7015063   575.878606   749.11913
...
```

```r
plot( fit )
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-1.png)

Now we use cross-validation to find the smoothing parameter and increase the time resolution: 


```r
( fit <- mlskygrid( tree , tau = NULL, tau_lower = .1, tau_upper = 20 , ncpu = 6, res = 100) )
```

```
## mlskygrid fit
## 	Smoothing parameter tau = 15.4640465609943 
## 
## Estimated Ne(t): 
##     Time before most recent sample        2.5%        MLE       97.5%
## 1                      -42.6682131    7.374377  715.40143 69402.36325
## 2                      -42.2415310    7.006665  537.46775 41228.11208
## 3                      -41.8148488    6.624695  403.75582 24607.73790
## 4                      -41.3881667    6.232627  303.30144 14759.70974
## 5                      -40.9614846    5.825396  227.49005  8883.81187
## 6                      -40.5348025    5.432208  171.11011  5389.82899
## 7                      -40.1081203    5.044602  128.84388  3290.79412
## 8                      -39.6814382    4.677240   97.38374  2027.60435
## 9                      -39.2547561    4.344917   74.14939  1265.41705
## 10                     -38.8280739    4.050746   56.97793   801.45367
...
```

```r
plot(fit) 
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-1.png)

Cross-validation indicates a higher smoothing parameter is optimal. 
