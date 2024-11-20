
<!-- README.md is generated from README.Rmd. Please edit that file -->

# povmap: Extension to the package emdi <img src="man/figures/logo.png" align="right" height="139" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/SSA-Statistical-Team-Projects/povmap/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/SSA-Statistical-Team-Projects/povmap/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Overview

The R package povmap is designed to facilitate the production of small
area estimates of means and poverty headcount rates. It adds several new
features to the [emdi](https://CRAN.R-project.org/package=emdi) package.
These include new options for:

- incorporating survey weights,
- ex-post benchmarking of estimates,
- two additional outcome variable transformations,
- and several new convenience functions to assist with reporting
  results.

## Installation

``` r
## To install the package from CRAN, run the following in R: 

install.packages("povmap")
```

## Development Version

``` r
## To get a bug fix or to use a feature from the development version, you can install 
## the development version of povmap from GitHub. Sometimes, povmap maybe unavailable in CRAN (although this is typically ## unlikely)

devtools::install_github("SSA-Statistical-Team-Projects/povmap")
## alternatively (remove the comment for the line of code below),
##remotes::install_github("SSA-Statistical-Team-Projects/povmap")
```

## Usage

``` r
library(povmap)

### estimate a unit level model with sample and population weights
ebp_model <- ebp(fixed = eqIncome ~ gender + eqsize + cash +
                 self_empl + unempl_ben + age_ben + surv_ben + sick_ben +
                 dis_ben + rent + fam_allow + house_allow + cap_inv +
                 tax_adj,
                 pop_data = eusilcA_pop, 
                 pop_domains = "district",
                 smp_data = eusilcA_smp, 
                 smp_domains = "district",
                 na.rm = TRUE, 
                 weights = "weight",
                 pop_weights = "hhsize", 
                 MSE = TRUE, 
                 weights_type = "nlme",
                 B = 2, 
                 L = 2)

summary(ebp_model)
#> Empirical Best Prediction
#> 
#> Call:
#>  ebp(fixed = eqIncome ~ gender + eqsize + cash + self_empl + unempl_ben + 
#>     age_ben + surv_ben + sick_ben + dis_ben + rent + fam_allow + 
#>     house_allow + cap_inv + tax_adj, pop_data = eusilcA_pop, 
#>     pop_domains = "district", smp_data = eusilcA_smp, smp_domains = "district", 
#>     L = 2, MSE = TRUE, B = 2, na.rm = TRUE, weights = "weight", 
#>     pop_weights = "hhsize", weights_type = "nlme")
#> 
#> Out-of-sample domains:  24 
#> In-sample domains:  70 
#> 
#> Sample sizes:
#> Units in sample:  1945 
#> Units in population:  25000 
#>                    Min. 1st Qu. Median      Mean 3rd Qu. Max.
#> Sample_domains       14    17.0   22.5  27.78571   29.00  200
#> Population_domains    5   126.5  181.5 265.95745  265.75 5857
#> 
#> Explanatory measures for the mixed model:
#>  Marginal_R2 Conditional_R2
#>     0.673562      0.7251496
#> 
#> Residual diagnostics for the mixed model:
#>                Skewness Kurtosis Shapiro_W    Shapiro_p
#> Error         0.6572179 9.757903 0.9521260 9.944428e-25
#> Random_effect 0.4185469 2.852405 0.9801002 3.295740e-01
#> 
#> ICC:  0.1580318 
#> 
#> Transformation:
#>  Transformation Method Optimal_lambda Shift_parameter
#>         box.cox   reml      0.6046901               0

### showcasing the actual poverty map
load_shapeaustria()  ## reading in the shapefile

## the new `map_plot` now requires sf shapefile objects
## rather than the older `sp` (it's dependencies are being phased out)
map_plot(object = ebp_model, 
         MSE = FALSE, 
         CV = FALSE, 
         map_obj = shape_austria_dis, 
         indicator = c("Head_Count"), 
         map_dom_id = "PB")
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

``` r


###### ---------------------------------------------------- #######
# compute direct estimates 
direct_est <- direct(y = "eqIncome", 
                     smp_data = eusilcA_smp,
                     smp_domains = "district", 
                     weights = "weight",
                     var = TRUE, 
                     B = 2)

### some general pre-estimation statistics
ebp_reportdescriptives(model = ebp_model, 
                       direct = direct_est,
                       smp_data = eusilcA_smp, 
                       weights = "weight",
                       pop_weights = "hhsize", 
                       CV_level = "state",
                       pop_data = eusilcA_pop, 
                       pop_domains = "district")
#> $cv_table
#>                    indicator       ebp_cv  direct_cv
#> 1 CV for Area: Lower Austria 0.0005623766 0.03030733
#> 2    CV for Area: Vorarlberg 0.0096720887 0.16370318
#> 3 CV for Area: Upper Austria 0.3875383775 7.39427881
#> 4        CV for Area: Styria 0.0053397075 0.31487144
#> 5      CV for Area: Salzburg 0.0768486572 0.45859197
#> 6         CV for Area: Tyrol 0.0192320593 0.68041843
#> 7     CV for Area: Carinthia 0.0154801556 0.25072764
#> 8    CV for Area: Burgenland 0.0209849881 0.87719916
#> 9        CV for Area: Vienna 0.0269278962 0.17973145
#> 
#> $basicinfo_df
#>                indicator census survey
#> 1        Number of Units  25000   1945
#> 2      Number of Regions      9      9
#> 3 Number of Target Areas     94     70
#> 
#> $poverty_df
#>               indicator        model        survey
#> 1 National Poverty Rate     0.141604     0.1528866
#> 2 National Poverty Line 10924.320000 10924.3200000


### compare the means between survey and census prior to model estimation
variables <- c("gender", "eqsize", "cash", "self_empl",
               "unempl_ben", "age_ben", "surv_ben",
               "sick_ben", "dis_ben", "rent", "fam_allow",
               "house_allow", "cap_inv", "tax_adj")

ebp_test_means(varlist = variables,
               pop_data = eusilcA_pop,
               smp_data = eusilcA_smp,
               weights = "weight")
#>       variable    smp_means    pop_means          diff    pvalue
#> 1      age_ben  5478.609860  5304.056507 -1.745534e+02 0.9898989
#> 2      cap_inv   487.265886   456.851454 -3.041443e+01 0.9944134
#> 3         cash 12712.233319 12706.449297 -5.784022e+00 0.9997571
#> 4      dis_ben   523.440415   542.325741  1.888533e+01 0.9965279
#> 5       eqsize     1.610035     1.603376 -6.659056e-03 0.9937560
#> 6    fam_allow  1617.487932  1624.957604  7.469672e+00 0.9987536
#> 7       gender     1.371915     1.390840  1.892518e-02 0.9780223
#> 8  house_allow    61.819822    65.413716  3.593893e+00 0.9947154
#> 9         rent   608.004341   770.357579  1.623532e+02 0.9856212
#> 10   self_empl  2039.125009  1966.965080 -7.215993e+01 0.9949896
#> 11    sick_ben    61.836570    67.450918  5.614348e+00 0.9957620
#> 12    surv_ben    61.917778    88.952360  2.703458e+01 0.9838119
#> 13     tax_adj   -76.558371   -82.539963 -5.981592e+00 0.9978321
#> 14  unempl_ben   463.508414   429.382134 -3.412628e+01 0.9903023

### report EBP model coefficients with significance levels
ebp_reportcoef_table(ebp_model, 4)
#>        Variable       coeff std_error
#> 1   (Intercept) 406.4743***   10.6335
#> 2  genderfemale      7.8761    5.0849
#> 3        eqsize      -32***    4.8715
#> 4          cash   0.0125***   0.00025
#> 5     self_empl   0.0104***   0.00031
#> 6    unempl_ben   0.0087***    0.0012
#> 7       age_ben   0.0129***   0.00031
#> 8      surv_ben   0.0128***    0.0026
#> 9      sick_ben   0.0110***    0.0033
#> 10      dis_ben   0.0137***    0.0008
#> 11         rent   0.0083***    0.0004
#> 12    fam_allow    -0.00006    0.0009
#> 13  house_allow    0.0156**    0.0064
#> 14      cap_inv   0.0088***    0.0007
#> 15      tax_adj  -0.0054***    0.0014


### compute the CVs (several CV types are estimated including
### Horowitz Thompson CVs, CVs accounting for the design effect,
### bootstraped CVs i.e. calibrated or naive)
head(ebp_compute_cv(model = ebp_model, calibvar = "gender"))
#>               Domain Direct_Head_Count EBP_Head_Count HT_Head_Count_CV
#> 1          Amstetten         0.2727273    0.272976680        0.3157348
#> 2              Baden         0.0250000    0.041176471        0.9484184
#> 3            Bludenz         0.4705882    0.392251816        0.3344887
#> 4     Braunau am Inn         0.3793103    0.383753501        0.2855876
#> 5            Bregenz         0.0000000    0.006823821              NaN
#> 6 Bruck-Mürzzuschlag         0.0000000    0.020138889              NaN
#>   CB_Head_Count_CV DesignEffect_CV EBP_Head_Count_CV
#> 1        0.1963673       0.2027130        0.05618261
#> 2        1.9158379       1.8528466        0.50857165
#> 3        0.1512468       0.1524379        0.37168619
#> 4        0.1700101       0.1747679        0.11225965
#> 5              NaN             NaN        3.58601508
#> 6              NaN             NaN        2.56021421

### report the model fit and normality assumptions
ebp_normalityfit(model = ebp_model)
#>          indicator     value
#> 1     rsq_marginal 0.6735620
#> 2  rsq_conditional 0.7251496
#> 3 epsilon_skewness 0.6572179
#> 4 epsilon_kurtosis 9.7579034
#> 5  random_skewness 0.4185469
#> 6  random_kurtosis 2.8524049
```

## Getting help

If you encounter a clear bug, please file an issue with a minimal
reproducible example on
[GitHub](https://github.com/SSA-Statistical-Team-Projects/povmap/issues)
