#P-value Weighting R Package

##Introduction
This R package contains p-value weighting methods for multiple hypothesis testing. These statistical methods are used for improving power in multiple testing via the use of prior information. The iGWAS method is provided for applications of p-value weighting in Genome-Wide Association Studies.

Some of these methods were were developed by the authors in the following paper:
*Optimal Multiple Testing Under a Gaussian Prior on the Effect Sizes* by Dobriban, Fortney, Kim, Owen:  http://arxiv.org/abs/1504.02935

## An Example
Suppose we want to find the significant effects in a large pool of candidates, by performing  multiple testing with the p-values `P_current`. We have some prior information about the size of each effect size in the form of prior test statistics `t1` with estimated variances `sigma`. The prior effects `t1` with standard errors `sigma` are our prior guesses for the current effects.

To use this information for improving power in multiple hypothesis testing, we can give each hypothesis a different weight. For instance, we can compute the Bayes p-value weights `w`, and use them to weight the current p-values `P_weighted` as follows: 

```{r}
w <- bayes_weights(t1,sigma,q)$w
P_weighted <- P_current/w
```

Finally, we perform weighted Bonferroni multiple testing, controlling the Family-Wise Error Rate, in the usual way: 
```{r}
P_w_adjusted <- p.adjust(P_weighted,"bonferroni")
```

Please see the vignette or the help files for examples and a description of the methods.

## Installation instructions

To install from GitHub, make sure that Hadley Wickham's devtools (https://github.com/hadley/devtools) is installed, then run:

```{r}
devtools::install_github("dobriban/pweight")
```

The package is also available from CRAN at:  https://cran.r-project.org/web/packages/pweight/.
To install from CRAN within R, type

```{r}
install.packages("pweight")
```

## Main components

The core of the package consists of the following p-value weighting methods, whose details are described in the documentation:

* Spjotvoll weights: `spjotvoll_weights()` 
* Exponential weights: `exponential_weights()` 
* Bayes weights: `bayes_weights()`

The iGWAS method, `iGWAS()`, is provided for applications of p-value weighting in Genome-Wide Association Studies.

Please contact the authors if you have any questions. 
