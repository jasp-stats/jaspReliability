Frequentist Unidimensional Reliability Analysis
===

The frequentist unidimensional reliability analysis allows the user to test the scale's ability to consistently measure a unidimensional construct. In other words the analysis indicates the amount of error captured in the mesaurement.

## Input
---
- All columns (variables/items) of the dataset 

### Variables Box
- Variables: All variables of interest that are ordinally or metrically scaled

### Scale Statistics
- Confidence interval: default is 95%
- Coefficient omega: The same as McDonald's omega (for unidimensional data, based on the single-factor model). Note the total test variance in the denominator of the reliability equation is the model implied total variance, that is, the summed model implied covariance matrix.
- Coefficient alpha: The same as Cronbach's alpha (for binary items the coefficient equals KR20)
- Guttman's lambda 2
- Split-half coefficient: Correlates the sum scores of two test-halves. By default the variables are split into odd and even numbered items in order or appearance in the variables window. If another split is desired the variables just need to be reordered.
- Average interitem correlation
- Mean:
	- of the sum scores of participants
	- of the mean scores of participants
- Variance:
	- of the sum scores of participants
	- of the mean scores of participants
- Standard deviation: 
	- of the sum scores of participants
	- of the mean scores of participants
	
### Individual Item Statistics
- Coefficient omega
- Coefficient alpha
- Guttman's lambda 2
- Split-half coefficient
- Item-rest correlation: The item-rest correlation indicates how the item correlates with the rest of the items
- Mean of each item
- Variance of each item
- Standard deviation of each item

## Reverse-Scaled Items
- This allows the user to select reverse-scaled items that need to be recoded.

## Advanced Options
### Confidence Intervals
- Reliability coefficients: Control over the confidence intervals for the reliability coefficients
	- Analytic: Interval is created using the analytically derived standard error in a Wald-type CI, that is, normal-theory based. The standard errors are from van der Ark (2024)
	- Bootstrapped: Interval is created through bootstrapping and is percentile type:
		- No. of bootstrap samples: The number of times bootstrapped data sets are created and statistics are calculated.
		- Non-parametric bootstrap: The bootstrapped data sets are created by resampling with replacement 
		- Parametric bootstrap: The bootstrapped data sets are created by repeatedly sampling from a multivariate Normal with the original data as parameters
- Variance and SD: Confidence intervals for the variance and standard deviation may be based on the chi-square distribution or non-parametric (see van der Ark, 2024)
- The CI for the mean is always analytic and thus Wald-type.

### Missing values
 - Exclude cases pairwise: For the data covariance matrix, each covariance will be computed using all cases with valid data for the corresponding variables. Sample sizes may therefore vary across the covariances.
- Exclude cases listwise: Each row in the data set with one or more missing values will be deleted. Subsequently, the analysis continues with a data set with reduced observations.

### Samples
- Disable the saving of bootstrap samples: In case you want to save space for your output file, you can check this box. Beware that this will also lead to a loss in speed for the analysis. This happens because some samples inside the reliability module are precomputed and stored, so that the analysis can move forward in a much faster way. However, this also results in an increased size of the output object, and if you decide to save your analysis the resulting file will contain these samples. If you decide to run the analysis with a large number of bootstrap samples you might want to check that box if you do not want an increased file size for your output. 

### Coefficient omega estimation: 
- CFA: The single factor model is fit in a confirmatory factor analysis. 
	- Single Factor Model Fit: common fit indices of a CFA
- PFA: The single factor model is fit in a principal factor analysis. 
- Standardized factor loadings: Display a table with the standardized loadings of the single-factor model
	
### Coefficients: 
- Unstandardized
- Standardized

### Repeatability
When bootstrapping is involved, set a seed, so that the background calculations in R yield equal results for equal seeds.



## Output 
--- 
### Tables
#### Frequentist Scale Reliability Statistics: 
- Coefficient
- Estimate: point estimate value
- Std. error and CI: Analytic or bootstrapped

#### Frequentist Individual Item Reliability Statistics: 
- The first column contains all the variables included in the analysis. 
- If item dropped: reliability coefficients if a particular item is dropped
- May also show: 
	- the correlation of each item with the rest of the scale
	- the mean, variance, and standard deviation of each item

#### Fit Measures of Single Factor Model Fit
- Chi-Square: test statistic of model fit
- df: the degrees of freedom of the test statistic
- p.value: indicates how likely it is to obtain a test statistic at least as extreme as the one observed under the H0 (the model fits perfectly); values closer to 0 indicate poor fit
- RMSEA: root mean square error of approximation, a value close to 0 is desired
- Lower and Upper CI RMSEA: the lower and upper bound of the confidence interval of the RMSEA value
- SRMR: standardized root mean square residual, a value close to 0 is desired 

#### Standardized loadings of the Single-Factor Model:
- Contains the standardized single-factor loadings

## References
---
- van der Ark, A. (2024). Standard Errors for Reliability Coefficients. [Manuscript under revision]
- Bonett, D. G., & Wright, T. A. (2015). Cronbach's alpha reliability: Interval estimation, hypothesis testing, and sample size planning. *Journal of Organizational Behavior, 36*(1), 3-15. https://doi.org/10.1002/job.1960
- Cronbach, L. J. (1951). Coefficient alpha and the internal structure of tests. *Psychometrika, 16*(3), 297–334. https://doi.org/10.1007/BF02310555
- Guttman, L. (1945). A basis for analyzing test-retest reliability. *Psychometrika, 10*(4), 255–282. https://doi.org/10.1007/BF02288892
- McDonald, R. P. (2013). *Test theory: A unified treatment*. New York, NJ, US: Psychology Press. https://doi.org/10.4324/9781410601087
- Rencher, A. C.  (2002). *Methods of multivariate analysis*. New York, NY, USA: John Wiley & Sons, Inc.  https://doi.org/10.1002/0471271357
- Woodhouse, B., & Jackson, P. H. (1977). Lower bounds for the reliability of the total score on a test composed of non-homogeneous items:  II: A search procedure to locate the greatest lower bound. *Psychometrika, 42*(4), 579–591. https://doi.org/10.1007/bf02295980

## R Packages
---
- Bayesrel
- MASS
- psych

## Example 
- For an example go to `Open` --> `Data Library` --> `13. Reliability` --> `Fear of Statistics`. 
