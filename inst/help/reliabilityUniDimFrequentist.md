Frequentist Unidimensional Reliability Analysis
===

The frequentist unidimensional reliability analysis allows the user to test the scale's ability to consistently measure a unidimensional construct. In other words the analysis indicates the amount of error captured in the mesaurement.

## Input
- All columns (variables/items) of the dataset 

### Variables Box
- Variables: All variables of interest that are ordinally or metrically scaled

### Scale Statistics
- Confidence interval: default is 95%
- McDonald's omega
- Cronbach's alpha
- Guttman's lambda 2
- Guttman's lambda 6
- Greatest lower bound
- Average interitem correlation
- Mean:
	- of the sum scores of participants
	- of the mean scores of participants
- Standard deviation: 
	- of the sum scores of participants
	- of the mean scores of participants
	
### Individual Item Statistics
- McDonald's omega
- Guttman's lambda 2
- Guttman's lambda 6
- Greatest lower bound
- Item-rest correlation: The item-rest correlation indicates how the item correlates with the rest of the items
- Mean of the items
- Standard deviation of the items

## Reverse-Scaled Items
- This allows the user to select reverse-scaled items that need to be recoded.

## Advanced Options
### Missing Values
 - Exclude cases pairwise: For the data covariance matrix, each covariance will be computed using all cases with valid data for the corresponding variables. Sample sizes may therefore vary across the covariances.
- Exclude cases listwise: Each row in the data set with one or more missing values will be deleted. Subsequently, the analysis continues with a data set with reduced observations.

### Bootstrap 
The number of times bootstrapped data sets are created and statistics are calculated. The bootstrapped intervals are percentile type.
- Non-parametric bootstrap: The bootstrapped data sets are created by resampling with replacement 
- Parametric bootstrap: The bootstrapped data sets are created by repeatedly sampling from a multivariate Normal with the original data as parameters
	
### McDonald's omega Estimation: 
- CFA: The single factor model is fit in a confirmatory factor analysis. 
	- Single Factor Model Fit: common fit indices of a CFA
	Interval: 
		- Analytic interval: Wald-CI from the CFA
		- Bootstrapped interval: CI from bootstrapped CFA
- PFA: The single factor model is fit in a principal factor analysis. 
	
 ### Cronbach's alpha Estimation: 
- Unstandardized: Cronbach's alpha is calculated from the data covariance matrix
- Standardized: Cronbach's alpha is calculated from the data correlation matrix
- Interval: 
	- Analytic interval: The interval is calculated by means of the method from (Bonett & Wright, 2015)
	- Bootstrapped interval: The interval is calculated by means of the percentile interval based on alpha's bootstrapped sample
		
### Repeatability
When bootstrapping is involved, set a seed, so that the background calculations in R yield equal results for equal seeds


## Output 
--- 
### Tables
#### Frequentist Scale Reliability Statistics: 
- Point estimate: 
- `...`% CI:
  - lower bound: The lower bound of the confidence interval. 
  - upper bound: The upper bound of the confidence interval. 
- McDonald's omega: by default the method to obtain the point estimate is a CFA, and to obtain the interval it is the Wald-CI
- Cronbach's alpha: by default the point estimate is unstandardized, the interval is by default analytic
- Guttman's lambda 2: the confidence interval is bootstrapped
- Guttman's lambda 6: the confidence interval is bootstrapped
- Greatest lower bound: the confidence interval is bootstrapped

#### Frequentist Individual Item Reliability Statistics: 
- The first column contains all the variables included in the analysis. 
- If item dropped: reliability coefficients if a particular item is dropped

#### Fit Measures of Single Factor Model Fit
- Chi-Square: test statistic of model fit
- df: the degrees of freedom of the test statistic
- p.value: indicates how likely it is to obtain a test statistic at least as extreme as the one observed under the H0 (the model fits perfectly); values closer to 0 indicate poor fit
- RMSEA: root mean square error of approximation, a value close to 0 is desired
- Lower and Upper CI RMSEA: the lower and upper bound of the confidence interval of the RMSEA value
- SRMR: standardized root mean square residual, a value close to 0 is desired 

## References
-------

- Bonett, D. G., & Wright, T. A. (2015). Cronbach's alpha reliability: Interval estimation, hypothesis testing, and sample size planning. *Journal of Organizational Behavior, 36*(1), 3-15. doi: 10.1002/job.1960
- Cronbach, L. J. (1951). Coefficient alpha and the internal structure of tests. *Psychometrika, 16*(3), 297–334. doi: 10.1007/BF02310555
- Guttman, L. (1945). A basis for analyzing test-retest reliability. *Psychometrika, 10*(4), 255–282. doi: 10.1007/BF02288892
- McDonald, R. P. (2013). *Test theory: A unified treatment*. New York, NJ, US: Psychology Press. doi: 10.4324/9781410601087
- Rencher, A. C.  (2002). *Methods of multivariate analysis*. New York, NY, USA: John Wiley & Sons, Inc.  doi:  10.1002/0471271357
- Woodhouse, B., & Jackson, P. H. (1977). Lower bounds for the reliability of the total score on a test composed of non-homogeneous items:  II: A search procedure to locate the greatest lower bound. *Psychometrika, 42*(4), 579–591. doi: 10.1007/bf02295980

## R Packages
---
- Bayesrel
- MASS
- psych

## Example 
- For an example go to `Open` --> `Data Library` --> `Descriptives` --> `Fear of Statistics`. 
