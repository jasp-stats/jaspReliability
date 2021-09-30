Bayesian Unidimensional Reliability Analysis
===

The frequentist unidimensional reliability analysis allows the user to test the scale's ability to consistently measure a unidimensional construct. In other words the analysis indicates the amount of error captured in the mesaurement.

## Input
---
- All columns (variables/items) of the dataset 

### Variables Box
- Variables: All variables of interest that are ordinally or metrically scaled

### Scale Statistics
- Credible interval: default is 95%
- McDonald's omega
	- Posterior predictive check: display a graphical check for the fit of the single factor model, in other words for the unidimensionality of the data
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
- If item dropped plot: display posterior densities of the reliability of the remaining items when a particular item is dropped
	- Order items:
		- Order items by mean: The densities are ordered so that the deletion of the item correpsonding to the topmost density brings the biggest change to the mean of the posterior distribution compared to the original distribution.
		- Order items by KL-divergence: Same principle as for the mean, except that the densities are now compared to the original distribution not by the mean but by the Kullback-Leibler divergence
		- Order items by KS-distance: Same principle as above, except that the densities are now compared to the original distribution by the Kolmogorov-Smirnov distance
- Item-rest correlation: The item-rest correlation indicates how the item correlates with the rest of the items
- Mean of the items
- Standard deviation of the items

### Plot posteriors
Display the posterior densities of the reliability coeffcients
- Fix range to 0-1: fix the x-axis of the plot to the interval [0, 1]
- Display priors: display the prior distributions of the coefficients

### Probability for
- Give the probability that the reliability coefficient is larger or smaller than some value or lies in a particular interval
- Shade posterior region in plot: the region in the posterior plot that corresponds to the interval specified above is shaded

## Convergence
### MCMC parameters
- No.samples: how many samples to draw from the posterior distribution
- No.burnin samples: how many of those samples to throw away after the start of the sampling
- Thinning: which indices of the posterior samples to use, e.g., each data point, every second, every third...
- No.chains: how many chains to obtain, meaning how often to start sampling form the posterior distribution

### Diagnostics
- R-hat: also called potential scale reduction factor
- Traceplots: display a traceplot per coefficient

### Repeatability
Since sampling from the posterior distribution is subjected to random processes, one can set a seed so that the background calculations in R yield equal results for equal seeds

## Reverse-Scaled Items
- This allows the user to select reverse-scaled items that need to be recoded.

## Missing Data Handling
### Missing Values
 - Exclude cases pairwise: The missing data is treated as unknown parameters in the posterior sampling process and is thereby imputed simultaneously 
 - Exclude cases listwise: Each row in the data set with one or more missing values will be deleted. Subsequently, the analysis continues with a data set with reduced observations.

## Output 
--- 
### Tables
#### Bayesian Scale Reliability Statistics: 
- Posterior mean: the mean of the posterior distribution 
- `...`% CI:
  - lower bound: The lower bound of the credible interval. 
  - upper bound: The upper bound of the credible interval. 
- R-hat: indicates the factor by how much the between chain variance would be reduced when the sampling would continue infinitely long, should be close to 1 and not larger than 1.1

#### Bayesian Individual Item Reliability Statistics:
- The first column contains all the variables included in the analysis. 
- If item dropped: 
	- mean of the posterior distribution of the coefficient for the remaining items
	- Lower `...` credible interval
	- Upper `...` credible interval

#### Probability that Reliability Statistic is Larger Than...
- Statistic: reliability coefficient 
- Probability: 
	- Prior: prior probability that the coefficient is larger than `...` and smaller than `...`
	- Posterior: posterior probability that the coefficient is larger than `...` and smaller than `...`

### Plots
#### Posterior Plots
- x-axis: reliability values, can be fixed to [0, 1]
- y-axis: density values
- black line: posterior density curve
- horizontal bar: limits of the `...`% credible interval
- optional: dotted line: prior density curve

#### If Item Dropped Posterior Plots:
- x-axis: reliability values
- y-axis: the items (unordered by default, meaning in order of appearance)
- display of the posterior densities of a reliability analysis with the remaining items when a particular item is dropped

#### Posterior Predictive Check Omega:
- x-axis: eigenvalue number
- y-axis: eigenvalue magnitude
- black dots: eigenvalues of the data covariance matrix
- gray bars: 95% credible interval of the posterior distribution of the eigenvalues of the model implied covariance matrix
- Omega is an adequate measure of reliability when the single factor model fits approximately, specifically, when all the eigenvalues of the sample (black dots) lie within the range of model implied eigenvalues (gray bars)

#### Convergence Traceplot:
- x-axis: data points in posterior sample chain
- y-axis: reliability values
- different colors represent different chains

## References
---
- Cronbach, L. J. (1951). Coefficient alpha and the internal structure of tests. *Psychometrika, 16*(3), 297–334. doi: 10.1007/BF02310555
- Gelman, A., & Rubin, D. B. (1992). Inference from iterative simulation using multiple sequences. *Statistical Science, 7*(4), 457–472. doi:10.1214/ss/1177011136
- Guttman, L. (1945). A basis for analyzing test-retest reliability. *Psychometrika, 10*(4), 255–282. doi: 10.1007/BF02288892
- Kolmogorov, A. N. (1933). Sulla determinazione empirica di une legge di distribuzion. *Instituto Italiano degli Attuari, Giornale, 4*, 83–91.
- Kullback, S., & Leibler, R. A. (1951). On information and sufficiency. *The Annals of Mathematical Statistics, 22*(1), 79–86. doi:10.1214/aoms/1177729694
- McDonald, R. P. (2013). *Test theory: A unified treatment*. New York, NJ, US: Psychology Press. doi: 10.4324/9781410601087
- Pfadt, J. M., van den Bergh, D., Sijtsma, K., Moshagen, M., & Wagenmakers, E.-J. (in press). Bayesian estimation of single-test reliability coefficients. *Multivariate Behavioral Research*. Preprint available at https://psyarxiv.com/exg2y
- Smirnov, N. (1939).  On the estimation of the discrepancy between empirical curves of distribution for two independent samples. *Bulletin Mathématique l’Université Moscou, 2*, 3–6.
- Woodhouse, B., & Jackson, P. H. (1977). Lower bounds for the reliability of the total score on a test composed of non-homogeneous items:  II: A search procedure to locate the greatest lower bound. *Psychometrika, 42*(4), 579–591. doi: 10.1007/bf02295980

## R Packages
---
- Bayesrel
- coda
- ggplot2
- ggridges
- LaplacesDemon
- MASS

## Example 

