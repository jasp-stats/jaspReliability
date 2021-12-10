Bayesian Unidimensional Reliability Analysis
===

The Bayesian unidimensional reliability analysis allows the user to test the scale's ability to consistently measure a unidimensional construct. In other words the analysis indicates the amount of error captured in the measurement.

## Input
---
- All columns (variables/items) of the data set 

### Variables Box
- Variables: All variables of interest that are ordinally or metrically scaled

### Scale Statistics
- Credible interval: default is 95%
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

The CTT-coefficients alpha, lambda 2, lambda 6, and the glb are computed from the data covariance matrix.
Coefficient omega is computed from the centered data matrix.
	
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

### Samples
- Disable the saving of posterior MCMC samples:
In case you want to save space for your output file, you can check this box. Beware that this will also lead to a loss in speed for the analysis. This happens because some samples inside the reliability module are precomputed and stored, so that the analysis can move forward in a much faster way. However, this also results in an increased size of the output object, and if you decide to save your analysis the resulting file will contain these samples. If you decide to run the analysis with a large number of iterations you might want to check that box if you do not want an increased file size for your output. 

## Priors
### CTT-Coefficients (α, λ2, λ6, glb)
The prior distributions for alpha, lambda2, lambda6, the glb, and the average inter-item correlation are induced by the prior distribution on the covariance matrix, which, by default, is an inverse Wishart distribution with the identity matrix as a scaling matrix and the number of items k as the degrees of freedom. 

- Scale: Precision values of the diagonal of the scaling matrix of the inverse Wishart distribution
- Df: Degrees of freedom of the inverse Wishart distribution

### McDonald's ω residual variances
The prior distribution on McDonald’s omega is induced by the prior distributions on the single-factor model parameters, which are: a normal distribution centered on zero for the factor loadings and scores; an inverse gamma distribution with shape=2 and scale=1 for the residuals; and for the variance of the latent variables an inverse Wishart distribution with the number of items k as a scaling matrix (scalar, since it is of dimension one) and k+2 as the degrees of freedom.

- Shape: of inverse gamma prior distribution on the residual variances
- Scale: of inverse gamma prior distribution on the residual variances
- Mean: of normal prior distribution on the factor loadings


## Reverse-Scaled Items
This allows the user to select reverse-scaled items that need to be recoded.

## Advanced Options
### Missing Values
 - Bayesian imputation: The missing data are treated as unknown parameters in the posterior sampling process.
 The missing values are sampled conditional on the remaining data and the sampled model parameters. This way we obtain a posterior distribution for each missing value.
 - Exclude cases listwise: Each row in the data set with one or more missing values will be deleted. Subsequently, the analysis continues with a data set with reduced observations.
 
### McDonald's omega Estimation
- Posterior predictive check: Display a graphical check for the fit of the single factor model, in other words for the unidimensionality of the data. Answers the question if the posterior single-factor model predicts data similar to the original data set?
- Fit measures: Display a table with Bayesian fit indices for the single factor model: 
Bayesian LR, B-RMSEA, B-CFI, B-TLI
  - p(RMSEA < .08): probability that the B-RMSEA is smaller than a cutoff
  - p(CFI/TLI < .90): probability that the B-CFI or the B-TLI are larger than a cutoff
- Standardized factor loadings: Display a table with the standardized loadings of the single-factor model


#### Coefficients
- Unstandardized (the default)
- Standardized: The use of standardized coefficients is not undisputed in psychometric literature. 
Keep in mind, that the covariance matrix, given the item scales are equal, is more informative about the measurement. 
See:

Carl F. Falk & Victoria Savalei (2011) The relationship between unstandardized and standardized alpha, true reliability, and the underlying measurement model, *Journal of Personality Assessment, 93*(5), 445-453. https://doi.org/10.1080/00223891.2011.594129

Hayashi, K. and Kamata, A. 2005. A note on the estimator of the alpha coefficient for standardized variables under normality. *Psychometrika, 70*, 579–586.

Sun, W., Chou, C. P., Stacy, A. W., Ma, H., Unger, J. and Gallaher, P. 2007. SAS and SPSS macros to calculate standardized Cronbach's alpha using the upper bound of the phi coefficient for dichotomous items. *Behavior Research Methods, 39*, 71–81.

Moss, J. (2020). Please avoid the standardized alpha and the ordinal alpha. https://doi.org/10.31234/osf.io/nvg5d

Warrens, M.J. Some relationships between Cronbach’s alpha and the Spearman-Brown formula. *Journal of Classification, 32*, 127–137 (2015). https://doi.org/10.1007/s00357-015-9168-0

#### Posterior point estimate
- Mean (the default)
- Median

## Output 
--- 
### Tables
#### Bayesian Scale Reliability Statistics: 
- Posterior mean or median: of the posterior distribution 
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

#### Probability that Reliability Statistic is Larger Than...:
- Statistic: reliability coefficient 
- Probability: 
	- Prior: prior probability that the coefficient is larger than `...` and smaller than `...`
	- Posterior: posterior probability that the coefficient is larger than `...` and smaller than `...`

#### Fit Measures of the Single-Factor Model
- Point estimate is either mean or median (see ''posterior point estimate''), except for LR, which is always the mean
- Relative to cutoff: probability RMSEA and SRMR smaller than their specified cutoff, CFI and TLI larger than their cutoff
  - B-LR: Bayesian likelihood ratio of tested model and saturated model 
  - B-RMSEA: Bayesian RMSEA adapted from Garnier-Villareal and Jorgensen (2020), method equal to blavaan fit indices; 
  uses devM method and degrees of freedom from dic
  - B-CFI, B-TLI: Bayesian CFI and TLI adapted from Garnier-Villareal and Jorgensen (2020), equal to blavFitIndices,
  requires a null-model, here a model with only variances of the observed variables; uses devM method and degrees of freedom from dic

#### Standardized loadings of the Single-Factor Model:
- Mean or median of standardized single-factor loadings

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
- Cronbach, L. J. (1951). Coefficient alpha and the internal structure of tests. *Psychometrika, 16*(3), 297–334. https://doi.org/10.1007/BF02310555
- Garnier-Villarreal M., & Jorgensen T. D. (2020). Adapting fit indices for Bayesian structural equation modeling: Comparison to maximum likelihood. *Psychological Methods. 25*(1), 46-70. https://doi.org/10.1037/met0000224.
- Gelman, A., & Rubin, D. B. (1992). Inference from iterative simulation using multiple sequences. *Statistical Science, 7*(4), 457–472. https://doi.org/10.1214/ss/1177011136
- Guttman, L. (1945). A basis for analyzing test-retest reliability. *Psychometrika, 10*(4), 255–282. https://doi.org/10.1007/BF02288892
- Kolmogorov, A. N. (1933). Sulla determinazione empirica di une legge di distribuzion. *Instituto Italiano degli Attuari, Giornale, 4*, 83–91.
- Kullback, S., & Leibler, R. A. (1951). On information and sufficiency. *The Annals of Mathematical Statistics, 22*(1), 79–86. https://doi.org/10.1214/aoms/1177729694
- McDonald, R. P. (2013). *Test theory: A unified treatment*. New York, NJ, US: Psychology Press. https://doi.org/10.4324/9781410601087
- Pfadt, J. M., van den Bergh, D., Sijtsma, K., Moshagen, M., & Wagenmakers, E.-J. (in press). Bayesian estimation of single-test reliability coefficients. *Multivariate Behavioral Research*. Preprint available at https://psyarxiv.com/exg2y
- Smirnov, N. (1939).  On the estimation of the discrepancy between empirical curves of distribution for two independent samples. *Bulletin Mathématique l’Université Moscou, 2*, 3–6.
- Woodhouse, B., & Jackson, P. H. (1977). Lower bounds for the reliability of the total score on a test composed of non-homogeneous items:  II: A search procedure to locate the greatest lower bound. *Psychometrika, 42*(4), 579–591. https://doi.org/10.1007/bf02295980

## R Packages
---
- Bayesrel
- coda
- ggplot2
- ggridges
- LaplacesDemon
- MASS

## Example 
Go to: `Open` --> `Data Library` --> `13. Reliability` --> `ASRM - Mania Scale`.
