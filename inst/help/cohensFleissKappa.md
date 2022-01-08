Cohen's and Fleiss' kappa
==========================

Cohen's and Fleiss' kappa are coefficients that measure the agreement between raters on a nominal or ordinal scale (Cohen, 1960; Cohen, 1968; Fleiss, 1971). Cohen's kappa is limited to measure the agreement between two raters. Fleiss' kappa measures the agreement between two or more raters.

*Rater* can refer to different judges, tests or other forms of rating here and in the section below.

### Input
-------

- Variables: The variables / columns to compute kappa for. Each variable corresponds to one rater, with different rows corresponding to different subjects being rated. 
The maximum number of raters for Cohen's kappa is 2. If more than 2 variables/raters are entered, all possible pairs of the variables/raters will be calculated.

- Cohen's kappa: Whether Cohen's kappa is calculated.
  - Unweighted or weighted: The unweighted kappa treats all disagreements equally, whereas the weighted kappa takes degrees of disagreement into account (Cohen, 1968). To calculate degrees of disagreement, a meaningful order of the rating categories is necessary. Therefore, weighted kappa can only be calculated for ordinal ratings. 
  
- Fleiss' kappa: Whether Fleiss' kappa is calculated.

- Confidence Interval: Whether a confidence interval should be reported for Cohen's and Fleiss' kappa and the size of the interval.

### Output
-------

#### Cohen's kappa (table)
This table shows Cohen's kappa for all possible pairs of raters and the average kappa. If selected, a confidence interval for the estimate will be reported.

Landis and Koch (1977) suggest the following guideline for the interpretation of Cohen's kappa:
- Less than 0: poor agreement
- Between 0.01 and 0.20: slight agreement
- Between 0.21 and 0.40: fair agreement
- Between 0.41 and 0.60: moderate agreement
- Between 0.61 and 0.80: substanital agreement
- Between 0.81 and 1: Almost perfect agreement


#### Fleiss' kappa (table)
This table shows Fleiss' kappa for the overall agreement and per rating category. If selected, a confidence interval for the estimate will be reported.

Landis and Koch (1977) suggest the following guideline for the interpretation of Fleiss' kappa:
- Less than 0: poor agreement
- Between 0.01 and 0.20: slight agreement
- Between 0.21 and 0.40: fair agreement
- Between 0.41 and 0.60: moderate agreement
- Between 0.61 and 0.80: substanital agreement
- Between 0.81 and 1: Almost perfect agreement

### References
-------
- Cohen, J. (1960). A coefficient of agreement for nominal scales. Educational and psychological measurement, 20(1), 37-46.
- Cohen, J. (1968). Weighted kappa: nominal scale agreement provision for scaled disagreement or partial credit. Psychological bulletin, 70(4), 213.
- Fleiss, J. L. (1971). Measuring nominal scale agreement among many raters. Psychological bulletin, 76(5), 378.
- Landis, J. R., & Koch, G. G. (1977). The measurement of observer agreement for categorical data. biometrics, 159-174.

### R Packages
---
- psych
