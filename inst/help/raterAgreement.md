Rater Agreement
==========================

Cohen's and Fleiss' kappa and Krippendorff's alpha are coefficients that measure the agreement between raters on a nominal or ordinal scale (Cohen, 1960; Cohen, 1968; Krippendorff, 1970; Fleiss, 1971).
 Cohen's kappa is limited to measure the agreement between two raters. Fleiss' kappa and Krippendorff's alpha measure the agreement between two or more raters.

*Rater* can refer to different judges, tests or other forms of rating here and in the section below.

### Input
-------

- Variables: The variables / columns to compute kappa for. Each variable corresponds to one rater, with different rows corresponding to different subjects being rated. 
Note that this data format is different from the common data format for Krippendorff's alpha, which has raters in rows and units of measurements (subjects) in columns. 
For Cohen's kappa the maximum number of raters  is 2. If more than 2 variables/raters are entered, all possible pairs of the variables/raters will be calculated. 

- Cohen's kappa: Whether Cohen's kappa is calculated.
  - Unweighted or weighted: The unweighted kappa treats all disagreements equally, whereas the weighted kappa takes degrees of disagreement into account (Cohen, 1968). To calculate degrees of disagreement, a meaningful order of the rating categories is necessary. Therefore, weighted kappa can only be calculated for ordinal ratings. 
  
- Fleiss' kappa: Whether Fleiss' kappa is calculated.

- Krippendorff's alpha: Whether Krippendorff's alpha is calculated.
  - Method: Krippendorff's alpha is calculated differently depending on the level of measurement of the variables. This option allows to select the method corresponding to the correct level of measurement. The available levels are nominal, ordinal, interval, and ratio.

- Confidence Interval: Whether a confidence interval should be reported for Cohen's and Fleiss' kappa and Krippendorff's alpha and the width of the interval.

### Output
-------

#### Cohen's kappa (table)
This table shows Cohen's kappa for all possible pairs of raters and the average kappa. If selected, a confidence interval for the kappa estimate will be reported.

Landis and Koch (1977) suggest the following guideline for the interpretation of Cohen's kappa:
- Less than 0: poor agreement
- Between 0.01 and 0.20: slight agreement
- Between 0.21 and 0.40: fair agreement
- Between 0.41 and 0.60: moderate agreement
- Between 0.61 and 0.80: substanital agreement
- Between 0.81 and 1: Almost perfect agreement


#### Fleiss' kappa (table)
This table shows Fleiss' kappa for the overall agreement and per rating category. If selected, a confidence interval for the kappa estimate will be reported.

Landis and Koch (1977) suggest the following guideline for the interpretation of Fleiss' kappa:
- Less than 0: poor agreement
- Between 0.01 and 0.20: slight agreement
- Between 0.21 and 0.40: fair agreement
- Between 0.41 and 0.60: moderate agreement
- Between 0.61 and 0.80: substanital agreement
- Between 0.81 and 1: Almost perfect agreement

#### Krippendorff's alpha (table)
This table shows Krippendorff's alpha for the overall agreement. If selected, a confidence interval for the alpha estimate will be reported.

Krippendorff (2004) suggests the following guideline for the interpretation of Krippendorff's alpha:
- Less than 0.66: unacceptable agreement
- Between 0.66 and 0.80: tentatively acceptable agreement
- Between 0.81 and .99: acceptable agreement
- 1: perfect agreement

### References
-------
- Cohen, J. (1960). A coefficient of agreement for nominal scales. *Educational and Psychological Measurement, 20*(1), 37-46. https://doi.org/10.1177/001316446002000104
- Cohen, J. (1968). Weighted kappa: nominal scale agreement provision for scaled disagreement or partial credit. *Psychological Bulletin, 70*(4), 213. https://doi.org/10.1037/h0026256
- Fleiss, J. L. (1971). Measuring nominal scale agreement among many raters. *Psychological Bulletin, 76*(5), 378. https://doi.org/10.1037/h0031619
- Krippendorff, K (1970). Estimating the reliability, systematic error, and random error of interval data. *Educational and Psychological Measurement, 30* (1), 61–70. https://doi.org/10.1177/001316447003000105
- Krippendorff, K. (2004). *Content analysis: An introduction to its methodology*. Thousand Oaks, California: Sage.
- Landis, J. R., & Koch, G. G. (1977). The measurement of observer agreement for categorical data. *Biometrics*, 159-174. https://doi.org/10.2307/2529310

### R Packages
---
- psych
- irr
