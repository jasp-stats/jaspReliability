Intraclass Correlation
==========================

The intraclass correlation (ICC) is a correlation coefficient, that assesses the consistency between measures of the same class (Field, 2012). It is a common measure of inter-rater reliability and consistency.

Shrout & Fleiss (1979) distinguish between 6 different ways of calculating the ICC, based on different designs. These are identified by different values for A and B in *ICC(A, B)*. *A* is dependent on how raters are selected and subjects are rated and *B* on whether or not ratings are averaged in the end.

*Rater* can refer to different judges, tests or other forms of rating here and in the section below.

### Input
-------

- Variables: The variables / columns to compute the ICC for. Each variable corresponds to one rater, with different rows corresponding to different subjects being rated.
- Each subject is rated by...
  - a different rater (randomly selected): Each subject is rated by a different rater. Raters are selected randomly. Corresponds to *ICC(1, B)* in accordance with the conventions by Shrout & Fleiss (1979).
  - the same set of randomly selected raters/tests: A random sample of *k* raters rate all subjects. Corresponds to *ICC(2, B)* in accordance with the conventions by Shrout & Fleiss (1979).
  - the same fixed set of raters/tests: A fixed set of raters rate all subjects. These are all the raters of interest and there is no generalization to a wider population. Corresponds to *ICC(3, B)* in accordance with the conventions by Shrout & Fleiss (1979).
- Ratings are averaged: Are ratings by all raters averaged in the end? This strongly impacts the ICC coefficient. If yes, it corresponds to an *ICC(A, k)* in accordance with the conventions by Shrout & Fleiss (1979), if not it corresponds to an *ICC(A, 1)*.
- Confidence Interval: Whether a confidence interval should be reported for the ICC and the size of the interval.
- Bland-Altman plot (and table):
  - Displays a scatterplot of the means and differences between the two measurement variables. Additionaly, a table is created displaying the values regarding the mean difference and its respective limits.
    - Plot components: The variables included in the analysis.
    - Measurement Pair: At least two variables have to be selected (pair)
    - Confidence Interval: Coverage of the confidence intervals in percentages. By default, the confidence interval is set to 95%. This can be changed into the desired percentage. The confidence interval is drawn around the mean difference and its upper and lower limits. Calculations for confidence intervals are based on Bland & Altman (1986).
	- Shading: Highlights the confidence bounds around the mean difference and its limits.
  

### Output
-------

#### Intraclass Correlation (table)
This table shows the type of ICC that was computed and the value of the ICC coefficient's estimate. If selected, a confidence interval for the estimate will be reported.

Cicchetti (1994) provides the following guideline for interpretation of the ICC coefficient:
- Less than 0.40: poor
- Between 0.40 and 0.59: fair
- Between 0.60 and 0.74: good
- Between 0.75 and 1.00: excellent

Koo and Li (2016) provide a newer, slightly more conservative guideline for interpretation of the ICC coefficient:

- below 0.50: poor
- between 0.50 and 0.75: moderate
- between 0.75 and 0.90: good
- above 0.90: excellent

Bland-Altman plot (and table):
- Displays a scatterplot of the means and differences between the two measurement variables. The x-axis represents the means of the two measurements, the y-axis represents the differences between the two measurements.
- By default, three lines (dashed) are added to the plot, indicating the mean difference between the two variables and its upper and lower limits. Limits are based on the critical difference which is calculated using the standard deviation 1.96.
- By enabling "Confidence Interval", the x% confidence intervals for the mean difference, the upper limit, and the lower limit are added to the plot (dashed). Enabling "Shading" highlights these bounds.
- By default, a table containing the above information is created. The table displays the point value, lower and upper confidence interval numbers for the mean difference and its limits, respectively.

### References
-------
- Field, A. P., Miles, J., & Field, Z. (2012). *Discovering statistics using R*. Sage.
- Revelle, W. (2021) psych: Procedures for Personality and Psychological Research, Northwestern University, Evanston, Illinois, USA, https://CRAN.R-project.org/package=psych Version = 2.1.9.
- Shrout, P. E., & Fleiss, J. L. (1979). Intraclass correlations: uses in assessing rater reliability. *Psychological Bulletin, 86*(2), 420. https://doi.org/10.1037/0033-2909.86.2.420
- Cicchetti, D. V. (1994). Guidelines, criteria, and rules of thumb for evaluating normed and standardized assessment instruments in psychology. *Psychological Assessment, 6*(4), 284. https://doi.org/10.1037/1040-3590.6.4.284
- Koo, T. K., & Li, M. Y. (2016). A guideline of selecting and reporting intraclass correlation coefficients for reliability research. *Journal of Chiropractic Medicine, 15*(2), 155-163. https://doi.org/10.1016/j.jcm.2016.02.012
- Bland, J. M., & Altman, D. (1986). Statistical methods for assessing agreement between two methods of clinical measurement. *The lancet, 327*(8476), 307-310. https://doi.org/10.1016/S0140-6736(86)90837-8

### R Packages
---
- psych
