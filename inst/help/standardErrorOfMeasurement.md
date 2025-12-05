
Standard Error of Measurement
===
The Standard Error of Measurement (sem) quantifies the precision of a measurement. One can go the traditional way and only estimate a single sem value for a test, or estimate sem values for each sum score of a test, that is, 

## Input
---
#### Variables Box
- Variables have to be ordinally or nominally scaled

#### Methods
- Split-test methods: based on CTT
  - Thorndike: Method that splits the test into two halves to estimate the error variances from the variances of the differences. The default split is one part contains the odd numbered, and the other part the even numbered items. If another split is desired rearrange the order of variables in the window.
  - Feldt: Method that splits the test into multiple halves with the same goal as the Thorndike method
    - Number of splits: How many splits to use for the Feldt method, can only be a divisor of the number of items. By default items are split into similar parts, for instance, with 9 items and 3 splits, item 1, 4, 7 are in part 1, item 2, 5, 8 are in part 2, and item 3, 6, 9 are in part 3.
  - Mollenkopf-Feldt: Method that splits in multiple test parts and predicts the differences with a polynomial regression
    - Number of splits: How many splits to apply, can only be a divisor of the number of items. Splits follow same logic as Feldt method.
    - Degree of polynomial: For instance, 3 means a regression like Y=X+X^2+X^3

- ANOVA: Method based on a repeated measures ANOVA, simplified by the Emons (2023) to use the ICC 3k
- IRT: Method based on item response theory, for dichomotously scored items this is the 2-parameter logistic model (2PLM), for polytomously scored items it is the graded response model (GRM). Both models assume a single underlying latent variable.

- Binomial methods: Based on the idea that the item scores follow a binomial distribution
  - Lord: Method that only requires number of correct and incorrect items
  - Keats: Corrects the Lord method for supposed bias and uses a reliability coefficient in the process
  - Lord generalized: Essentially the Lord method for multiple test parts
    - Number of splits: How many splits to apply, can only be a divisor of the number of items

### Options
- Sum score table: Displays a table with the sum scores and confidence intervals around them that differ per sem method. The CIs are normal-theory CIs that apply the sem in their calculation.
  - CI
- User defined reliability: Users can specify a custom reliability value to be used in the calculation of the score-unrelated sem and the calculation of the Keats method
- Minimum number of observations per score group: The sem methods except the Mollenkopf-Feldt and the IRT method require a minimum number of observations for a sum score to properly estimate the sem, if that number is not reached the sum score is merged with the subsequent sum score(s) until there are enough observations in the merged score group. 
- Hide sem table

### Plots
- Histogram of sum score counts: Displays a histogram with the count (y-axis) of each sum score (x-axis)
- Plot per method: Displays a point plot per method of the sems (y-axis) for each sum score (x-axis). Only for the IRT method the plot is a line.
- Combined plot: Displays the sems (y-axis) for all methods in one plot.
- Sum score plots: Displays the sum scores(x-axis) and the CIs around them per method. The y-axis denotes the hypothetical true score that would be covered by 95% of the CIs if they were repeatedly constructed in the same way.
  - CI
  - Display cutoff score: Displays a horizontal cutoff line

---
## Output

- Standard error of measurement: The table contains the sem values for the different sum scores of the different sem methods. For some methods (Thorndike, Feldt, ANOVA) the sem for several sum scores may be summarized in one value since the counts for those scores do not reach the minimum size which is 10 by default. 

- Sum score CI table: The table contains the confidence intervals (CI) around each sum score differing for each method. Specifically, the 95% CI for a given sum score contains the range of values that we assume cover the true score in 95% of times when repeatedly constructed in the same way. A CI is constructed using the sem and the z-score from a standard normal distribution, making the CIs normal-theory based intervals. 

- Histogram of counts per sum score group

- Plots: The point plots show the sem values for the sum scores for each method. For the IRT method the plot is a line-plot since IRT provides a continuous sum score. The combined plot combines all methods in one plot. 

- Sum Score CI Plots: For each sem method the plot displays the sum scores and the 95% confidence intervals around them. A red line indicates a potential cutoff. The CI around a sum score contains the corresponding true score 95% of times if repeatedly constructed in the same way. The IRT method shows a line with a ribbon since IRT provides a continous sum score.

---
## References
- Emons, W. H. M. (2023). Methods for estimating conditional standard errors of measurement and some critical reflections. In L. A. Van Der Ark, W. H. M. Emons, & R. R. Meijer (Eds.), *Essays on Contemporary Psychometrics* (pp. 195–216). Springer International Publishing. https://doi.org/10.1007/978-3-031-10370-4_11
- Feldt, L. S., & Qualls, A. L. (1996). Estimation of measurement error variance at specific score levels. *Journal of Educational Measurement, 33*(2), 141–156. https://doi.org/10.1111/j.1745-3984.1996.tb00486.x
- Feldt, L. S., Steffen, M., & Gupta, N. C. (1985). A comparison of five methods for estimating the standard error of measurement at specific score levels. *Applied Psychological Measurement, 9*(4), 351–361. https://doi.org/10.1177/014662168500900402
- Keats, J. A. (1957). Estimation of error variances of test scores. *Psychometrika, 22*(1), 29–41. https://doi.org/10.1007/BF02289207
- Keulers, E. H. H., & Hurks, P. P. M. (2021). Psychometric properties of a new ADHD screening questionnaire: Parent report on the (potential) underlying explanation of inattention in their school-aged children. *Child Neuropsychology, 27*(8), 1117–1132. https://doi.org/10.1080/09297049.2021.1937975
- Lord, F. M. (1955). Estimating test reliability. *Educational and Psychological Measurement, 15*(4), 325–336. https://doi.org/10.1177/001316445501500401
- Lord, F. M. (1965). A strong true-score theory, with applications. *Psychometrika, 30*(3), 239–270. https://doi.org/10.1007/BF02289490
- Samejima, F. (1968). Estimation of latent ability using a response pattern of graded scores. *Psychometrika, 34*(Suppl 1), 1–97. https://doi.org/10.1007/bf03372160
- Thorndike, R. L. (1951). Reliability. In E. F. Lindquist (Ed.), *Educational measurement* (pp. 560–620). American Council on Education.
