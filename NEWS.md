# jaspReliability Changelog

> **HOW TO READ AND UPDATE THIS CHANGELOG:**
> 
> This document follows a modified [Keep a Changelog](https://keepachangelog.com/) format adapted for the R/JASP ecosystem. Releases are listed in reverse chronological order (newest first).
> As an example see [jaspModuleTemplate](https://github.com/jasp-stats/jaspModuleTemplate/blob/master/NEWS.md)
> * **Adding New Changes (For Contributors):** All new commits should be logged at the very top of the file under the `# jaspModuleTemplate (development version)` header. Place your bullet point under the appropriate category (`## Added`, `## Fixed`, etc.). 
> * **Issue References:** Please reference the relevant GitHub Issue (if any) at the end of your line (e.g., `([Issue #19](https://github.com/jasp-stats/jaspModuleTemplate/issues/19)`). 
> * **Format Categories:** >   * **Added:** New template features, QML examples, or build tools.
>   * **Changed:** Updates to default configurations, boilerplate code, or dependencies. 
>   * **Fixed:** Bug fixes in the build pipeline, R wrappers, or QML layouts.
>   * **Deprecated / Removed:** Outdated template components or legacy code.


---

# jaspReliability (development version)

## Fixed
* Rater Agreement: declared factor level order is now authoritative for weighted Cohen's kappa and ordinal Krippendorff's alpha even when the labels look numeric (e.g. levels "1","3","2" declaring low < medium < high); numeric sorting is only used for columns with no declared order. Contradictory numeric-looking declared orders across raters now show the common-scale error instead of being silently sorted.
* Rater Agreement: pairwise weighted Cohen's kappa now uses the full declared common scale for every rater pair, so a pair that never observes an interior category still gets that category's correct distance/weight instead of a rescaled subset.
* Rater Agreement: Kendall's W (table and bootstrap CI) now rejects raters-in-rows data whose merged ordinal scale is ambiguous, matching the existing policy for weighted Cohen's kappa and ordinal Krippendorff's alpha; previously the result depended on the order of the subject columns.


---

# jaspReliability 0.97.1

## Added
* Rater Agreement: added Kendall's W coefficient with bootstrap CI support.
* Rater Agreement: added the F test for Kendall's W alongside the chi-square test ([Issue #2151](https://github.com/jasp-stats/jasp-issues/issues/2151)).

## Fixed
* Rater Agreement: declared ordinal level order is now respected (Kendall's W, Krippendorff's alpha, weighted Cohen's kappa), including with raters in rows and when raters use different subsets of the categories; validation runs on the analyzed (complete/pairwise) data; degenerate inputs (all-missing, constant ratings, single category, one rater) show clear errors instead of failing or reporting invalid results.
* Rater Agreement: Krippendorff's alpha is validated on the ratings that actually enter a coincidence (subjects rated at least twice); previously such data could report a non-estimable coefficient together with a bootstrap confidence interval.
* Rater Agreement: Fleiss' kappa standard errors are computed from the Fleiss, Nee & Landis (1979) formulas instead of being reconstructed from rounded output.

## Changed
* Rater Agreement: grouped bootstrap samples, CI level, and seed into an "Advanced Options" section; removed pre-checked defaults for all coefficients; added placeholder table prompting users to select a coefficient when variables are assigned.
* Rater Agreement: the tie correction for Kendall's W is enabled by default, and bootstrap confidence intervals are only available for the tie-corrected coefficient, because resampling with replacement introduces ties.
* Rater Agreement: weighted Cohen's kappa and ordinal Krippendorff's alpha require the raters to share one ordinal scale. When the declared category orders are contradictory or do not determine a unique order, an error is shown instead of an order-dependent result.
* Rater Agreement: Cohen's kappa validates and reports each rater pair on its own pairwise-complete data. Pairs with fewer than three jointly rated subjects, without variation, or with a non-estimable coefficient are excluded from the average kappa and listed in a note. Pair rows are ordered V1-V2, V1-V3, V1-V4, ...

---

# jaspModuleTemplate 0.2.0
## Added
* Added NEWS.md
* Added workflow to remind users to update their `NEWS.md`.
* Added workflow to auto-bump version when user does not do so.

---

# jaspModuleTemplate 0.1.0

## Added
* Initial examples to showcase JASP module development

## Changed
* Use best practices for checking input ([Issue #19](https://github.com/jasp-stats/jaspModuleTemplate/issues/19)).
* The main results table now defaults to displaying 95% Confidence Intervals for effect sizes.

## Fixed
* Remove deprecated dependencies from qml files ([Issue #14](https://github.com/jasp-stats/jaspModuleTemplate/issues/14)).
