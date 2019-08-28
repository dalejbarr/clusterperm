#' Visual-World Eyetracking Data from Kronmüller & Barr (2007)
#'
#' A dataset containing visual-world eyetracking data from 56 subjects over 35 time bins (from -200 to 1500 ms relative to expression onset). The data are from a 2x2x2 experiment, with factors Speaker (Same, Different), Precedent (Break, Maintain), and Load (Yes, No).  For details, see the citation below.
#'
#' @format A data frame with 15680 rows and 6 variables:
#' \describe{
#'   \item{SubjID}{Unique identifier for each subject}
#'   \item{Speaker}{Speaker condition (Same, Diff)}
#'   \item{Precedent}{Precedent (Break, Maintain)}
#'   \item{Load}{Cognitive Load (Yes, No)}
#'   \item{bin}{Time bin (-200 to 1500 ms in steps of 50 ms)}
#'   \item{TAS}{Target Advantage Score}
#' }
#' @source Kronmüller & Barr (2007), Perspective-free pragmatics: Broken precedents and the recovery-from-preemption hypothesis. Journal of Memory and Language, 56, 436-455.
"kb07bins"
