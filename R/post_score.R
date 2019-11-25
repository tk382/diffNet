#' Returns new score function that was corrected for small sample size
#'
#' @param score score statistic before correction
#' @param coef cubic coefficients from cubic_coeff_c function
#'
#' @export
post_score = function(score, coef){
  roots = polyroot(c(-score, coef))
  return(Re(roots)[abs(Im(roots)) < 1e-6][1])
}
