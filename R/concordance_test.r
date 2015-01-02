
#' Test the difference of two data sets based on their concordance
#' 
#' @param x1 the data set to compare to a reference (normalizing) data set
#' @param x2 the reference data set
#' @param alternative "two.sided", "less", or "greater. Default is "two.sided".
#' @param model "Beta", "F", "Cauchy", or "all". The default is "all".
#' A non-parametric test is coming.
#' @param conf.level confidence level for the interval
#' @param eigen_basis should \code{x1} and \code{x2} be projected onto 
#' \code{x2}'s basis? Default is \code{FALSE}.
#' @export
concordance_diff_test = function(x1, x2, alternative="two.sided", model="all",
  conf.level=0.95, eigen_basis=FALSE) {
  if (model == "all") {
    ret = list( 
      beta_test=concordance_diff_test_beta(x1, x2, alternative, eigen_basis),
      F_test=concordance_diff_test_F(x1, x2, alternative, eigen_basis),
      cauchy_test=concordance_diff_test(x1, x2, alternative, eigen_basis))
    class(ret) = "htests"
  } else {
    test = switch(model,
                  "beta", concordance_diff_test_beta,
                  "F", concordance_diff_test_F,
                  "Cauchy", concordance_diff_test_Cauchy)
    ret = test(x1, x2, alternative, conf.level, eigen_basis)
  }
  ret 
}

concordance_diff_test_beta = function(x1, x2, alternative, conf.level,
                                      eigen_basis) {
  mc = match.call()
  
  # Get the data_name.
  data_name = paste(deparse(mc$x1), "and", deparse(mc$x2))

  # Get the estimates for the parameters of the beta.
  estimate = c(nrow(x2)/nrow(x1), nrow(x1)/2, (nrow(x2)-nrow(x1))/2)
  names(estimate) = c("coefficient", "shape1", "shape2")
  if (eigen_basis) {
    method = "Concordance Test (Eigen Basis, Beta Distribution)"
  } else {
    method = "Concordance Test (Beta Distribution)"
  }
  parameter=c("Beta", ifelse(eigen_basis, "Eigen Basis", "Trace"))
  names(parameter) = c("Distribution", "Concordance Type")
 
  # Get the concordance and the beta-scaled concordance.
  conc = concordance(x1, x2, eigen_basis)
  statistic=conc
  names(statistic) = "concordance"
  conc_normalized = conc * nrow(x1) / nrow(x2)
  
  # Get the p-value and confidence interval.
  if (alternative == "two.sided") {
    p.value = min(pbeta(conc_normalized, estimate['shape1'], estimate['shape2'],
                        lower.tail=FALSE),
                  pbeta(conc_normalized, estimate['shape1'], estimate['shape2'],
                        lower.tail=TRUE))
  } else if (alternative == "less") {
    p.value=pbeta(conc_normalized, estimate['shape1'], estimate['shape2'], 
                  lower.tail=TRUE)
  } else if (alternative == "greater") {
    p.value=pbeta(conc_normalized, estimate['shape1'], estimate['shape2'], 
                  lower.tail=FALSE)
  } else{
    stop('The alternative parameter must be "two.sided", "less", or "greater".')
  }

  # Get the confidence interval.
  a = estimate['shape1']
  b = estimate['shape2']
  alpha = (1-conf.level)/2
  conf_int = c(qbeta(alpha, estimate['shape1'], estimate['shape2'], 
                     lower.tail=TRUE) * nrow(x2) / nrow(x1),
               qbeta(alpha, estimate['shape1'], estimate['shape2'], 
                     lower.tail=FALSE) * nrow(x2) / nrow(x1))

  # Get the null.value
  null.value = 1
  names(null.value) = "concordance"
  ret = list(
    statistic=statistic,
    parameter=c(),
    p.value=p.value,
    conf.int=conf_int,
    estimate=estimate,
    null.value=null.value,
    alternative=alternative,
    method=method, 
    data.name=data_name
  )
  class(ret) = "htest"
  ret
}

concordance_diff_test_F = function(x1, x2, alternative, conf.level, 
                                   eigen_basis) {
  stop("Not implemented yet.")
}
concordance_diff_test_Cauchy = function(x1, x2, alternative, conf.level,
                                        eigen_basis) {
  stop("Not implemented yet.")
}

