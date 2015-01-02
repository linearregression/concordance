
#' Calculate the concordance between two data sets.
#'
#' The \code{concordance} function takes two data sets (either matrices
#' or data.frames) and calculates the concordance between the two.
#'
#' @param x1 the data set to compare to a reference (normalizing) data set
#' @param x2 the reference data set
#' @param eigen_basis should \code{x1} and \code{x2} be projected onto 
#' \code{x2}'s basis? Default is \code{FALSE}.
#' @return the concordance value between the two data sets
#' @examples
#' concordance(iris[iris$Species=="setosa",c("Sepal.Length", "Sepal.Width")],
#'             iris[iris$Species=="virginica",c("Sepal.Length","Sepal.Width")]
#' @export
concordance = function(x1, x2, eigen_basis=FALSE) {
  UseMethod("concordance", x1)
}

#' @export
concordance.matrix = function(x1, x2, eigen_basis=FALSE) {
  if (!(inherits(x2, "matrix") || inherits(x2, "Matrix")))
    stop("Second argument must be a matrix (or Matrix)")
  x1_scatter = crossprod(x1)
  x2_scatter = crossprod(x2)
  ret = 0
  if (eigen_basis) {
    # Project x1 onto x2's eigen basis.
    e = eigen(x2_scatter)
    ret = nrow(x2) / nrow(x1) / ncol(x1) * 
      sum(diag(t(e$vectors) %*% x1_scatter %*% e$vectors)/ e$values)
  } else {
    ret = nrow(x2) / nrow(x1) /ncol(x1) * 
      sum(diag(x1_scatter %*% solve(x2_scatter)))
  }
  ret
}

#' @export
concordance.Matrix = function(x1, x2, eigen_basis=FALSE) {
  concordance.matrix(x1, x2, eigen_basis)
}

#' @export
concordance.data.frame = function(x1, x2, eigen_basis=FALSE) {
  if (!setequal(names(x1), names(x2)))
    stop("Arguments do not have the same column names.")
  form = as.formula(paste("~", paste(names(x1), collapse="+")))
  x1m = model.matrix(form, x1)
  x2m = model.matrix(form, x2)
  concordance(x1m, x2m, eigen_basis)
}

