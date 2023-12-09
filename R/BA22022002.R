#' @importFrom stats cov sd
#' @useDynLib BA22022002, .registration = TRUE
NULL


#' @title Pearson correlation
#' @description Calculate pearson correlation between xs and ys
#' @param xs the first vector
#' @param ys the second vector
#' @return Pearson correlation between xs and ys
#' @examples
#' \dontrun{
#' xs <- runif(100)
#' ys <- runif(100)
#' res <- pearson.r(xs, ys)
#' }
#' @export
pearson.r <- function(xs, ys) cov(xs, ys) / sd(xs) / sd(ys)



#' @title Spearman correlation
#' @description Calculate spearman correlation between xs and ys
#' @param xs the first vector
#' @param ys the second vector
#' @return Spearman correlation between xs and ys
#' @examples
#' \dontrun{
#' xs <- runif(100)
#' ys <- runif(100)
#' res <- spearman.r(xs, ys)
#' }
#' @export
spearman.r <- function(xs, ys) {
  xs.rank <- rankdata(xs)
  ys.rank <- rankdata(ys)
  pearson.r(xs.rank, ys.rank)
}