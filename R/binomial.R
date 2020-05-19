#' Calculate binomial distribution P-values between all columns in a matrix
#'
#' Calculate P-values from the binomial distribution and take their negative
#' logarithms as an indicator of coexpression, as proposed by Mohammadi et al.
#'
#' @param mat a matrix of data, with samples in rows and features in columns
#' @return the binomial P-values between non-zero/missing values in each pair of
#'   columns
#'
#' @importFrom slam as.simple_triplet_matrix crossprod_simple_triplet_matrix
#' @importFrom Matrix rowMeans rowSums
#' @export
binomial = function(mat) {
  present = !is.na(mat) & mat > 0
  size = nrow(mat)
  q = present %>% as.simple_triplet_matrix() %>% crossprod_simple_triplet_matrix()
  prob = t(colMeans(present)) %>% as.simple_triplet_matrix() %>% crossprod_simple_triplet_matrix()
  cor = pbinom(q, size, prob, lower.tail = F)
  dimnames(cor) = list(colnames(mat), colnames(mat))
  -log10(cor)
}
