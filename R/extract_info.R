#'Subset the file information
#'
#'Subset the file to get summary statistics information
#'
#' @param file_to_rd the summary statistics file
#'
#' @return a list contains \code{beta}, \code{se}, \code{n};
#' or the \code{expo.beta}, \code{expo.se}, \code{expo.n} of exposure trait in exposure ascertainment data set
#' the \code{out.beta}, \code{out.se}, \code{out.n} of outcome trait in exposure ascertainment data set
#'
#' @export

rd_info <- function(file_to_rd){
  file_to_rd <- dplyr::filter(file_to_rd, file_to_rd[, 5]>0)
  x.beta <- file_to_rd[, 3] %>% unlist() %>% unname();
  x.se <- file_to_rd[, 5] %>% unlist() %>% unname();
  n_x <- length(x.beta);
  if(length(file_to_rd) < 6){
    return(res <- list(beta = x.beta,
                       se = x.se,
                       n = n_x))
  }

  file_to_rd <- dplyr::filter(file_to_rd, file_to_rd[, 8]>0)
  y.beta <- file_to_rd[, 6] %>% unlist() %>% unname();
  y.se <- file_to_rd[, 8] %>% unlist() %>% unname();
  n_y <- length(y.beta)
  return(res <- list(expo.beta = x.beta,
                     expo.se = x.se,
                     expo.n = n_x,
                     out.beta = y.beta,
                     out.se = y.se,
                     out.n = n_y))
}
