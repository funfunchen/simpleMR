#'Subset the file information
#'
#'Subset the file to get summary statistics information
#'
#' @param file_to_rd the summary statistics file
#'
#' @return a list contains \code{beta}, \code{se};
#' or the \code{expo.beta}, \code{expo.se} of exposure trait in exposure ascertainment data set
#' the \code{out.beta}, \code{out.se}  of outcome trait in exposure ascertainment data set
#'
#' @export

rd_info <- function(file_to_rd){
  x.beta <- file_to_rd[, 3] %>% unlist() %>% unname();
  x.se <- file_to_rd[, 5] %>% unlist() %>% unname();
  if(length(file_to_rd) < 6){
    return(res <- list(beta = x.beta,
                       se = x.se));
  }
  y.beta <- file_to_rd[, 6] %>% unlist() %>% unname();
  y.se <- file_to_rd[, 8] %>% unlist() %>% unname();
  return(res <- list(expo.beta = x.beta,
                     expo.se = x.se,
                     out.beta = y.beta,
                     out.se = y.se));
}
