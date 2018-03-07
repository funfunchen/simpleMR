#'Inverse-variance weighted method for MR using the input data set
#'
#'Calculating the causal inference estimate(beta), se and p value using IVW
#'
#' @param data.set the input summary statistic data set
#'
#' @return a list with beta, se and p
#'
#' @export
#'
mr_ivw <- function(data.set){
  data.info <- rd_info(data.set);
  expo.beta <- data.info$expo.beta;
  out.se <- data.info$out.se;
  out.beta <- data.info$out.beta

  if(length(expo.beta) <2 ) {
    return(lsit(IVW.beta = NA,
                IVW.se = NA,
                IVW.p = NA))
  }

  mod.fit <- lm(out.beta ~ expo.beta.xx - 1, weights = out.se^-2)
  IVW.beta <- summary(mod.fit)$coef[1];
  IVW.se <- summary(mod.fit)$coef[2];
  IVW.p <- summary(mod.fit)$coef[4];
  return(lsit(IVW.beta = IVW.beta,
              IVW.se = IVW.se,
              IVW.p = IVW.p))
}
