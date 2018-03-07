#' Fisher’s Z-transformation
#'
#' Compute the Fisher’s Z-transformation of each data set
#'
#' x vectors of effect sizes on Traint X
#'
#'@param  y vectors of effect sizes on Traint Y
#'
#'@return  A list with  correlation coefficient \code{rho} and  Fisher z-transformation \code{z}
#'
#'
z_trans <- function(x, y){

  rho<- cor(x, y, method = "spearman", use = "pairwise.complete.obs");
  z <- 0.5 * (log(1+rho) - log(1-rho));
  output <- list(rho=rho, z=z);
  return(output);

}


#' approximate likelyhood function
#'
#' approximate likelihood for the two correlation coefficients
#'
#' @param   x_hat Fisher z-transformation of trait X data set
#' @param  x  causal coefficients of trait X data set
#' @param  n_x length of trait X data set
#' @param  y_hat Fisher z-transformation of trait Y data set
#' @param  y causal coefficients of trait Y data set
#' @param  n_y length of trait Y data set
#'
#' @return likelyhood
#'
#'
lh_M <- function(x_hat, x, n_x, y_hat, y, n_y){

  dnorm(x_hat, x, sqrt(1/(n_x - 3)))*dnorm(y_hat, y, sqrt(1/(n_y - 3)));

};



#' comparing bi-direction causal models using AIC
#'
#' calculate the Akaike information criterion (AIC) for four causal inference models;
#' implementation of the Pickrell's method
#'
#' @param  X.set Trait X ascertainment data set
#' @param  Y.set Trait Y ascertainment data set
#'
#' @return  relative likelihood of the best non-causal model compared to the best causal model.
#'
#' @export
bi_mr_aic <- function(X.set, Y.set){
  n_x <- rd_info(X.set)$expo.n;
  n_y <- rd_info(Y.set)$out.n;
  if(n_x <= 3 | n_y <= 3){
    warning("Not enough data for AIC caculation");
    return(NA);
  }
  # X.set ascertainment data
  beta.xx <- rd_info(X.set)$expo.beta;
  beta.xy <- rd_info(X.set)$out.beta;
  sebeta.xx <- rd_info(X.set)$expo.se;
  sebeta.xy <- rd_info(X.set)$out.se;

  # Y.set ascertainment data
  beta.yx <- rd_info(Y.set)$out.beta;
  beta.yy <- rd_info(Y.set)$expo.beta;
  sebeta.yx <- rd_info(Y.set)$out.beta;
  sebeta.yy <- rd_info(Y.set)$expo.beta;

  #transformation
  z_x <- z_trans(beta.xx, beta.xy);
  z_x_h <- z_x$z;
  x_rho <- z_x$rho;
  z_y <- z_trans(beta.yx, beta.yy);
  z_y_h <- z_y$z;
  y_rho <- z_y$rho;

  ### calculate AIC-r
  lh_M1 <- lh_M(z_x_h, z_x_h, n_x, z_y_h, 0, n_y); #M1: if trait X causes Y, then we estimate ZX and set z_y_h = 0.
  lh_M2 <- lh_M(z_x_h, 0, n_x, z_y_h, z_y_h, n_y); # M2: If trait Y causes X, then we estimate ZY and set z_x_h = 0.
  lh_M3 <- lh_M(z_x_h, 0, n_x, z_y_h, 0, n_y); # M3: If there are no relationships between the traits, then z_x_h = z_y_h = 0.
  AIC_score_M1 <- 2-2*log(lh_M1);
  AIC_score_M2 <- 2-2*log(lh_M2);
  AIC_score_M3 <- 0-2*log(lh_M3);
  lh_fun_M4 <- function(z, x_h, y_h, n.x, n.y){
    -1*dnorm(x_h, z, sqrt(1/(n.x - 3)))*dnorm(y_h, z, sqrt(1/(n.y - 3)))
  };
  z_m4 <- nlminb(0.1, lh_fun_M4, x_h=z_x_h, y_h=z_y_h, n.x=n_x, n.y=n_y, gradient = NULL)$par;
  lh_M4 <- lh_M(z_x_h, z_m4, n_x, z_y_h, z_m4, n_y); #M4: If the correlation does not depend on how the variants were ascertained, z_x_h = z_y_h  .
  AIC_score_M4 <- 2-2*log(lh_M4);
  r <- exp((min(AIC_score_M1, AIC_score_M2) - min(AIC_score_M3, AIC_score_M4))/2);
}
