#' Fisher’s Z-transformation
#' 
#' Compute the Fisher’s Z-transformation of each data set
#' 
#' @param x vectors of effect sizes on Traint X 
#' @param y vectors of effect sizes on Traint Y
#' 
#' @return A list with  correlation coefficient \code{rho} and  Fisher z-transformation \code{z} 
#' 
#' @export
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
#' @param  x_hat Fisher z-transformation of trait X data set
#' @param x  causal coefficients of trait X data set
#' @param n_x length of trait X data set
#' @param y_hat Fisher z-transformation of trait Y data set
#' @param y causal coefficients of trait Y data set
#' @param n_y length of trait Y data set
#' 
#' @return   likelyhood
#' 
#' @export
lh_M <- function(x_hat, x, n_x, y_hat, y, n_y){
  
  dnorm(x_hat, x, sqrt(1/(n_x - 3)))*dnorm(y_hat, y, sqrt(1/(n_y - 3)));
  
};



#' comparing causal models using AIC
#' 
#' calculate the Akaike information criterion (AIC) for four causal inference models
#' 
#'  @param X.set Trait X ascertainment data set
#'  @param Y.set Trait X ascertainment data set
#'  
#'  @return relative likelihood of the best non-causal model compared to the best causal model.
#'    
pickrell_aic <- function(X.set, Y.set){
  n_x <- nrow(X.set);
  n_y <- nrow(Y.set);
  if(n_x <= 3 | n_y <= 3){
    warning("Not enough data for AIC caculation");
    return(NA);
  } else{
    message("Analysing ACI between '", X.set$expo.id[1], "' with '", Y.set$expo.id[1], "'")
  }
  X.set <- X.set[, c(2:5)]
  colnames(X.set) <- c("X_x_beta", "X_x_se", "X_y_beta", "X_y_se")
  
  Y.set <- Y.set[, c(2:5)]
  colnames(Y.set) <- c("Y_x_beta", "Y_x_se", "Y_y_beta", "Y_y_se")
  
  # X.set ascertainment data
  beta.xx <- X.set$X_x_beta;
  beta.xy <- X.set$X_y_beta;
  sebeta.xx <- X.set$X_x_se;
  sebeta.xy <- X.set$X_y_se;
  
  # Y.set ascertainment data
  beta.yx <- Y.set$Y_x_beta;
  beta.yy <- Y.set$Y_y_beta;
  sebeta.yx <- Y.set$Y_x_se;
  sebeta.yy <- Y.set$Y_y_se;
  
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