#' Implementation for all Mendelian randomization tests
#'
#' @param dat Harmonised exposure and outcome data. Output from \code{harmonise_exposure_outcome}
#' @param parameters Parameters to be used for various MR methods. Default is output from \code{dafault_param}.
#' @param method_list List of methods to use in analysis. See \code{mr_method_list()} for details.
#'
#' @export
#' @return List with the following elements:
#'         mr: Table of MR results
#'         extra: Table of extra results
mr <- function(dat, parameters=default_parameters(), method_list=subset(mr_method_list(), use_by_default==TRUE)$obj){
  
  mr_tab <- plyr::ddply(dat, c("expo.id", "out.id"), function(x) {
    if(nrow(x) == 0)
    {
      message("No SNPs available for MR analysis of '", x$expo.id[1], "' on '", x$out.id[1], "'")
      return(NULL)
    } else {
      message("Analysing '", x$expo.id[1], "' on '", x$out.id[1], "'")
    }
    
    res <- lapply(method_list, function(meth)
    {
      get(meth)(x$beta.exposure, x$beta.outcome, x$se.exposure, x$se.outcome, parameters)
    })
    
    methl <- mr_method_list()
    
    mr_tab <- data.frame(
      # outcome = x$out.id[1],
      # exposure = x$expo.id[1],
      method = methl$name[match(method_list, methl$obj)],
      nsnp = sapply(res, function(x) x$nsnp),
      b = sapply(res, function(x) x$b),
      se = sapply(res, function(x) x$se),
      pval = sapply(res, function(x) x$pval)
    )
    
    mr_tab <- subset(mr_tab, !(is.na(b) & is.na(se) & is.na(pval)))
    return(mr_tab)
  })
  return(mr_tab)
}

mr_method_list <- function(){
  a <- list(
    list(
      obj="mr_ivw",
      name="Inverse variance weighted",
      Description="",
      use_by_default=TRUE,
      heterogeneity_test=TRUE
    ),
    list(
      obj="mr_egger_regression",
      name="MR Egger",
      Description="",
      use_by_default=TRUE,
      heterogeneity_test=TRUE
    ),
    list(
      obj="mr_wald_ratio",
      name="Wald ratio",
      PubmedID="",
      Description="",
      use_by_default=TRUE,
      heterogeneity_test=FALSE
    ),
    list(
      obj="mr_simple_median",
      name="Simple median",
      PubmedID="",
      Description="",
      use_by_default=FALSE,
      heterogeneity_test=FALSE
    ),
    list(
      obj="mr_weighted_median",
      name="Weighted median",
      PubmedID="",
      Description="",
      use_by_default=TRUE,
      heterogeneity_test=TRUE
    ),		
    list(
      obj="mr_weighted_mode",
      name="Weighted mode",
      PubmedID="",
      Description="",
      use_by_default=TRUE,
      heterogeneity_test=TRUE
    ),
    list(
      obj="mr_weighted_mode_nome",
      name="Weighted mode (NOME)",
      PubmedID="",
      Description="",
      use_by_default=FALSE,
      heterogeneity_test=FALSE
    ),
    list(
      obj="mr_ivw_base",
      name="Inverse variance weighted(mrbase)",
      PubmedID="",
      Description="",
      use_by_default=TRUE,
      heterogeneity_test=TRUE
    )
  )
  a <- lapply(a, as.data.frame)
  a <- plyr::rbind.fill(a)
  a <- as.data.frame(lapply(a, as.character), stringsAsFactors=FALSE)
  a$heterogeneity_test <- as.logical(a$heterogeneity_test)
  a$use_by_default <- as.logical(a$use_by_default)
  return(a)
}

default_parameters <- function()
{
  list(
    nboot = 1000,
    Cov = 0,
    penk = 20,
    phi = 1,
    alpha = 0.05,
    Qthresh = 0.05,
    over.dispersion = TRUE,
    loss.function = "huber"
  )
}

mr_ivw <- function(b_exp, b_out, se_exp, se_out, parameters=default_parameters())
{
  if(sum(!is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) & !is.na(se_out)) < 2)
    return(list(b=NA, se=NA, pval=NA, nsnp=NA))
  
  mod <- lm(b_out ~ -1 + b_exp, weights = 1/se_out^2)
  ivw.res <- summary(mod)
  b <- ivw.res$coef["b_exp","Estimate"] # res$coef[1,1]
  se <- ivw.res$coef["b_exp","Std. Error"] # res$coef[1,2]
  pval <- ivw.res$coef["b_exp", "Pr(>|t|)"] # res$coef[1,4]
  Q_df <- length(b_exp) - 1
  Q <- ivw.res$sigma^2 * Q_df
  Q_pval <- pchisq(Q, Q_df, low=FALSE)
  cooks_d <- cooks.distance(mod)
  # from formula phi =  Q/DF rearranged to to Q = phi*DF, where phi is sigma^2
  # Q.ivw<-sum((1/(se_out/b_exp)^2)*(b_out/b_exp-ivw.reg.beta)^2)
  return(list(b = b, se = se, pval = pval, nsnp=length(b_exp), Q = Q, Q_df = Q_df, Q_pval = Q_pval, 
              cooks_d = cooks_d))
}

mr_egger_regression <- function(b_exp, b_out, se_exp, se_out, parameters){
  
  stopifnot(length(b_exp) == length(b_out))
  stopifnot(length(se_exp) == length(se_out))
  stopifnot(length(b_exp) == length(se_out))
  
  # print(b_exp)
  
  nulllist <- list(
    b = NA,
    se = NA,
    pval = NA,
    nsnp = NA,
    b_i = NA,
    se_i = NA,
    pval_i = NA,
    Q = NA,
    Q_df = NA,
    Q_pval = NA,
    mod = NA,
    smod = NA,
    dat = NA
  )
  if(sum(!is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) & !is.na(se_out)) < 3)
  {
    return(nulllist)
  }
  
  mod <- lm(b_out ~ b_exp, weights=1/se_out^2)
  smod <- summary(mod)
  if(nrow(coefficients(smod)) > 1)
  {
    b <- coefficients(smod)[2,1]
    se <- coefficients(smod)[2,2]
    pval <- coefficients(smod)[2,4]
    b_i <- coefficients(smod)[1,1]
    se_i <- coefficients(smod)[1,2] 
    pval_i <- coefficients(smod)[2,4]
    
    Q <- smod$sigma^2 * (length(b_exp) - 2)
    Q_df <- length(b_exp) - 2
    Q_pval <- pchisq(Q, Q_df, low=FALSE)
  } else {
    warning("Collinearities in MR Egger, try LD pruning the exposure variables.")
    return(nulllist)
  }
  return(list(b = b, se = se, pval = pval, nsnp = length(b_exp), 
              b_i = b_i, se_i = se_i, pval_i = pval_i, Q = Q, Q_df = Q_df, Q_pval = Q_pval, mod = smod))
}

mr_wald_ratio <- function(b_exp, b_out, se_exp, se_out, parameters)
{
  if(length(b_exp) > 1)
  {
    return(list(b=NA, se=NA, pval=NA, nsnp=NA))
  }
  b <- b_out / b_exp
  se <- se_out / abs(b_exp)
  # sqrt((segd^2/gp^2) + (gd^2/gp^4)*segp^2 - 2*(gd/gp^3)) #full delta method with cov set to 0
  pval <- pnorm(abs(b) / se, lower.tail=FALSE) * 2
  return(list(b=b, se=se, pval=pval, nsnp=1))
}


## Perform bootstraps 
weighted_median_bootstrap <- function(b_exp, b_out, se_exp, se_out, weights, nboot)
{
  med <- rep(0, nboot)
  for(i in 1:nboot){
    b_exp.boot = rnorm(length(b_exp), mean=b_exp, sd=se_exp)
    b_out.boot = rnorm(length(b_out), mean=b_out, sd=se_out)
    betaIV.boot = b_out.boot/b_exp.boot
    med[i] = weighted_median(betaIV.boot, weights)
  }
  return(sd(med))
}

weighted_median <- function(b_iv, weights)
{
  betaIV.order <- b_iv[order(b_iv)]
  weights.order <- weights[order(b_iv)]
  weights.sum <- cumsum(weights.order)-0.5*weights.order
  weights.sum <- weights.sum/sum(weights.order)
  below <- max(which(weights.sum<0.5))
  b = betaIV.order[below] + (betaIV.order[below+1]-betaIV.order[below])*
    (0.5-weights.sum[below])/(weights.sum[below+1]-weights.sum[below])
  return(b)
}

mr_simple_median <- function(b_exp, b_out, se_exp, se_out, parameters=default_parameters())
{
  if(sum(!is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) & !is.na(se_out)) < 3)
    return(list(b=NA, se=NA, pval=NA, nsnp=NA))
  
  b_iv <- b_out / b_exp
  b <- weighted_median(b_iv, rep(1/length(b_exp), length(b_exp)))
  se <- weighted_median_bootstrap(b_exp, b_out, se_exp, se_out, rep(1/length(b_exp), length(b_exp)), parameters$nboot)
  pval <- 2 * pnorm(abs(b/se), low=FALSE)
  return(list(b=b, se=se, pval=pval, nsnp=length(b_exp)))
}

mr_weighted_median <- function(b_exp, b_out, se_exp, se_out, parameters=default_parameters())
{
  if(sum(!is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) & !is.na(se_out)) < 3)
    return(list(b=NA, se=NA, pval=NA, nsnp=NA))
  
  b_iv <- b_out / b_exp
  VBj <- ((se_out)^2)/(b_exp)^2 + (b_out^2)*((se_exp^2))/(b_exp)^4
  b <- weighted_median(b_iv, 1 / VBj)
  se <- weighted_median_bootstrap(b_exp, b_out, se_exp, se_out, 1 / VBj, parameters$nboot)
  pval <- 2 * pnorm(abs(b/se), low=FALSE)
  return(list(b=b, se=se, pval=pval, Q=NA, Q_df=NA, Q_pval=NA, nsnp=length(b_exp)))
}

### MR Base ivw implement
mr_ivw_base <- function(b_exp, b_out, se_exp, se_out, parameters=default_parameters())
{
  if(sum(!is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) & !is.na(se_out)) < 2)
    return(list(b=NA, se=NA, pval=NA, nsnp=NA))
  
  ivw.res <- summary(lm(b_out ~ -1 + b_exp, weights = 1/se_out^2))
  b <- ivw.res$coef["b_exp","Estimate"]
  se <- ivw.res$coef["b_exp","Std. Error"]/min(1,ivw.res$sigma) #sigma is the residual standard error
  pval <- 2 * pnorm(abs(b/se), low=FALSE)
  Q_df <- length(b_exp) - 1
  Q <- ivw.res$sigma^2 * Q_df
  Q_pval <- pchisq(Q, Q_df, low=FALSE)
  # from formula phi =  Q/DF rearranged to to Q = phi*DF, where phi is sigma^2
  # Q.ivw<-sum((1/(se_out/b_exp)^2)*(b_out/b_exp-ivw.reg.beta)^2)
  return(list(b = b, se = se, pval = pval, nsnp=length(b_exp), Q = Q, Q_df = Q_df, Q_pval = Q_pval))
}
