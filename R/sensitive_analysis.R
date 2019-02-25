mr_leaveoneout <- function(dat, parameters=default_parameters(), method=mr_ivw)
{
  if(!"samplesize.outcome" %in% names(dat))
  {
    dat$samplesize.outcome <- NA
  }
  
  stopifnot("out.id" %in% names(dat))
  stopifnot("expo.id" %in% names(dat))
  stopifnot("beta.exposure" %in% names(dat))
  stopifnot("beta.outcome" %in% names(dat))
  stopifnot("se.exposure" %in% names(dat))
  stopifnot("se.outcome" %in% names(dat))
  
  
  res <- plyr::ddply(dat, c("expo.id", "out.id"), function(x)
  {
    nsnp <- nrow(x)
    if(nsnp == 0)
    {
      x <- x[1,]
      d <- data.frame(
        SNP = "All",
        b = NA,
        se = NA,
        p = NA,
        samplesize = NA,
        outcome = x$out.id[1],
        exposure = x$expo.id[1]
      )
      return(d)
    }
    if(nsnp > 2)
    {
      l <- lapply(1:nsnp, function(i)
      {
        with(x, method(beta.exposure[-i], beta.outcome[-i], se.exposure[-i], se.outcome[-i], parameters))
      })
      l[[nsnp+1]] <- with(x, method(beta.exposure, beta.outcome, se.exposure, se.outcome, parameters))
      d <- data.frame(
        SNP = c(as.character(x$SNP), "All"),
        b = sapply(l, function(y) y$b),
        se = sapply(l, function(y) y$se),
        p = sapply(l, function(y) y$pval),
        samplesize = x$samplesize.outcome[1]
      )
      d$outcome <- x$out.id[1]
      d$exposure <- x$expo.id[1]
      
    } else {
      a <- with(x, method(beta.exposure, beta.outcome, se.exposure, se.outcome, parameters))
      d <- data.frame(
        SNP = "All",
        b = a$b,
        se = a$se,
        p = a$pval,
        samplesize = x$samplesize.outcome[1]
      )
      d$outcome <- x$out.id[1]
      d$exposure <- x$expo.id[1]
    }
    return(d)
  })
  res <- subset(res, select=c(exposure, outcome, samplesize, SNP, b, se, p))
  return(res)
}



mr_singlesnp <- function(dat, parameters=default_parameters(), single_method="mr_wald_ratio", all_method=c("mr_ivw", "mr_egger_regression"))
{
  
  if(!"samplesize.outcome" %in% names(dat))
  {
    dat$samplesize.outcome <- NA
  }
  
  stopifnot("out.id" %in% names(dat))
  stopifnot("expo.id" %in% names(dat))
  stopifnot("beta.exposure" %in% names(dat))
  stopifnot("beta.outcome" %in% names(dat))
  stopifnot("se.exposure" %in% names(dat))
  stopifnot("se.outcome" %in% names(dat))
  
  res <- plyr::ddply(dat, c("expo.id", "out.id"), function(x)
  {
    nsnp <- nrow(x)
    if(nsnp == 0)
    {
      x <- x[1,]
      d <- data.frame(
        SNP = "No available data",
        b = NA,
        se = NA,
        p = NA,
        samplesize = NA,
        outcome = x$out.id[1],
        exposure = x$expo.id[1]
      )
      return(d)
    }
    l <- lapply(1:nsnp, function(i)
    {
      with(x, get(single_method)(beta.exposure[i], beta.outcome[i], se.exposure[i], se.outcome[i], parameters))
    })
    nom <- c()
    for(i in 1:length(all_method))
    {
      l[[nsnp+i]] <- with(x, get(all_method[i])(beta.exposure, beta.outcome, se.exposure, se.outcome, parameters))
      
      nom <- c(nom, paste0("All - ", subset(mr_method_list(), obj==all_method[i])$name))
    }
    
    d <- data.frame(
      SNP = c(as.character(x$SNP), nom),
      b = sapply(l, function(y) y$b),
      se = sapply(l, function(y) y$se),
      p = sapply(l, function(y) y$pval),
      samplesize = x$samplesize.outcome[1]
    )
    d$outcome <- x$out.id[1]
    d$exposure <- x$expo.id[1]
    return(d)
  })
  res <- subset(res, select=c(exposure, outcome, samplesize, SNP, b, se, p))
  return(res)
}


mr_cook_distance <- function(dat, parameters=default_parameters(), method=mr_ivw)
{
  if(!"samplesize.outcome" %in% names(dat))
  {
    dat$samplesize.outcome <- NA
  }
  
  stopifnot("out.id" %in% names(dat))
  stopifnot("expo.id" %in% names(dat))
  stopifnot("beta.exposure" %in% names(dat))
  stopifnot("beta.outcome" %in% names(dat))
  stopifnot("se.exposure" %in% names(dat))
  stopifnot("se.outcome" %in% names(dat))
  
  nsnp <- nrow(dat)
  if(nsnp == 0)
  {
    d <- NA
    return(d)
  }
  else {
    a <- with(dat, method(beta.exposure, beta.outcome, se.exposure, se.outcome, parameters))
    d <- data.frame(
      SNP = as.character(dat$SNP),
      distance = a$cooks_d) 
    return(d)
  }
}

