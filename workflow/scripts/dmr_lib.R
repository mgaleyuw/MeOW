 coordTLong <- function(df, test.sample, var.equal, correction, minPositions=100, bootstraps=25, population="", paired=TRUE){
  robust.CpGs <- unique(df %>% group_by(cpg) %>% filter(n() > minPositions) %>% pull(cpg))
  small.df <- df %>% filter(cpg %in% robust.CpGs)
  
  results.df <- small.df %>% group_by(cpg) %>% dplyr::group_modify( ~ runTTest(.x, test.sample, var.equal, bootstraps, population, paired=paired,  length(robust.CpGs)))
  
  #results.df["bhpvalue"] <- p.adjust(results.df %>% pull(pvalue), method = correction)
  return(results.df)
}

runTTest <- function(df, test.sample, var.equal, bootstraps=25, population="", paired=TRUE, correction="BH"){
  if(is.null(match.call()$population)){
    control.samples <- setdiff(colnames(df), c("pos", "chr", test.sample, "cpg", "gene"))
  } else {
    control.samples <- population
  }
  n.positions <- dim(df)[1]
  # control.mat stores the reference to each sample (column) that should be used to pull methylation from.
  control.mat <- matrix(sample(1:length(control.samples), size=n.positions*bootstraps, replace=TRUE), nrow=bootstraps)
  corrective.add <- seq(0,(n.positions*length(control.samples)-1), length(control.samples))
  corrective.control.mat <- t(control.mat) + corrective.add
  meth.mat <- t(df[,control.samples])
  
  #get string of test values that will remain stable
  testvals <- df %>% pull(all_of(test.sample))
  
  # make a container for the results that is #stats X #bootstraps
  results.mat <- matrix(0,3,bootstraps)
  
  # run test for each bootstrap column ( meth.mat[corrective.control.mat[,i]]) and store p, cohen, and diff
  for(i in 1:bootstraps){
    cont.vec <- meth.mat[corrective.control.mat[,i]]
   
    tdata <- t.test(x=testvals, y=cont.vec, var.equal=var.equal, paired=paired)
    results.mat[1,i] <- tdata$p.value
    #print(tdata$parameter)
    pooled.var <- sd(c(cont.vec, testvals))
    mean.diff <- tdata$estimate[1]
    results.mat[2,i] <- mean.diff/pooled.var
    results.mat[3,i] <- mean.diff
  } 
  meanResults <- rowMeans(results.mat)
  corrected.p <- p.adjust(meanResults[1], method="BH", n=correction*bootstraps)
  #print(results.mat)
  return(tibble(pvalue=corrected.p, p_sd=sd(results.mat[1,]), effect.size=meanResults[2], effect_sd=sd(results.mat[2,]), diff=meanResults[3], diff_sd=sd(results.mat[3,]), n=n.positions))
}

runBetaSub <- function(df, test.sample, bootstraps=25, subsample=100, fit.model=FALSE){
  control.samples <- setdiff(colnames(df), c("pos", "chr", test.sample,  "cpg", "gene"))
  n.positions <- dim(df)[1]

  # control.mat stores the reference to each sample (column) that should be used to pull methylation from.
  control.mat <- matrix(sample(1:length(control.samples), size=n.positions*bootstraps, replace=TRUE), nrow=bootstraps)
  corrective.add <- seq(0,(n.positions*length(control.samples)-1), length(control.samples))
  corrective.control.mat <- t(control.mat) + corrective.add
  meth.mat <- t(df[,control.samples])
  
  #get string of test values that will remain stable
  testvals <- df %>% pull(all_of(test.sample))
  
  # make a container for the results that is #stats X #bootstraps
  results.mat <- matrix(0,3,bootstraps)
  
  # run test for each bootstrap column ( meth.mat[corrective.control.mat[,i]]) and store p, cohen, and diff
  for(i in 1:bootstraps){
    rand.positions <- sample(1:n.positions, subsample, FALSE)
    cont.vec <- meth.mat[corrective.control.mat[,i]][rand.positions]
    test.vec <- testvals[rand.positions]
    testmean <- mean(test.vec)
    meancont <- mean(cont.vec)

    if(fit.model){
      fitdist
      suppressWarnings(fit <- fitdist(cont.vec, "beta", method="mge", gof="KS", lower=c(0,0), silent=TRUE))
      betastat <- pbeta(testmean, fit$estimate[1], fit$estimate[2])
    } else {
      fit <- beta_mom(cont.vec)
      #print(fit)
      betastat <- pbeta(testmean, fit$alpha, fit$beta)
    }
    
    meandiff <- testmean - meancont
    
    if(meandiff > 0){
      results.mat[1,i] <- 1-betastat
    } else {
      results.mat[1,i] <- betastat
    }
    pooled.var <- sd(c(cont.vec, test.vec))
    results.mat[2,i] <- meandiff/pooled.var
    results.mat[3,i] <- meandiff
  } 
  meanResults <- rowMeans(results.mat)
  #corrected.p <- p.adjust(meanResults[1], method="BH", n=correction*bootstraps)
  #print(results.mat)
  return(tibble(pvalue=meanResults[1], p_sd=sd(results.mat[1,]), effect.size=meanResults[2], effect_sd=sd(results.mat[2,]), diff=meanResults[3], diff_sd=sd(results.mat[3,]), n=n.positions))
}

coordBetaLong <- function(df, test.sample, correction="BH", minPositions=100, bootstraps=25, exclusions="", fit.model=FALSE){
  defaultW <- getOption("warn")
  options(warn = -1)
  robust.CpGs <- unique(df %>% group_by(cpg) %>% filter(n() > minPositions) %>% pull(cpg))
  small.df <- df %>% filter(cpg %in% robust.CpGs) %>% dplyr::select(-any_of(exclusions))
  
  results.df <- small.df %>% group_by(cpg) %>% dplyr::group_modify( ~ runBetaSub(.x, test.sample, bootstraps=bootstraps, fit.model=fit.model, subsample=minPositions))
  
  results.df["bhpvalue"] <- p.adjust(results.df %>% pull(pvalue), method = correction)
  options(warn=defaultW)
  return(results.df)
}

beta_mom <- function(x) {
  m_x <- mean(x, na.rm = TRUE)
  s_x <- sd(x, na.rm = TRUE)

  alpha <- m_x*((m_x*(1 - m_x)/s_x^2) - 1)
  beta <- (1 - m_x)*((m_x*(1 - m_x)/s_x^2) - 1)

  return(list(alpha = alpha, beta = beta))
}

