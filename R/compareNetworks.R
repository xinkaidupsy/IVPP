# compare networks
compareNetworks <- function(true,est, directed = FALSE){
  cor0 <- function(x,y,...){
    if (sum(!is.na(x)) < 2 || sum(!is.na(y)) < 2 || sd(x,na.rm=TRUE)==0 | sd(y,na.rm=TRUE) == 0){
      return(0)
    } else {
      return(cor(x,y,...))
    }
  }

  bias <- function(x,y) mean(abs(x-y),na.rm=TRUE)

  if (is.matrix(true) & is.matrix(est)){
    if (directed){
      real <- c(true)
      est <- c(est)
    } else {
      real <- true[upper.tri(true,diag=FALSE)]
      est <- est[upper.tri(est,diag=FALSE)]
    }
  } else {
    real <- true
  }

  # True positives:
  TruePos <- sum(est != 0 &  real != 0)

  # False pos:
  FalsePos <- sum(est != 0 & real == 0)

  # True Neg:
  TrueNeg <- sum(est == 0 & real == 0)

  # False Neg:
  FalseNeg <- sum(est == 0 & real != 0)

  out <- list()

  ### Sensitivity:
  out$sensitivity <- TruePos / (TruePos + FalseNeg)

  # Specificity:
  out$specificity <- TrueNeg / (TrueNeg + FalsePos)

  # Correlation:
  out$correlation <- cor0(est,real)

  out$bias <- bias(est,real)

  out$truePos <- TruePos
  out$falsePos <- FalsePos
  out$trueNeg <- TrueNeg
  out$falseNeg <- FalseNeg

  out <- as.data.table(out, keep.rownames = TRUE)

  return(out)
}
