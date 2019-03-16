ASRV <- function(genoTypeData, phenoTypeStatus, nDropOut = 1000, nReplicates = 1000000, isCheckParam = TRUE)
{
  ###################### checking function parameters ################# 
  if (isCheckParam) {
    if (all(phenoTypeStatus == 0 | phenoTypeStatus == 1) != TRUE)
      stop('All values of vector phenoTypeStatus must be 0 or 1.')
    if (all(genoTypeData == 0 | genoTypeData == 1 | genoTypeData == 2) != TRUE)
      stop('All values of vector genoTypeData must be coded as 0, 1, 2.')
    if (nDropOut >= nReplicates)
      stop('Invalid parameters(nDropOut < nReplicates):\n','\tburn = ', 
           nDropOut, ', nReplicates = ', nReplicates)
  }
  nPTS <- length(phenoTypeStatus)
  ncolGTD <- ncol(genoTypeData)
  nPTS1 <- nPTS- 1
  tPTS <- t(phenoTypeStatus)
  ###################### Contingency Table ################# 
  coTbl <- array(c(0), dim = c(2, nPTS))
  cnt <- 0
  cpGTD <- genoTypeData
  while (nrow(genoTypeData) > 0) {
    cnt <- cnt + 1
    for (i in 1:nrow(genoTypeData)) {
      diffGTD <- genoTypeData[1, ] - genoTypeData[i, ]
      if (all(diffGTD == 0)) { 
        nd <- i - sum(coTbl[, cnt])
        cpGTD <- cpGTD[-nd, ]
        dim(cpGTD) <- c(nPTS1 - sum(coTbl), ncolGTD)
        ty <- tPTS[nd]
        tPTS <- tPTS[-nd]
        if (ty == 1) 
          coTbl[1, cnt] <- coTbl[1, cnt] + 1
        else 
          coTbl[2, cnt] <- coTbl[2, cnt] + 1
      }     
    }
    genoTypeData <- cpGTD
  }
  if (sum(coTbl) ==  nPTS1) {
    cnt <- cnt + 1
    if (tPTS == 1) 
      coTbl[1, cnt] <- 1
    else 
      coTbl[2, cnt] <- 1
  }
  cutCoTbl <- coTbl[, 1:cnt]
  ###################### Markov Chain Monte Carlo #################
  if (!is.matrix(cutCoTbl))
    stop('MCMC:\nNot enough sample data to estimate')
  rsum <- rowSums(cutCoTbl)
  csum <- colSums(cutCoTbl)
  expectedFreqs <- rsum %*% t(csum)
  expectedFreqs <- expectedFreqs / nPTS
  tbl4CalcStat <- (cutCoTbl - expectedFreqs) * (cutCoTbl - expectedFreqs) / expectedFreqs
  statVal <- sum(tbl4CalcStat)
  cpStatVal <- statVal
  sigCnt <- c(0)
  for (i in 1:(nReplicates + nDropOut)) { 
    z <- sample(1:cnt, 2)
    a <- sort(z)[1]
    b <- sort(z)[2]
    e <- sample(c(-1, 1), 1, replace = TRUE, prob = c(0.5, 0.5))
    c <- cutCoTbl[1, a] + e
    d <- cutCoTbl[2, a] - e
    f <- cutCoTbl[1, b] - e
    g <- cutCoTbl[2, b] + e
    if ((c >= 0) && (d >= 0) && ( f>= 0) && (g >= 0)) {
      u <- runif(1, 0, 1)
      u0 <- ifelse(e == 1, 
                   cutCoTbl[1, b] * cutCoTbl[2, a] / (c * g), 
                   cutCoTbl[1, a] * cutCoTbl[2, b] / (d * f))
      if (u <= u0) { 
        cutCoTbl[1, a] <- c
        cutCoTbl[2, a] <- d
        cutCoTbl[1, b] <- f
        cutCoTbl[2, b] <- g
        tbl4CalcStat <- (cutCoTbl - expectedFreqs) * (cutCoTbl - expectedFreqs) / expectedFreqs
        cpStatVal <- sum(tbl4CalcStat)
      }
    }
    if ((cpStatVal >= statVal) && (i > nDropOut)) 
      sigCnt <- sigCnt + 1
  }
  return (sigCnt / nReplicates)
}
