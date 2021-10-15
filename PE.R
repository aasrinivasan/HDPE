HD_PE = function(Y1, Y2, n1, n2, alpha, type = "F"){
  p = ncol(Y1)
  print(c(n1,n2,p))
  
  maxTest = maximum_type_test(Y1, Y2, n1, n2, p)
  
  quadTest = cqtest(Y1, Y2)
  quadTest = list(stat = quadTest[1], pvalue = quadTest[2])
  
  ### Add a small amount to pvalues to make them not machine rounded to 0
  pvals = c(maxTest$pvalue, quadTest$pvalue)
  if(pvals[1] == 0){pvals[1] = pvals[1] + .Machine$double.eps}
  if(pvals[2] == 0){pvals[2] = pvals[2] + .Machine$double.eps}

  if(type == "F"){
    combin_test = fishersMethod(pvals[1], pvals[2])
  }
  
  if(type == "C"){
    # Cauchy Combination
    combin_test = cauchyMethod(pvals[1], pvals[2])
  }
  
  pv = combin_test$pvalue
  decision = comparePvalues(pv, alpha)
  stat = combin_test$stat
  res = list(alternative = decision, statistic = stat, p.value = pv)
  return(res)
}

maximum_type_test = function(Y1, Y2, n1, n2, p){
  n = n1*n2/(n1+n2)
  
  meanY1 = colMeans(Y1)
  meanY2 = colMeans(Y2)
  
  meanDiff = (meanY1 - meanY2)^2
  
  # calculate Gamma
  cols = 1:p
  gamma = lapply(X = cols, FUN = calcGamma, n1 = n1, n2 = n2, Y1 = Y1, Y2 = Y2)
  gamma = unlist(gamma)
  
  # difference ratio
  ratio = meanDiff/gamma
  maxRatio = max(ratio)
  Mn = (n1*n2/(n1+n2)) * maxRatio
  
  gumbelStat = Mn-2*log(p)+log(log(p))
  
  pv =  1-exp(-1/sqrt(pi)*exp(-(Mn-(2*log(p)-log(log(p))))/2))
  return(list(stat = Mn, pvalue = pv))
}

calcGamma = function(x,n1,n2,Y1,Y2){
  mY1 = colMeans(Y1)
  mY2 = colMeans(Y2)
  colMean1 = mY1[x]
  colMean2 = mY2[x]
  
  diff1 = sum((Y1[,x] - rep(colMean1, length(Y1[,x])))^2)
  diff2 = sum((Y2[,x] - rep(colMean2, length(Y2[,x])))^2)
  
  sumDiff = diff1+diff2
  gammajj = as.numeric(sumDiff/(n1+n2))
  
  return(gammajj)
}
