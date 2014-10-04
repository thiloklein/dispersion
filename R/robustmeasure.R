# ----------------------------------------------------------------------------
# R-code (www.r-project.org/) for dispersion statistics for small samples
#
# Copyright (c) 2014 Paul Kattuman
#
# This library is distributed under the terms of the GNU Public License (GPL)
# for full details see the file LICENSE
#
# ----------------------------------------------------------------------------

#' @title Dispersion statistics for small samples (general)
#'
#' @description This program generates sampling distributions of CoV, Gini, Entropy(H), 
#' GE, HHI, and P(S more equal than s), P(S less equal than s) from draws from Pareto
#' distribution, for various alpha values.
#' Also generates data for the power curve for all the indices CoV, Gini, Entropy(H), GE, 
#' HHI and P(S more equal than s), P(S less equal than s) for the left end of the samp dist 
#' for chosen size of the test: ie. prop from the true dist to the left of the 5# (1#) cutoff 
#' from the null hypothesized  distribution for the right end of the samp dist for chosen 
#' size of the test: ie. prop from the true dist to the right of the 95# (99#) cutoff from 
#' the null hypothesized dist.
#'
#' @param shape Gamma distribution parameter for shape, for D(1), shape=1: uniform distribution on the simplex
#' @param s The observed sample (not in share vector form).
#' @param n bootstrap sample size , typically 3
#' @param bss number of bootstrap samples
#' @param lb lower bound 5 or .05 or 10 or .1
#' @param ub upper bound 95 or .95 or 90 or .9
#' @param l1 the number of samples drawn for prior - S~D(lambda,) to calculate by simulation
#' @export
#' @return
#' \code{robustmeasure} returns samp dists of different measures for different alphas and
#' power test proportions at the right of the right side cutoff and left of the left side 
#' cutoff  of the samp dist for different measures.
#' @author Paul Kattuman, Thilo Klein, Jun Ma
#' @keywords dispersion statistics
#' @references Kattuman, P.A. and Ma, J. (2014). Dispersion Statistics for
#' Small Samples. \emph{Working Paper}, Cambridge Judge Business School.
#' @examples
#' ## share vector of length 3
#' robustmeasure(shape=1,s=c(1,3,7),n=3,bss=2000,lb=0.1,ub=0.9,l1=1000)
#' 
#' ## share vector of length larger 3
#' robustmeasure(shape=1,s=c(1,2,3,4,5,6,7,8,9,10),n=3,bss=2000,lb=0.1,ub=0.9,l1=100)
#' 
#' \dontrun{
#' ## matrix of share vectors
#' data(libor)
#' robustmeasure(s=libor, shape=1, l1=100, n=3, bss=2000, lb=10, ub=90)
#' }
robustmeasure <- function(s,shape,l1,n,bss,lb,ub){
  
  if(is.null(dim(s))){  ## if share vector then run core function
    core(s=s, shape=shape, l1=l1, n=n, bss=bss, lb=lb, ub=ub)
  } 
  else{ ## if matrix of share vectors, then iterate core function over number of rows
    s = as.matrix(s)
    r <- sapply(c("mean","sd","lower","upper"),function(x) NULL)
    
    ## create progress bar
    total <- dim(s)[1]
    pb <- txtProgressBar(min = 0, max = total, style = 3)
    
    for(i in 1:dim(s)[1]){
      h = core(s=s[i,], shape=1, l1=1000, n=3, bss=2000, lb=10, ub=90)
      
      for(j in 1:4){
        r[[j]] = rbind(r[[j]], h[j,])
      }
      
      setTxtProgressBar(pb, i)
      
    }
    
    ## close progress bar
    close(pb)
    
    ## add rownames to results
    for(j in 1:4){
      rownames(r[[j]]) = rownames(s)
    }
    
    ## return list of results
    return(r)
  }
}


core <- function(s,shape,l1,n,bss,lb,ub){
  
  if(lb > 1){
    lb <- lb/100
  }
  if(ub > 1){
    ub <- ub/100
  }
  
  a <- 1     # "a" indexes the col number
  # pb <- 1  # count
  # ps <- 1  # counter for the full resampled  means
  
  # finding P(s more equal than S), and P(s less equal than S) 
  
  # resample bss times from the data vector
  x <- matrix(NA, nrow=bss, ncol=n)
  for(i in 1:bss){ # for the number of bootstrap samples
    b <- s[sample(length(s))]   # randomly permute the smaples of size l, and choose the first n elements (bootstrap smaple size)
    x[i,] <- b[1:n]             # x is the matrix containing bss no of bootstrap samples of size n.    
  }                  

  meanx <- t(apply(x,1,mean))
  vasx  <- t(apply(x,1,var))  # variance with n-1
  y     <- t(apply(x,1,sort)) # orginal bootstrap sample (renamed) and sorted for Gini calculation
  meany <- t(apply(y,1,mean))
  
  ## convert data vector resamples to share vectors
  x      <- t(apply(x,1,sort))  # share vectors 
  xs     <- t(apply(x,1,cumsum))
  xsct   <- t(sapply(1:bss, function(i) x[i,]/xs[i,n])) # bsxfun(@rdivide, x',xs(:,3)')
  x      <- xsct
  xst    <- t(apply(xs, 1, function(i) i/i[n])) # bsxfun(@rdivide, xs',xs(:,3)')
  sx     <- xst
  sx[,n] <- 1   # ensuring final cumulative sum element is 1.
  
  # use simulation if the required resample size is > 3, or if shape
  # parameter > 1 otherwise, use the geometric area 
  # if the required resample size = 3, and shape
  # parameter = 1 use geometric area
  if(shape==1 && n==3){
    # Area  Approach, applied to bss resamples from the data
    ######################################################
    
    numwlp <- matrix(NA, nrow=bss, ncol=a)
    numwhp <- matrix(NA, nrow=bss, ncol=a)
    
    for(c in 1:bss){
      a1 = x[c,1]
      a2 = x[c,2]
      a3 = x[c,3]
      
      d1 = sqrt((a2-a3)^2*2)
      d2 = sqrt((a2-a1)^2*2)
      numwlp[c,a] = ((d1+2*d2)^2-3*d2^2)/2
      
      if(a2 == a3){
        numwhp[c,a] = (3*((d2+a1*2*sqrt(2))^2-d2^2)-6*a1^2)/2
      }
      else if( 2*(2*sqrt(2)*a1+d2) >= sqrt(2) ){
        numwhp[c,a] = (3*((d2+a1*2*sqrt(2))^2-d2^2)-3*(2*(2*sqrt(2)*a1+d2)-sqrt(2))^2)/2  
      }
      else{
        numwhp[c,a] = 3*((d2+a1*2*sqrt(2))^2-d2^2)/2
      }
    }  
    
    ######################################################
  }
  else{
    # Simulation where area is not used
    ###################################################
    
    #generate l1 random draws from D(lambda), and translate to share vectors
    scale  = 1
    z      = matrix( rgamma(n=l1*n, shape=shape, scale=scale), nrow=l1, ncol=n) # gamrnd(shape,scale,l1,n) 
    zco    = t(apply(z,1,sort))
    zs     = t(apply(zco,1,cumsum)) # cumsum(zco,2)
    zsct   = t(sapply(1:n, function(i) zs[,i]/zs[,n])) # bsxfun(@rdivide, zs',zs(:,n)')
    zscosc = t(zsct)
    zscosc[,n] = 1  # only used for prec and succ, as comparing lorenz curves of sample with DD LC
    
    numwlp <- matrix(NA, nrow=bss, ncol=a)
    numwhp <- matrix(NA, nrow=bss, ncol=a)
    
    for(c in 1:bss){
    
      msx = matrix(NA, nrow=l1, ncol=n)
      for(m1 in 1:l1){
        for(m2 in 1:n){
          msx[m1,m2] = sx[c,m2]
        }  
      }
      
      mcheck = matrix(NA, nrow=l1, ncol=1)
      for(m1 in 1:l1){
        mcheck[m1,1] = n 
      }
      
      indexwl     = zscosc >= msx
      indexwh     = zscosc <= msx
      indexwls    = apply(indexwl,1,sum) # sum(indexwl,2)
      indexwhs    = apply(indexwh,1,sum) # sum(indexwh,2)
      indexwlsc   = indexwls == mcheck
      indexwhsc   = indexwhs == mcheck
      numwlp[c,a] = sum(indexwlsc)/l1   #the probability of P(S more equal than s)  clue less - towards eq line: Schur-convex
      numwhp[c,a] = sum(indexwhsc)/l1   #the probability of P(S less equal than s)  clue less - towards max ineq : Schur-concave
    }
    #######################################################
  }
  
  numwlp = sort(numwlp)
  numwhp = sort(numwhp)
  
  #Output collected.
  P <- matrix(NA, nrow=4, ncol=2)
  P[1,1] = mean(numwlp[1:bss])  #Mean of bss samples P(S more equal than s)
  P[2,1] = sd(numwlp[1:bss])   #SD of bss samples P(S more equal than s)
  P[3,1] = numwlp[round(bss*lb)]
  P[4,1] = numwlp[round(bss*ub)]
  
  P[1,2] = mean(numwhp[1:bss])  #Mean of bss samples P(S less equal than s)  
  P[2,2] = sd(numwhp[1:bss])   #SD of bss samples P(S less equal than s)
  P[3,2] = numwhp[round(bss*lb)]
  P[4,2] = numwhp[round(bss*ub)]
  
  colnames(P) <- c("P(s>=S)","P(s<=S)")
  rownames(P) <- c("mean","sd","lower","upper")
  return(P)
  
  #xlswrite('SMes.xls',numwlp)   # output: bss rows
  #xlswrite('SLes.xls',numwhp)   # output: bss rows  
  
}
