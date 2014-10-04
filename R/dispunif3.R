# ----------------------------------------------------------------------------
# R-code (www.r-project.org/) for dispersion statistics for small samples
#
# Copyright (c) 2014 Paul Kattuman
#
# This library is distributed under the terms of the GNU Public License (GPL)
# for full details see the file LICENSE
#
# ----------------------------------------------------------------------------

#' @title Dispersion statistics for small samples
#'
#' @description This program generates dispersion measures P(S more equal than s) and P(S less
#' equal than s) for data vectors of size n. The approach is to resample from the 
#' data vector, choosing 3 elements randomly, and find the defined measures 
#' with S drawn from the uniform distribution on the unit 2-simplex, with bss resamples.
#'
#' @param s The observed sample; (not in share vector form) read in from an external file
#' @param lb lower bound for confidence interval
#' @param ub upper bound for confidence interval
#' @param bss number of bootstrap samples
#' @export
#' @return
#' \code{dispunif3} returns a matrix with two columns for dispersion measures
#' P(s >= S) and P(s <= S).
#' @author Paul Kattuman, Thilo Klein, Jun Ma
#' @keywords dispersion statistics
#' @references Kattuman, P.A. and Ma, J. (2014). Dispersion Statistics for
#' Small Samples. \emph{Working Paper}, Cambridge Judge Business School.
#' @examples
#' ## share vector of length 3
#' dispunif3(s=c(1,3,7),lb=0.1,ub=0.9,bss=2000)
#' 
#' ## share vector of length larger 3
#' dispunif3(s=c(1,2,3,4,5,6,7,8,9,10),lb=0.1,ub=0.9,bss=2000)
dispunif3 <- function(s=c(1,3,7),lb=0.1,ub=0.9,bss=2000){

    ## Measures: P(s more equal than S), and P(s less equal than S) 
 
    if(length(s) > 3){
      
        # Resample bss times from the data vector
        x <- matrix(NA, nrow=bss, ncol=3)
        for(i in 1:bss){ # for the number of bootstrap resamples
            b     <- s[sample(length(s))]  # randomly permute smaple, and choose the first 3 elements 
            x[i,] <- b[1:3]                # x: matrix containing bss no of bootstrap samples of size 3.
        }                    

        # Convert data vector resamples to share vectors
        x       <- t(apply(x,1,sort))  # share vectors 
        xs      <- t(apply(x,1,cumsum))
        xsct    <- t(sapply(1:bss, function(i) x[i,]/xs[i,3]))  # bsxfun(@rdivide, x',xs(:,3)')
        x       <- xsct     # double(xsct')
        xst     <- t(apply(xs, 1, function(i) i/i[3])) # bsxfun(@rdivide, xs',xs(:,3)')
        sx      <- xst      # double(xst')
        #sx[,3]  <- 1        # final cumulative sum element is 1.

        # Area  Approach, applied to bss resamples of size 3 from the data
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        a <- 1  # a  indexes the col. number

        numwlp <- matrix(NA, nrow=bss, ncol=a)
        numwhp <- matrix(NA, nrow=bss, ncol=a)
        for(c in 1:bss){
            a1  <- x[c,1]
            a2  <- x[c,2]
            a3  <- x[c,3]

            d1  <- sqrt((a2-a3)^2*2)
            d2  <- sqrt((a2-a1)^2*2)
            numwlp[c,a]  <- ((d1+2*d2)^2-3*d2^2)/2

            if(a2 == a3){
                numwhp[c,a] <- (3*((d2+a1*2*sqrt(2))^2-d2^2)-6*a1^2)/2
            }
            else if( 2*(2*sqrt(2)*a1+d2) >= sqrt(2) ){
                numwhp[c,a] <- (3*((d2+a1*2*sqrt(2))^2-d2^2)-3*(2*(2*sqrt(2)*a1+d2)-sqrt(2))^2)/2
            }
            else{
                numwhp[c,a] <- 3*((d2+a1*2*sqrt(2))^2-d2^2)/2
            }
        }

        # Output
        numwlp <- sort(numwlp)
        numwhp <- sort(numwhp)

        #Output 
        P <- matrix(NA, nrow=4, ncol=2)
        
        # P(s more equal than S)
        P[1,1] <- mean(numwlp[1:bss])   # Mean of bss samples P(S more equal than s)
        P[2,1] <- sd(numwlp[1:bss])    # SD of bss samples P(S more equal than s)
        P[3,1] <- numwlp[round(bss*lb)] # Lower confidence band for P(S more equal than s)
        P[4,1] <- numwlp[round(bss*ub)] # Upper confidence band for P(S more equal than s)

        # P(s less equal than S) 
        P[1,2] <- mean(numwhp[1:bss])   # Mean of bss samples P(S less equal than s)  
        P[2,2] <- sd(numwhp[1:bss])    # SD of bss samples P(S less equal than s)
        P[3,2] <- numwhp[round(bss*lb)] # Lower confidence band for P(S less equal than s)
        P[4,2] <- numwhp[round(bss*ub)] # Upper confidence band for P(S less equal than s)

        #return(list(wlp.x=density(numwlp)$x, wlp.y=density(numwlp)$y
        #  , whp.x=density(numwhp)$x, whp.y=density(numwhp)$y, P=P))
        
        colnames(P) <- c("P(s>=S)","P(s<=S)")
        rownames(P) <- c("mean","sd","lower","upper")
        return(P)
    
    } else if(length(s) == 3){

        x      <- s
        x      <- sort(x)  # share vectors 
        xs     <- cumsum(x)
        xsct   <- x/xs[3]  # bsxfun(@rdivide, x',xs(:,3)')
        x      <- xsct     # double(xsct')
        xst    <- xs/xs[3] # bsxfun(@rdivide, xs',xs(:,3)')
        sx     <- xst      # double(xst')
        #sx[3]  <- 1        # final cumulative sum element is 1.

        #Area  Approach, applied to the one sample of size 3 
        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        a1 <- x[1]
        a2 <- x[2]
        a3 <- x[3]

        d1 <- sqrt((a2-a3)^2*2)
        d2 <- sqrt((a2-a1)^2*2)
        numwlp <- ((d1+2*d2)^2-3*d2^2)/2

        if(a2 == a3){
            numwhp <- (3*((d2+a1*2*sqrt(2))^2-d2^2)-6*a1^2)/2
        } else if( 2*(2*sqrt(2)*a1+d2) >= sqrt(2) ){
            numwhp <- (3*((d2+a1*2*sqrt(2))^2-d2^2)-3*(2*(2*sqrt(2)*a1+d2)-sqrt(2))^2)/2
        } else{
            numwhp <- 3*((d2+a1*2*sqrt(2))^2-d2^2)/2
        }


        #Output 
        P <- matrix(NA, nrow=1, ncol=2)
      
        # P(s more equal than S)
        P[1,1] <- numwlp   # Mean of bss samples P(S more equal than s)
        # P(s less equal than S) 
        P[1,2] <- numwhp   # Mean of bss samples P(S less equal than s)  

        #return(list(wlp=numwlp, whp=numwhp, P=P))
        colnames(P) <- c("P(s>=S)","P(s<=S)")
        rownames(P) <- "mean"
        return(P)

    }
  
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    # xlswrite('Smeqs.xls',numwlp)     # output: all bss rows of
    # xlswrite('Sleqs.xls',numwhp)     # output: bss rows 

}


