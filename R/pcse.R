#
# delia@caltech.edu
# pcse.R
# function to estimate pcses
#

pcse <- function(object, groupN, groupT, pairwise=FALSE){

 mc <- match.call()

 # Input Checks

 check <- class(object)
 if ("lm" %in% check == FALSE){
   stop("Formula object must be of class 'lm'.")
 }

 check <- length(groupN) == length(groupT)
 if (check == FALSE){
   stop("Length of groupT and groupN be of equal length.")
 }

 check <- length(groupN) == dim(eval(object$call$data))[1]
 if (check == FALSE){
   stop("Length of groupN and groupT must equal nrows of using data.")
 }
 
 check <- is.na(groupN)
 if (any(check == TRUE)){
   stop("There must not be any missing values in the CS groupN!")
 }
 check <- is.na(groupT)
 if (any(check == TRUE)){
   stop("There must not be any missing values in the TS groupT!")
 }

 N     <- length(na.omit(unique(groupN)))
 T     <- length(na.omit(unique(groupT)))
 check <- N*T >= dim(eval(object$call$data))[1]
 if (check == FALSE){
   stop("There cannot be more than N*T rows in the using data!")
 }

 # Make factors numeric for easier manipulation.

 if ("factor" %in% class(groupN)){
   groupN <- as.numeric(groupN)
 }
 
 if ("factor" %in% class(groupT)){
   groupT <- as.numeric(groupT)
 }
 
 # Check for balanced data.

 units <- unique(groupN)
 units <- na.omit(units)
 N     <- length(units)
 time  <- unique(groupT)
 time  <- na.omit(time)
 T     <- length(time)

 check <- 0
 for (i in 1:N){
   if (sum(groupN  == units[i]) == T){
     check <- check + 1
   }
   else{
     check <- check
   }
 }
 
 flag <- ifelse(check == N, TRUE, FALSE)
 

 # END CHECKS

 # Data manipulation.

 using       <- data.frame(groupN = groupN, groupT = groupT)
 using       <- na.omit(using)
 using$resid <- resid(object)
 using       <- data.frame(using, model.matrix(object))

 ord <- order(using$groupN, using$groupT)
 using <- using[ord, ]

 units <- unique(using$groupN)
 N     <- length(units)
 time  <- unique(using$groupT)
 T     <- length(time)

 # calculate avg number of obs per panel.
 avgN <- dim(using)[1] / N

 # get largest balanced subset of data.

 brows <- c()
 for (i in 1:T){
   br    <- which(using$groupT == time[i])
   check <- length(br) == N
   if (check == TRUE){
     brows <- c(brows, br)
   }
 }

 balanced <- using[brows, ]
 ord      <- order(balanced$groupN, balanced$groupT)
 balanced <- balanced[ord, ]

 Bunits <- unique(balanced$groupN)
 BN     <- length(Bunits)
 Btime  <- unique(balanced$groupT)
 BT     <- length(Btime)
 
 # Get rectangular data (if not already).
 rect <- using

 # Fill in missing Time and Cross-Section Units.

 rect$groupN <- as.numeric(rect$groupN)
 Runits      <- unique(rect$groupN)
 Runits      <- na.omit(Runits)
 
 # what are we missing in each panel?

 missN <- matrix(NA, N, 2)
 for (i in 1:N){
 if (sum(rect$groupN == Runits[i]) != T){
     missN[i, 1] <- Runits[i]
     missN[i, 2] <- T - (sum(rect$groupN == Runits[i]))
   }
 }
 missN <- na.omit(missN)
 missT <- c()
 tmp   <- c()


 if (dim(missN)[1] != 0 & dim(missN)[2] != 0){
   for (i in 1:dim(missN)[1]){
     tt <- time %in% rect$groupT[rect$groupN == missN[i, 1]]
     missT <- c(missT, time[!tt])
   }
   for (i in 1:dim(missN)[1]){
     tmp <- c(tmp, rep(missN[i, 1], missN[i, 2]))
   }
   missN <- tmp
   M     <- length(missN)

   # Fill in using data with NAs to get N*T rows.

   R <- dim(rect)[1]
   C <- dim(rect)[2]
   if (R != N*T){
     for (i in (R+1):(N*T)){
       rect[i, ] <- rep(NA, C)
     }
   }
   rect[c((R+1):(R+M)), 1] <- missN
   rect[c((R+1):(R+M)), 2] <- missT
 }

 ord    <- order(rect$groupN, rect$groupT)
 rect   <- rect[ord, ]
 Runits <- unique(rect$groupN)
 Runits <- na.omit(Runits)
 RN     <- length(Runits)
 Rtime  <- unique(rect$groupT)
 Rtime  <- na.omit(Rtime)
 RT     <- length(Rtime)
    
 # ESTIMATION
 
 if (flag == TRUE){
   # estimate Sigma.hat using whole data frame and get middle matrix.

   e <- using$resid
   E <- matrix(NA, N, T)
   
   for (i in 1:N){
     E[i, ] <- e[c(((i - 1) * T + 1):(i * T))]
   }
   E <- t(E)

   Sigma.hat <- crossprod(E)/T
   X         <- as.matrix(using[, 4:dim(using)[2]])
   omega     <- kronecker(Sigma.hat, diag(1, T))
   middle    <- t(X) %*% omega %*% X
   nobs      <- length(e)
   dataX     <- X
 }

 # END flag == TRUE

 
 if (flag == FALSE){
   if (pairwise == FALSE){

     e <- balanced$resid
     E <- matrix(NA, BN, BT)
     
     for (i in 1:BN){
       E[i, ] <- e[c(((i - 1) * BT + 1):(i * BT))]
     }
     E <- t(E)

     Sigma.hat <- crossprod(E) / BT

     if (avgN/2 > BT){
       print(paste("Caution! The number of CS observations per panel,", BT,
                   ", used to compute the vcov matrix is less than half the",
                   "average number of obs per panel in the original data.",
                   "You should consider using pairwise selection."))
     }

     X           <- as.matrix(rect[ , 4:dim(rect)[2]])
     X[is.na(X)] <- 0
     omega       <- kronecker(Sigma.hat, diag(1, T))
     middle      <- t(X) %*% omega %*% X
     nobs        <- length(resid(object))
   }

   # END Pairwise == FALSE
   
   if (pairwise == TRUE){

     # The next section of code follows gauss procedure from Franzeze (1996).
     # Get vector of 1/0 for valid row of obs or not.
     
     V <- rect[ , 4:dim(rect)[2]]

     for (i in 1:dim(V)[1]){
       for (j in 1:dim(V)[2]){
         if (is.na(V[i, j])){
           V[i, j] <- 0
         }
         else {
           V[i, j] <- 1
         }
       }
     }
     valid <- apply(V, 1, prod)
     nobs  <- sum(valid)

     # Reshape valid and E to RTxRN

     e           <- rect$resid
     e[is.na(e)] <- 0
     E           <- matrix(NA, RN, RT)
     for (i in 1:RN){
       E[i, ] <- e[c(((i - 1) * RT + 1):(i * RT))]
     }
     E <- t(E)

     V <- matrix(NA, RN, RT)
     for (i in 1:RN){
       V[i, ] <- valid[c(((i - 1) * RT + 1):(i * RT))]
     }
     V <- t(V)
   
     numer <- crossprod(E)
     denom <- crossprod(V)
     check <- is.na(denom)
     if (sum(check) != 0){
       stop("Error! A CS-unit exists without any obs or without any obs in
             common with another CS-unit. You must remove that unit from the
             data passed to pcse().")
     }

     # Element by element division to get right T in denom.

     Sigma.hat <- numer/denom
     X           <- as.matrix(rect[ , 4:dim(rect)[2]])
     X[is.na(X)] <- 0
     omega       <- kronecker(Sigma.hat, diag(1, T))
     middle      <- t(X) %*% omega %*% X
   }
   # End pairwise == TRUE
 }
 
 # End FLAG == FALSE
 
 # estimate X'X with whole X.
 # estimate vcov with middle and X'X
 # calculate results.


 XX     <- t(X) %*% X
 XXinv  <- solve(XX)
 vcov   <- XXinv %*% middle %*% XXinv
 pcse   <- sqrt(diag(vcov))
   
 b      <- summary(object)$coef[ , 1]
 tstats <- b/pcse
 df     <- nobs - ncol(X)
 pval   <- 2*pt(abs(tstats), df, lower.tail=FALSE)
 
 res <- list(vcov=vcov, pcse=pcse, b=b, tstats=tstats, df=df, pval=pval,
             pairwise=pairwise, nobs=nobs, nmiss=(N*T)-nobs, call=mc)

 class(res) <- "pcse"
 return(res)
}
 
# End pcse.
