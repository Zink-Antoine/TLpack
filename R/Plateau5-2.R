##
#' Plateau5_2
#'
#' MC Analysis TL from linear regression using plateau samplerS
#'
#' @inheritParams Slice1
#' @param Ti [numeric] (**with default**) temperature minimum for the plateau
#' @param Tf [numeric] (**with default**) temperature maximum for the plateau
#' @param inv [logical] (**with default**) TRUE to calculate the inverse estimator.
#' Default is FALSE: i.e. the direct estimator
#' @param n.burnin [numeric] (**with default**)	number of iterations to discard at the beginning (burn in).
#'  Default is n.iter/2, that is, discarding the first half of the simulations.
#' @param n.thin [numeric] (**with default**) thinning rate. Must be a positive integer. Set n.thin > 1 to save memory and computation time if n.iter is large.
#' Default is max(1, floor((n.iter-n.burnin) / 10)) which will only thin if there are at least 20 simulations.
#'
#' @import Slice
#' @import coda
#'
#' @return an (r × 6)-matrix,
#' \tabular{lll}{
#'  **column** \tab **Type** \tab **Description**\cr
#'  `alpha` \tab  `numeric` \tab intercept\cr
#'  `beta` \tab `numeric` \tab slope \cr
#'  `sigma2` \tab  `numeric` \tab variance\cr
#'  `T1` \tab `numeric` \tab lower limit of the temperature range \cr
#'  `T2` \tab `numeric` \tab upper limit of the temperature range \cr
#'  `De` \tab `numeric` \tab 'True' Dose value \cr
#' }
#'
#' @export
#'
#' @examples
#' ##load data
#' if(dev.cur()!=1) dev.off()
#' data(multiTL, envir = environment())
#' Dose<-multiTL$Dose
#' df.T<-multiTL$df.T
#' df.y<-multiTL$df.y
#' n.iter<-multiTL$n.iter
#' test<-Plateau5(Dose,df.T,df.y,Ti=50,Tf=400,n.iter=n.iter)
#' #
#' if(dev.cur()!=1) dev.off()
#' data(TLpan, envir = environment())
#' Dose<-c(0,0,80,80,80,160,160,160)
#' df.T<-matrix(rep(seq(26,500),8),475,8)
#' df.y<-TL.Pan[,1:8]
#' Pan<-Plateau5(Dose,df.T,df.y,n.iter=100,inv=TRUE)
#'
#'
Plateau5_2<-
  function (Dose,df.T,df.y, Ti=200,Tf=500,
            n.iter,inv=FALSE,n.burnin=n.iter/2,n.thin=max(1,floor(n.iter-n.burnin)/10))
    {

    #parameters
    n<-length(Dose)
    n.y<-seq(1,n)
    Rmx<-apply(df.T,2,max)
    Lmin<-apply(df.T,2,min)
    if(Tf>Rmx[[1]]) return (print("Tf out of range"))

    #variables: initial values
    mat <- matrix(ncol=6, nrow=n.iter)

    alpha<- 1
    beta<- 1
    sigma2<- 1
    xn<-0
    mat[1, ] <- c(alpha,beta,sigma2, Ti,Tf,xn)

    for (i in 2:n.iter) {

      #Temperature range calculation
      T2<-round(runif(1,mat[i-1,4],Tf))
      T1<-round(runif(1,Ti,T2))
      if (T1==T2) {T1<-T1-1}

      #sampling's end

      m<-which(df.T[,1]<T2+0.1&df.T[,1]>T1-0.1)
      if (inv) {
        Y<-Dose
        X<-apply(df.y[m,n.y],2,mean)
      }
      else{
        X<-Dose
        Y<-apply(df.y[m,n.y],2,mean)
      }

      p<-runif(n.y,0,1)
      lambda<-rbinom(n.y,1,p)
      Y<-(1-lambda)*Y+lambda*(alpha+beta*X)

      Sxx<-sum((X-mean(X))^2)
      Sxy<-sum((X-mean(X))*(Y-mean(Y)))
      beta<-rnorm(1,mean=(Sxy/Sxx),sd=sqrt(sigma2/Sxx))

      S1<-mean(Y)-beta*mean(X)
      S2<-(1/n+mean(X)^2/Sxx)
      alpha<-rnorm(1,mean=S1,sd=sqrt(sigma2*S2))

      SE<-(Y-alpha-beta*X)^2
      SSE<-sum(SE)
      tau<-rgamma(1,shape=((n-2)/2),rate=(SSE/2))
      sigma2<-1/tau

      De<-rnorm(1,-alpha/beta,sqrt(sigma2))

      mat[i,]<-c(alpha,beta,sigma2,T1,T2,De)
    }
    colnames(mat)<-c("intercept","x","sigma2","Temperature1","Temperature2","natural dose")
    mat<-mat[seq(n.burnin+1,n.iter,n.thin),]
    mcmc(mat,start=n.burnin+1,end=n.iter,thin=n.thin)
  }

