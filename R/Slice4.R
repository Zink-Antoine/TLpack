#' Slice4
#'
#' MC Analysis TL (slice +gibbs) following gibbs4.R Hickey (2006)
#' Bayesian approach without expert judgment
#' additional scalar variable T (Slice sampler)
#'
#'
#' @inheritParams Slice1
#' @param yn [numeric] (**with default**) the reference luminescence signal. For MAAD TL, it is extrapolated to yn=0.
#' @param k [numeric] (**required**) the number of parameters in the group
#' involving sigma
#' @param n.burnin [numeric] (**with default**)	number of iterations to discard at the beginning (burn in).
#'  Default is n.iter/2, that is, discarding the first half of the simulations.
#' @param n.thin [numeric] (**with default**) thinning rate. Must be a positive integer. Set n.thin > 1 to save memory and computation time if n.iter is large.
#' Default is max(1, floor((n.iter-n.burnin) / 10)) which will only thin if there are at least 20 simulations.
#' @param threshold [numeric] (**with default**) measured point below or equal to the threshold are rejected.
#' Default is the minimum measured value.
#'
#' @import Slice
#'
#' @return an (r × 5)-matrix,
#'  \tabular{lll}{
#'  **column** \tab **Type** \tab **Description**\cr
#'  `alpha` \tab  `numeric` \tab intercept\cr
#'  `beta` \tab `numeric` \tab slope \cr
#'  `sigma2` \tab  `numeric` \tab variance\cr
#'  `T` \tab `numeric` \tab Temperature \cr
#'  `De` \tab `numeric` \tab 'True' Dose value \cr
#' }
#'
#' @references Gibbs sampler: Hickey, G. L. 2006. « The Linear Calibration Problem: A Bayesian Analysis ». PhD Thesis, PhD dissertation, University of Durham. 1–148. http://www.dur.ac.uk/g.l.hickey/dissertation.pdf.
#' @references chapter 6.3.6 and Appendix G.7
#' @references Slice sampler: Neal, R. 2003. « Slice Sampling ». Annals of Statistics 31 (3): 705‑67.
#'
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
#' test<-Slice4(Dose,df.T,df.y,k=1,n.iter)
#'
#'
Slice4<-
function (Dose,df.T,df.y,yn=0, k=1,
          n.iter,n.burnin=n.iter/2,n.thin=max(1,floor(n.iter-n.burnin)/10),threshold=min(df.y[df.y>0])) {

  n<-length(Dose)
  df.x<-Dose
  n.x<-seq(1,n)
  Rmx<-apply(df.T,2,max)
  Lmin<-apply(df.T,2,min)

  #variables: initial values
  mat <- matrix(ncol=5, nrow=n.iter)
  mcInit<-list()
  for (j in 1:ncol(df.y)){
    repeat{
      mcInit[[j]]<-Slice_Init(df.T[,j],df.y[,j])
      Tj<-mcInit[[j]]$x0
      foo_x<-mcInit[[j]]$foo_x
      if (foo_x(Tj)>threshold) break
    }
  }

  alpha<- 1
  beta<- 1
  sigma2<- 1
  T<-mcInit[[1]]$x0
  De<-0
  mat[1, ] <- c(alpha,beta,sigma2, T,De)

for (i in 2:n.iter) {
  #Temperature calculation using Slice sampler
  repeat{
    run<-Slice_Run(T,mcInit[[1]]$foo_x,mcInit[[1]]$foo_y,mcInit[[1]]$hist_y,Rmx=Rmx[[1]])
    print(i)
    if (mcInit[[1]]$foo_x(run[[1]])>threshold) break

  }
  T<-run[[1]]
  L<-run[[2]]
  R<-run[[3]]
  y0<-run[[4]]
  for (j in 2:ncol(df.y)){
    foo_x<-mcInit[[j]]$foo_x
    if (y0>foo_x(T)){
      Sol.hat<-Shrink(foo_x,T,y0,L,R,R,Rmx=Rmx[[j]],Lmin=Lmin[[j]]) #shrinkage
      T<-Sol.hat[[1]]
      L<-Sol.hat[[2]]
      R<-Sol.hat[[3]]
    }
    print(j)
  }
#slice's end

m<-which(df.T[,1]==round(T))

sxx<-sum((df.x-mean(df.x))^2)
S1<-sum(df.y[m,n.x]-beta*df.x)+yn-(beta*De)
alpha<-rnorm(1,mean=(S1/(n+1)),sd=sqrt(sigma2/(n+1)))
S2<-sum(df.x*(df.y[m,n.x]-alpha))+(yn*De)-(alpha*De)
S4<-sum((df.x)^2)+(De^2)
beta<-rnorm(1,mean=(S2/S4),sd=sqrt(sigma2/S4))
S3<-sum((df.y[m,n.x]-alpha-beta*(df.x))^2)+((yn-alpha-beta*De)^2)
tau<-rgamma(1,shape=((n+k)/2),rate=(S3/2))
sigma2<-1/tau
u<-((n+1)*sxx)+(n*((De-mean(df.x))^2))
t<-rgamma(1,shape=0.5,rate=u)
mu1<-((2*sigma2*n*t*(mean(df.x)))+(beta*(yn-alpha)))/(2*sigma2*n*t+beta^2)
var1<-sigma2/(2*sigma2*n*t+beta^2)
De<-rnorm(1,mu1,sqrt(var1))

mat[i,]<-c(alpha,beta,sigma2,T,De)

}
colnames(mat)<-c("intercept","x","sigma2","Temperature","natural dose")
mat<-mat[seq(n.burnin+1,n.iter,n.thin),]
mcmc(mat,start=n.burnin+1,end=n.iter,thin=n.thin)
}
