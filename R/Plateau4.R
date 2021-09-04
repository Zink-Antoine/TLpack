#' Plateau4
#'
#' MC Analysis TL (Plateau +gibbs) following gibbs4.R Hickey (2006)
#' Bayesian approach without expert judgement
#' additional scalar variable T (Plateau sampler)
#'
#'
#' @inheritParams Slice1
#' @param Ti [numeric] (**with default**) temperature minimum for the plateau
#' @param Tf [numeric] (**with default**) temperature maximum for the plateau
#' @param k [numeric] (**required**) the number of parameters in the group
#' involving sigma
#' @param yn [numeric] (**with default**) intensity value of the 'true' (natural) Dose (in extrapolation, yn = 0)
#' @param n.burnin [numeric] (**with default**)	number of iterations to discard at the beginning (burn in).
#'  Default is n.iter/2, that is, discarding the first half of the simulations.
#' @param n.thin [numeric] (**with default**) thinning rate. Must be a positive integer. Set n.thin > 1 to save memory and computation time if n.iter is large.
#' Default is max(1, floor((n.iter-n.burnin) / 10)) which will only thin if there are at least 20 simulations.
#'
#'
#' @return an (r × 5)-matrix (mcmc class),
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
#' test<-Plateau4(Dose,df.T,df.y,k=1,n.iter=n.iter)
#'
#'
Plateau4<-
function (Dose,df.T,df.y, Ti=200,Tf=500, k=1,yn=0,
          n.iter,n.burnin=n.iter/2,n.thin=max(1,floor(n.iter-n.burnin)/10)) {

mat <- matrix(ncol=6, nrow=n.iter)
n<-length(Dose)
df.x<-Dose
n.x<-seq(1,n)

alpha<- 1
beta<- 1
sigma2<- 1
De<-0
mat[1, ] <- c(alpha,beta,sigma2, Ti,Tf,De)

for (i in 2:n.iter) {
  #Temperature range calculation
  T2<-round(runif(1,mat[i-1,4],Tf))
  T1<-round(runif(1,Ti,T2))
  if (T1==T2) {T1<-T1-1}

  #sampling's end

m<-which(df.T[,1]<round(T)+0.1&df.T[,1]>round(T)-0.1)

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

mat[i,]<-c(alpha,beta,sigma2,T1,T2,De)

}
colnames(mat)<-c("intercept","x","sigma2","Temperature1","Temperature2","natural dose")
mat<-mat[seq(n.burnin+1,n.iter,n.thin),]
mcmc(mat,start=n.burnin+1,end=n.iter,thin=n.thin)
}
