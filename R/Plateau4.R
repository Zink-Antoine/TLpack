#' Plateau4
#'
#' MC Analysis TL (Plateau +gibbs) following gibbs4.R Hickey (2006)
#' Bayesian approach without expert judgement
#' additional scalar variable T (Plateau sampler)
#'
#'
#' @inheritParams Slice1
#' @inheritParams Plateau5
#' @inheritParams Slice4
#' @param inv [logical] (**with default value**) inverse (TRUE) or direct (FALSE) regression (default value FALSE)
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
          n.iter,inv=FALSE,n.burnin=n.iter/2,n.thin=max(1,floor(n.iter-n.burnin)/10)) {


n<-length(Dose)
#df.x<-Dose
n.x<-seq(1,n)
Rmx<-apply(df.T,2,max)
Lmin<-apply(df.T,2,min)
if(Tf>Rmx[[1]]) return (print("Tf out of range"))

#variables: initial values
mat <- matrix(ncol=6, nrow=n.iter)

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

m<-which(df.T[,1]<T2+0.1&df.T[,1]>T1-0.1)
if (inv) {
  Y<-Dose
  X<-apply(df.y[m,n.x],2,mean)
}
else{
  X<-Dose
  Y<-apply(df.y[m,n.x],2,mean)
}

sxx<-sum((X-mean(X))^2)
S1<-sum(Y-beta*X)+yn-(beta*X0)
alpha<-rnorm(1,mean=(S1/(n+1)),sd=sqrt(sigma2/(n+1)))
S2<-sum(X*(Y-alpha))+(yn*X0)-(alpha*X0)
S4<-sum(X^2)+(X0^2)
beta<-rnorm(1,mean=(S2/S4),sd=sqrt(sigma2/S4))
S3<-sum((Y-alpha-beta*X)^2)+((yn-alpha-beta*X0)^2)
tau<-rgamma(1,shape=((n+k)/2),rate=(S3/2))
sigma2<-1/tau
u<-((n+1)*sxx)+(n*((X0-mean(X))^2))
t<-rgamma(1,shape=0.5,rate=u)
mu1<-((2*sigma2*n*t*(mean(X)))+(beta*(yn-alpha)))/(2*sigma2*n*t+beta^2)
var1<-sigma2/(2*sigma2*n*t+beta^2)
X0<-rnorm(1,mu1,sqrt(var1))

mat[i,]<-c(alpha,beta,sqrt(var1),T1,T2,-X0)

}
colnames(mat)<-c("intercept","x","std dev","Temperature1","Temperature2","natural dose")
mat<-mat[seq(n.burnin+1,n.iter,n.thin),]
mcmc(mat,start=n.burnin+1,end=n.iter,thin=n.thin)
}
