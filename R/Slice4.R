#' Slice4
#'
#' MC Analysis TL (slice +gibbs) following gibbs4.R Hickey (2006)
#' Bayesian approach without expert judgement
#' additional scalar variable T (Slice sampler)
#'
#'
#' @inheritParams Slice1
#' @param k [numeric] (**required**) the number of parameters in the group
#' involving sigma
#'
#' @import Slice
#'
#' @return an (r × 5)-matrix,
#' De 'True' dose value
#' sigma2 variance
#' alpha intercept
#' beta slope
#' T Temperature
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
#' attach(multiTL)
#' test<-Slice4(Dose,df.T,df.y,k=1,n.iter)
#' detach(multiTL)
#'
#'
Slice4<-
function (Dose,df.T,df.y, k=1, n.iter) {

mcInit<-list()
for (j in 1:ncol(df.y)){
	mcInit[[j]]<-Slice_Init(df.T[,j],df.y[,j])
}

mat <- matrix(ncol=5, nrow=n.iter)
De<- 1
sigma2<- 1
alpha<- 1
beta<- 1
T<-mcInit[[1]]$x0
mat[1, ] <- c(De, sigma2,alpha,beta, T)

n<-length(Dose)
df.x<-Dose
n.x<-seq(1,n)
Rmx<-max(df.T)
Lmin<-min(df.T)

for (i in 2:n.iter) {
#Temperature calculation using Slice sampler
  run<-Slice_Run(T,mcInit[[1]]$foo_x,mcInit[[1]]$foo_y,mcInit[[1]]$hist_y,Rmx=Rmx[[1]])
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
	}
#slice's end

m<-which(df.T[,1]<round(T)+0.1&df.T[,1]>round(T)-0.1)

sxx<-sum((df.x-mean(df.x))^2)
S1<-sum(df.y[m,n.x]-beta*df.x)+y0-(beta*De)
alpha<-rnorm(1,mean=(S1/(n+1)),sd=sqrt(sigma2/(n+1)))
S2<-sum(df.x*(df.y[m,n.x]-alpha))+(y0*De)-(alpha*De)
S4<-sum((df.x)^2)+(De^2)
beta<-rnorm(1,mean=(S2/S4),sd=sqrt(sigma2/S4))
S3<-sum((df.y[m,n.x]-alpha-beta*(df.x))^2)+((y0-alpha-beta*De)^2)
tau<-rgamma(1,shape=((n+k)/2),rate=(S3/2))
sigma2<-1/tau
u<-((n+1)*sxx)+(n*((De-mean(df.x))^2))
t<-rgamma(1,shape=0.5,rate=u)
mu1<-((2*sigma2*n*t*(mean(df.x)))+(beta*(y0-alpha)))/(2*sigma2*n*t+beta^2)
var1<-sigma2/(2*sigma2*n*t+beta^2)
De<-rnorm(1,mu1,sqrt(var1))

mat[i,]<-c(De,sigma2,alpha,beta,T)

}
mat
}
