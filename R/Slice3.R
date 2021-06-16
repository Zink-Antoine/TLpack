#' Slice3
#'
#' MC Analysis TL (slice +gibbs) following gibbs3.R Hickey (2006)
#'
#' additional scalar variable T (Slice sampler)
#'
#' @inheritParams Slice1
#' @param mu_alpha [numeric] (**required**) Prior of the intercept's average
#' @param var_alpha [numeric] (**required**) Prior of the intercept's variance
#' @param mu_beta [numeric] (**required**) Prior of the slope's average
#' @param var_beta [numeric] (**required**) Prior of the slope's variance
#' @param y0 [numeric] (**with default**) Prior of the luminescence for the true dose (Default y0=0, extrapolation) **need ???**
#' @param a [numeric] (**required**) the lower end points of the experimental calibration range
#' @param b [numeric] (**required**) the upper end points of the experimental calibration range
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
#' @references chapter 5 and Appendix G.6
#' @references Slice sampler: Neal, R. 2003. « Slice Sampling ». Annals of Statistics 31 (3): 705‑67.
#'
#' @export
#'
#' @examples
#' ##load data
#' if(dev.cur()!=1) dev.off()
#' data(multiTL, envir = environment())
#' mu_0= -105
#' var_0=50
#' mu_alpha=6.3e-3
#' var_alpha=9e-8
#' mu_beta=5e5
#' var_beta=2e-12
#' y0=0
#' a=6.8274
#' b=1.3819
#' attach(multiTL)
#' test<-Slice3(Dose,df.T,df.y,mu_0,var_0,mu_alpha,var_alpha, mu_beta, var_beta,y0,a,b,n.iter)
#' detach(multiTL)
#'
Slice3<-
function (Dose,df.T,df.y, mu_0, var_0, mu_alpha, var_alpha, mu_beta, var_beta, y0=0,a,b, n.iter) {

n<-length(Dose)
df.x<-Dose
n.x<-seq(1,n)
Rmx<-max(df.T)
Lmin<-min(df.T)

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

for (i in 2:n.iter) {
#Temperature calculation using Slice sampler
	run<-Slice_Run(T,mcInit[[1]]$foo_x,mcInit[[1]]$foo_y,mcInit[[1]]$hist_y,Rmx=Rmx[[1]])
	T<-run[[1]]
	L<-run[[2]]
	R<-run[[3]]
	ys<-run[[4]]
	for (j in 2:ncol(df.y)){
		foo_x<-mcInit[[j]]$foo_x
		if (ys>foo_x(T)){
			Sol.hat<-Shrink(foo_x,T,ys,L,R,R,Rmx=Rmx[[j]],Lmin=Lmin[[j]]) #shrinkage
			T<-Sol.hat[[1]]
			L<-Sol.hat[[2]]
			R<-Sol.hat[[3]]
		}
	}
#slice's end

m<-which(df.T[,1]<round(T)+0.1&df.T[,1]>round(T)-0.1)

S3<-sum((df.y[m,n.x]-alpha-beta*df.x)^2)+(y0-alpha-beta*De)^2
tau<-rgamma(1, shape=(n+a+1)/2, rate=(b+S3)/2)
mean1<-((tau*beta*(y0-alpha))+(mu_0*(1/var_0)))/(((beta^2)*tau)+(1/var_0))
var1<-1/(((beta^2)*tau)+(1/var_0))
mean2<-((tau*(sum(df.y[m,n.x])-(beta*sum(df.x))+y0-(beta*De)))+((1/var_alpha)*mu_alpha))/((tau*(n+1))+(1/var_alpha))
var2<- 1/((tau*(n+1))+(1/var_alpha))
S4<-sum((df.x)^2)+(De^2)
mean3<-((tau*(sum((df.x)*(df.y[m,n.x]))-(alpha*sum(df.x))+(y0*De)-(alpha*De)))+((1/var_beta)*mu_beta))/((tau*S4)+(1/var_beta))
var3<-1/((tau*S4)+(1/var_beta))
beta<-rnorm(1, mean=mean3, sd=sqrt(var3))
alpha<-rnorm(1, mean=mean2, sd=sqrt(var2))
De<-rnorm(1, mean=mean1, sd=sqrt(var1))
sigma2<-1/tau

mat[i,]<-c(De,sigma2,alpha,beta,T)

}
mat
}


