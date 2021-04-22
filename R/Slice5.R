##
#' Slice5
#'
#'
#'
#' @inheritParams Slice1
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
#' @return an (r × 4)-matrix,
#' alpha intercept
#' beta slope
#' sigma2 variance
#' T Temperature
#'
#' @export
#'
Slice5<-
function (Dose,df.T,df.y, n.iter,inv=FALSE,n.burnin=n.iter/2,n.thin=max(1,floor(n.iter-n.burnin)/10)) {

mat <- matrix(ncol=4, nrow=n.iter)
T<-300
alpha<- 1
beta<- 1
sigma2<- 1

mat[1, ] <- c(alpha,beta,sigma2, T)
Err<-function(y){var(y)*(length(y)-1)}
x.factor<-factor(Dose)
n<-length(Dose)
n.y<-seq(1,n)

mcInit<-list()
for (j in 1:ncol(df.y)){
	mcInit[[j]]<-Slice_Init(df.T[,j],df.y[,j])
}

for (i in 2:n.iter) {
  #Temperature calculation using Slice sample
	run<-Slice_Run(T,mcInit[[1]]$foo_x,mcInit[[1]]$foo_y,mcInit[[1]]$hist_y)
	T<-run[[1]]
	L<-run[[2]]
	R<-run[[3]]
	y0<-run[[4]]
	for (j in 2:ncol(df.y)){
		foo_x<-mcInit[[j]]$foo_x
		foo_y<-mcInit[[j]]$foo_y
		hist_y<-mcInit[[j]]$hist_y
		if (y0>foo_x(T)){
			Sol.hat<-Shrink(foo_x,T,y0,L,R) #shrinkage
			T<-Sol.hat[[1]]
			L<-Sol.hat[[2]]
			R<-Sol.hat[[3]]
			if (y0>foo_x(T)){
					run<-Slice_Run(T,foo_x,foo_y,hist_y)
					T<-run[[1]]
					L<-run[[2]]
					R<-run[[3]]
					y0<-run[[4]]
					}
		}
	}
#slice's end

m<-which(df.T[,1]<round(T)+0.1&df.T[,1]>round(T)-0.1)
if (inv) {
	Y<-Dose
	X<-df.y[m,n.y]
}
else{
	X<-Dose
	Y<-df.y[m,n.y]
}

Sxx<-sum((X-mean(X))^2)
Sxy<-sum((X-mean(X))*(Y-mean(Y)))
beta<-rnorm(1,mean=(Sxy/Sxx),sd=sqrt(sigma2/Sxx))

S1<-mean(Y)-(Sxy/Sxx)*mean(X)
S2<-(1/n+mean(X)^2/Sxx)
alpha<-rnorm(1,mean=S1,sd=sqrt(sigma2*S2))

SSE<-sum(tapply(df.y[m,n.y],x.factor,Err))
tau<-rgamma(1,shape=((n-2)/2),rate=(SSE/2))
sigma2<-1/tau

mat[i,]<-c(alpha,beta,sigma2,T)
}
colnames(mat)<-c("intercept","x","sigma2","Temperature")
mat<-mat[seq(n.burnin,n.iter,n.thin),]
as.mcmc(mat)
}

