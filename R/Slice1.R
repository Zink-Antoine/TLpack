#####
#' Slice1
#'
#' MC Analysis TL (slice +gibbs) following gibbs1.R Hickey (2006)
#' extrapolation y0=0 (gibbs1)
#' additional scalar variable T (Slice sampler)
#'
#' @param Dose [numeric] (**required**) set of the irradiation doses
#' @param df.T [matrix] (**required**) Luminescence data: the temperatures
#' @param df.y [matrix] (**required**) Luminescence data; the luminescence signal
#' @param mu_b [numeric] (**required**) Prior of the slope's average
#' @param mu_0  [numeric] (**required**) Prior of the true dose's average
#' @param var_b [numeric] (**required**) Prior of the slope's variance
#' @param var_0 [numeric] (**required**) Prior of the true dose's variance
#' @param var_y [numeric] (**required**) Prior of the intensity's variance
#' @param n.iter [numeric] (**required**) the number of iteration
#'
#' @import Slice
#' @import stats
#'
#' @return an (n.iter x 3)-matrix
#' x0 'True' dose value
#' b slope
#' T Temperature
#'
#' @references Gibbs sampler: Hickey, G. L. 2006. « The Linear Calibration Problem: A Bayesian Analysis ». PhD Thesis, PhD dissertation, University of Durham. 1–148. http://www.dur.ac.uk/g.l.hickey/dissertation.pdf.
#' @references chapter 3 and Appendix G.1
#' @references Slice sampler: Neal, R. 2003. « Slice Sampling ». Annals of Statistics 31 (3): 705‑67.
#'
#' @export
#'
#' @examples
#' ##load data
#' data(multiTL, envir = environment())
#' attach(multiTL)
#' Slice1(Dose,df.T,df.y,mu_b,mu_0,var_b,var_0,var_y,n.iter)
#' detach(multiTL)
Slice1<-
function (Dose,df.T,df.y, mu_b, mu_0, var_b, var_0, var_y, n.iter) {

mat <- matrix(ncol=3, nrow=n.iter)
x0 <- 1
b <- 1
T<-300

mat[1, ] <- c(x0, b, T)
data.x<-seq(1,length(Dose))
data.T <- df.T
data.y<-df.y

mcInit<-list()
for (j in 1:ncol(data.y)){
	mcInit[[j]]<-Slice_Init(data.T[,j],data.y[,j])
}

for (i in 2:n.iter) {
#Temperature calculation using Slice sampler
	run<-Slice_Run(T,mcInit[[1]]$foo_x,mcInit[[1]]$foo_y,mcInit[[1]]$hist_y)
	T<-run[[1]]
	L<<-run[[2]]
	R<<-run[[3]]
	y0<<-run[[4]]
	for (j in 2:ncol(data.y)){
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
#Slice's end
#print(T)
m<-which(data.T[,1]<round(T)+0.1&data.T[,1]>round(T)-0.1)

sum_xy <- sum((Dose)*(data.y[m,data.x]))
sum_xx <- sum((Dose)^2) + x0^2
mean1 <- (mu_0/var_0)/((b*b/var_y)+(1/var_0))
var1 <- 1/((b*b/var_y)+(1/var_0))
mean2 <- ((sum_xy/var_y)+(mu_b/var_b))/((sum_xx/var_y)+(1/var_b))
var2 <- 1/((sum_xx/var_y)+(1/var_b))
x0 <- rnorm(1, mean1, sqrt(var1))
b <- rnorm(1, mean2, sqrt(var2))
mat[i, ] <- c(x0, b,T)
}
mat
}

