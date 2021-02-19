############Stepping out################
#' Step.Out
#'
#' Slice sampler: stepping out procedure
#' Neal(2003) fig.3
#'
#' @param foo [function] (**required**) function proportional to the density
#' @param x0 [numeric] (**required**) the current point
#' @param y [numeric] (**required**) the vertical level defining the slice
#' @param w [numeric] (**required**) estimate of the typical size of a slice
#' @param m [numeric] (**required**) integer limiting the size of a slice to mw
#' @param Rmx [numeric] (**required**) bound value to limit edge effects (added)
#'
#' @import stats
#'
#' @return
#' @export
#'
#'
Step.Out<-function(foo,x0,y,w,m,Rmx=500){

	if(y>foo(x0)) return(print("x0 invalid"))

	U<-runif(1,0,1)
	L<-x0-w*U
	R<-L+w
	V<-runif(1,0,1)
	J<-floor(m*V)
	K<-(m-1)-J
	if (L<30){L<-30}

	while(J>0 && y<foo(L)){
		L<-L-w
		J<-J-1
		if (L<30){L<-30}
		}

	while(K>0 && y<foo(R)&& R<Rmx){
		R<-R+w
		K<-K-1
		if (R>Rmx) {R<-Rmx} #control line added to limit edge effects
		}

	return(c(L,R,J,K))
}

