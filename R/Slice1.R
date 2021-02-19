#####Analyse MC TL (slice +gibbs) d'après gibbs1.R Hickey (2006)
Slice1<-
function (Dose,df.T,df.y, mu_b, mu_0, var_b, var_0, var_y, n.iter) {
#extrapolation y0=0
#nouvelle inconnue T scalaire

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
	mcInit[[j]]<-Slice.Init(data.T[,j],data.y[,j])
}

for (i in 2:n.iter) {
#calcul par slice de la température T
	run<-Slice.Run(mcInit[[1]]$foo.x,mcInit[[1]]$foo.y,mcInit[[1]]$hist.y,data.T)
	T<-run[1]
	L<<-run[2]
	R<<-run[3]
	y0<<-run[4]
	#print(c(i,T))
	for (j in 2:ncol(data.y)){
		foo.x<-mcInit[[j]]$foo.x
		foo.y<-mcInit[[j]]$foo.y
		hist.y<-mcInit[[j]]$hist.y
		if (y0>foo.x(T)){
			Sol.hat<-Shrink(foo.x,T,y0,L,R) #shrinkage
			T<-Sol.hat[1]
			L<-Sol.hat[2]
			R<-Sol.hat[3]
			#print(c(i,T))
			if (y0>foo.x(T)){					
					run<-Slice.Run(foo.x,foo.y,hist.y,data.T)
					T<-run[1]
					L<-run[2]
					R<-run[3]
					y0<-run[4]
					}
		}
	}
#fin du slice
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

