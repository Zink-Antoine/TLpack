
#####Analyse MC TL (slice +gibbs) d'après regression lineaire
Slice5<-
function (Dose,df.T,df.y, n.iter,inv=FALSE,burnin=1,thin=1) {
#nouvelle inconnue T scalaire

require(coda)

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
	mcInit[[j]]<-Slice.Init(df.T[,j],df.y[,j])
}

for (i in 2:n.iter) {
#calcul par slice de la température T
	run<-Slice.Run(T,mcInit[[1]]$foo.x,mcInit[[1]]$foo.y,mcInit[[1]]$hist.y,df.T)
	T<-run[1]
	L<-run[2]
	R<-run[3]
	y0<-run[4]
	for (j in 2:ncol(df.y)){
		foo.x<-mcInit[[j]]$foo.x
		foo.y<-mcInit[[j]]$foo.y
		hist.y<-mcInit[[j]]$hist.y
		if (y0>foo.x(T)){
			Sol.hat<-Shrink(foo.x,T,y0,L,R) #shrinkage
			T<-Sol.hat[1]
			L<-Sol.hat[2]
			R<-Sol.hat[3]
			if (y0>foo.x(T)){					
					run<-Slice.Run(T,foo.x,foo.y,hist.y,df.T)
					T<-run[1]
					L<-run[2]
					R<-run[3]
					y0<-run[4]
					}
		}
	}
#fin du slice


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
mat<-mat[seq(burnin,n.iter,thin),]
as.mcmc(mat)
}


