#####Analyse MC TL (slice +gibbs) d'après gibbs4.R Hickey (2006)
Slice4<-
function (Dose,df.T,df.y, k=1, n.iter) {
#nouvelle inconnue T scalaire

mat <- matrix(ncol=5, nrow=n.iter)
T<-300
alpha<- 1
beta<- 1
sigma2<- 1
x0<-1
t<- 1


mat[1, ] <- c(x0, sigma2,alpha,beta, T)
n<-length(Dose)
df.x<-Dose
n.x<-seq(1,n)

mcInit<-list()
for (j in 1:ncol(df.y)){
	mcInit[[j]]<-Slice.Init(df.T[,j],df.y[,j])
}

for (i in 2:n.iter) {

#calcul par slice de la température T
	run<-Slice.Run(mcInit[[1]]$foo.x,mcInit[[1]]$foo.y,mcInit[[1]]$hist.y,df.T)
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
					run<-Slice.Run(foo.x,foo.y,hist.y,df.T)
					T<-run[1]
					L<-run[2]
					R<-run[3]
					y0<-run[4]
					}
		}
	}
#fin du slice

m<-which(df.T[,1]<round(T)+0.1&df.T[,1]>round(T)-0.1)

sxx<-sum((df.x-mean(df.x))^2)
S1<-sum(df.y[m,n.x]-beta*df.x)+y0-(beta*x0)
alpha<-rnorm(1,mean=(S1/(n+1)),sd=sqrt(sigma2/(n+1)))
S2<-sum(df.x*(df.y[m,n.x]-alpha))+(y0*x0)-(alpha*x0)
S4<-sum((df.x)^2)+(x0^2)
beta<-rnorm(1,mean=(S2/S4),sd=sqrt(sigma2/S4))
S3<-sum((df.y[m,n.x]-alpha-beta*(df.x))^2)+((y0-alpha-beta*x0)^2)
tau<-rgamma(1,shape=((n+k)/2),rate=(S3/2))
sigma2<-1/tau
u<-((n+1)*sxx)+(n*((x0-mean(df.x))^2))
t<-rgamma(1,shape=0.5,rate=u)
mu1<-((2*sigma2*n*t*(mean(df.x)))+(beta*(y0-alpha)))/(2*sigma2*n*t+beta^2)
var1<-sigma2/(2*sigma2*n*t+beta^2)
x0<-rnorm(1,mu1,sqrt(var1))

mat[i,]<-c(x0,sigma2,alpha,beta,T)

}
mat
}
