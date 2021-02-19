#####d'après gibbs3.R Hickey (2006) appliqué aux TL
Slice3<-
function (Dose,df.T,df.y, mu_0, var_0, mu_alpha, var_alpha, mu_beta, var_beta, y0=0,a,b, n.iter) {
#nouvelle inconnue T scalaire

mat <- matrix(ncol=5, nrow=n.iter)
T<-300
alpha<- 1
beta<- 1
sigma2<- 1
x0<- 1

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
	ys<-run[4]
	for (j in 2:ncol(df.y)){
		foo.x<-mcInit[[j]]$foo.x
		foo.y<-mcInit[[j]]$foo.y
		hist.y<-mcInit[[j]]$hist.y
		if (ys>foo.x(T)){
			Sol.hat<-Shrink(foo.x,T,ys,L,R) #shrinkage
			T<-Sol.hat[1]
			L<-Sol.hat[2]
			R<-Sol.hat[3]
			if (ys>foo.x(T)){					
					run<-Slice.Run(foo.x,foo.y,hist.y,df.T)
					T<-run[1]
					L<-run[2]
					R<-run[3]
					ys<-run[4]
					}
		}
	}
#fin du slice

m<-which(df.T[,1]<round(T)+0.1&df.T[,1]>round(T)-0.1)

S3<-sum((df.y[m,n.x]-alpha-beta*df.x)^2)+(y0-alpha-beta*x0)^2
tau<-rgamma(1, shape=(n+a+1)/2, rate=(b+S3)/2)
mean1<-((tau*beta*(y0-alpha))+(mu_0*(1/var_0)))/(((beta^2)*tau)+(1/var_0))
var1<-1/(((beta^2)*tau)+(1/var_0))
mean2<-((tau*(sum(df.y[m,n.x])-(beta*sum(df.x))+y0-(beta*x0)))+((1/var_alpha)*mu_alpha))/((tau*(n+1))+(1/var_alpha))
var2<- 1/((tau*(n+1))+(1/var_alpha))
S4<-sum((df.x)^2)+(x0^2)
mean3<-((tau*(sum((df.x)*(df.y[m,n.x]))-(alpha*sum(df.x))+(y0*x0)-(alpha*x0)))+((1/var_beta)*mu_beta))/((tau*S4)+(1/var_beta))
var3<-1/((tau*S4)+(1/var_beta))
beta<-rnorm(1, mean=mean3, sd=sqrt(var3))
alpha<-rnorm(1, mean=mean2, sd=sqrt(var2))
x0<-rnorm(1, mean=mean1, sd=sqrt(var1))
sigma2<-1/tau

mat[i,]<-c(x0,sigma2,alpha,beta,T)

}
mat
}

