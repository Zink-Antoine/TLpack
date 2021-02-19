#### cycle de base pour slice sampling #######################

Slice.Init<-function(x,y){
	#fonction affichant les données mesurées (x,y) ety établissant les fonctions correspondantes y(x) et p(y)

	if (dev.cur()==1) plot(x,y) else points(x,y)#affichage des données (x,y)
	foo.x<-approxfun(x,y) #transformation des données en une fonction y(x)
	lines(x,foo.x(x),col=3)#affichage de la fonction y(x)

	hist.y<-hist(y,breaks=100,plot=FALSE) #établissement d'une courbe (histogramme) de densité des y
	foo.y<-approxfun(hist.y$mids,hist.y$density) #transformation en fonction p(y)
	return(list(hist.y=hist.y,foo.y=foo.y,foo.x=foo.x))
}

Slice.Run<-function(x0,foo.x,foo.y,hist.y,x,w=20,m=20){
	#fonction déterminant le slice et le nouveau point x1
	repeat{
	y0<-Fn3(runif(1),hist.y$mids,foo.y)#tirage aléatoire d'un y0 de la distribution p(y)
	if(y0<foo.x(x0)) break
	#x0<-Fn3(runif(1),x,foo.x)#tirage aléatoire d'un x0 de la distribution y(x)
	#print(c(x0,y0,foo.x(x0)))
	}

	Sol<-Step.Out(foo.x,x0,y0,w,m) #stepping-out
	lines(Sol[1:2],c(y0,y0),col=2)
	Sol.hat<-Shrink(foo.x,x0,y0,Sol[1],Sol[2]) #shrinkage
	lines(Sol.hat[2:3],c(y0,y0),col=6)

	return(c(Sol.hat,y0)) #données finales x1, L-hat, R-hat, y0
}

###### Cumulative distribution functions #####
###CDF pour une fonction foo
Fn1<-function(x,foo){
CDF<- CDF0<-0
for (i in x) {
CDF0<-CDF0+foo(i)
CDF<-c(CDF,CDF0)
}
return(CDF/CDF0)
}

###probabilité pour x de densité foo
Fn2<-function(p,x,foo){
CDF<-CDF0<-pCDF<-0
for (i in x) {
CDF0<-CDF0+foo(i)
CDF<-c(CDF,CDF0)
if (i<p) pCDF<-CDF0
}
return(pCDF/CDF0)
}


####quantile q pour la distri foo
Fn3<-function(q,x,foo){
CDF0<-0
i<-min(x)
norm<-sum(foo(x))
step<-(max(x)-min(x))/length(x)
while (i<max(x) && CDF0<q) {
CDF0<-CDF0+foo(i)/norm
i<-i+step
}
return(i)
}


