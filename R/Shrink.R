##############shrinkage########################
Shrink<-function(foo,x0,y0,L,R,m=100,Rmx=475,Lmin=200){
#shrinkage procedure
#d'après Neal(2003) Fig.5

	L.hat<-L
	R.hat<-R

	V<-runif(1,0,1)#ajout
	J<-floor(m*V) #ajout


	while (J>=0){
		U<-runif(1,0,1)
		x1<-L.hat+U*(R.hat - L.hat)
		if (x1>Rmx) {x1<-Rmx}#ligne de contrôle rajouté effet de bord
		if (x1<Lmin) {x1<-Lmin}#ligne de contrôle rajouté effet de bord
		if (y0<foo(x1)) break
		if (x1<x0) {L.hat<-x1}
		else {R.hat<-x1}
		J<-J-1
	}

return(c(x1,L.hat,R.hat))
}


