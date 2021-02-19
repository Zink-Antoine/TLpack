############Stepping out################
Step.Out<-function(foo,x0,y,w,m,Rmx=500){
#stepping- out procedure
#Neal (2003) fig.3

	if(y>foo(x0)) return(print("x0 non valable"))

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
		if (R>Rmx) {R<-Rmx}#ligne de contrôle rajouté effet de bord
		}

	return(c(L,R,J,K))
}

