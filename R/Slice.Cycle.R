######### mc procedure plusieures courbes ###################

Slice.Cycle<-function(x,y,n.iter=1000,burn.in=n.iter/2){
#sequentiel
	mcInit<-list()
	for (j in 1:ncol(y)){
		mcInit[[j]]<-Slice.Init(x[,j],y[,j])
	}

	attach(mcInit[[1]])
	mcSlice<-alist(x1=0,L=0,R=0)
	for (i in 1:n.iter){
		run<-Slice.Run(foo.x,foo.y,hist.y,x)
		x1<-run[1]
		L<-run[2]
		R<-run[3]
		y0<-run[4]
		for (j in 2:ncol(y)){
			detach()
			attach(mcInit[[j]])
			if (y0>foo.x(x1)){
				Sol.hat<-Shrink(foo.x,x1,y0,L,R) #shrinkage
				x1<-Sol.hat[1]
				L<-Sol.hat[2]
				R<-Sol.hat[3]
				if (y0>foo(x1){					
					run<-Slice.Run(foo.x,foo.y,hist.y,x)
					x1<-run[1]
					L<-run[2]
					R<-run[3]
					y0<-run[4]
					}
			}
		}

		mcSlice$x1[i]<-x1
		mcSlice$L[i]<-L
		mcSlice$R[i]<-R
		mcSlice$iter[i]<-i
	}
	detach()
	return(mcSlice)
}

