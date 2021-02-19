######### mc procedure###################

Slice.mc<-function(x,y,n.iter=1000,burn.in=n.iter/2){
	mcInit<-Slice.Init(x,y)
	attach(mcInit)
	mcSlice<-alist(x1=0,L=0,R=0)
	for (i in 1:n.iter){
	run<-Slice.Run(foo.x,foo.y,hist.y,x)
	mcSlice$x1[i]<-run[1]
	mcSlice$L[i]<-run[2]
	mcSlice$R[i]<-run[3]
	}
	detach(mcInit)
	return(mcSlice)
}

