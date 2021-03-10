#' TL.plot
#'
#' plot for TL measurements
#'
#' @param file [Risoe.BINfileData] (**required**) the BIN/BINX file(s)
#' @param nomFile [character] (**with default**) name of the BIN/BINX file
# #' @param File [list] (**required**) an output obtained by [OSLpack::ReadFile] function, i.e a list object incluing [Risoe.BINfileData] objects and the names of the corresponding BIN/BINX files
#' @param ech [numeric] (**with default**) the sample number
#' @param Doseb0 [numeric] (**with default**)  the reference of beta irradiation (in seconds)
#' @param Dosea0 [numeric] (**with default**)  the reference of alpha irradiation (in seconds)
#' @param supra [logical] (**with default**) TRUE if supra measurements are inclued. Based on the number of data files, the value is corrected
#' @param norm [logical] (**with default**) TRUE if the measurements are normalized.
#' @param plateau [numeric],[list] (**with default**) the plateau range for normalization
#' @param Temp [numeric],[list] (**with default**) the temperature range
#' @param NumInv [character] (**with default**) the inventory number
#' @param b.check [numeric],[list] (**with default**) list of beta curves to plot
#' @param a.check [numeric],[list] (**with default**) list of alpha curves to plot
#' @param sup.check [numeric],[list] (**with default**) list of supralinarity curves to plot
#'
#' @return plot figures
#'
#' @importFrom graphics axis box lines par plot.default text title
#'
#' @export
#'

'TL.plot'<-
function(file,nomFile="",ech=1,Doseb0=90,Dosea0=90,supra=FALSE,norm=FALSE,plateau=seq(200,500),Temp=seq(26,599),NumInv="",b.check=rep(1,9),a.check=rep(1,4),sup.check=rep(1,9)) {

	alpha<-Dosea0>1
	plateau.corr<-mapply(Canal,Temp=plateau,file=list(file[[1]]))
	if (supra) 	plateau.corr.sup<-mapply(Canal,Temp=plateau,file=list(file[[2]]))

	L<-Lum(file,ech=ech,Doseb0=Doseb0,Dosea0=Dosea0,alpha=alpha,supra=supra,Temp=Temp)

	par(mfrow=c(2,2),mar=c(3,3,2,2))

	Ib<-L$b
	if (alpha) Ia<-L$a
	if (supra) Isup<-L$bsup
	if (norm){
		for (i in 1:3){
			for (j in 1:3){
				Ib[,j,i]<-Ib[,j,i]/sum(L$n[plateau.corr,,(i-1)*3+j])
				if (supra) {Isup[,j,i]<-Isup[,j,i]/sum(L$nsup[plateau.corr.sup,,(i-1)*3+j])}
			}
		}
		if (alpha)
		{
			for (i in 1:2){
				for (j in 1:2){
					Ia[,j,i]<-Ia[,j,i]/sum(L$n[plateau.corr,,(i-1)*2+9+j])
				}
			}
		Ia[NaN]<-0
		}
		Ib[NaN]<-0
		if (supra) Isup[NaN]<-0
	}

	plot(Temp,xlim=c(20,500),ylim=c(0,max(Ib[seq(175,450),,])*1.25),col=1,type="l",axes=FALSE)
	axis(1,tcl=-0.2,padj=-2,cex.axis=0.7)
	axis(2,tcl=-0.2,padj=2,cex.axis=0.7)
	box()
	title(xlab=paste("Temperature (\u00B0","C)",sep=""),ylab="Luminescence (ua)",cex.lab=0.7,line=1)
	text(50,max(Ib[seq(200,450),,])*1.05,"Irradiation beta",cex=1,adj=0)
	for (i in 1:3){
		for (j in 1:3){
	   lines(Temp,Ib[,j,i],col=j*b.check[3*(i-1)+j])
		}
	}

	if (alpha)
	{
		plot(Temp,xlim=c(20,500),ylim=c(0,max(Ia[seq(175,450),,])*1.25),col=1,type="l",axes=FALSE)
		axis(1,tcl=-0.2,padj=-2,cex.axis=0.7)
		axis(2,tcl=-0.2,padj=2,cex.axis=0.7)
		box()
		title(xlab=paste("Temperature (\u00B0","C)",sep=""),ylab="Luminescence (ua)",cex.lab=0.7,line=1)
		text(50,max(Ia[seq(200,450),,])*1.05,"Irradiation alpha",cex=1,adj=0)
		for (i in 1:3){
			lines(Temp,Ib[,1,i],col=1*b.check(i))
		}
		for (i in 1:2){
			lines(Temp,Ia[,1,i],col=2*a.check[i])
			lines(Temp,Ia[,2,i],col=3*a.check[2+i])
		}
	}

	if (supra)
	{
		plot(Temp,xlim=c(20,500),ylim=c(0,max(Isup[seq(175,450),,])*1.25),col=1,type="l",axes=FALSE)
		axis(1,tcl=-0.2,padj=-2,cex.axis=0.7)
		axis(2,tcl=-0.2,padj=2,cex.axis=0.7)
		box()
		title(xlab=paste("Temperature (\u00B0","C)",sep=""),ylab="Luminescence (ua)",cex.lab=0.7,line=1)
		text(50,max(Isup[seq(200,450),,])*1.05,paste("Supralin\u00e9","arit\u00e9",sep=""),cex=1,adj=0)
		for (j in 1:3){
			for (i in 1:3){
				lines(Temp,Isup[,j,i],col=i*sup.check[3*(i-1)+j])
			}
		}
	}

	plot.default(c(0,100), c(0,100), type="n", axes=FALSE,ylab="", xlab="")

	position<-c(-5,30,70)
	couleur<-c("black","red","green")
	beta.text<-c("beta : Nat","Nat+10.4Gy","Nat+21,0Gy")
	alpha.text<-c("alpha : Nat",paste("Nat+10.6\u00b5","m-2",sep=""),paste("Nat+21.1\u00b5","m-2",sep=""))
	supra.text<-c("supra : 10.4Gy","     21.0Gy","31.4Gy")

	text(50,90,paste("n? inv. ",NumInv))
	text(50,80,"Thermoluminescence")
	text(50,70,paste("1re chauffe : fichier ",nomFile[[1]][1]))
	text(position,60,beta.text,col=couleur,pos=4,cex=0.9)
	if (alpha) text(position,50,alpha.text,col=couleur,pos=4,cex=0.9)
	if (supra) {text(50,40,paste("2nde chauffe : fichier ",nomFile[[1]][2]))
	text(position,30,supra.text,col=couleur,pos=4,cex=0.9)}

}
