#' Lum
#' extracts TL measurements
#'
#' @inheritParams TL.plot
#'
#' @param Doseb0 [numeric] (**with default**)  the reference of beta irradiation (in seconds)
#' @param Dosea0 [numeric] (**with default**)  the reference of alpha irradiation (in seconds)
#' @param alpha [logical] (**with default**) TRUE if alpha measurements are inclued. Based on the number of data file, the value is corrected
#' @param supra [logical] (**with default**) TRUE if supra measurements are inclued. Based on the number of data file, the value is corrected
#' @param TypLum [character],[list] (**with default**) luminescence type "TL" = signal; "BG" = background
#'
#' @return a list object
#' @return $alpha
#' @return $supra
#' @return $b
#' glowcurve corrected from background for beta irradiation
#' @return $a
#' glowcurve corrected from background for alpha irradiation
#' @return $n
#' glowcurve corrected from background for normalisation
#' @return $Bb
#' uncorrected glowcurve for beta irradiation
#' @return $Ba
#' uncorrected glowcurve for alpha irradiation
#' @return $Bn
#' uncorrected glowcurve for normalisation
#' @return $bsup
#' glowcurve corrected from background for supralinearity irradiation
#' @return $nsup
#' glowcurve corrected from background for supralinearity nomalisation
#' @return $Bbsup
#' uncorrected glowcurve for supralinearity irradiation
#' @return $Bnsup
#' uncorrected glowcurve for supralinearity nomalisation
#' @return $CanalTemp
#' Temperature channel. channel = temperature for Risoe (SystemID=88); channel = Temperature-25 for Lexsyg reader (SYSTEMID=0)
#'
#' @export
#'
'Lum' <-
function(file,ech=1,Doseb0=90,Dosea0=90,alpha=TRUE,supra=TRUE,TypLum=c("TL","BG"),Temp=seq(26,599))
	{

	D1<-Raw.Data(file,ech,1)

	if (D1$L==36) alpha<-FALSE

	if (supra) {
		if (length(file)==2)	D2<-Raw.Data(file,ech,2)
		else supra<-FALSE
		}

	B1<-array(unlist(D1$Brut),dim=c(D1$T,D1$L))#pur matrix to calculate
	Net1<-B1[,D1$TL]-B1[,D1$BG]
	Net1[Net1<0]<-0

	#no preheat
	cycle0<-2 #TL+BG
	nbd<-9+alpha*4 #9 discs (with alpha, 13)
	cycleTLb<-cycle0*3*3 #(TL,BG)*{nat,nat+b,nat+2*b}*3
	if (alpha)cycleTLa<-cycle0*2*2 #(TL,BG)*2*{nat+a,nat+2*a}
	cycleTLn<-cycle0*nbd #(TL,BG)*nbd discs

	Dosen<-c(Doseb0)
	Doseb<-c(0,Doseb0,2*Doseb0)
	if (alpha) Dosea<-c(Dosea0,2*Dosea0)

	Disk<-seq(1,nbd)

	CanalTemp<-mapply(Canal,Temp=Temp,file=list(file[[1]]))

	if (D1$L==36){
		Brutb<-array(D1$Brut[seq(1,18)],dim=c(cycle0,3,3),dimnames=list(TypLum,Doseb,seq(1,3)))
		Brutn<-array(D1$Brut[seq(19,36)],dim=c(cycle0,1,9),dimnames=list(TypLum,Dosen,Disk))

		Netb<-array(Net1[CanalTemp,seq(1,9)],dim=c(length(Temp),3,3),dimnames=list(Temp,Doseb,seq(1,3)))
		Netn<-array(Net1[CanalTemp,seq(10,18)],dim=c(length(Temp),1,9),dimnames=list(Temp,Dosen,Disk))
		}
	else {
		Brutb<-array(D1$Brut[seq(1,18)],dim=c(cycle0,3,3),dimnames=list(TypLum,Doseb,seq(1,3)))
		Bruta<-array(D1$Brut[c(19,20,23,24,21,22,25,26)],dim=c(cycle0,2,2),dimnames=list(TypLum,Dosea,seq(4,5)))
		Brutn<-array(D1$Brut[c(seq(27,44),45,46,49,50,47,48,51,52)],dim=c(cycle0,1,13),dimnames=list(TypLum,Dosen,Disk))

		Netb<-array(Net1[CanalTemp,seq(1,9)],dim=c(length(Temp),3,3),dimnames=list(Temp,Doseb,seq(1,3)))
		Neta<-array(Net1[CanalTemp,c(10,12,11,13)],dim=c(length(Temp),2,2),dimnames=list(Temp,Dosea,seq(4,5)))
		Netn<-array(Net1[CanalTemp,c(seq(14,22),23,25,24,26)],dim=c(length(Temp),1,13),dimnames=list(Temp,Dosen,Disk))
		}

	if (supra==TRUE){
	  CanalTemp.sup<-mapply(Canal,Temp=Temp,file=list(file[[2]]))

		B2<-array(unlist(D2$Brut),dim=c(D2$T,D2$L))#matrix need to calculate

		Net2<-B2[,D2$TL]-B2[,D2$BG]
		Net2[Net2<0]<-0

		#no monitored preheat
		cycle0<-2 #TL+BG
		nbd2<-9 #9 discs
		cycleTL<-cycle0*3*3 #(TL,BG)*{b,2*b,3*b}*3

		Doseb2<-c(Doseb0,2*Doseb0,3*Doseb0)

		Disk2<-seq(1,nbd2)

		D2$Brutb<-array(D2$Brut[seq(1,18)],dim=c(cycle0,3,3),dimnames=list(TypLum,seq(1,3),Doseb2))
		D2$Brutn<-array(D2$Brut[seq(19,36)],dim=c(cycle0,1,9),dimnames=list(TypLum,Dosen,Disk2))

		Net2b<-array(Net2[CanalTemp.sup,seq(1,9)],dim=c(length(Temp),3,3),dimnames=list(Temp,seq(1,3),Doseb2))
		Net2n<-array(Net2[CanalTemp.sup,seq(10,18)],dim=c(length(Temp),1,9),dimnames=list(Temp,Dosen,Disk2))

		}

	if (!alpha){Neta<-"N/A";Bruta<-"N/A"}
	if (!supra){Net2b<-"N/A";Net2n<-"N/A";D2<-list(Brutb="N/A",Brutn="N/A");CanalTemp.sup<-"N/A"}

	Lum<-list(alpha=alpha,supra=supra,b=Netb,a=Neta,n=Netn,Bb=Brutb,Ba=Bruta,Bn=Brutn,bsup=Net2b,nsup=Net2n,Bbsup=D2$Brutb,Bnsup=D2$Brutn,CanalTemp=CanalTemp,CanalTempsup=CanalTemp.sup)
	return(Lum)
	}

##########################

#' Raw.Data
#'
#' extract only the useful data, removing the pre-annealing measurements (if need)
#'
#' @param file [Risoe.BINfileData] (**required**) the BIN/BINX file
#' @param ech [numeric] (**with default**) the sample number
#' @param n.chauf [numeric](**required**) heat experiment number (1 - first heat; 2 - supralinearity)
#'
#' @return alist object with the following elements
#' @return $Brut the data
#' @return $L number of measurements without the furnace pre-annealing
#' @return $T number of measuring points (1 point per degree)
#' @return $TL Thermoluminescence numbering
#' @return $BG background numbering
#'
#' @importFrom OSLpack ExtractFile
#'
#' @noRd
#'
Raw.Data<-function(file,ech,n.chauf) {
  file.sel<-OSLpack::ExtractFile(file=file,n_file=n.chauf)[[n.chauf]]
  L<-length(file.sel)   #number of measurements
  corr<-L%%36
  L<-L-corr   #number of measurements without the furnace pre-annealing
  T<-unique(file.sel@METADATA$NPOINT)  # number of measuring points (1 point per degree)
  if (L==72)	L<-36#two samples
  Brut<-file.sel@DATA[seq(1,L)+corr]  #Raw data from heating n.chauf after removal of the furnace pre-annealing
  if (ech==2){
    Brut<-file.sel@DATA[seq(L,2*L)+corr]
  }
  TL<-seq(1,L-1,2) # TL numbering
  BG<-seq(2,L,2) # BG numbering

  RD<-list(Brut=Brut,L=L,T=T,TL=TL,BG=BG)
  return(RD)
}


###########################################
#' Canal
#'
#' temperature channel
#'
#' @param Temp [vector]  (required)
#' @param file [Risoe.BINfileData] (**required**) single BIN/BINX file
#'
#' @return Canal channels corresponding to the temperatures.
#'
#' @noRd
#'
Canal<-function(Temp,file){
  if(unique(file@METADATA$SYSTEMID)==0)CanalTemp<-Temp-unique(file@METADATA$LOW)
  if(unique(file@METADATA$SYSTEMID)==88)CanalTemp<-Temp
  CanalTemp
}


