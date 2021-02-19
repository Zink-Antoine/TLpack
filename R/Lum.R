############################# fonction Lum ####################################
`Lum` <-
function(file,ech=1,Doseb0=90,Dosea0=90,alpha=TRUE,supra=TRUE,TypLum=c("TL","BG"),Temp=seq(26,599)) #détermination la table des valeurs normalisées des OSL ech= numéro échantillon, OSL= type OSL
	{

	D1<-Raw.Data(file,ech,1)

	if (supra) {
		if (length(file)==2)	D2<-Raw.Data(file,ech,2)
		else supra<-FALSE
		}
	if(file[[1]]@METADATA$SYSTEMID[1]==0)Temp.corr<-Temp-file[[1]]@METADATA$LOW[1]
	if(file[[1]]@METADATA$SYSTEMID[1]==88)Temp.corr<-Temp

	B1<-array(unlist(D1$Brut),dim=c(D1$T,D1$L))#pur matrice pour calcul

	Net1<-B1[,D1$TL]-B1[,D1$BG]
	Net1[Net1<0]<-0

	#pas de préchauffe enregistrée
	cycle0<-2 #TL+BG
	nbd<-9+alpha*4 #9 disques, 13 si alpha
	cycleTLb<-cycle0*3*3 #(TL,BG)*{nat,nat+b,nat+2*b}*3 
	if (alpha)cycleTLa<-cycle0*2*2 #(TL,BG)*2*{nat+a,nat+2*a} 
	cycleTLn<-cycle0*nbd #(TL,BG)*13 disques

	Dosen<-c(Doseb0)
	Doseb<-c(0,Doseb0,2*Doseb0)
	if (alpha) Dosea<-c(Dosea0,2*Dosea0)

	Disk<-seq(1,nbd)

	if (D1$L==36){
		alpha<-FALSE

		Brutb<-array(D1$Brut[seq(1,18)],dim=c(cycle0,3,3),dimnames=list(TypLum,Doseb,seq(1,3)))
		Brutn<-array(D1$Brut[seq(19,36)],dim=c(cycle0,1,9),dimnames=list(TypLum,Dosen,Disk))

		Netb<-array(Net1[Temp.corr,seq(1,9)],dim=c(length(Temp),3,3),dimnames=list(Temp,Doseb,seq(1,3)))
		Netn<-array(Net1[Temp.corr,seq(10,18)],dim=c(length(Temp),1,9),dimnames=list(Temp,Dosen,Disk))
		}
	else {
		Brutb<-array(D1$Brut[seq(1,18)],dim=c(cycle0,3,3),dimnames=list(TypLum,Doseb,seq(1,3)))
		Bruta<-array(D1$Brut[c(19,20,23,24,21,22,25,26)],dim=c(cycle0,2,2),dimnames=list(TypLum,Dosea,seq(4,5)))
		Brutn<-array(D1$Brut[c(seq(27,44),45,46,49,50,47,48,51,52)],dim=c(cycle0,1,13),dimnames=list(TypLum,Dosen,Disk))

		Netb<-array(Net1[Temp.corr,seq(1,9)],dim=c(length(Temp),3,3),dimnames=list(Temp,Doseb,seq(1,3)))
		Neta<-array(Net1[Temp.corr,c(10,12,11,13)],dim=c(length(Temp),2,2),dimnames=list(Temp,Dosea,seq(4,5)))
		Netn<-array(Net1[Temp.corr,c(seq(14,22),23,25,24,26)],dim=c(length(Temp),1,13),dimnames=list(Temp,Dosen,Disk))
		}

	if (supra==TRUE){
		if(file[[2]]@METADATA$SYSTEMID[1]==0)Temp.corr<-Temp-file[[2]]@METADATA$LOW[1]
		if(file[[2]]@METADATA$SYSTEMID[1]==88)Temp.corr<-Temp

		B2<-array(unlist(D2$Brut),dim=c(D2$T,D2$L))#pur matrice pour calcul

		Net2<-B2[,D2$TL]-B2[,D2$BG]
		Net2[Net2<0]<-0

		#pas de préchauffe enregistrée
		cycle0<-2 #TL+BG
		nbd2<-9 #9 disques
		cycleTL<-cycle0*3*3 #(TL,BG)*{b,2*b,3*b}*3 

		Doseb2<-c(Doseb0,2*Doseb0,3*Doseb0)

		Disk2<-seq(1,nbd2)

		D2$Brutb<-array(D2$Brut[seq(1,18)],dim=c(cycle0,3,3),dimnames=list(TypLum,seq(1,3),Doseb2))
		D2$Brutn<-array(D2$Brut[seq(19,36)],dim=c(cycle0,1,9),dimnames=list(TypLum,Dosen,Disk2))

		Net2b<-array(Net2[Temp.corr,seq(1,9)],dim=c(length(Temp),3,3),dimnames=list(Temp,seq(1,3),Doseb2))
		Net2n<-array(Net2[Temp.corr,seq(10,18)],dim=c(length(Temp),1,9),dimnames=list(Temp,Dosen,Disk2))

		}

	if (!alpha){Neta<-"N/A";Bruta<-"N/A"}
	if (!supra){Net2b<-"N/A";Net2n<-"N/A";D2<-list(Brutb="N/A",Brutn="N/A")}

	Lum<-list(alpha=alpha,supra=supra,b=Netb,a=Neta,n=Netn,Bb=Brutb,Ba=Bruta,Bn=Brutn,bsup=Net2b,nsup=Net2n,Bbsup=D2$Brutb,Bnsup=D2$Brutn,Temp.corr=Temp.corr)
	return(Lum)
	}

##########################
Raw.Data<-function(file,ech,n.chauf) {
	file.sel<-RW.File(file,n.chauf,file[[n.chauf]]@METADATA$SEL)
	L<-length(file.sel)  #nombre de mesures
	corr<-L%%36
	L<-L-corr  #nombre de mesures corrigé des préchauffes du four
	T<-file.sel@METADATA$NPOINT[1] #nombre de point de mesure (1pt par degré)
	if (L==72)	L<-36#deux échantillons
	Brut<-file.sel@DATA[seq(1,L)+corr] #données brutes issue de la chauffe n° n.chauf après élimination des préchauffes du four
	if (ech==2){		
			Brut<-file.sel@DATA[seq(L,2*L)+corr]
			}
	TL<-seq(1,L-1,2) # numérotation des TL
	BG<-seq(2,L,2) #numérotation des BF

	RD<-list(Brut=Brut,L=L,T=T,TL=TL,BG=BG)	
	return(RD)
}


#############reconstitution#########################
#reconstitue un fichier théorique à partir 
#d’un fichier contenant une séquence inusuelle

RW.File<-function(file,nfile,sequence){
file.n<-file[[nfile]]
file.n@DATA<-file.n@DATA[sequence]
file.n@METADATA<-file.n@METADATA[sequence,]
file.n@.RESERVED<-file.n@.RESERVED[sequence]
return(file.n)
}

