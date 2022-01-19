#' TLpack
#'
#'  Various function to analyze and plot TL measurement using MCMC code based on Gibbs and Slice samplers.
#'
#' @docType package
#' @name <TLpack>
#' @import OSLpack
#' @importFrom graphics axis box lines par plot.default text title

NULL

#' TL data of Pan figurine
#'
#' An example of actual data from a TL measurement
#'
#' This data set gives the TL growth curves for Pan figurine Br3459
#'
#' @docType data
#' @keywords datasets
#' @name TL.Pan
#' @usage data(TLpan)
#' @format A matrix 475 x 8 measurement points.
#' @source L2609193.binx TL grain fin  Pan Br4359 C2RMF78067
#'

NULL

#' Supralinearity data of Pan figurine
#'
#' An example of actual data from a supralinearity measurement
#'
#' This data set gives the supralinearity for Pan figurine Br3459
#'
#' @docType data
#' @keywords datasets
#' @name TL.sup.Pan
#' @usage data(TLsuppan)
#' @format A matrix 475 x 9 measurement points.
#' @source 2601201.binx TL supra Tanagra C2RMF77756 2-Base; Pan Br4359 C2RMF78067
#'

NULL


#' TL data of a gouge
#'
#' An example of actual data from a TL measurement
#'
#' This data set gives the TL growth curves for gouge
#'
#' @docType data
#' @keywords datasets
#' @name TL.Gouge
#' @usage data(TLgouge)
#' @format A matrix 475 x 8 measurement points.
#' @source 2305181.binx TL grain fin gouge à douille
#'

NULL

#' TL data of an Osiris figurine
#'
#' An example of actual data from a TL measurement
#'
#' This data set gives the TL growth curves for Osiris FGM0873
#'
#' @docType data
#' @keywords datasets
#' @name TL.Ser3
#' @usage data(TLser3)
#' @format A matrix 475 x 9 measurement points.
#' @source Serapeum E33273 /FGM0873 1510202.binx
#'

NULL

#' TL data of an etruscan vase
#'
#' An example of actual data from a TL measurement
#'
#' This data set gives the TL growth curves for Etruscan vase D171
#'
#' @docType data
#' @keywords datasets
#' @name TLetru
#' @usage data(TLetru)
#' @format A Risoe.BINfileData.object
#' @source 1601141.bin D171
#'

NULL

#' TL data of an etruscan vase (2)
#'
#' An example of actual data from a TL measurement
#'
#' This data set gives the TL growth curves for Etruscan vase D173.
#'
#'to find how to get Dose, df.T, df.y from Risoe.BINFile file, see example section
#'
#'
#' @docType data
#' @keywords datasets
#' @name TLetru2
#' @usage data(TLetru2)
#' @format A Risoe.BINfileData.object
#' @source 2703141.bin D173
#'
#' @examples
#' data(TLetru2)
#' table<-Lum(TLetru2,Doseb=180,alpha=FALSE,supra=FALSE)
#'
#' B<-table$b
#' N<-table$n
#'
#' table.norm<-table$b
#'
#' for (j in 1:3){
#'   for (i in 1:3)
#'   {
#'     table.norm[,i,j]<-(B[,i,j])/sum(N[seq(355,405),1,(j-1)*3+i])  #norm sur 380-430°C
#'   }
#' }
#'
#' table.données<-cbind(table.norm[,,1],table.norm[,,2],table.norm[,,3])
#' ii<-c(1,4,7,2,5,8,3,6,9)
#' table.données<-table.données[seq(1,475),ii]
#'
#' Dose<-c(0,0,0,180,180,180,360,360,360)
#' df.T<-matrix(rep(seq(26,500),9),475,9)
#' df.y<-table.données[,1:9]

NULL
