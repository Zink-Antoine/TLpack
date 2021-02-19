
###### présentation de quantité final (output) MC########
######            GUM suppl1                    ##########

Output.Summary <- 
function (output,p){
require(data.table) #pour la fonction shift
output.sort<-sort(output,decreasing = FALSE)
#test croissance absolue
m<-length(output.sort)
	while(length(m)>0){
	output.diff<-output.sort-shift(output.sort,n=1)
	m<-which(output.diff==0) #indice où l'écart avec le nombre précédent est nul
	output.sort[m]<-output.sort[m]+1e-128 #inclu une perturbation minime
}

#calcul coverage
M<-length(output)
q<-p*M

#symmétrique
#r<-(M-q)/2
#low<-r
#high<-r+q
#cov<-output.sort[high]-output.sort[low]

#shortest
cov.star<-M
for (r in 1:(M-q)){
low<-r
high<-r+q
cov<-output.sort[high]-output.sort[low]
if(cov<cov.star) {
cov.star<-cov
low.star<-low
high.star<-high
}

}

return(list(output.sort=output.sort,esp.sort=mean(output.sort),dev.sort=sd(output.sort),cov.star=cov.star,low.cov=output.sort[low.star],high.cov=output.sort[high.star]))
}

