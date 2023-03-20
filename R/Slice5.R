##
#' Slice5
#'
#' MC Analysis TL (slice +gibbs) from linear regression
#'
#'
#'
#' @inheritParams Slice1
#' @param n.burnin [numeric] (**with default**)	number of iterations to discard at the beginning (burn in).
#'  Default is n.iter/2, that is, discarding the first half of the simulations.
#' @param n.thin [numeric] (**with default**) thinning rate. Must be a positive integer. Set n.thin > 1 to save memory and computation time if n.iter is large.
#' Default is max(1, floor((n.iter-n.burnin) / 10)) which will only thin if there are at least 20 simulations.
#' @param threshold [numeric] (**with default**) measured point below or equal to the threshold are rejected.
#' Default is the minimum measured value.
#'
#' @import Slice
#' @import coda
#'
#' @return an mcmc object,
#' \tabular{lll}{
#'  **column** \tab **Type** \tab **Description**\cr
#'  `intercept` \tab  `numeric` \tab intercept\cr
#'  `x` \tab `numeric` \tab slope \cr
#'  `std.dev` \tab  `numeric` \tab standard deviation\cr
#'  `Temperature` \tab `numeric` \tab Temperature \cr
#'  `natural dose` \tab `numeric` \tab estimated Dose value \cr
#' }
#'
#' @export
#'
#' @examples
#' ##load data
#' if(dev.cur()!=1) dev.off()
#' data(multiTL, envir = environment())
#' Dose<-multiTL$Dose
#' df.T<-multiTL$df.T
#' df.y<-multiTL$df.y
#' n.iter<-multiTL$n.iter
#' test<-Slice5(Dose,df.T,df.y,n.iter=n.iter)
#' #
#' if(dev.cur()!=1) dev.off()
#' data(TLpan, envir = environment())
#' Dose<-c(0,0,80,80,80,160,160,160)
#' df.T<-matrix(rep(seq(26,500),8),475,8)
#' df.y<-TL.Pan[,1:8]
#' Pan<-Slice5(Dose,df.T,df.y,n.iter=100)
#' #
#' \dontrun{
#' data(TLetru)
#' table<-Lum(TLetru,Doseb0=360,Dosea=0,alpha=FALSE,supra=FALSE)
#' B<-table$b
#' N<-table$n
#' table.norm<-B
#'  for (j in 1:3){
#'    for (i in 1:3)
#'      {
#'        table.norm[,i,j]<-(B[,i,j])/sum(N[seq(350,400),1,(j-1)*3+i])
#'        }
#'    }
#'  table.data<-cbind(table.norm[,,1],table.norm[,,2],table.norm[,,3])
#'  ii<-c(1,4,7,2,5,8,3,6,9)
#'  table.data<-table.data[seq(1,475),ii]
#'  Dose<-as.numeric(colnames(table.data))
#'  df.T<-matrix(rep(seq(26,500),9),475,9)
#'  df.y<-table.data[,1:9]

#'  plot(seq(26,500),table.data[,8])
#'  Slice5(Dose,df.T,df.y,n.iter=10)
#'  }
#'
Slice5<-
function (Dose,df.T,df.y,
            n.iter,n.burnin=n.iter/2,n.thin=max(1,floor(n.iter-n.burnin)/10),
              threshold=min(df.y[df.y>0])) {

#parameters
#x.factor<-factor(Dose)
n<-length(Dose)
n.y<-seq(1,n)
Rmx<-apply(df.T,2,max)
Lmin<-apply(df.T,2,min)


#variables: initial values
mat <- matrix(ncol=5, nrow=n.iter)
mcInit<-list()
for (j in ncol(df.y):1){
  repeat{
    mcInit[[j]]<-Slice_Init(df.T[,j],df.y[,j])
    Tj<-mcInit[[j]]$x0
    foo_x<-mcInit[[j]]$foo_x
    if (foo_x(Tj)>threshold) break
  }
}

alpha_n<-alpha<- 1
beta_n<-beta<- 1
var<- 1
Sxx<-sum((Dose-mean(Dose))^2)
var1_n<-var1<- (var/beta^2)*(1/n+(mean(Dose)+(-alpha/beta)^2)/Sxx)
T_n<-T<-mcInit[[1]]$x0
ch<-mcInit[[1]]$hist_y$mids
threshold<-min(ch)+(max(ch)-min(ch))/length(ch)
De_n<-De<-0
mat[1, ] <- c(alpha,beta,var, T,De)
print(threshold)

for (i in 2:n.iter) {
  #Temperature calculation using Slice sampler
  repeat{
    if(mcInit[[1]]$foo_x(T)<threshold){T<-mcInit[[1]]$x0}
    run<-Slice_Run(T,mcInit[[1]]$foo_x,mcInit[[1]]$foo_y,mcInit[[1]]$hist_y,Rmx=Rmx[[1]])
# print(i)
#   print(mcInit[[1]]$foo_x(run[[1]])>threshold)
#
    if (mcInit[[1]]$foo_x(run[[1]])>threshold) break

  }
  T<-run[[1]]
  L<-run[[2]]
  R<-run[[3]]
  y0<-run[[4]]
  for (j in 2:ncol(df.y)){
    foo_x<-mcInit[[j]]$foo_x
    if (y0>foo_x(T)){
      Sol.hat<-Shrink(foo_x,T,y0,L,R,R,Rmx=Rmx[[j]],Lmin=Lmin[[j]]) #shrinkage
      T<-Sol.hat[[1]]
      L<-Sol.hat[[2]]
      R<-Sol.hat[[3]]
    }
    #print(j)
  }
  #slice's end

  m<-which(df.T[,1]==round(T))

	X<-Dose
	Y<-df.y[m,n.y]

  Sxx<-sum((X-mean(X))^2)
  Sxy<-sum((X-mean(X))*(Y-mean(Y)))
  beta<-rnorm(1,mean=(Sxy/Sxx),sd=sqrt(var/Sxx))

  S1<-mean(Y)-beta*mean(X)
  S2<-(1/n+mean(X)^2/Sxx)
  alpha<-rnorm(1,mean=S1,sd=sqrt(var*S2))

  SE<-(Y-alpha-beta*X)^2
  SSE<-sum(SE)
  var<-SSE/(n-2)

  var1<- (var/beta^2)*(1/n+(mean(X)+(-alpha/beta)^2)/Sxx)

  De<-rnorm(1,alpha/beta,sqrt(var1))

  #slope condition
  u<-rexp(1,rate = 15000/(mean(Dose)))
  if(u<abs((De-De_n)/(T-T_n))){
    var1<-var1_n;
    De<-De_n
    beta<-beta_n
    alpha<-alpha_n
    T<-T_n
  }

  alpha_n<-alpha
  beta_n<-beta
  De_n<-De
  var1_n<-var1
  T_n<-T


  mat[i,]<-c(alpha,beta,sqrt(var1),round(T),De)
}
colnames(mat)<-c("intercept","x","std dev","Temperature","natural dose")
mat<-mat[seq(n.burnin+1,n.iter,n.thin),]
mcmc(mat,start=n.burnin+1,end=n.iter,thin=n.thin)
}

