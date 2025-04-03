#' multiTL_fun
#'
#' multiTL example from RLumModel
#'
#' @inheritParams RLumModel::model_LuminescenceSignals
#'
#' @param dose [numeric][list] (**with default**) irradiation dose in the sequence
#' @param distri_scatter [character] (**with default**) a distribution function to scatter the glow curves.
#' @param distri_noise [character] (**with default**) a distribution function to noise the glow curve.
#' @param SeqType [character] (**with default**) the type of measurement TL or OSL
#'
#' @import RLumModel
#' @import Luminescence
#'
#' @return  a dataset containing multiples TL calculated with RLumModel
#'
#' @export
#'
#' @examples
#'
#' \dontrun{
#' multiTL_fun()
#' }
#'

'multiTL_fun'<-
  function(dose=c(20,20,20,40,40,40,60,60,60),model = "Bailey2001",
           distri_scatter="runif(1,min=0.8,max=1.2)",
           distri_noise="rnorm(n,mean=1,sd=0.1)",
           SeqType="TL=c(20,400,5)"){

temp<-strsplit(SeqType,"=")[[1]]
Type<-paste(temp[1])
Param<-paste(temp[2])

model.output <- lapply(dose, function(x){
  sequence <- list(IRR = c(20, x, 0.1)
                  )


  eval(substitute(
    var <- val,
    list(var=str2lang(paste("sequence$",Type,sep="")),val=str2lang(Param)
    )))

  data <- model_LuminescenceSignals(
    sequence = sequence,
    model,
    plot = FALSE,
    verbose = FALSE,
    simulate_sample_history = TRUE )
  return(get_RLum(data, recordType = paste(Type,"$",sep=""), drop = FALSE))
})

##combine output curves
TL_curve.merged <- merge_RLum(model.output)
get_TL_curve.merged<-get_RLum(TL_curve.merged)

n.pt<-length(get_TL_curve.merged[[1]]@data [,1])
n.irr<-length(dose)

##plot
plot_RLum( object = TL_curve.merged,
           xlab = paste("Temperature (\u00B0","C)",sep=""),
           ylab = "TL signal [a.u.]",
           main = "TL signal with various dose",
           legend.text = paste("dose", dose, "Gy"),
           combine = TRUE)

y<-x<-array(dim=c(n.pt,n.irr))

scatter<-function(distri_scatter){
  eval(parse(text=distri_scatter))
}

noise<-function(distri_noise){
  dnoise<-parse(text=distri_noise)
  eval(dnoise, c(list(n=n.pt)))
}

## scattering and/or noising glow curves
for (i in 1:n.irr){
  x[,i]<-round(get_TL_curve.merged[[i]]@data [,1],1)
  y[,i]<-get_TL_curve.merged[[i]]@data [,2]*scatter(distri_scatter)*noise(distri_noise)#
}


##plot
plot(x[,c(7,8,9)],y[,c(7,8,9)],
     col=4,
     type="l",
     xlab = paste("Temperature (\u00B0","C)",sep=""),
     ylab = "TL signal [a.u.]",
     main = "TL signal with various dose"
)
lines(x[,c(4,5,6)],y[,c(4,5,6)],col=3)
lines(x[,c(1,2,3)],y[,c(1,2,3)],col=2)
legend("topright",
       legend = paste("dose",unique(dose)," Gy"),col=c(2,3,4),lty=1)

multiTL<-list(Dose=dose-dose[1], df.T=x, df.y=y, n.iter=10 )

}
