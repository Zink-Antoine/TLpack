#' multiTL_fun
#'
#' multiTL example from RLumModel
#'
#' @param irradiation_dose [numeric][list] (**with default**) irradiation dose
#' @param model [character] (**with default**) model as RLumModel
#' @param distri_scatter [character] (**with default**) a distribution function to scatter the glow curves.
#' @param par_scatter [numeric][list] (**with default**) scatter parameters associated to scatter distribution.
#' @param distri_noise [character] (**with default**) a distribution function to noise the glow curve.
#' @param par_noise [numeric][list] (**with default**) noise parameter associated to the noise distribution.
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
  function(irradiation_dose=c(20,20,20,40,40,40,60,60,60),model = "Bailey2001",
           distri_scatter="runif(1,s,t)",par_scatter=list(s=0.8,t=1.2),
           distri_noise="runif(n,a,b)",par_noise=list(a=0.9,b=1.1)){
#The simulations were performed at 20 s, considered as the natural irradiation, and at 40 and 80 s,
#corresponding respectively to Nat +20 s and Nat + 40 s.


model.output <- lapply(irradiation_dose, function(x){
  sequence <- list(IRR = c(20, x, 0.1),
                   #PH = c(220, 10, 5),
                   TL=c(20,400,5))
  data <- model_LuminescenceSignals(
    sequence = sequence,
    model,
    plot = FALSE,
    verbose = FALSE,
    simulate_sample_history = TRUE )
  return(get_RLum(data, recordType = "TL$", drop = FALSE))
})

##combine output curves
TL_curve.merged <- merge_RLum(model.output)
get_TL_curve.merged<-get_RLum(TL_curve.merged)

n.pt<-length(get_TL_curve.merged[[1]]@data [,1])
n.irr<-length(irradiation_dose)

##plot
plot_RLum( object = TL_curve.merged,
           xlab = paste("Temperature (\u00B0","C)",sep=""),
           ylab = "TL signal [a.u.]",
           main = "TL signal with various dose",
           legend.text = paste("dose", irradiation_dose, "Gy"),
           combine = TRUE)

y<-x<-array(dim=c(n.pt,n.irr))

scatter<-function(distri_scatter,par_scatter){
  dscatter<-parse(text=distri_scatter)
  eval(dscatter, par_scatter)
}

noise<-function(distri_noise,par_noise){
  dnoise<-parse(text=distri_noise)
  eval(dnoise, c(list(n=n.pt),par_noise))
}

## scattering and/or noising glow curves
for (i in 1:n.irr){
  x[,i]<-round(get_TL_curve.merged[[i]]@data [,1],1)
  y[,i]<-get_TL_curve.merged[[i]]@data [,2]*scatter(distri_scatter,par_scatter)*noise(distri_noise,par_noise)#
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
       legend = paste("dose",unique(irradiation_dose)," Gy"),col=c(2,3,4),lty=1)

multiTL<-list( Dose=irradiation_dose-irradiation_dose[1], df.T=x, df.y=y, n.iter=10 )

}
