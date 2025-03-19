#' multiTL example from RLumModel
# version sous forme de fonction
#'
#' @param irradiation_dose [list] (**required**) irradiation dose
#'
#' @return  a dataset containing multiples TL calculated with RLumModel
#'
#' @importFrom graphics axis box lines par plot.default text title
#'
#' @import RLumModel
#' @import Luminescence
#'
#' @examples
#'
#' \dontrun{
#'
#' }
#'
require(RLumModel)
require(Luminescence)
'multiTL_fun'<-
  function(irradiation_dose=c(200,200,200,400,400,400,600,600,600)){
#The simulations were performed at 200 sβ, considered as the natural irradiation, and at 400 and 800 sβ,
#corresponding respectively to Nat +200 sβ and Nat + 400 sβ.


model.output <- lapply(irradiation_dose, function(x){
  sequence <- list(IRR = c(20, x, 0.1),
                   #PH = c(220, 10, 5),
                   TL=c(20,400,5))
  data <- model_LuminescenceSignals(
    sequence = sequence,
    model = "Bailey2001",
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
           xlab = "Temperature [°C]",
           ylab = "TL signal [a.u.]",
           main = "TL signal with various dose",
           legend.text = paste("dose", irradiation_dose, "Gy"),
           combine = TRUE)

y<-x<-array(dim=c(n.pt,n.irr))

for (i in 1:n.irr){
  x[,i]<-round(get_TL_curve.merged[[i]]@data [,1],1)
  y[,i]<-get_TL_curve.merged[[i]]@data [,2]*runif(1,0.8,1.2)
}

##plot
plot(x[,c(7,8,9)],y[,c(7,8,9)],
     col=4,
     type="l",
     xlab = "Temperature [°C]",
     ylab = "TL signal [a.u.]",
     main = "TL signal with various dose"
)
lines(x[,c(4,5,6)],y[,c(4,5,6)],col=3)
lines(x[,c(1,2,3)],y[,c(1,2,3)],col=2)
legend("topright",
       legend = paste("dose",unique(irradiation_dose)," Gy"),col=c(2,3,4),lty=1)

multiTL<-list( Dose=c(0,0,0,200,200,200,400,400,400), df.T=x, df.y=y, n.iter=10 )

}
