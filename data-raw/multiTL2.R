#multiTL_2 example from RLumModel
require(RLumModel)
require(Luminescence)

irradiation_dose <- c(0,0,0,20,20,20,40,40,40)
model.output <- lapply(irradiation_dose,
                       function(x){
                         sequence <- list(IRR = c(20, x, 0.1),
                                          #PH = c(220, 10, 5),
                                          TL=c(20,400,5))
                         data <- model_LuminescenceSignals(
                           sequence = sequence,
                           model = "Bailey2001",
                           plot = FALSE,
                           verbose = FALSE,
                          )
                         return(get_RLum(data, recordType = "TL$", drop = FALSE))
                       })


##combine output curves
TL_curve.merged <- merge_RLum(model.output)

n.pt<-length(TL_curve.merged[1]$data[,1])
n.irr<-length(irradiation_dose)


##plot
plot_RLum(
  object = TL_curve.merged,
  xlab = "Temperature [°C]",
  ylab = "TL signal [a.u.]",
  main = "TL signal with various dose",
  legend.text = paste("dose", irradiation_dose, "Gy"),
  combine = TRUE)
##

y<-x<-array(dim=c(n.pt,n.irr))

for (i in 1:n.irr){
  x[,i]<-round(TL_curve.merged[i]$data[,1],1)
  y[,i]<-TL_curve.merged[i]$data[,2]*runif(1,0.75,1.25)
}

plot(x[,c(7,8,9)],y[,c(7,8,9)],col=4,type="l")
lines(x[,c(4,5,6)],y[,c(4,5,6)],col=3)
lines(x[,c(1,2,3)],y[,c(1,2,3)],col=2)

multiTL2<-list(
  Dose=c(0,0,0,20,20,20,40,40,40),
  df.T=x,
  df.y=y,
  mu_b=2500,
  mu_0= 0,
  var_b=2500,
  var_0=10,
  var_y=10,
  n.iter=10
)


usethis::use_data(multiTL2, overwrite = TRUE)
