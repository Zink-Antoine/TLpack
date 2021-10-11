#multiTL example from RLumModel
require(RLumModel)
require(Luminescence)

irradiation_dose <- seq(from = 0,to = 100,by = 20)
model.output <- lapply(irradiation_dose,
                       function(x){
                         sequence <- list(IRR = c(20, x, 1),
                                          #PH = c(220, 10, 5),
                                          TL=c(20,400,5))
                         data <- model_LuminescenceSignals(
                           sequence = sequence,
                           model = "Bailey2001",
                           plot = FALSE,
                           verbose = FALSE)
                         return(get_RLum(data, recordType = "TL$", drop = FALSE))
                       })
##combine output curves
TL_curve.merged <- merge_RLum(model.output)
##plot
plot_RLum(
  object = TL_curve.merged,
  xlab = "Temperature [°C]",
  ylab = "TL signal [a.u.]",
  main = "TL signal with various dose",
  legend.text = paste("dose", irradiation_dose, "Gy"),
  combine = TRUE)
##
n.pt<-length(TL_curve.merged[1]$data[,1])
n.irr<-length(irradiation_dose)
y<-x<-array(dim=c(n.pt,n.irr))
for (i in 1:n.irr){
  x[,i]<-TL_curve.merged[i]$data[,1]
  y[,i]<-TL_curve.merged[i]$data[,2]
}

multiTL<-list(
  Dose=seq(20,120,20),
  df.T=x,
  df.y=y,
  mu_b=2500,
  mu_0= 0,
  var_b=2500,
  var_0=10,
  var_y=10,
  n.iter=10
)


usethis::use_data(multiTL, overwrite = TRUE)
