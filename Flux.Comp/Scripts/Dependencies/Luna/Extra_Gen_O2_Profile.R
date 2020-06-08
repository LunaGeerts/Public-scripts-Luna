################################################################################
##                                           
##  Content - Auxillary functions
##  What does this script do? - Put noise on generated data (noise)
##                            - Calculate microbial consumption rate associated with a certain flux and depth
##  
##  Written by: Luna Geerts
##  Contact: Luna.Geerts@hotmail.com 
##  
##  version 0.1
##
################################################################################










noise<-function(x,sd=1){
  return(x+rnorm(x,mean=0,sd=sd)) }

R_Flux_Convert<- function(Flux,       #[ mmol m^-2 day^-1]
                          depth_m){   # meter
  temp_flux   <- Flux/depth_m         # To convert [mmol m^-2 day^-1] to [mmol m^-3 day^-1]
  return     (temp_flux*0.001*365.25) # [mmol m^-3 day^-1] to [micromol cm^-3 year^-1] (0.001 = 1000 micromol/ 1 mmol * 1 m^3/1000 000 cm^3)
}
