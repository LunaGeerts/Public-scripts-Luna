################################################################################
##                                           
##  Content - plot flux.comp full output
##  Written by: Luna Geerts
##  Contact: Luna.Geerts@hotmail.com 
##  
##  Script takes full output from Flux.comp and makes a plot of it
##
################################################################################

#plot_IDs<-names(Real_R$Full_Output_FLIPPER)

#Filepathplot= "Plots"
#Plot_type = "png"

source("./Scripts/Dependencies/Dependency.R")
plot_full_output<- function(Data,IDs, Filepathplot,Plot_type="png",closeplot=TRUE,Gen.O2=FALSE){

Real_R_Full<- Data
  
for (i in IDs ){


try({
  x11(  width  = 100,height = 60)
current.dev <- dev.cur()
plot.FLIPPER(Real_R_Full$Full_Output_FLIPPER[[i]])


input.Depth <- Real_R_Full$Full_Output_FLIPPER[[i]]$input$user.input$x
input.C     <- Real_R_Full$Full_Output_FLIPPER[[i]]$input$user.input$C
not.real.fit<- data.frame(x= Real_R_Full$Full_Output_PROFILE[[i]]$Depth *10^-2, #To correct back to meters
                          C= Real_R_Full$Full_Output_PROFILE[[i]]$Concentration)

input.prod     <-  data.frame(Prod=Real_R_Full$Full_Output_PROFILE[[i]]$Production,
                              depth=Real_R_Full$Full_Output_PROFILE[[i]]$Depth*10^-2) #To correct back to meters
input.prod     <-  input.prod[!duplicated(input.prod$Prod),]


PROFILE.flux   <- as.numeric(Real_R_Full$Results$Fluxes$Discrete.Berg[Real_R_Full$Results$Fluxes$ID == as.numeric(i)])

plot.discrete (input.Depth,
               input.C,
               not.real.fit,
               input.prod,
               R.int=PROFILE.flux)

if (Gen.O2){
plot.continuous(depth = Real_R_Full$Full_Output_FLIPPER[[i]]$input$continuous.input$x ,
                conc  = Real_R_Full$Full_Output_FLIPPER[[i]]$input$continuous.input$C,
                modelfit = Real_R_Full$Full_Output_TRUE[[i]],
                R.int = NULL,TRUE_FLUX = TRUE)
}
name<-paste("/plots ID",i,".",Plot_type,sep = "")
#Still make a function where one can specify where the plots should be written to for now its a junk folder to not
#overwrite past plots and if no folder is given that it should just keep the plot open for easy comparison

path<-paste(Filepathplot ,sep="")
pathing<-paste(path,name,sep="")

savePlot(filename=pathing,type = Plot_type)
if (closeplot){graphics.off() }
  
})
}
}
