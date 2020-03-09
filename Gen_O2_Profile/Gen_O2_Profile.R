################################################################################
##                                           
##  Content - Gen.O2.profile
##  What does this script do? - Generate oxygen profiles in ocean sediments
##  given a flux or oxygen concentration  
##  
##  Written by: Luna Geerts
##  Contact: Luna.Geerts@hotmail.com 
##  
##  version 1.0
##
################################################################################
#NOTES
################################################################################

#the script will generate an oxygen profile given a depth in meters and the amount 
#of data points you want for this given length. Initially this script was written to compare the outputs of
#Flux.comp. Because of this the script only accepts one environmental parameter for the entire sediment
# (same for porosity).
# It is a very simple model thus only diffusive transport and consumption of oxygen (following a monod distribution)
# is taken into account.  


#The output can be directly plugged in into flux.comp !

#############################################################
# libraries needed
#############################################################

require(ReacTran)
require(marelac)

#############################################################
# UNITS
#############################################################

#C       [mmol m^-3]
#x.cor   [micrometers]
#por     [-]

#############################################################
# ASSUMPTIONS
#############################################################

#No lateral transport occurs -> thus it can be approximated by an 1D model
#There is no production 
#No bioturbation hence only diffusion causes transport of O2 in the water collumn to the sediment.
#We are modeling to achieve steady state, hence our flux.comp also assumes that the sediment is in steady state
#might be interesting though to compare model results when this is not the case.
#we assume constant porosity (and thus also turtousity ) over the sediment depth since we cannot input different porosities
#in our function flux.comp :(
#we also assume all oxygen is used up at depth "L"
#we assume the volume is constant and the area of the block of sediment we are considering is equal through the entire depth


#since diffusion is the only mode of transport the flux in is equal to ficks law of diffusion for that depth "x"
#giving us -> dJ/dx 

#since there is no production and no flux out the flux in for the whole block of sediment considered must equal
#the consumption if we reached steady state since

#dC/dt= 0 = dj/dx - consumption  


#############################################################
# limitations
#############################################################

#Consumption is according to a monod 






#############################################################
# Make simple O2 models to test my data on
#############################################################





Gen.O2.Profile<- function(N,                       # #datapoints
                          L,                       # length in METERS
                          por,                     # Porosity
                          Flux.top=NULL,           #  mmol m^-2 day^-1
                          O2.ow= NULL,             # concentration overlying water
                          S  = 35,                 # Salinity no units  
                          P  = 1.013253,           # Pressure in bars
                          TC = 25,                 # Temperature in degrees Celsius
                          Ks = NULL,               # Half saturation constant (monod)
                          R.O2 = NULL              # Oxygen consumption 
                          ) {

if(is.null(O2.ow)) {water.given<-FALSE } else {water.given<-TRUE}
  if(is.null(Flux.top&O2.ow))  {warning(paste("Neither flux or concentration water given, will assume water concentration of 250 mmol m^-3 "))}
if(!is.null(Flux.top)&!is.null(O2.ow)) {warning(paste("Both flux and concentration specified will continue with flux"))}
  

if(!is.null(Flux.top))  {Flux.top<-sqrt(Flux.top^2)} #so that the input variable is always positive


O2.model<- function(time,state,parms){
  
  with(as.list(parms),{ 
    
#we make a variable start O2 that will contain all oxygen concentrations at every depth interval dx.mid
#this will be used as input for our transport model, we basically have to start from somewhere
    
    O2             <- state[1:N] 
    
    if(!is.null(Flux.top)){
    tran.O2        <- tran.1D(C=O2,
                          C.down = 0, #since all oxygen is used up at the end of our modelled block aka no flux out
                          flux.up = Flux.top ,
                            D     = Grid.Ds,
                            VF    = Grid.por,
                            dx    = Grid)$dC
    }else if(!is.null( O2.ow)){    
      tran.O2        <- tran.1D(C=O2,
                                C.down = 0, #since all oxygen is used up at the end of our modelled block aka no flux out
                                C.up  = Flux.top ,
                                D     = Grid.Ds,
                                VF    = Grid.por,
                                dx    = Grid)$dC
    }  else {warning(paste("check line 97")) }
  #using excercise 4 as a baseline and using the monod kinetic expression
 #end of with expression  
Consumption.O2 <- - R.O2*( O2 / (O2+Ks))    
reaction.O2    <- Consumption.O2 
    
    dC.O2.dt       <- tran.O2+reaction.O2

    
    return(list(dC.O2.dt=dC.O2.dt, production =  Consumption.O2))
  }  )}
                    


# parameters

N   <- N       # essentially the amount of data points you have
L   <- L        # length 1 cm in meters
S   <- S           # if functions here or something preferably its input for diffcoeff
TC  <- TC
P   <- P
por <- por 
tort<- 1-2*log(por)    #same way we find tortuosity in flux.comp
D   <-    diffcoeff(S=S,t=TC,P=P,species="O2")[[1]]*3600*24#meters^2 d^-1
Ds  <-    D/tort


#diffcoeff in for flux has is in units m^2 s^-1, here it needs to be in micrometer^2 per year
#but diffcoeff computes m^2 /s so conversion to micrometers^2 yr^-1
#*3600*24*365.25 for a year
#*(10^6)^2=10^12  for m^2->micrometer^2 

conversion<- 1000 #to go from cm^3 to meter^3
#and go from micromol to millimol divide by 1000

if(is.null(Ks))                 {Ks   <- 0.005*conversion}# [ mmol m^-3] 
if(is.null(O2.ow)){O2.ow   <- 0.25*conversion }# [ mmol m^-3] 
if(is.null(R.O2))               {R.O2 <- 500 *conversion/365.25}# [ mmol m^-3] 


parm<-c(
O2.ow   =    O2.ow,   # [ mmol m^-3]
R.O2    =    R.O2,    # [ mmol m^-3 day^-1]
Ks      =     Ks )    # [ mmol m^-3] 


# Grid


Grid     <- setup.grid.1D ( x.up  = 0   , N    = N   , L=L)
Grid.por <- setup.prop.1D ( value = por , grid = Grid)
Grid.Ds  <- setup.prop.1D ( value = Ds  , grid = Grid)






# Steady state

#however in flux.comp we are interested in the steady state of the profile, this is one of the assumptiosn however
#Nontheless it can be interesting to test how each model performs against violation of this assumption

state<-rep(0,N)

out <- steady.1D(y=state, func=O2.model, parms=parm, nspec=1)

steady.state.reached <- attributes(out)$steady
if (steady.state.reached) {O2.SS <- out$y} else stop
if (steady.state.reached) {Prod.SS <- out$production} else stop


#in units mmol per liter or micromol per square cm


if (water.given)
  {
True.flux<-with(as.list(parm),{ 
  tran.O2        <- tran.1D(C      = O2.SS,
                            C.up   = O2.ow ,  #Concentration oxygen in water collumn
                            C.down = 0,
                            D      = Grid.Ds,
                            VF     = Grid.por,
                            dx     = Grid)
     return( round(-tran.O2$flux.up,4)) })
} else {True.flux<- -Flux.top}

Gen.oxygenprofile<-data.frame(ID     = 1,
                              x.cor  = Grid$x.mid*10^6,
                              C      = O2.SS,
                              Por    = por,
                              T      = TC,
                              S      = S,
                              P      = P,
                              Author = "Luna.gen profile",
                              R.O2= R.O2,
                              True.Flux=True.flux,
                              Production= Prod.SS )
return(Gen.oxygenprofile)

}

