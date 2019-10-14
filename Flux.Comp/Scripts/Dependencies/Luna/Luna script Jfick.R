################################################################################
##                                           
##  Content - Luna script Jfick
##  What does this script do? - Calculate jfick using a slope between two dots
##  Written by: Luna Geerts
##  Contact: Luna.Geerts@hotmail.com 
##  
##  version 1.6
##
################################################################################
##  VERSION HISTORY               
##----------------------------------------------------------------------------                    
##  Created script 9/9/2019

##  10/9/2019 - Made it so it does not matter anymore if you select a point deeper in
##  the sediment first or one higher in the sediment. aka order of point selection does
## not matter anymore

##  13/9/2019 - Made code more user friendly giving prompts when certain values that need to be defined are missing
##  eg porosity lacking gives a prompt and still calculates the slope regardless but with no flux estimate
##  Function is more customisable allowing multiple Porosities and diffusioncoefficients since conditions might change
##  depending on depth/place
##  Point selection on plots for slope calculation is more streamlined. No need to press escape to go to the next plot
##  Only one profile can be calculated or multiple, if multiple are present the dataframe should contain a collumn 
##  with ID names
##  BUG fixes:"There needs to be at least 2 points between the two points selected" -> should be fixed!


## 16/9/2019 - Finalised the code !

## 17/9/2019 - fixed bug where ID's not starting from 1 would not work

## 20/9/2019 - Added option to save plots in a folder

## 30/9/2019 - Removed the need to attach df.int as a result if code stops working midway and df.int is still attached
#              No error should pop up anymore, nor should there be an issue if there is an existing value called C or x


## ADDITIONAL NOTES / SUGGESTIONS
##----------------------------------------------------------------------------                    

## Perhaps an iterative function could be written that calculates a slope for every given 
## interval of point this could maybe reduce the subjective nature of the jfick flux estimation? 


## KNOWN BUGS
##----------------------------------------------------------------------------
## If the input terms for the function Jfick are not quoted the function cannot find
## the collumn name since instead it will give the values that are associated with said
## name.

## Secondly if the dataframe in "Dataframename" contains a collumn that has the same
## name than "Depth.name" or "Oxygen", the function does not work

## If in point selection the user would click on the plot somewhere in the top left and next bottem right
## then the script gets stuck.
## This is because when the user clicks on the plot, the x coordinates get used to determine which point they
## likely wanted to select and the closest point to this x coordinate gets used further for the slope calculations.
## I included an if statement that essentially checks whether the user clicked first on a 
## point on the left hand of the oxygen profile(which is deeper in the sediment thus a higher y value)
## a dot right from this first clicked point or vice versa. If not for this the user would always have
## to click the deepest point in the sediment first and
## then the second point would have to be one higher in the sediment.
## The if statement essentially checks whether:
## the y coordinates of the first clicked point are smaller or bigger than the second point.
## this works fine if the user clicks somewhat near the line plotted. 
## Issues arise however when the user clicks for example on the top left of the profile and thene the bottem right. 
## See code for more details


################################################################################
# Set working directory and load packages
################################################################################

#setwd("C:/Users/install/Dropbox/Thesis/Raw data")

#install.packages("DescTools")
library(DescTools)
#install.packages("marelac")
library(marelac)




################################################################################
#Input needed for Jfick / manual
################################################################################

#The following parameters should be given to the function Jfick
# 1 -First, Dataframename, the name of the dataframe containing the oxygenprofile(s)
#    that fluxes need to be calculated for.

# 2 -ID.check is a check if multiple oxygen profiles within the same dataframe have to be considered
#    eg, if one dataframe contains multiple oxygen profiles these should be seperated first via a collumn "ID"
#    This can be done via the function "assign ID" which considers also ID's of other datasets.

# 3 -Depth.name includes the name of the collumn with depths (these can be zero corrected or not)
#    See function "zerodepthcor" for more info

# 4 -Oxygen includes the name of the collumn with the oxygen concentrations
#    UNITS: Micromol per liter

#    The following variables can be given as a collumn name or a constant, when a collumn name is given then each 
#    of the following parameters will use the constant given at their respective ID

# 5 -Tort can be given two options depending on how tortuosity should be calculated
#    0 (and default) Tortuosity = ( 1-log(porosity)) 
#    1               Tortuosity = (1-2*log(porosity))

# 6 - Porosity can be a fixed value or should include the name of the collumn containing the 
#     porosities in "dataframename".

# 7 - A diffusioncoefficient can be given if this is desired, if not then the temperature and species need to be given
#     See ?diffcoeff for more info

# 8 - Species input, species needed for diffusioncoefficient calculation needs to go here, if left blank O2 is assumed
#     are at standard conditions see ?diffcoeff for more info

# 9- Temperature, input of temperature for diffusioncoefficient calculation if left as a constant or collumn name
#     blank, 20 degrees will be assumed. See ?diffcoeff for more info.

# 10 & 11- Pressure and salinity if left blank are assumed to be 35 (UNITS?) and 1.013253 in bar 
#          See ?diffcoeff for more info. Can be filled in as a collumn name or constant

# 12- FLIPPER is for integration with the function "FLIPPER" which can calculate gradient fluxes as well, to ensure the 
#     same points are selected in both Jfick will use the same points the user indicated in FLIPPER. The object the results of
#     FLIPPER are written to should be put here

# 13- Plot.save.path is an optional input, if you want you can specify an existing folder to save all your files in
#     example Plot.save.path="./Save here/" would save all the plots in a folder "save here" located in your working directory

################################################################################
#Jfick code + Jfickflux.correction code
################################################################################

#This first function we call upon later in our jfick code, however if for whatever reason you can't calculate fluxes
#Easily in a given dataset but can calculate it's slopes, or a correction needs to be done, you can use this function
#If slopes are known. Other than that it should not be used by the user

Jfickflux.correction<-function(Slope,Porosity=NA,
                               Diffusioncoeff=NA,
                               Tort=NULL,
                               species.input=NA,
                               Temperature = NA,
                               Pressure = NA,
                               Salinity = NA
                               ){
  
#If no salinity or pressure is given use these values  
  if  (is.na(Salinity)) {Salinity=35}
  if  (is.na(Pressure)) {Pressure=1.013253}
  
  
  
#A whole lot of if statements, to give different prompts depending on if the user gave a diffusion coefficient or not
#and if temperature or species were given, if not the function still continues but will make certain assumptions
  if  (  !is.na(Temperature) 
         & !is.na(species.input)
         & !is.na(Diffusioncoeff)) {print("The diffusion coefficient given will be used, Temperature and species variables won't be considered")
    Diffusioncoeff=Diffusioncoeff
    
  }else if (   !is.na(Temperature) 
               &!is.na(species.input)
               & is.na(Diffusioncoeff)) {Diffusioncoeff<-diffcoeff(t=Temperature, species = species.input,S=Salinity,P=Pressure)[[1]]
               
  }else if (is.na(Temperature) &
            is.na(Diffusioncoeff)&
            !is.na(species.input))   {print("No temperature is given, 20 degrees temperature will be assumed")
    Temperature    <- 20 
    Diffusioncoeff <- diffcoeff(t=Temperature, species = species.input)[[1]]
    
  }else if (is.na(species.input) &
            is.na(Diffusioncoeff)&
            !is.na(Temperature))   {print("No species is given, Oxygen will be assumed") 
    Diffusioncoeff<-diffcoeff(t=Temperature, species = "O2" )[[1]]
  }else if (is.na(species.input) &
            is.na(Temperature)&
            is.na(Diffusioncoeff))    {print("No parameters to calculate diffusioncoefficient were given") 
    print("And no diffusioncoefficient was given")
    print("Hence standard conditions at 20 degrees for oxygen was used, see '?diffcoeff'")
    Diffusioncoeff<-diffcoeff(t=20,species="O2")
    
    
  }else {print("Something went wrong with wrong during diffusioncoefficient determination, make sure either")
    print("either a diffusion coefficient or temperature and the species are given.")}
  
 
#unit adjustements
  
  eenheden<-86400/(10^-6)
  por<-as.vector(Porosity)
  
   if (is.null(Tort)) {  t<-(1-log(por)) 
  } else  if (Tort==1)       {  t<-(1-2*log(por))
  } else {stop("Please enter a valid number for tortuosity or leave blank") } 
  
  result<-Slope*por*(1/t)*Diffusioncoeff*eenheden
  return(result)
}

Jfick     <-function(Dataframename,
                     ID.check = NULL, 
                     Depth.name = NA,
                     Oxygen = NA,
                     Tort = NULL,
                     Por.cte = NA,
                     Porosity = NA,
                     Diffusioncoeff = NA,
                     species.input = NA,
                     Temperature = NA,
                     Pressure = NA,
                     Salinity = NA,
                     FLIPPER = NULL,
                     Plot.save.path = NULL) {

  
#Added these if statements in case someone forgets about entering either depth or oxygen in brackets since without these
  #The code wont function

  
  if (is.character(Depth.name)==FALSE){
    stop("Enter Depth.name in brackets")}
  
  if (is.character(Oxygen)==FALSE){
    stop("Enter Oxygen in brackets")
    }
  
  
#Here we create some variables we need later on in our function
  
  f                  <-NULL
  
  ID.total=1
 
  
#This part checks if the user wants multiple profiles (within the same dataset seperated by an ID) checked   
  if (is.null(ID.check)){
              ID       = 1
              ID.total = 1
   
    #Instructions should only be printed if FLIPPER is not present, since when it is present no points should be selected          
    if (is.null(FLIPPER)){          
    print("Only one profile will be estimated, Write 'ID.check=1'for multiple ID comparison ")
    print("Click on two spots on the curve,inbetween these points a slope will be calculated")}
    
  } else if (ID.check==1) {
    ID.total           <- sort(unique(Dataframename$ID))
    print(paste("Estimation of profiles:",length(ID.total)))
    if (is.null(FLIPPER)){ 
    print("Click on two spots on the curve, inbetween these points a slope will be calculated.")
    print("After selecting two points the next plot will be created")}
    
  } else{
    print("Please leave ID.check open or write 1")
    stop()
    }
  
#Here we create variables we need later in our for loop
  
  
  Coordinatelist     <-1:length(ID.total)
  Coordinatelist     <-as.data.frame(Coordinatelist)
  
  mat                <-1:length(ID.total)
  mat                <-as.data.frame(mat)
  mat$ID             <-ID.total
  
  mat$por            <-NA
  mat$Diffusioncoeff <-NA
  mat$Tort           <-NULL
  mat$species.input  <-NA
  mat$Temperature    <-NA
  mat$Pressure       <-NA
  mat$Salinity       <-NA
  
  
  
  for        (i in ID.total) {
    
#The get function will search for a variable with the given name, these are extracted and put in a temporary dataframe
  
    df.int<-NULL
    df.int<-subset (Dataframename, ID == i , select = c(get(Depth.name), get((Oxygen))))
    
    colnames(df.int) <- c("x","C")
    

    
#using our temporary dataframe we will construct the plot on which the user selects two points for slope calculation
#However we need to construct two linear regressions since we do not want the slope where concentration is in function of
#Depth but the other way around (reversed to what the plot makes you believe)
    plot(df.int$x~df.int$C, ylim=rev(range(df.int$x)), 
         main = paste("Oxygen profile from ID",i),
         xlab = "", 
         ylab = bquote("Depth"~mu~"meter"))
    mtext(side=1, line=2.5, bquote(mu~mol ~O[2]~l^-1), cex=1.0)
    
    
    abline (1,0)
    
    points(df.int$C[df.int$x==0], df.int$x[df.int$x==0], col="green",
           cex = 2, pch = 21, bg = "red")
    if (!is.null(FLIPPER)){
      
      
      if(isTRUE(FLIPPER$method=="all")){
        
      FLIP<-as.data.frame(FLIPPER$output$gradient.output$fit)
        
      loc   <- list(y = c (FLIP$x[1] , FLIP$x [length(FLIP$x)] ),
                    x = c (FLIP$C[1] , FLIP$C [length(FLIP$C)] ))
      loc$y <- loc$y*10^6}
      
      if(isTRUE(FLIPPER$method=="gradient")){
        
        FLIP<-as.data.frame(FLIPPER$output$fit)
        
        loc  <- list(y = c (FLIP$x[1] , FLIP$x [length(FLIP$x)] ),
                     x = c (FLIP$C[1] , FLIP$C [length(FLIP$C)] ))
        loc$y<- loc$y*10^6}
      
    } else {loc       <- locator(n=2)}
    
    
    
#This "if" part I included so it does not matter anymore if someone would select a dot that is at a
#lower depth or at a higher depth, however it still gets stuck if the user clicks on a point first that is higher 
#on the y axis but has an x value coresponding with it lower than the second selected point. This is due to my if
# function looking at the closest point on the graph for evaluation but rather where the user clicked on the point
#In the lm part of the function no point exists where point1[2] is then > point2[2] which is the reason for the error.
    
    if (loc$y[2]>loc$y[1]) {
      
      temploc<-loc$y
      loc$y[1]<-temploc[2]
      loc$y[2]<-temploc[1]
      
      templocx<-loc$x
      loc$x[1]<-templocx[2]
      loc$x[2]<-templocx[1]
    }
    
    
    
    Depth2    <- c((Closest(df.int$C,loc$x[1])),
                   (Closest(df.int$C,loc$x[2])))

# I use the x values since these are further apart and facilitates easier selection
# Here point1 is a vector with on place 1 the y coordinates and on place two the x coordinates

    
    point1     <- NULL
    point2     <- NULL
    point1     <- c(Depth2[1],df.int$x[df.int$C==Depth2[1]]) 
    point2     <- c(Depth2[2],df.int$x[df.int$C==Depth2[2]])
    
    
    points(point1[1], point1[2], col="blue", cex=1.5, pch=21, bg="blue")
    points(point2[1], point2[2], col="blue", cex=1.5, pch=21, bg="blue")
    
# Here I make a regression model that only considers the values between the points selected.
# But we need to make two regressions, one for the plot and one for the flux estimate since for flux we need the 
# slope of how much oxygen changes when depth changes by one micrometer
    
    visualfit <-  lm(x~C,
                     data =df.int[which(df.int$x <= point1[2] & df.int$x >= point2[2]),])
    
    abline(visualfit,col="red",lty=2,lwd=2)
   
#we do not want to plot the legend if FLIPPER is to be compared (within the same plot or in R), if you still want to save
    #the plot then the legend will be plotted regardless
  if (is.null(FLIPPER)){
    legend("bottomright",legend=c("Depth at 0","Points used for slope calculation","Slope")
           ,pch=20 ,lwd=c(0,0,2) ,col=c(2,4,2) ,lty=c(0,0,2))}
    
    fit <- lm(C~x,
              data =df.int[which(df.int$x <= point1[2] & df.int$x >= point2[2]),])
    
  
  
    summary(fit)
    slope  <- (summary(fit))$coefficients[2,1]


    
    mat$Slope[mat$ID==i]      <- slope  
 
  
    
 
   
     mat$Jfick.flux[mat$ID==i] <- NA
    
    
#This part of the code is to include collumn names in the variables given. For example
# Temperature= "T" if this is the case then for every oxygen flux the Temperature will be used for that ID 
# this option is also possible for datasets without an ID numbering but that have their constants in a seperate collumn
# If no collumn is given then simply the constant is used, if neither is supplied then nothing is filled in and 
# assumptions will be made

      
      if(is.character(Porosity)){
        mat$por[mat$ID==i]<-get(Porosity,Dataframename)[Dataframename$ID==i][1]}
      if(is.numeric(Porosity)){
        mat$por[mat$ID==i]<-Porosity
        }
      
      if(is.character(Diffusioncoeff)){
      mat$Diffusioncoeff[mat$ID==i]  <- get(Diffusioncoeff,Dataframename)[Dataframename$ID==i][1]}
      if(is.numeric(Diffusioncoeff)){
        mat$Diffusioncoeff[mat$ID==i]<- Diffusioncoeff  
      }
    
      if(is.character(Tort)){
      mat$Tort[mat$ID==i]            <- get(Tort,Dataframename)[Dataframename$ID==i][1]}
      if(is.numeric(Tort)){
        mat$Tort[mat$ID==i]          <- Tort
      }
      
      if(is.character(species.input)){
        mat$species.input<-species.input
      }
    
      if (is.character(Temperature)){
      mat$Temperature[mat$ID==i]    <- get(Temperature,Dataframename)[Dataframename$ID==i][1]}
      if (is.numeric(Temperature)){
        mat$Temperature[mat$ID==i]  <- Temperature
      }
    
    
      if (is.character(Pressure)){
      mat$Pressure[mat$ID==i]    <- get(Pressure,Dataframename)[Dataframename$ID==i][1]}
      if (is.numeric(Pressure)){
        mat$Pressure[mat$ID==i]  <- Pressure
      }
    
      if (is.character(Salinity)){
      mat$Salinity[mat$ID==i]    <- get(Salinity,Dataframename)[Dataframename$ID==i]}
      if (is.numeric(Salinity)){
        mat$Salinity[mat$ID==i]  <- Salinity
      }
    
#in case Porosity wasnt given or known, slopes can still be calculated but no diffusion can be given
      if (is.na(mat$por[1]))      {print("Slopes calculated but no flux calculated since no porosity was given")}
    
    
    Coordinatelist$ID                                   <- ID.total
    Coordinatelist$X1.coordin.Cor[Coordinatelist$ID==i] <- Closest(df.int$C,loc$x[1])[1]
    Coordinatelist$X2.coordin.Cor[Coordinatelist$ID==i] <- Closest(df.int$C,loc$x[2])[1]
    Coordinatelist$Y1.coordin.Cor[Coordinatelist$ID==i] <- Closest(df.int$x,loc$y[1])[1]
    Coordinatelist$Y2.coordin.Cor[Coordinatelist$ID==i] <- Closest(df.int$x,loc$y[2])[1]
    
#This part is to print it on the plot.
    
    
    
    
        temp.result <- Jfickflux.correction(
        Slope          = mat$Slope[mat$ID==i],
        Porosity       = mat$por[mat$ID==i],
        Diffusioncoeff = mat$Diffusioncoeff[mat$ID==i],
        Tort           = mat$Tort[mat$ID==i],
        species.input  = mat$species.input[mat$ID==i],
        Temperature    = mat$Temperature[mat$ID==i],
        Pressure       = mat$Pressure[mat$ID==i],
        Salinity       = mat$Salinity[mat$ID==i])
        value.plot    <- paste("Flux / R.int= ",round(temp.result,digits=3))
        
        text(0,1000,value.plot,  pos = 4)
        
if (!is.null(Plot.save.path)) { 
          
          path<-file.path(Plot.save.path,paste("Oxygen profile from ID",i,".png",sep = ""))
          
          png(filename = path )
          
          plot(x~C, ylim=rev(range(x)), 
               main = paste("Oxygen profile from ID",i),
               xlab = bquote(mu~mol ~O[2]~l^-1), 
               ylab = bquote("Depth"~mu~"meter"))
          
          abline (1,0)
          
          points(C[x==0], x[x==0], col="green",
                 cex = 2, pch = 21, bg = "red")
          points(point1[1], point1[2], col="blue", cex=1.5, pch=21, bg="blue")
          points(point2[1], point2[2], col="blue", cex=1.5, pch=21, bg="blue")
          abline(visualfit,col="red",lty=2,lwd=2)
          
          legend("bottomright",legend=c("Depth at 0","Points used for slope calculation","Slope")
                 ,pch=20 ,lwd=c(0,0,2) ,col=c(2,4,2) ,lty=c(0,0,2))
          
          
          text(0,1000,value.plot,  pos = 4)
          
       
        }  
           
        try(detach(df.int),silent=TRUE)
    
    
    
    
  }
  
 for (i in ID.total) {
   mat$Jfick.flux[mat$ID==i] <- Jfickflux.correction(
                            Slope          = mat$Slope[mat$ID==i],
                            Porosity       = mat$por[mat$ID==i],
                            Diffusioncoeff = mat$Diffusioncoeff[mat$ID==i],
                            Tort           = mat$Tort[mat$ID==i],
                            species.input  = mat$species.input[mat$ID==i],
                            Temperature    = mat$Temperature[mat$ID==i],
                            Pressure       = mat$Pressure[mat$ID==i],
                            Salinity       = mat$Salinity[mat$ID==i])
}
  
#we only wish to retain ID,porosity, Slope and the Jfick.flux
  
  mat<-mat[,c("ID","por","Slope","Jfick.flux")]
  Coordinatelist<-Coordinatelist[,-1]
 
  
  return(list(Results     = mat,
              Coordinates = Coordinatelist)) 
  }


#returns a list with two dataframes, one containing the slopes and Jfick fluxes for every oxygen profile.
#Second dataframe contains the coordinates of every point clicked on.
  






################################################################################
#Example of script running
################################################################################

#The loaded file must include a collumn name named ID

#load("./Data/Test data for scripts/TestdataJfick.rda")
#df1<-df10[[1]]
#df2<-df1[df1$ID==2|df1$ID==3|df1$ID==1,]


#f<-NULL
#f<-Jfick(df2,
#     ID.check= 1, 
#    Depth.name="depth.cor",
#   Porosity = 0.8,
#  Oxygen="C",
# Temperature = 20 ,
# species.input = "O2")



#does work, this is because the variables need to be passed as names
#Jfick(df2,"depth.cor","C",0.85,diffcoeff(t=20,species = "O2")[[1]])


# EXAMPLE OF THE BUG
#colnames(df2)[1]<-"Depth.name"

#also does not work since dataframe df2 contains a collumn with the name "Depth.name" or "Oxygen"
#Jfick(df2,"depth.cor","C",0.85,diffcoeff(t=20,species = "O2")[[1]])


