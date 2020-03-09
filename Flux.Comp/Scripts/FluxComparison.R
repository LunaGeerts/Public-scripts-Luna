################################################################################
##                                           
##  Content - Flux comparison
##  What does this script do? - Calculate flux according to the "gradient" method (listed as 'LunaGradient' and 
##  FlipperGradient), the discrete method (As "FlipperDiscrete" and Discrete.Berg) and the continuous method
##  (FLipperSavGolay)
##  
##  Written by: Luna Geerts
##  Contact: Luna.Geerts@hotmail.com 
##  
##  version 0.2
##
################################################################################
#NOTES
################################################################################
#Due to how the script works only one environmental value is considered, what this means is that if you have different
#temperatures at different depths then the script won't run, this is mainly a limitation of FLIPPER and PROFILE
#These two methods only allow calculations with a constant Diffusivity (for which you need temperature salinity etc)
#See "?diffcoeff" for more info


################################################################################
#NEW IN V2
################################################################################
#You can now specify if you want to save the plot files as a png or EPS 
#Also Flux.Comp now works properly with Gen_O2_Profile and can create an additional plot with the true flux and
#consumption given that the appropriate consumptions and fluxes are given!


#NOTA BINNING



################################################################################
#Bugs
################################################################################

#If a varaiable T is named or used we get an error related to PROFILE since there we say intern=T and if we assign
#a value to T... yeaaahhh


################################################################################
#TO DO
################################################################################
#error handling for when multiple enviornmental parameters are filled in

################################################################################
#Input needed for Flux.Comp / manual
################################################################################
# A dataframe with oxygenprofiles with collumn:
# ID ,where seperated from one another via ID's.
# A collumn "x.cor", if no correction of depth has to take place simply name the depth collumn "x.cor" 
# UNITS: MICROMETERS
# C with the oxygen concentration UNITS:mmol m^-3

#this is the bare minimum, it is advised you also provide a collumn with:
#T , temperature (only first value is considered)
#Por , porosity   (only first value is considered)

#Further in IDs.to.analyse a vector of ID's can be specified, for these ID's flux comp will be done.
#if left blank all ID's will be analysed
#Filepathplot, the file path where plots should be saved, if no path is specified it will be saved in your working directory
#Filepathresults, the file path where results should be saved, if no path is specified it will be saved in your working directory
#if Filepathplot=FALSE then no matrix will be saved
#PROFILE.files, if set to true, keep output files from PROFILE and writes them to your current working directory



Flux.CompV2<-function(Data.frame.name,
                    IDs.to.analyse=NULL,
                    Filepathplot=NULL,
                    Filepathresults=NULL,
                    PROFILE.files=TRUE,
                    Species="O2", 
                    Gen_O2_int= TRUE ){


#Get the dependencies out of the way  
source("./Scripts/Dependencies/Dependency.R")


#First to check which conditions are met in Flux.comp  


if (is.null(IDs.to.analyse)) {IDs.to.analyse=unique(Data.frame.name$ID)}    
if (is.null(Filepathplot))   {Filepathplot=getwd()}
if (is.null(Filepathresults)){Filepathresults=paste( getwd(),
                                                     "/Results ID",
                                                     min(Data.frame.name$ID),
                                                     "-",max(Data.frame.name$ID), sep="")}      
  
  
  df.temp  <-as.data.frame(Data.frame.name)

#Now to give an error if any of the environmental parameters per ID have (which stretch over the whole ID)

#I have the feeling this can be made much more efficient but not sure how...  
for (i in IDs.to.analyse) {
  
  
  if( length(  unique( df.temp$TC  [df.temp$ID==i] ))>1) {stop(paste("Temperature in ID",i,"contains multiple values")) }
  if( length(  unique( df.temp$S  [df.temp$ID==i] ))>1) {stop(paste("Salinity in ID",i,"contains multiple values")) }
  if( length(  unique( df.temp$P  [df.temp$ID==i] ))>1) {stop(paste("Pressure in ID",i,"contains multiple values")) }
  if( length(  unique( df.temp$z  [df.temp$ID==i] ))>1) {stop(paste("Z in ID",i,"contains multiple values")) }
  if( length(  unique( df.temp$Dmol  [df.temp$ID==i] ))>1) {stop(paste("Dmol in ID",i,"contains multiple values")) }
  if( length(  unique( df.temp$Ds  [df.temp$ID==i] ))>1) {stop(paste("Ds in ID",i,"contains multiple values")) }
  
  
  }

#Creation of matrix that stores results
  
  mat<-matrix(NA,ncol=7,nrow = length(IDs.to.analyse))
  mat<-as.data.frame(mat) 
  names<-c("ID","LunaGradient","FlipperGradient","FlipperDiscrete","FlipperSavGolay","Author","Discrete.Berg")
  colnames(mat)<-(names)
  mat$ID<-IDs.to.analyse

  Zones.comp<-matrix(NA,ncol=3,nrow = length(IDs.to.analyse))
  Zones.comp<-as.data.frame(Zones.comp)  
  
  colnames(Zones.comp)<-c("ID","Flipper.zones","Berg.zones") 
  Zones.comp$ID<-IDs.to.analyse
  
  
for (i in IDs.to.analyse) {
  
  print(paste("Profile ID:",i))
  
  
  #First we need to transform the data in such a way it can be used by FLIPPER
  
  ForFLIPPER<-NULL
  ForFLIPPER$env.parms <- list(TC    = df.temp$T[df.temp$ID==i][1],
                               S     = df.temp$S[df.temp$ID==i][1],
                               P     = df.temp$P[df.temp$ID==i][1],
                               z     = df.temp$z[df.temp$ID==i][1],
                               Dmol  = df.temp$Dmol[df.temp$ID==i][1],
                               Ds    = df.temp$Ds[df.temp$ID==i][1],
                               ID    = df.temp$ID[df.temp$ID==i][1])
  ForFLIPPER$species   <- list("O2")    
  
  ForFLIPPER$input     <- list(x     = (df.temp$x.cor[df.temp$ID==i]*10^-6), #To correct micrometers back to meters
                               C     = df.temp$C[df.temp$ID==i],
                               por   = df.temp$Por[df.temp$ID==i],
                               ID    = df.temp$ID[df.temp$ID==i])
  
  #I want a better solution for this part so I dont need a minimum restriction but it has to do for now
  if (length(ForFLIPPER$input$x)>100) {print(paste("Data trimmed from",length(ForFLIPPER$input$x),"datapoints"))
                                       ForFLIPPER$input<-L.BinV2(ForFLIPPER$input)
                                       print(paste("till",length(ForFLIPPER$input$x),"datapoints"))}
  
  #The reason I use "try" here is that with these big functions they tend to give an error along the way for whatever
  #reason so this way the whole process of my for loop is not stopped 

  testgradient <- try(FLIPPER.func( input    = as.data.frame(ForFLIPPER$input),
                                    tort.dep = 1,
                                    species  = ForFLIPPER$species[[1]],
                                    method   = "all",
                                    env.parms = ForFLIPPER$env.parms,
                                    discrete.parms=list(LBC="conc.down",C.down=ForFLIPPER$input$C[length(ForFLIPPER$input$C)]),
                                    continuous.parms=list(optimal.window.size="interactive"),full.output = TRUE))


  x11(  width  = 100,height = 60)
  current.dev <- dev.cur()
  try(plot.FLIPPER(testgradient))
 
  ForJfick<-ForFLIPPER$input
  ForJfick$x<-ForJfick$x*10^6 #we need to adjust meters to micrometers for my own function
  
  
  ResultLuna <-Jfick(as.data.frame(ForJfick),
                     ID.check = 1,
                     Depth.name = "x",
                     Oxygen = "C",
                     Porosity = "por",
                     Diffusioncoeff = ForFLIPPER$env.parms$Ds,
                     species.input  = as.character(ForFLIPPER$species),
                     Temperature = ForFLIPPER$env.parms$TC,
                     Pressure = ForFLIPPER$env.parms$P,
                     Tort = 1,
                     Salinity = ForFLIPPER$env.parms$S,
                     FLIPPER= testgradient,
                     Plot.save.path =)
  
  mat$LunaGradient[mat$ID==i]<-ResultLuna$Results$Jfick.flux
  mat$FlipperGradient[mat$ID==i]<-testgradient$output$gradient.output$R.int

  if (isTRUE(testgradient$output$discrete=="FAILED")){
    
    mat$FlipperDiscrete     [mat$ID==i] <- testgradient$output$discrete
  }else {mat$FlipperDiscrete     [mat$ID==i] <- testgradient$output$discrete$R.int}
  
  
  if (isTRUE(testgradient$output$continuous=="FAILED")){
    
    mat$FlipperSavGolay     [mat$ID==i] <- testgradient$output$continuous
  }else{ mat$FlipperSavGolay     [mat$ID==i] <- testgradient$output$continuous$R.int}
  
  
  mat$Author              [mat$ID==i] <- as.character(df.temp$Author [df.temp$ID==i][1])
  
  
  #NOW TO INTEGRATE PROFILE FROM BERG IN THIS..... SOMEHOW
  
  Ds.O2          <- testgradient$input$user.input$Ds[1] # [m2 d^-1] , however does this not only work if Ds is the same?
  
  Ds.O2.Profile  <- 10000*Ds.O2/(3600*24) # [cm2 s-1] 
  
  
  
  #=============================================================================
  # Input of Data FOR PROFILE
  #=============================================================================
  
  
  ForPROFILE   <-   as.data.frame(ForFLIPPER$input)
  
  
  colnames(ForPROFILE)<-c("Depth","Signal","Por","ID")
  
  
  
  ForPROFILE$Timepoint       <- 1
  ForPROFILE$Ox.OLW          <- T
  ForPROFILE$Profile.Number  <- 1 
  ForPROFILE$Sensor.Name     <- c("O2")
  
  ForPROFILE                 <- ForPROFILE[!ForPROFILE$Depth<=-0.3,] # Cut off but why?
  ForPROFILE                 <- ForPROFILE[!ForPROFILE$Depth>=0.5,] #REMOVE THIS PART IF CODE RUNS 
  
  #=============================================================================
  # Rate via PROFILE.exe script  
  #=============================================================================
  
  try(
    {
      
      
      # Calculate via PROFILE
      {#Selection of data
        
        
        Ndata      <- nrow(as.data.frame(ForFLIPPER$input))
        
        # Construction of data frame with profile data                                                 
        
        ForPROFILE.raw  <- data.frame(
          X    = ForFLIPPER$input$x*10^2, #to correct meters to cm
          FI   = 1, #This value (porosity) is not technically 1 but we need a way to compare
          
          #Since we cannot access the code of berg but we need to start of with the same assumptions regarding tortuosity
          #berg took Ds(corrected for tortuosity) as D / porosity, by taking por as 1 we keep "Ds" af it were D, now by filling in
          # D (the same one we use in the calculations for FLIPPER etc) we can use the same D here as in FLIPPER
          #only limitation is that we cannot change porosity over depth
          
          
          DB   = rep(0,Ndata),
          ALFA = rep(0,Ndata),
          C    = ForFLIPPER$input$C)
        
        ##Instead of trimming first I will trim afterwards
        
        ForPROFILE<- ForPROFILE.raw[ForPROFILE.raw$X >= 0,]
        sel       <- which(ForPROFILE.raw$X >= 0)
        head(ForPROFILE)
        
        
        #Copied for the most part of your function but adjusted the names to fit my uses   
        interpolate.L <- function(Profile, h){
          
          # Artificially increasing resolution to stepsize h
          Depth             <- seq(from = Profile$X[1],
                                   to   = Profile$X[length(Profile$X)],
                                   by   = h)
          
          Concentration     <- approx(x=Profile$X,
                                      y=Profile$C,
                                      xout=Depth,
                                      rule=2)$y
          
          # Profile.i(nterpolated)
          Profile.i <- data.frame(Depth          = Depth, 
                                  Concentration  = Concentration,
                                  Profile.Number = 1)
          
          return(Profile.i)  
          
        }
        
        ForPROFILE.i  <- interpolate.L(ForPROFILE, h=0.005) # nmol cm-3
        head(ForPROFILE.i)
        
        
        # Construction of data frame with meta data
        ## Don't really get this part, is this the part above the actual data that has to be filled in?
        N <- 13
        meta  <- data.frame(value=rep(0,N),description=rep("",N),stringsAsFactors=FALSE)
        
        meta[1,1] <- min(ForPROFILE.raw$X[sel])
        meta[1,2] <- "Depth at top of calculation domain" 
        meta[2,1] <- max(ForPROFILE.raw$X[sel])
        meta[2,2] <- "Depth at bottom of calculation domain"
        meta[3,1] <- 7
        meta[3,2] <- "Max number of equally spaced zones in interpretation (1 to 12)"
        meta[4,1] <- 3
        meta[4,2] <- "Type of boundary conditions (1:t=C b=C, 2:t=C t=F, 3:b=C b=F 4:t=C b=F 5:t=F b=C)"
        #meta[5,1] <- Ox.profile$Concentration[sel][1]
        meta[5,1] <- 0.0
        meta[5,2] <- "First boundary condition"
        meta[6,1] <- 0.0
        meta[6,2] <- "Second boundary condition"
        meta[7,1] <- Ds.O2.Profile #THIS WE HAVE TO ADJUST FOR TORTUOSITY AND POROSITY AND THEN PUT POROSITY=1
        meta[7,2] <- "Diffusivity in water (D)"    
        meta[8,1] <- 1
        meta[8,2] <- "Expression for sediment diffusivity (Ds) (1: Ds=FI*D, 2: Ds=FI**2*D, 3: Ds=D/(1+3*(1-FI))"
        meta[9,1] <-  ForPROFILE.i$Concentration[1]
        meta[9,2] <- "Concentration in water column (C0)"
        meta[10,1] <- -1.0E+20
        meta[10,2] <- "Minimum for production rate"                                
        meta[11,1] <- 0.0
        meta[11,2] <- "Maximum for production rate"
        meta[12,1] <- 0.001
        meta[12,2] <- "Maximum deviation (in %) when accepting a calculated minimum"
        meta[13,1] <- 0.01
        meta[13,2] <- "Level of significance in the F statistics"        
        
        head(meta)
        
        # Write all data in suitable input file format for PROFILE.exe
        
        out.file    <- paste("ID",i,".inp",sep="")
        result.file <- paste("Profile results ID",i,".txt",sep="")
        write("O2 simulation header",file=out.file)
        write.table(x=meta,file=out.file,append=TRUE,row.names = FALSE,col.names=FALSE,sep = " ")                                                   
        write.table(x=ForPROFILE.raw,file=out.file,append=TRUE,row.names = FALSE)                                                   

        # Execute PROFILE.exe
        Fstart         <- 5
        Finput         <- Fstart 
        PROFILE.output <- system("./Scripts/Dependencies/PROFILE.exe",input=c(out.file,result.file,Fstart,Finput),intern=T)
        PROFILE.output
        x              <- getPROFILEpar(PROFILE.output)
        
        Fstat          <- unlist(x[1],use.names=F)
        
        #Only continue if Fstat is not a character aka error message, if it is, it should write Fstat in the spot of 
        #the flux in table mat
        if(!is.character(Fstat)){        
          counter        <- 0
          
          while(Fstat != Finput){
            Finput         <- Fstat
            PROFILE.output <- system("./Scripts/Dependencies/PROFILE.exe",input=c(out.file,result.file,Fstart,Finput),intern=T)
            
            x              <- getPROFILEpar(PROFILE.output)
            Fstat          <- unlist(x[1], use.names=F)
            counter        <- counter + 1
            print(counter)
            if(counter > 5){break()}
          }
          
          Top.flux       <- unlist(x[2], use.names=F)
          Top.flux
          # Read the output file from PROFILE.exe
          
          file.name               <- result.file
          out                     <- read.table(file.name,header=TRUE)
          colnames(out)[c(1,5,7)] <- c("Depth", "Concentration", "Production")
          out$Production          <- out$Production*3600*24 # [nmol cm-3 s-1] -> [mmol m-3 d-1]
          
          N.out        <- length(out$Production)
          PROFILE.cons <- -sum(1E-2*diff(out$Depth)*0.5*(out$Production[2:N.out]+out$Production[1:(N.out-1)]))
          
          PROFILE.flux <- Top.flux*1E-2*3600*24 # [nmol cm-2 s-1] -> [mmol m-2 d-1]
          
          
          PROFILE <- list(Top.flux=PROFILE.flux, Total.Consumption = PROFILE.cons, Fstat=Fstat, Profile = out[,c(1,5,7)])
        }else {PROFILE.flux <- Fstat }
        
        mat$Discrete.Berg[mat$ID==i]<-PROFILE.flux
        
        
        
      } 
      
      
    })
  
  
  
  if(isTRUE(testgradient$output$discrete=="FAILED")) 
  { Zones.comp$Flipper.zones [Zones.comp$ID==i] <- "FAILED"
  }else 
  {Zones.comp$Flipper.zones [Zones.comp$ID==i] <-  length(testgradient$output$discrete$R.vol$Prod)}  
  
  #Plotting BERGS results
  
  
  #same here, only continue if flux is an actual number  
  if(!is.character(PROFILE.flux)){
    
    
    input.Depth    <-  ForPROFILE.raw$X*(10^-2) #To correct back to meters
    input.C        <-  ForPROFILE.raw$C
    input.prod     <-  data.frame(Prod=out$Production,depth=out$Depth*(10^-2)) #To correct back to meters
    input.prod     <-  input.prod[!duplicated(input.prod$Prod),]
    
    not.real.fit<-data.frame(x=out$Depth*10^-2, #To correct back to meters
                             C=out$Concentration)
    
    
    try( {plot.discrete (input.Depth,
                         input.C,
                         not.real.fit,
                         input.prod,
                         R.int=PROFILE.flux)
      
      
      
      Zones.comp$Berg.zones    [Zones.comp$ID==i] <-  PROFILE$Fstat} )
    
    
  }
######################################################################  
  #PART THAT INCLUDES THE ACTUAL TRUE FLUX AS WELL AS CONSUMPTION
######################################################################  
  
  if(Gen_O2_int) {
  mat$True.Flux[mat$ID==i]<-df.temp$True.Flux[df.temp$ID==i][1]
  
  if(Gen_O2_int){
    

    
    df.R                     <- data.frame           (x     = (df.temp$x.cor[df.temp$ID==i]*10^-6), #To correct micrometers back to meters
                                                     C     = df.temp$C[df.temp$ID==i],
                                                     por   = df.temp$Por[df.temp$ID==i],
                                                     ID    = df.temp$ID[df.temp$ID==i],
                                                     Prod  = df.temp$Production[df.temp$ID==i], #I put in my model CONSUMPTIOn but for plotting
                                                     #we need PRODUCTION hence the negative!
                                                     True.Flux = df.temp$True.Flux[df.temp$ID==i]) 
 
    
    approx_Prod              <- approx(x = df.R$x , y = df.R$Prod , xout =   testgradient$input$continuous.input$x)$y
    approx_flux              <- approx(x = df.R$x , y = df.R$True.Flux , xout =   testgradient$input$continuous.input$x)$y
    
        
    True_flux_data           <- data.frame           (x        = testgradient$input$continuous.input$x,
                                                     C         = testgradient$input$continuous.input$C,
                                                     True.Flux = approx_flux,
                                                     Prod      = approx_Prod)
    
    
  

plot.true <- function(depth, conc, modelfit, R.int=NULL, y.limits = NULL, 
                              prod.limits = NULL, flux.limits = NULL, conc.limits = NULL){
    
    flux         <- modelfit$True.Flux
    prod         <- modelfit$Prod
    model.depth  <- modelfit$x
    
    if(is.null(y.limits))    y.limits <- c(max(depth, na.rm=T)*1.25,min(depth, na.rm = T)*0.75)
    
    if(is.null(prod.limits)) prod.limits <- c(range(c(prod*1.25,prod*0.75)))
    if(prod.limits[1]>0) prod.limits[1] <- 0
    if(prod.limits[2]<0) prod.limits[2] <- 0
    
    if(is.null(conc.limits)) conc.limits <- c(range(c(conc*1.25,conc*0.75)))
    if(is.null(flux.limits)) flux.limits <- c(range(c(flux*1.25,flux*0.75)))
    
    
    par(new=F)
    
    plot(x=conc, y=depth, ylim=y.limits, pch=21, cex=1, bg=gray(level=0.2), xlab="",ylab="", axes=F,
         xlim = conc.limits)
    lines(x=modelfit[,2], y=modelfit[,1], lwd=2, lty=1, col="red")
    
    axis(1, cex.axis=1.2, lwd=1.5, pos=par()$yaxp[1])
    #abline(h=par()$yaxp[1])
    
    mtext(side=1, line=2.5, "Concentration", cex=1.5)
    axis(2, cex.axis=1.2, lwd=1.5, pos=conc.limits[1])
    mtext(side=2, line=1.5, "Depth", cex=1.5)
    
    abline(h=0, lty=1)  
    
    par(new=T)
    
    plot(x=prod, y=model.depth, ylim=y.limits, lwd=2, lty=2, 
         xlab="",ylab="", axes=F, xlim = prod.limits, type="l")
    
    
    
    #abline(v=0)
    
    axis(3, cex.axis=1.2, lwd=1.5)
    # abline(h=par()$yaxp[2])
    mtext(text="Production", side=3, line=2.5, cex=1.5, lwd=1.5)
    
    if(!is.null(R.int)){
      text(y=max(depth), x=mean(prod.limits), adj=c(0,0),
           paste("TRUE Flux \nR.int =",round(R.int,3)), cex=1.5)}
    
    
    
  }


try  (    plot.true(depth = testgradient$input$user.input$x,
          conc  = testgradient$input$user.input$C,
          modelfit = True_flux_data,
          R.int = True_flux_data$True.Flux))
  }
  
  }
  
  
  if(is.character(Filepathplot)){
  name<-paste("/plots ID",i,".png",sep = "")
  #Still make a function where one can specify where the plots should be written to for now its a junk folder to not
  #overwrite past plots and if no folder is given that it should just keep the plot open for easy comparison
  
  path<-paste(Filepathplot ,sep="")
  pathing<-paste(path,name,sep="")
  
  savePlot(filename=pathing,type="png")
  
  }
  ##These two are to remove clutter, otherwise you end up with a dozen of files in your working directory   
  if(PROFILE.files){  unlink(out.file)
                      unlink(result.file)
                      unlink("PROFILE.txt")   }
   
  
  Results <- list(Fluxes=mat,Zone.Comparison=Zones.comp)
  dev.off(current.dev)
  
#what this does is basically that when Filepathplot is left open then filepathplot is null and it gets a character
#if a character or path is given the next condition is also true, if however filepathplot =FALSE then none of the
#if arguments are satisfied since is.null(FALSE)-> FALSE and thus no character is written and as a result
# is.character is also false
  if (is.character(Filepathplot))  save(Results,file=paste(Filepathresults,".Rdata",sep=""))
   
}  
  return(Results)
}
