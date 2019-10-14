###############################################
##                                           ##
##  Functions to import and handle MPdata    ##
##  Written by: Laurine Burdorf              ##
##  laurine.burdorf@nioz.nl                  ##
##                                           ##
###############################################

# Linux graph problem (Ubuntu machine)
Linux.graphing <- function(){
  require(grDevices)
  X11.options(type="nbcairo") 
}

###############################################
# Functions to import data                    #
###############################################

# Function getMPdata
# Reads Excel file or txt file from Unisense or Pyroscience
# Puts extra column in with Profile Number

# Output is a list with x levels, where x is number of sensors
getMPdata <- function(filename, 
                      system       = c("Unisense"))
{
  if (system == "Unisense") {
  # Load file  
  Excel              <- loadWorkbook(filename)
  Profiles           <- readWorksheet(Excel, sheet = "Profiles", startRow = 2)
  
  ## Loop to assign a profile number to each row
  
  Profiles$Profile.Number <- 1
  k <- 1
  for (i in 1: length( Profiles[,1]) )
  {
    if ( is.na( Profiles[i,1] ) == T) k <- k + 1 # if na is encountered in first column, k 1 up
    Profiles$Profile.Number[i] <- k              # use k-value for the profile number
  }
  
  Profiles           <- Profiles[!is.na(Profiles[,1]),]
  Unisense.names     <- c("Time", "Depth", "Concentration", "Signal")  
  number.of.sensors  <- (ncol(Profiles)-1)/4
  colnames(Profiles) <- c(rep(Unisense.names,number.of.sensors),"Profile.Number")
  
  Profiles.l       <- list()
  for (i in 1:number.of.sensors){
    Profiles.l[[i]] <- Profiles[,c((i*4-3):(i*4),ncol(Profiles))]
  }

  return(Profiles.l)
  }
  
  if (system == "Pyroscience") {
    Profiles <- read.delim2(filename, header=F)
    
    breaks   <- as.numeric(rownames(Profiles[Profiles$V2=="(?m)",]))
    if(nrow(Profiles) %in% breaks){breaks <- breaks[!breaks==nrow(Profiles)]}
    begin    <- c(breaks+1)
    end      <- c(breaks-5)
    rm.zero  <- c(end[2:length(end)],nrow(Profiles))-begin
    if(5 %in% rm.zero){begin <- begin[!rm.zero==5]}

    
    Profiles$Profile.number <- 1
    
    fill.vector <- c(begin, nrow(Profiles))
    start <- fill.vector[1]
    k <- 1
    for(i in 2:length(fill.vector)){
      Profiles$Profile.number[start:fill.vector[i]] <- k
      start <- fill.vector[i]
      k     <- k +1
    }
    
    for (i in length(end):1){
      Profiles <- Profiles[-c(end[i]:(end[i]+5)),]
    }
    
    df <- Profiles[,c(2,5,6)]
    
    Profiles[,c(2,5,6)] <- apply(apply(df, 2, gsub, patt=",", replace="."), 2, as.numeric)
    
    
    Pyroscience.names <- c("Time","Depth","Concentration",
                           "Signal", "Profile.Number")

    
    Profiles <- Profiles[!is.na(Profiles[,2])==T,]
    
    
    Profiles.l <- list()
    Profiles.l[[1]] <- Profiles[,c(1,2,3,5,ncol(Profiles))]
    colnames(Profiles.l[[1]]) <- Pyroscience.names 
    Profiles.l[[2]] <- Profiles[,c(1,2,4,6,ncol(Profiles))]
    colnames(Profiles.l[[2]]) <- Pyroscience.names 
    
   
    
    return(Profiles.l)
  }
}


getMPcaldata <- function(filename){
  Excel                  <- loadWorkbook(filename)
  cal.data.int           <- readWorksheet(Excel, sheet = "Calibration Data")[,c(7,8,11,12)]
  colnames(cal.data.int) <- c("Concentration","mV","Concentration.2","mV.2")
  cal.data               <- list(H2S=data.frame(Concentration=cal.data.int[,1],
                                               mV=cal.data.int[,2]),
                                 pH=data.frame(Concentration=cal.data.int[is.na(cal.data.int[,3])==F,3],
                                              mV=cal.data.int[is.na(cal.data.int[,3])==F,4]))
  return(cal.data)
  
}

###############################################
# Functions to plot data                      #
###############################################
# Function MP.plot
# Plots profiles of a list (new window per sensor) or of a dataframe
# Normally will read Signal for x-axis, can be changed
MP.plot.1 <- function(data, xax="Signal",xmin=NULL,xmax=NULL,ymin=NULL,ymax=NULL){
  if(is.null(xmin)==T) xmin <- min(data[,xax])
  if(is.null(xmax)==T) xmax <- max(data[,xax])
  if(is.null(ymax)==T) ymax <- max(data$Depth)
  if(is.null(ymin)==T) ymin <- min(data$Depth)
  plot(x=data[,xax], y=data$Depth, 
       ylim=c(ymax,ymin), 
       xlim=c(xmin,xmax), type='p',
       ylab=c("Depth (um)"), xlab=xax, 
       main=paste("Profile Number",data$Profile.Number[1],sep=" "))
  lines(x=data[,xax], y=data$Depth)
  abline(h=0)
}

MP.plot <- function(all.profiles, xax="Signal"){
  if (class(all.profiles)=="list"){
    number.of.sensors <- length(all.profiles)
    for (i in 1:number.of.sensors){
      Profile.Numbers   <- levels(as.factor(all.profiles[[i]][,"Profile.Number"]))
      sensor.data       <- all.profiles[[i]]
      xmin              <- min(sensor.data[,xax])
      xmax              <- max(sensor.data[,xax])
      ymin              <- min(sensor.data$Depth)
      ymax              <- max(sensor.data$Depth)
      maxplot.rows      <- 4
      plotrows          <- c(ceiling(length(Profile.Numbers)/4))
      windows()
      par(mfrow=c(plotrows,4))
      for (j in 1:length(Profile.Numbers)){
        data <- sensor.data[sensor.data$Profile.Number == Profile.Numbers[j],]
        MP.plot.1(data, xax=xax, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)  }
      
    }
  }
  if (class(all.profiles)=="data.frame"){
    if(is.null(all.profiles$Sensor.Name)==F){
      number.of.sensors <- length(levels(as.factor(all.profiles$Sensor.Name)))
      sensor.names      <- levels(as.factor(all.profiles$Sensor.Name))}
    if(is.null(all.profiles$Sensor.Name)==T){number.of.sensors <- 1}
    for (i in 1:number.of.sensors){
      if(is.null(all.profiles$Sensor.Name)==F) sensor.data <- subset(all.profiles, Sensor.Name==sensor.names[i])
      if(is.null(all.profiles$Sensor.Name)==T) sensor.data <- all.profiles
      
      windows()
      Profile.Numbers   <- levels(as.factor(sensor.data$Profile.Number))
      One               <- length(Profile.Numbers)==1
      xmin              <- min(sensor.data[,xax])
      xmax              <- max(sensor.data[,xax])
      ymin              <- min(sensor.data$Depth)
      ymax              <- max(sensor.data$Depth)
      maxplot.rows      <- 4
      plotrows          <- c(ceiling(length(Profile.Numbers)/4))
      par(mfrow=c(plotrows,4))
      if (One == T){par(mfrow=c(1,1))}
      for (j in 1:length(Profile.Numbers)){
        data <- sensor.data[sensor.data$Profile.Number == Profile.Numbers[j],]
        MP.plot.1(data, xax=xax, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
      }
      
    }   
  }
}

identifyMark <- function(x, y=NULL, n=length(x), col="red", pch=19, ...)
{
  xy <- xy.coords(x, y); x <- xy$x; y <- xy$y
  sel <- rep(FALSE, length(x)); res <- integer(0)
  while(sum(sel) < n) {
    ans <- identify(x[!sel], y[!sel], n=1, plot=TRUE, col=col,...)
    if(!length(ans)) break
    ans <- which(!sel)[ans]
    points(x[ans], y[ans], col = col, pch=pch)
    sel[ans] <- TRUE
    res <- c(res, ans)
  }
  res
}

plot.and.ask.for.OK <- function(x,y)
{
  satisfied <- FALSE
  
  xdot <- vector(length=2)
  xdot[1] <- min(x)+0.8*(max(x)-min(x))
  xdot[2] <- xdot[1]
  
  ydot <- vector(length=2)
  ydot[1] <- min(y)+ 0.8*(max(y)-min(y))
  ydot[2] <- min(y)+ 0.9*(max(y)-min(y))
  
  points(x=xdot[1], y=ydot[1],pch=21,col="blue")
  points(x=xdot[2], y=ydot[2],pch=21,col="red")
  text(x=xdot[1],y=ydot[1],labels=c("OK"),col="blue", pos=4)
  text(x=xdot[2],y=ydot[2],labels=c("REDO"),col="red", pos=4)
  mark <- identify(x=xdot, y=ydot, n = 1, plot = FALSE)
  if (mark == 1) points(x=xdot[mark], y=ydot[mark], pch=19, col="blue")
  if (mark == 2) points(x=xdot[mark], y=ydot[mark], pch=19, col="red")
  if (mark == 1) satisfied <- TRUE 
  
  return(satisfied)
}

click.SWI <- function(Profile, lim = NA, label = "Profile", suggest = NA)
{
  # Initialisation  
  
  if(is.na(lim[1]) == T) lim   <- sort(range(Profile$Depth),T)
  
  Depth <- Profile$Depth[Profile$Depth <= max(lim)]
  Value <- Profile$Signal[Profile$Depth <= max(lim)]
    
  # Loop to determine mean Signal (mV) value in overlying water
  satisfied <- FALSE
  while (!satisfied){
    
    plot(Value, Depth,
         ylab=(""), xlab="", ylim = lim, axes=FALSE, 
         frame.plot=FALSE, type ="p", col="red", pch=15)
    if(is.na(suggest[1]) == F) points(Value[suggest], Depth[suggest], col="black", pch=15)
    axis(pos=par()$xaxp[1], side=2,lwd=1)
    abline(h=0, lwd=1, lty=1)
    axis(pos=par()$yaxp[2], side=3,lwd=1)
    
    
    xdot <- ( max(Value) + min(Value) ) / 2
    ydot <- max(Depth)
    
    
    text(xdot,0.8*ydot,labels="Indicate the SWI")
    text(xdot,ydot,labels=label)
    mark <- identifyMark(x=Value, y=Depth, labels = c("SWI"),n = 1,col="blue",pch=19)
    
    satisfied <- plot.and.ask.for.OK(Value,Depth)
  }
  
  return(mark)
}

click.OPD <- function(Profile, lim = NA, label = "Profile", suggest = NA)
{
  # Initialisation  
  
  if(is.na(lim[1]) == T) lim   <- sort(range(Profile$Depth),T)
  
  Depth <- Profile$Depth[Profile$Depth <= max(lim)]
  Value <- Profile$Signal[Profile$Depth <= max(lim)]
  
  # Loop to determine mean Signal (mV) value in overlying water
  satisfied <- FALSE
  while (!satisfied){
    
    plot(Value, Depth,
         ylab=(""), xlab="", ylim = lim, axes=FALSE, 
         frame.plot=FALSE, type ="p", col="red", pch=15)
    if(is.na(suggest[1]) == F) points(Value[suggest], Depth[suggest], col="black", pch=15)
    axis(pos=par()$xaxp[1], side=2,lwd=1)
    abline(h=0, lwd=1, lty=1)
    axis(pos=par()$yaxp[2], side=3,lwd=1)
    
    
    xdot <- ( max(Value) + min(Value) ) / 2
    ydot <- max(Depth)
    
    
    text(xdot,0.8*ydot,labels="Indicate the OPD")
    text(xdot,ydot,labels=label)
    mark <- identifyMark(x=Value, y=Depth, labels = c("OPD"),n = 1,col="blue",pch=19)
    
    satisfied <- plot.and.ask.for.OK(Value,Depth)
  }
  
  return(Profile$Depth[mark])
}

click.mark <- function(Profile, lim = NA, label = "Profile", suggest = NA)
{
  # Initialisation  
  
  if(is.na(lim[1]) == T) lim   <- sort(range(Profile$Depth),T)
  
  Depth <- Profile$Depth[Profile$Depth <= max(lim)]
  Value <- Profile$Signal[Profile$Depth <= max(lim)]
  
  # Loop to determine mean Signal (mV) value in overlying water
  satisfied <- FALSE
  while (!satisfied){
    
    plot(Value, Depth,
         ylab=(""), xlab="", ylim = lim, axes=FALSE, 
         frame.plot=FALSE, type ="p", col="red", pch=15)
    if(is.na(suggest[1]) == F) points(Value[suggest], Depth[suggest], col="black", pch=15)
    axis(pos=par()$xaxp[1], side=2,lwd=1)
    abline(h=0, lwd=1, lty=1)
    axis(pos=par()$yaxp[2], side=3,lwd=1)
    
    
    xdot <- ( max(Value) + min(Value) ) / 2
    ydot <- max(Depth)
    
    
    text(xdot,0.8*ydot,labels="Indicate mark")
    text(xdot,ydot,labels=label)
    mark <- identifyMark(x=Value, y=Depth, labels = c("mark"),n = 1,col="blue",pch=19)
    
    satisfied <- plot.and.ask.for.OK(Value,Depth)
  }
  
  return(mark)
}

###############################################
# Functions to set zero                       #
###############################################

set.0.bulk <- function(data){
  windows()
  if (class(data) == "list"){
    number.of.sensors <- length(data)
    for (j in 1:number.of.sensors){
    
    Profile.Numbers <- levels(as.factor(data[[j]]$Profile.Number))
    
    for (i in 1:length(Profile.Numbers)){
      Profile        <- subset(data[[j]], Profile.Number == Profile.Numbers[i])
      mark           <- click.SWI(Profile)
      depth.offset   <- Profile$Depth[mark]
      Profile$Depth  <- Profile$Depth-Profile$Depth[mark]
      
      data[[j]][data[[j]]$Profile.Number==Profile.Numbers[i],]  <-Profile
    }
    }
  }
  if (class(data) == "data.frame"){
    
    if ("Sensor.Name" %in% colnames(data)==T){
      number.of.sensors <- length(levels(as.factor(data$Sensor.Name)))
      sensor.names      <-  levels(as.factor(data$Sensor.Name))
      for (j in 1:number.of.sensors){
        sensor.data      <- subset(data, Sensor.Name==sensor.names[j]) 
        Profile.Numbers  <- levels(as.factor(data$Profile.Number))       
        for (i in 1:length(Profile.Numbers)){
          Profile        <- subset(sensor.data, Profile.Number == Profile.Numbers[i])
          mark           <- click.SWI(Profile)
          depth.offset   <- Profile$Depth[mark]
          Profile$Depth  <- Profile$Depth-Profile$Depth[mark]
          sensor.data[sensor.data$Profile.Number==Profile.Numbers[i],] <- Profile
          } 
        data[data$Sensor.Name==sensor.names[j],] <- sensor.data
      }
    }
    
    if ("Sensor.Name" %in% colnames(data)==F){
    for (i in 1:length(Profile.Numbers)){
      Profile        <- subset(data, Profile.Number == Profile.Numbers[i])
      mark           <- click.SWI(Profile)
      depth.offset   <- Profile$Depth[mark]
      Profile$Depth  <- Profile$Depth-Profile$Depth[mark]
      
      data[data$Profile.Number==Profile.Numbers[i],]  <-Profile
    } 
    }
  }
  return(data)
  }

calibrate.O2 <- function(O2.object, O2_ow, label = "Profile", OLW.100=T, mV_OW=NULL)
{ windows() 
  # Initialisation
  Depth <- O2.object$Depth
  Value <- O2.object$Signal
  lim <- sort(range(Depth),T)
  if (OLW.100==T){
  # Loop to determine mean Signal (mV) value in overlying water
  satisfied <- FALSE
  while (!satisfied){
    
    plot(Value, Depth,
         ylab=(""), xlab="", ylim = lim, axes=FALSE, 
         frame.plot=FALSE, type ="p", col="red", pch=15, lwd=2)
    axis(pos=par()$xaxp[1], side=2,lwd=1)
    abline(h=0, lwd=1, lty=1)
    axis(pos=par()$yaxp[2], side=3,lwd=1)
    mtext(side=3, cex=0.8, text=expression(O[2]~(mu*mol~L^-1), sep=" "), line=1.5)
    
    xdot <- 0.5*max(Value)
    ydot <- max(Depth)
    
    
    text(xdot,0.8*ydot,labels="Indicate the two markers for the mean mV in overlying water ")
    text(xdot,ydot,labels=label)
    mark <- identifyMark(x=Value, y=Depth, labels = c("OW1"),n = 1,col="blue",pch=19)
    max.start <- mark
    mark <- identifyMark(x=Value, y=Depth, labels = c("OW2"),n = 1,col="blue",pch=19)
    max.stop <- mark
    
    n <- length(max.start:max.stop)
    mV_OW <- mean(Value[max.start:max.stop])
    lines(rep(mV_OW,n),Depth[max.start:max.stop],lwd="2",col="blue")
    
    
    satisfied <- plot.and.ask.for.OK(Value,Depth)
  }
  }
  
  # Loop to determine mean Signal (mV) value in anoxic sediment
  satisfied <- FALSE
  while (!satisfied){
    
    plot(Value, Depth,
         ylab=(""), xlab="", ylim = lim, axes=FALSE, 
         frame.plot=FALSE, type ="p", col="red", pch=15, lwd=2)
    axis(pos=par()$xaxp[1], side=2,lwd=1)
    abline(h=0, lwd=1, lty=1)
    axis(pos=par()$yaxp[2], side=3,lwd=1)
    mtext(side=3, cex=0.8, text=expression(O[2]~(mu*mol~L^-1), sep=" "), line=1.5)
    xdot <- 0.5*max(Value)
    ydot <- max(Depth)
    
    
    text(xdot,0.8*ydot,labels="Indicate the two markers for the mean mV in the anoxic zone")
    text(xdot,ydot,labels=label)
    
    mark <- identifyMark(x=Value, y=Depth, labels = c("AZ1"),n = 1,col="blue",pch=19,atpen=TRUE)
    max.start <- mark
    mark <- identifyMark(x=Value, y=Depth, labels = c("AZ2"),n = 1,col="blue",pch=19,atpen=TRUE)
    max.stop <- mark
    
    n <- length(max.start:max.stop)
    mV_AZ <- mean(Value[max.start:max.stop])
    lines(rep(mV_AZ,n),Depth[max.start:max.stop],lwd="2",col="blue")
    
    text(xdot,0.8*ydot,labels=label)
    satisfied <- plot.and.ask.for.OK(Value,Depth)
  }
  
  Conc <- (O2_ow)/(mV_OW-mV_AZ)*(O2.object$Signal-mV_AZ)
  
  
  
  return(list(C = Conc, mV_OW = mV_OW, mV_AZ = mV_AZ))
}

calibrate.O2.bulk <- function(data, O2_ow, OLW.100=T, mV_OW=NULL){
  
  # Empty and prepare concentration column
  data$Concentration <- NaN
  mV_AZ              <- NULL
  
  if (OLW.100 == T){
    Profile.Numbers <- levels(as.factor(data$Profile.Number))
    
    for (i in 1:length(Profile.Numbers)){
      Profile               <- subset(data, Profile.Number == Profile.Numbers[i])
      O2.calibration        <- calibrate.O2(O2.object=Profile, O2_ow = O2_ow)
      Profile$Concentration <- O2.calibration$C
      mV_OW                 <- c(mV_OW,O2.calibration$mV_OW)
      mV_AZ                 <- c(mV_AZ,O2.calibration$mV_AZ)
      data[data$Profile.Number==Profile.Numbers[i],]  <- Profile
      
    }    
  } 
 
  if (OLW.100 == F){
    Profile.Numbers <- levels(as.factor(data$Profile.Number))
    
    for (i in 1:length(Profile.Numbers)){
      Profile               <- subset(data, Profile.Number == Profile.Numbers[i])
      O2.calibration        <- calibrate.O2(Profile, O2_ow = O2_ow, OLW.100=F, mV_OW=mV_OW)
      Profile$Concentration <- O2.calibration$C
      
      data[data$Profile.Number==Profile.Numbers[i],]  <-Profile
    }    } 
    
  data[,3] <- as.numeric(data[,3])
  return(list(data=data, mV_OW = mV_OW, mV_AZ = mV_AZ))
}

interpolate <- function(Profile, h){
  
  # Artificially increasing resolution to stepsize h
  Depth             <- seq(from = Profile$Depth[1],
                           to   = Profile$Depth[length(Profile$Depth)],
                           by   = h)
  
  Concentration     <- approx(x=Profile$Depth,
                              y=Profile$Concentration,
                              xout=Depth,
                              rule=2)$y
  
  # Profile.i(nterpolated)
  Profile.i <- data.frame(Depth          = Depth, 
                          Concentration  = Concentration,
                          Profile.Number = Profile$Profile.Number[1])
  
  return(Profile.i)  
  
}

add.0 <- function(Profile.i, N.smooth){
  test.Concentration <- c(Profile.i$Concentration, rep(0,N.smooth) )
  test.depth         <- c(Profile.i$Depth, 
                          Profile.i$Depth[length(Profile.i$Depth)] +
                          seq(0,(N.smooth-1)*h.O2,length.out=N.smooth))
  Profile.i  <- data.frame(Depth=test.depth, Concentration=test.Concentration)
  return(Profile.i)
}

burst.list <- function(profile.list, sensor.names=NULL){
  
  if(class(profile.list) != 'list'){return(c("input is not a list"))}
  
  number.of.sensors      <- length(profile.list)
  if(is.null(sensor.names)==T){sensor.names <- seq(1,number.of.sensors)}
  
  profile.dataframe <- data.frame()
  for (i in 1:number.of.sensors){
    df                <- profile.list[[i]]
    df$Sensor.Name    <- sensor.names[i]
    profile.dataframe <- rbind(profile.dataframe,df)
    df                <- NULL
  }
  return(profile.dataframe)
}



