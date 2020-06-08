
plot.continuous <- function(depth, conc, modelfit, R.int=NULL, y.limits = NULL, 
                    prod.limits = NULL, flux.limits = NULL, conc.limits = NULL,
                    TRUE_FLUX=FALSE #Luna Added this argument 19/4/2020
                    ){

  flux         <- modelfit$J.tot
  prod         <- modelfit$R
  model.depth  <- modelfit$x
  
  if(is.null(y.limits))    y.limits <- c(max(depth, na.rm=T)*1.25,min(depth, na.rm = T)*0.75)
  
  if(is.null(prod.limits)) prod.limits <- c(range(c(prod*1.25,prod*0.75)))
  if(prod.limits[1]>0) prod.limits[1] <- 0
  if(prod.limits[2]<0) prod.limits[2] <- 0
  
  if(is.null(conc.limits)) conc.limits <- c(range(c(conc*1.25,conc*0.75)))
  if(is.null(flux.limits)) flux.limits <- c(range(c(flux*1.25,flux*0.75)))
  
  
  par(new=F)
  
  
  
  plot(x=prod, y=model.depth, ylim=y.limits, lwd=2, lty=2, 
       xlab="",ylab="", axes=F, xlim = prod.limits, type="n")
  axis(3, cex.axis=1.2, lwd=1.5)
  
  y1<- c(0,0,model.depth,max(model.depth))
  x <- c(0,min(prod),prod,0)
  

  
  polygon(x,y1, border=NA,  col=gray(level=0.95))
   
  lines(x=prod, y=model.depth, lwd=2, lty=2) 
  lines(x=seq(0,0,length.out = length(model.depth)), y=seq(0,max(model.depth),length.out = length(model.depth))
                                                         , lwd=1)
  
  par(new=T)
  
  
  plot(x=conc, y=depth, ylim=y.limits, pch=21, cex=1, bg=gray(level=0.2), xlab="",ylab="", axes=F,
       xlim = conc.limits)
  lines(x=modelfit[,2], y=modelfit[,1], lwd=2, lty=1, col="red")
  
  axis(1, cex.axis=1.2, lwd=1.5, pos=par()$yaxp[1])
  #abline(h=par()$yaxp[1])
  
  mtext(side=1, line=2, bquote(mmol ~O[2]~m^-3), cex=1.0)
  axis(2, cex.axis=1.2, lwd=1.5, pos=conc.limits[1])
  mtext(side=2, line=2, "Depth (meter)", cex=1.0)
  
 # abline(h=0, lty=1)  
  
  
   #abline(v=0)
  
  # abline(h=par()$yaxp[2])
  mtext(bquote(Production ~O[2]~m^-2~day^-1), side=3, line=2, cex=1, lwd=1.5)  
  abline(h=0, lty=1) 
  
  if(!is.null(R.int)){
    text(y=max(depth), x=mean(conc.limits),
         paste("FL_Con \nR.int =",round(R.int,3)), cex=1.5)}
  
  if(TRUE_FLUX){
    R.int <- flux[1]
      text(y=max(depth), x=mean(conc.limits),
          paste("True Flux =",round(R.int,3)), cex=1.5)}
  
  
  
}



plot.discrete <- function(depth, conc, modelfit, prod, R.int=NULL, 
                            conc.limits=NULL, prod.limits=NULL, y.limits=NULL){
  
  zones <- prod$depth
  prod  <- prod$Prod
  zones <- c(zones, max(depth))
  if(is.null(y.limits))    y.limits <- c(max(depth, na.rm=T)*1.25,min(depth, na.rm = T)*0.75)
  if(is.null(prod.limits)) prod.limits <- c(range(c(prod*1.25,prod*0.75)))
  if(prod.limits[1]>0) prod.limits[1] <- 0
  if(prod.limits[2]<0) prod.limits[2] <- 0
  if(is.null(conc.limits)) conc.limits <- c(range(c(conc*1.25,conc*0.75)))
  
  
  par(new=F)
  plot(x=prod, y=zones[1:(length(zones)-1)], type="n", ylim=y.limits, axes=F, ylab="", xlab="",
       xlim = prod.limits)
  rect(xleft=0, ybottom=zones[2:length(zones)], ytop=zones[1:(length(zones)-1)], xright=prod, 
       col=gray(level=0.95))
  
  #abline(v=0)
  
  axis(3, cex.axis=1.2, lwd=1.5)
  # abline(h=par()$yaxp[2])
  mtext(bquote(Production ~O[2]~m^-2~day^-1), side=3, line=2, cex=1, lwd=1.5)  
  
  par(new=T)
  
  plot(x=conc, y=depth, ylim=y.limits, pch=21, cex=1, bg=gray(level=0.2), xlab="",ylab="", axes=F,
       xlim = conc.limits)
  lines(x=modelfit[,2], y=modelfit[,1], lwd=2, lty=1, col="red")
  
  axis(1, cex.axis=1.2, lwd=1.5, pos=par()$yaxp[1])
  #abline(h=par()$yaxp[1])
  
  mtext(side=1, line=2, bquote(mmol ~O[2]~m^-3), cex=1.0)
  axis(2, cex.axis=1.2, lwd=1.5, pos=conc.limits[1])
  mtext(side=2, line=2, "Depth (meter)", cex=1.0)
  
  abline(h=0, lty=1)  
 
   if(!is.null(R.int)){
    text(y=max(depth), x=mean(conc.limits),
         paste("FL_Disc \nR.int =",round(R.int,3)), cex=1.5)}
  
 
}



plot.gradient <- function(depth, conc, modelfit, R.int=NULL, 
                          conc.limits=NULL, y.limits=NULL){
  
  
  if(is.null(y.limits))    y.limits <- c(max(depth, na.rm=T)*1.25,min(depth, na.rm = T)*0.75)
  if(is.null(conc.limits)) conc.limits <- c(range(c(conc*1.25,conc*0.75)))
  
  
  par(new=F)
  
  plot(x=conc, y=depth, ylim=y.limits, pch=21, cex=1, bg=gray(level=0.2), xlab="",ylab="", axes=F,
       xlim = conc.limits)
  
  
  lines(x=modelfit$x*modelfit$slope+modelfit$intercept, y=modelfit$x, lwd=2, lty=1, col="red")
  
  axis(1, cex.axis=1.2, lwd=1.5, pos=par()$yaxp[1])
  #abline(h=par()$yaxp[1])
  
  mtext(side=1, line=2, bquote(mmol ~O[2]~m^-3), cex=1.0)
  axis(2, cex.axis=1.2, lwd=1.5, pos=conc.limits[1])
  mtext(side=2, line=2, "Depth (meter)", cex=1.0)
  
  
  abline(h=0, lty=1)  
  
  if(!is.null(R.int)){
    text(y=max(depth), x=mean(conc.limits),
         paste("FL_Grad \nR.int =",round(R.int,3)), cex=1.5)}
  
  
}



plot.FLIPPER <- function(result){
  
  if (result$method=="continuous"){
    
    par(mfrow=c(1,1))
    depth    <- result$input$continuous.input[,"x"]
    conc     <- result$input$continuous.input[,"C"]
    modelfit <- result$output$overview[,c("x","C","J.tot","R.vol")]
    R.int    <- result$output$R.int
    
    plot.continuous(depth,conc,modelfit,R.int)
  }
  
  if (result$method=="discrete"){
    
    par(mfrow=c(1,1))
    depth     <- result$input[,"x"]
    conc      <- result$input[,"C"]
    modelfit  <- result$output$fit
    prod      <- result$output$R.vol
    R.int     <- result$output$R.int
    
    plot.discrete(depth,conc,modelfit,prod,R.int)
  }
  
  if (result$method=="gradient"){
    
    par(mfrow=c(1,1))
    depth     <- result$input[,"x"]
    conc      <- result$input[,"C"]
    modelfit  <- result$output$fit
    R.int     <- result$output$R.int
    
    plot.gradient(depth,conc,modelfit,R.int)
  }
  
  if (result$method=="all"){
    par(mfrow=c(1,5))
    
    # plot gradient
try({   
    depth     <- result$input$user.input[,"x"]
    conc      <- result$input$user.input[,"C"]
    modelfit  <- result$output$gradient.output$fit
    R.int     <- result$output$gradient.output$R.int
    
    plot.gradient(depth,conc,modelfit,R.int)
   })
    # plot discrete
try({  
    modelfit  <- result$output$discrete$fit
    prod      <- result$output$discrete$R.vol
    R.int     <- result$output$discrete$R.int
    
    plot.discrete(depth,conc,modelfit,prod,R.int)
  })    
    # Plot continuous
try({    modelfit <- result$output$continuous$overview[,c("x","C","J.tot","R.vol")]
    R.int    <- result$output$continuous$R.int
    
    plot.continuous(depth,conc,modelfit,R.int)
    
   
       })
    }
    
 }
