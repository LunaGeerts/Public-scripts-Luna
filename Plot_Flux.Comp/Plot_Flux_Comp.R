#########################################################
#Personal use scripts to quickly plot results of a list of compared fluxes
#########################################################



##########################################################
#Needed liberaries
##########################################################
require(reshape2)

##########################################################
#CODE
##########################################################

Error<- function(Result) {
  df.temp <- as.data.frame(Result)
  error.L    <- ((df.temp$LunaGradient-df.temp$True.Flux)^2)/nrow(Result)
  error.F.G  <- ((df.temp$FlipperGradient-df.temp$True.Flux)^2)/nrow(Result)
  error.F.D  <- ((df.temp$FlipperDiscrete-df.temp$True.Flux)^2)/nrow(Result)
  error.F.S  <- ((df.temp$FlipperSavGolay-df.temp$True.Flux)^2)/nrow(Result)
  error.B    <- ((df.temp$Discrete.Berg-df.temp$True.Flux)^2)/nrow(Result)
  error.true <- ((df.temp$True.Flux-df.temp$True.Flux)^2)/nrow(Result)
  
  method <- colnames(Result)[c(2:5,7,8)]
  d<-data.frame(Errors=c(error.L, error.F.G, error.F.D, error.F.S, error.B, error.true),
                Method=rep( method, each=length(df.temp[,1]) ) )
  return(d)
}
plot.dens    <- function(Results, ...){
  
  df.temp  <- as.data.frame(Results)          
  L.dens   <- density( df.temp$LunaGradient)
  FG.dens  <- density( df.temp$FlipperGradient)
  FD.dens  <- density( df.temp$FlipperDiscrete)
  FSG.dens <- density( df.temp$FlipperSavGolay)
  B.dens   <- density( df.temp$Discrete.Berg)
  True.dens<- density( df.temp$True.Flux)
  
  names    <-c("L.dens","FG.dens","FD.dens","FSG.dens","B.dens","True.dens")
  listing  <-list(L.dens,FG.dens,FD.dens,FSG.dens,B.dens,True.dens)
  names(listing)<- names
  x.range<- range ( sapply(listing, function(list){   return( c(min(list$x),max(list$x))) }  ) )
  y.range<- range ( sapply(listing, function(list){   return( c(min(list$y),max(list$y))) }  ) )
  
  plot(True.dens,main= ...,xlim=x.range,ylim = y.range,lwd=2)
  
  #sapply(listing[-(length(listing))],lines,col=2) #you probably can do this with a sapply but not sure how
  
  #for gradient
  lines(listing[["L.dens"]],lty=2,col=2,lwd=1.7)
  lines(listing[["FG.dens"]],lty=2,col=3,lwd=1.7)
  
  #for discrete
  lines(listing[["FD.dens"]],lty=3,col=2,lwd=1.7)
  lines(listing[["B.dens"]],lty=3,col=3,lwd=1.7)
  
  #Golay
  lines(listing[["FSG.dens"]],lty=4,col=4,lwd=1.7)
  
  legend("topright",legend = c("True flux","Luna.Grad","Flip.Grad","Flip.Disc","Berg.disc","Flip.sav"),
         lty = c(1,2,2,3,3,4),lwd=1.7,col = c(1,2,3,2,3,4) )
  
}

L.plot.func<- function(x,title.line= "line graph" ,
                       title.box= "Squared error of each method versus true flux",
                       title.dens="Density plot",
                       text.graph= ("Data generated with:\n
                                    Standard deviation of \n
                                    Depth of x meters \n
                                    x points spread over this depth (but trimmed)")) {
  
  #Creating the data needed later in plots
  
  df<-as.data.frame(x)
  
  error.df<-Error(x)
  methods<-c(levels(error.df$Method))
  
  
  df.temp<- melt(df,id.vars=c("ID","True.Flux"),measure.vars  = methods)
  df.temp.2<-df.temp[!df.temp$variable=="True.Flux",]
  
  df.temp.2$ID<-as.factor(df.temp.2$ID)
  
  
  
  xrange<-range(df.temp.2$True.Flux)
  yrange<-range(df.temp.2$value)
  
  x11(width=210,height=120)
  par(mfrow=c(2,2))
  
  plot(True.Flux~1 ,data=df.temp.2,type="n",xlim=rev(xrange),ylim=rev(yrange),lwd=3,
       xlab="True.Flux",ylab="predicted flux",main=title.line)
  abline(a=0,b=1,lwd=3)
  
  points(value~True.Flux ,data=df.temp.2[df.temp.2$variable=="LunaGradient",],type="p",xlim=rev(xrange),ylim=rev(yrange))
  abline(lm(value~True.Flux ,data=df.temp.2[df.temp.2$variable=="LunaGradient",]))
  
  points(value~True.Flux ,data=df.temp.2[df.temp.2$variable=="FlipperGradient",],type="p",col=2)
  abline(lm(value~True.Flux ,data=df.temp.2[df.temp.2$variable=="FlipperGradient",]),col=2)
  
  points(value~True.Flux ,data=df.temp.2[df.temp.2$variable=="FlipperDiscrete",],type="p",col=3)
  abline(lm(value~True.Flux ,data=df.temp.2[df.temp.2$variable=="FlipperDiscrete",]),col=3,lty=2)
  
  points(value~True.Flux ,data=df.temp.2[df.temp.2$variable=="Discrete.Berg",],type="p",col=4)
  abline(lm(value~True.Flux ,data=df.temp.2[df.temp.2$variable=="Discrete.Berg",]),col=4,lty=2)
  
  points(value~True.Flux ,data=df.temp.2[df.temp.2$variable=="FlipperSavGolay",],type="p",col=6)
  abline(lm(value~True.Flux ,data=df.temp.2[df.temp.2$variable=="FlipperSavGolay",]),col=6,lty=3)
  
  legend("bottomright",legend=c("Luna.Grad","Flip.Grad","Flip.Disc","Berg.disc","Flip.sav"),col=c(1,2,3,4,6),lty=c(1,1,2,2,3))
  
  
  
  
  boxplot(error.df$Errors~error.df$Method,ylab="Flux error",main=title.box)
  
  plot.dens(df,main=title.dens)
  
  
  par(mar = c(0,0,0,0))
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  text(x = 0.5, y = 0.5, paste(text.graph), 
       cex = 1.6, col = "black")
  par(mar = c(5, 4, 4, 2) + 0.1)
  par(mfrow=c(1,1))
  
  
  
  
  
  
  }  


