######################################################
#Example of Gen.O2.Profile
######################################################

source("./Gen_O2_Profile.R")


#For example we want 500 data points over a 10 cm (or 0.01 m depth) with a flux at the top of 2 and porosity of 0.8
y<-Gen.O2.Profile(N=500,L=0.01,Flux.top =2 ,por = 0.8,T=10)

#plotting the data
plot(x.cor~C,data=y,ylim=rev(range(x.cor)),main="example profile",xlab="mmol m^-3 O2",ylab="micrometer")

#example calculating multiple profiles in one line


N<-50
L<-0.01
por<-0.6
Flux<-c(1,2,3)


d<-mapply(Gen.O2.Profile,N,L,por,Flux,SIMPLIFY = FALSE)
d
#and plotting these
par(mfrow=c(1,length(d)))
lapply(d,function(x) {plot(x.cor~C,
                           data=x,
                           ylim=rev(range(x.cor)),
                           main="example profile",xlab="mmol m^-3 O2",ylab="micrometer",type="p")
  text(x=max(x$C),pos=4,labels=paste(x$True.Flux,"True flux"))
 } )
par(mfrow=c(1,1))


