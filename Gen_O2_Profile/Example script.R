######################################################
#Example of Gen.O2.Profile
######################################################

source("Gen_O2_Profile.R")

# 500 data points (N), length of 0.01 meter (L) and with a flux at the top of -2 mmol m^-2 d^-1, porosity of 0.8 and temp
# 10 celsius:

y<-Gen.O2.Profile(N=500,L=0.01,Flux.top =2 ,por = 0.8,T=10)

plot(x.cor~C,data=y,ylim=rev(range(x.cor)),main="example profile",xlab="mmol m^-3 O2",ylab="micrometer")

N<-50
L<-0.01
por<-0.6
Flux<-c(1,2,3)

#example calculating multiple profiles in one line

d<-mapply(Gen.O2.Profile,N,L,por,Flux,SIMPLIFY = FALSE)
d
#and plotting these
par(mfrow=c(1,length(d)))
lapply(d,function(x) {plot(x.cor~C,
                           data=x,
                           ylim=rev(range(x.cor)),
                           main="example profile",xlab="mmol m^-3 O2",ylab="micrometer",type="p")
  text(x=max(x$C),pos=4,labels=paste(x$True.Flux,"True flux"))
 })
par(mfrow=c(1,1))


