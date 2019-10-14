

Test<-read.table("./Example data.txt",header=T)

source("./Scripts/FluxComparison.R")



#Example

x<-Flux.Comp(Test,IDs.to.analyse = 18)

#Example running without Temperature input and with file locations specified. for series of IDs

Test2<-Test[,-3]
y<-Flux.Comp(Test2,IDs.to.analyse = 18:20,Filepathplot = "./Dummymap",Filepathresults = "./Dummymap/ID Testing")

args(Flux.Comp)
