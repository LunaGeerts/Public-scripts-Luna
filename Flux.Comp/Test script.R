

Test<-read.table("./Example data.txt",header=T)

source("./Scripts/FluxComparison.R")



#Example

Test_adj<- Test[Test$x.cor>=0,] #currently Continuous method can only handly positive values

x<-Flux.Comp(Test_adj,IDs.to.analyse = 18,Gen_O2_int = FALSE)

#Example running with multiple and full output
y<-Flux.Comp(Test_adj,IDs.to.analyse = 19:20,Gen_O2_int = FALSE,Filepathplot = "./Dummymap",Filepathresults = "./Dummymap/ID Testing",
             Full.Output = TRUE)


