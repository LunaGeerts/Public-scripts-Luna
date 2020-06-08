################################################################
#
#
#
#
#
#20/3/2020 changed line 133 to compute MEDIAN instead of mean, mean gives slightly OFF values




L.Bin <- function (input,Bin.size.K){

  input.f <- as.data.frame(input)
  
  For.Bin            <- matrix(data=NA, 
                        nrow = (length((input.f$x))/ Bin.size.K ),
                        ncol =  length(names(input.f)))
  
  colnames(For.Bin)  <- names(input.f)
  For.Bin            <- as.data.frame(For.Bin)  
  
for (i in 1:length(For.Bin$x)) {

  
#we need to keep track of where we are in the i loop. we want to form groups of size bin size K
#If we run the first loop then count = 0 then the mean gets taken from
# 0*3+1=1 (=start count) till (0*3+3) = 3 (=end count) and this gets written in spot 1 (first datapoint)
# the second loop becomes (i=2) :
# 1*3+1=4 (=start count) till (1*3+3) = 6 (=end count) and this gets written in spot 2
#third loop...
# 2*3+1=7 (=start count) till (2*3+3) = 9 (=end count) and this gets written in spot 2
#Each loop considers only 3 values and does not use these 3 values further  
  
  
count<-(i-1)
  
start.count   <-  count * Bin.size.K + 1
end.count     <-  count * Bin.size.K + Bin.size.K

mean.depth    <- mean (input.f$x[start.count:end.count])
mean.C        <- mean (input.f$C[start.count:end.count])
mean.por      <- mean (input.f$por[start.count:end.count])
ID            <- input.f$ID[1]   


 
  

For.Bin$x[i]<-mean.depth   
For.Bin$C[i]<-mean.C   
For.Bin$por <-mean.por
For.Bin$ID  <-ID
}
  return(as.list(For.Bin))
  } 




#For L.BinV2 to run one needs an input with collumns x,C,ID and porosity (same as ForFlipper input)
#A grouping percentage needs to be given:
#this is the minimum amount of change needed in two following datapoints to 
#NOT trigger grouping, a standard value of 1 is given, or in other words values need to be at least 1% apart in concentration
#to not get grouped
#third a percentage needs to be given that tells how many datapoints are allowed to be grouped, standard is 3%
#which means that of the total amount of data points only 3% can be grouped into one data point


L.BinV2<- function(input,
                   treshold.percentage=1 ,
                   Grouping.Percentage= 20){

input.f <- as.data.frame(input)

min.C<-min(input.f$C) #maybe filter those values with "negative" concentrations of course you need to then also 
#adjust for the variation in the positive direction otherwise you skew your residual error, no?

max.C<-max(input.f$C)

range.C<-max.C-min.C #is the range of concentrations
if (treshold.percentage>100) { stop("Please enter a valid percentage (below 100)") }
if (Grouping.Percentage>100) { stop("Please enter a valid percentage (below 100)") }

#different methods of grouping
#First we can group depending on the range and a percentage given, eg group all values that dont at least make 10% change
#in total range. Issue with this lies in data with high range and few points it overgroups
percentile<- ((treshold.percentage/100)*range.C) #the amount of change needed in C to trigger grouping

#percentile <- (length(input.f$C)/range.C)


input.f <- as.data.frame(input)

For.Bin            <- matrix(data=NA, 
                             nrow = (length((input.f$x) )),
                             ncol =  length(names(input.f)))

colnames(For.Bin)  <- names(input.f)
For.Bin            <- as.data.frame(For.Bin)  


j<-1

for (i in 1:length(input.f$x)) {
 
#we start off every iteration of i with the counter being 0 once again and the difference of C being 0
#this way the while loop gets initiated
  
  counter     <- 0
  diff.C      <- 0
  
while (diff.C < percentile) {

#The loop keeps going till either the difference in C is not smaller than the given percentile or when
# we grouped a given treshhold percentage of all data points
    
  counter     <- counter +1
  diff.C      <- input.f$C[j]-input.f$C[j+counter]
  
#We build in 2 breaks, one that makes it so if the counter becomes larger than the treshold % of the whole dataset it stops so 
# the biggest grouping that can happen is the treshold percentage amount of datapoints
# can be grouped at max, the second break is needed so the while does not continue if j+counter becomes bigger than
# the possible input length of variable x, if this is not in place, NA's get introduced in diff.C and the code breaks
  
  if (j+counter  >= length    ( input.f$x))                 {break()} 
  if (counter    == as.integer( length( input.f$C )/100*5)) {break()}
}
#when the while loop stops we will group those values together that contributed to either count =20 or 
#there where the difference in C was bigger than the percentile given
  

For.Bin$x[i]     <- median(input.f$x[j: (j+counter) ])  
For.Bin$C[i]     <- mean(input.f$C[j: (j+counter) ])  

#Trying with geometric mean which is given as the product of every data point
#and the square root of the amount of points considered


#Geom.l<-length(input.f$x[j: (j+counter) ])  
#For.Bin$x[i]  <- prod(input.f$x[j: (j+counter) ])^(1/Geom.l)
#For.Bin$C[i]  <- prod(input.f$C[j: (j+counter) ])^(1/Geom.l)
  
  
#This is so we start from the value that follows up the next set of values of our input in the next for loop
j  <- counter+j 



}


For.Bin$por <- input.f$por
For.Bin$ID  <- input.f$ID

For.Bin<-na.omit(For.Bin)

return(For.Bin)
}

#d1<-L.BinV2(ForFLIPPER$input)


#plot(x~C,data=d1)



