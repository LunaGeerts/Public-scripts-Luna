################################################################################
##                                           
##  Content - Luna script Zerodepthcor
##  What does this script do? - Adjust 0 depth of oxygen profile to the
##                              one selected on the graph of raw data
##  Written by: Luna Geerts
##  Contact: Luna.Geerts@hotmail.com 
##  
##  version 1.0
##
################################################################################
##  VERSION HISTORY               
##----------------------------------------------------------------------------                    
##  Created script 9/6/2019
##  Edited 9/8/2019


## KNOWN BUGS
##----------------------------------------------------------------------------
## Script does not deal well with ID numberings with gaps in them: eg 
## a series of 1 till 5 with no ID 2

################################################################################
# Set working directory and load packages
################################################################################

#setwd("C:/Users/install/Dropbox/Thesis/Raw data")

#install.packages("DescTools")
library(DescTools)





################################################################################
#Input needed for ZeroDepthcor / manual
################################################################################

#The depth units have to be in micrometer whereas a negative depth means 
#there was sampling above the sediment in the water column, a positive depth value means
#the datapoint was taken inside the sediment

#Dataframes should contain a collumn named ID for the script to work, if no such collumn exist
#Make one manually or use the function "assign.ID" from another script I wrote

#The data must be in a dataframe and 3 variables have to be included 
# 1--Name of the dataset to be considered
# 2--The second argument contains the collumn name for depth
# 3--The third argument is the name of the collumn that contains the oxygen data


################################################################################
#ZeroDepthcor function code
################################################################################

ZeroDepthcor<-function(dataframename,depthname,oxygen){
  
  nulldepth        <- NULL
  Coordinatelist   <- NULL
  dataframename$depth.cor     <-NULL
  attach(dataframename)
  
  
  dataframename$depth.cor<-NULL
  for (i in sort(unique(ID)) ) {
   
     #The title of the plot contains the ID number
    
    plot(depthname[ID==i] ~ oxygen[ID==i], ylim = rev(range(depthname[ID==i])), data = dataframename, main=i)
    
   
    cor                      <- locator(n=1)$y
    nulldepth                <- Closest(depthname[ID==i],cor[1])
    dataframename$x.cor[ID==i] <- depthname[ID==i]-nulldepth
    

    
    
    #Now we create the coordinate list for the points selected so that we know afterwards what depth we considered 0
    Coordinatelist$ID<- sort(unique(ID))
   
    
    #These are the depths I clicked on in the original plots to be stored for potential later use
    Coordinatelist$zerodepth[Coordinatelist$ID==i] <- nulldepth 
    Coordinatelist            <-as.data.frame(Coordinatelist)
  }
  
  
    detach(dataframename)
    output                    <-list(dataframename,Coordinatelist)
    return( output )
  
  
}



################################################################################
#Example of script running
################################################################################
#The loaded file must include a collumn name that makes each profile distinct from the other


#load("./testdataZeroDepthcor.rda")
#First we convert the depth to micrometers
#head(testdataZeroDepthcor)
#testdataZeroDepthcor$Depth<-testdataZeroDepthcor$Depth*1000
#df10<-ZeroDepthcor(testdataZeroDepthcor,Depth,C)

#save(df10,file="TestdataJfick.rda")

#We now have an object that contains the adjusted zero depths as well as the coordinates of the point selected
#if one wishes to only extract the adjusted dataframe then add "[[1]]" to the end of the
#function so only the first dataframe is stored

