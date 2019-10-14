################################################################################
##                                           
##  Content - Luna script assign ID
##  What does this script do? - Assign unique ID's to a new dataset
##                              to facillitate merging with another dataset
##  Written by: Luna Geerts
##  Contact: Luna.Geerts@hotmail.com 
##  
##  version 1.1
##
################################################################################
##  VERSION HISTORY               
##----------------------------------------------------------------------------                    
##  Created script 6/9/2019
##  11/9/2019 - changed that IDs now are passed on as numeric values rather than 
##  factors which caused issues in scripts such as zerodepthcor that relied on the ID
##  being numeric

##  KNOWN BUGS:
##----------------------------------------------------------------------------
##  When metadata with more than one row is loaded in, the following warning pops up

##----------------------------------------------------------------------------
##Warning message:
##In if (Meta[1] == 0) { :
##  the condition has length > 1 and only the first element will be used
##----------------------------------------------------------------------------

##The script still works fine and as intended regardless.

################################################################################
# Set working directory and load packages
################################################################################
#setwd("C:/Users/install/Dropbox/Thesis/Raw data")
#No packages needed! :)




################################################################################
#Input needed for assign.ID / manual
################################################################################

#The data must be in a dataframe and 3 variables have to be included 
# 1--Name of the collumn used for assigning unique ID's, This collumn will be  
#    used to as basis for unique ID assignment
# 2--The second argument contains the dataframe name
# 3--The third argument contains the metadata (if present), if no metadata 
#    is present enter a 0
#    The metadata used needs to have a collumn named ID and it needs to be 
#    preferably in the first collumn, also no ID 0 should exist, tl:dr as long
#    as no value "0" is present on the first collumn first row things should work fine



################################################################################
#assign.ID function code
################################################################################

assign.ID<-function(Identifier,dataframename,Meta) {
  
  #First I let the function check if there is already metadata present, 
  #if there is, it will start numbering The new ID's +1 the 
  #highest ID number of the previous metadata set
  #if no data is present then it will start counting from 1

  if (Meta[1]==0) {
    d=0
  }else{  
    d=max(Meta$ID)
  }
  
  #I attach the names of the dataframe so you dont have to write down the whole path of 
  #the identifier, only writing the name of the identifier is sufficient
  

  
  
  #Now I will assign an ID collumn based on the identifier name 
  attach(dataframename)
  dataframename$ID<-as.factor(Identifier)
  detach(dataframename)
  #If no metadata was available then d=0 and a sequence from d+1 till the length of all the levels is made
  #If metadata is present; eg 5 ID's then a sequence from 5+1(=6) is created till the length of the new to be introduced
  #dataset + the previous ID numbering present
  #if the old dataset had 5 id's and the new one has 3 unique Id's the 3 to be introduced ID's need to have a numbering from
  #6 till 8
  
  levels(dataframename$ID)<-seq(d+1,length(levels(dataframename$ID))+d)
  
  #This way we transform our factored numbers into actual numeric numbers
  attach(dataframename)
  dataframename$ID<-as.numeric(levels(ID))[ID]
  detach(dataframename)
  return(dataframename)
}


################################################################################
#Example of script running
################################################################################
#The loaded file must include a collumn name that makes each profile distinct from the other

#load("./dfO2.RData")
#If a file containing previous data already exists it must be loaded in as well, this is to ensure when merging
#the datasets later no two profiles got the same ID
#Example of how the script can generate a new ID list if no metadata is present

#assign.ID(Core.Name,plot.dataframe.O2,0)

#Example of how script can generate a new ID list when metadata was already present

#load("./metadata.Rda")
#assign.ID(Core.Name,plot.dataframe.O2,metadata2)



