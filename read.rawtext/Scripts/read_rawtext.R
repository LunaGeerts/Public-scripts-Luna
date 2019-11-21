################################################################################
##                                           
##  Content - script for read.rawtext
##  What does this script do? - Read, raw text data files and delete the comment section so R can handle these
##  
##  Written by: Luna Geerts
##  Contact: Luna.Geerts@hotmail.com 
##  
##  version 1.0
##
################################################################################
#NOTES
################################################################################

# I made this simple (and likely very buggy) script to read raw txt files such as those found on 
# PANGAEA, but these datafiles contain a comment section till the actual data , seperated by "*/"
# In order to more conveniently read in these datafiles without having to manually alter these data by hand
# thus I made this short script. read.rawtext is a wrapper function for read.delim
# (read.delim is a wrapper function itself of read.table)

# in File you need to specify the path to your data, in the argument skiptill a character that the function which
# has to search for, then all

#############################################################
# libraries needed
#############################################################

#None !

read.rawtext<-function(file,skiptill="*/", ...){

#First read the input text file, line by line and determine on which row we have an "*/"
line.temp   <-  readLines(file)
del.row   <-  which(line.temp=="*/")


#read the data frame skipping till the row containing "*/"
the.data  <- read.delim(file,skip=del.row)

return(the.data)

}


