# function to get parameters from screen output of PROFILE
getPROFILEpar <- function(PROFILE.output){

x     <- strsplit(PROFILE.output,"The F statistics suggest") # splits all lines on "The "
Fstat <- PROFILE.output[lapply(x, length)==2]                # Takes out the lines that had the string in the line
##Input Luna for error handling

if (identical(Fstat,character(0))) {
       Fstat    <- PROFILE.output[4]
       Top.flux <- PROFILE.output[4]
       
}else {Fstat <- unlist(strsplit(Fstat[[length(Fstat)]]," "))                    # Splits the line 
       Fstat <- as.numeric(Fstat[is.na(as.numeric(Fstat))==F])      # Takes out only number in line

x        <- strsplit(PROFILE.output,"Depth integration of production")
Top.flux <- PROFILE.output[lapply(x, length)==2]
Top.flux <- unlist(strsplit(Top.flux[[1]]," "))
Top.flux <- as.numeric(Top.flux[is.na(as.numeric(Top.flux))==F])}

return(list(F.stat=Fstat,Top.flux=Top.flux))
}
