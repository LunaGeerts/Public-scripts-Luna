# function to get parameters from screen output of PROFILE
getPROFILEpar <- function(PROFILE.output){

x     <- strsplit(PROFILE.output,"The F statistics suggest") # splits all lines on "The "
Fstat_org <- PROFILE.output[lapply(x, length)==2]                # Takes out the lines that had the string in the line

Zone_sug <- any(grepl("The F statistics have no suggestion. Choose the number of zones for further",x))
if(Zone_sug) {
              failedF<- "Failed"
   }else{     failedF<- "OK"
                
              }



if (identical(Fstat_org,character(0))) { ##Input Luna for error handling
       Fstat    <- PROFILE.output[4]
       Top.flux <- PROFILE.output[4]
       
}else {Fstat <- unlist(strsplit(Fstat_org[[length(Fstat_org)]]," "))                    # Splits the line 
       Fstat <- as.numeric(Fstat[is.na(as.numeric(Fstat))==F])      # Takes out only number in line

       Fstat3 <- unlist(strsplit(Fstat_org[[1]]," "))                    # Splits the line 
       Fstat3 <- as.numeric(Fstat3[is.na(as.numeric(Fstat3))==F])     
       
       
x        <- strsplit(PROFILE.output,"Depth integration of production")
Top.flux <- PROFILE.output[lapply(x, length)==2]
Top.flux <- unlist(strsplit(Top.flux[[1]]," "))
Top.flux <- as.numeric(Top.flux[is.na(as.numeric(Top.flux))==F])}

return(list(F.stat=Fstat,Top.flux=Top.flux,Fstat3=Fstat3,failedF))
}
