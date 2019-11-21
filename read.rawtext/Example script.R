#####################################
#Example script of read.rawtext
#####################################



source("./Scripts/read_rawtext.R")

df.temp<-read.rawtext("./Braeckman_2014.tab" )

#simple as that
head(df.temp)
