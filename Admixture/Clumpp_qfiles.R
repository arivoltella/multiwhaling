#Load pophelper package - if not installed use command: install.packages("pophelper")
library(pophelper)

#Create a list of "Q" files
#Note: $ at the end means that this is end of string. As "." is a special character in regular expressions we need to escape it using "\\."
Q.files <- list.files(pattern = ".Q") 

#Create input for CLUMPP - this will output a new direcory with the name "pop_KX" where "X" is the K value being summarised.
clumppExport(readQ(Q.files), exportpath=getwd())

#Exit R
quit()