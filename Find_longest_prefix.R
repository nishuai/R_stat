####find the longest prefix in a set of given strings
strings=c('AATTGGDDAF','AATEFDGD','AATTGCGDSFS', 'AATFGDERDDSS', 'AATDKFJLSDUROEIHDHFLSKFJSKDJFKLSDFHKLSAHFKARJKEQW', 'AAA')

Find_prefix=function(strings){
  
prefix=c()
###the longest possible prefix is 
###limited by the shorted string
max_len=min(nchar(strings))
  for (i in 1:max_len){
    key=substr(strings[1], i, i)
    for (k in 2:length(strings)){
      ###It jumps out the way and outputs the 
      ###prefix once the first mismatch is met
      if (substr(strings[k], i ,i) != key) return(prefix)
      }
    prefix=c(prefix, key)
  }
  return(prefix)
}

