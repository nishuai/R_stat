###find match in a list of sorted integers
###this function outputs the index of the match in 
###sorted vector, otherwise return -1

Find_target=function(sorted_vector, target){
    ####fill the vector to the nearest larger
    ####2-powered number+1, for example, 209->256+1
    Power_index=2^(ceiling(log2(length(sorted_vector))):1)
    sorted_vector=c(sorted_vector, rep(Inf, Power_index[1]-length(sorted_vector)+1))
    Power_index=c(Power_index[-1], 1)
    p=Power_index[1]
    if (sorted_vector[p]==target) return(p)
    for (i in 2:length(Power_index)){
        p=ifelse(sorted_vector[p]<target, p+Power_index[i], p-Power_index[i])
        if (sorted_vector[p]==target) return(p)
    }
    return(-1)
}

sorted_vector=sort(sample(1000, 444))
target=666
sum(sorted_vector==target)
sorted_vector[Find_target(sorted_vector, target)]
