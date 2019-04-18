###Use dynamic programming to determine if a 
###combined string is the interleaf of two strings
string1=c('A','G','C','G')
string2=c('T','T','C','G','T','A','C')
combined=c('A','T','G','C','T','C','G','T','A','C','G')

len1=length(string1)
len2=length(string2)
##intialize the matrix
solution=matrix(NA, len1+1, len2+1)
solution[1,1]=TRUE

for (i in 1:len1){
    solution[i+1,1] = string1[i]==combined[i] & solution[i,1]
}
for (k in 1:len2){
    solution[1,k+1] = string2[k]==combined[k] & solution[1,k]
}

#fill the matrix
for (i in 1:len1){
    for (k in 1:len2){
        solution[i+1, k+1]=(solution[i+1,k] & combined[i+k]== string2[k]) | 
            (solution[i, k+1] & combined[i+k] == string1[i])
    }
}
