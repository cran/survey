## Paley construction of Hadamard matrices
## Only implemented for GF(p), because it's
## not entirely straightforward to find
## representations of GF(p^m)

paley<-function(n, nmax=2*n){

  ## these are primes with p+1 a multiple of 4
  small.primes<-c(3, 7, 11, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83, 103, 107, 
                  127, 131, 139, 151, 163, 167, 179, 191, 199, 211, 223, 227, 239, 
                  251, 263, 271, 283, 307, 311, 331, 347, 359, 367, 379, 383, 419, 
                  431, 439, 443, 463, 467, 479, 487, 491, 499, 503, 523, 547, 563, 
                  571, 587, 599, 607, 619, 631, 643, 647, 659, 683, 691, 719, 727, 
                  739, 743, 751, 787, 811, 823, 827, 839, 859, 863, 883, 887, 907, 
                  911, 919, 947, 967, 971, 983, 991, 1019, 1031, 1039, 1051, 1063, 
                  1087, 1091, 1103, 1123, 1151, 1163, 1171, 1187, 1223, 1231, 1259, 
                  1279, 1283, 1291, 1303, 1307, 1319, 1327, 1367)


  if (n>max(small.primes)) return(NULL)
  p<-min(small.primes[small.primes>=n])
  if (p>nmax) return(NULL)
  
  m<-outer(0:(p-1) ,0:(p-1),"+") %% p
  
  res<-integer(1+floor((p-1)/2))
  res[1]<-0
  res[2]<-1
  for(i in 2:floor((p-1)/2))
    res[i+1]<- (i*i) %% p
  
  m[m %in% res]<-0
  m[m>0]<-1
  
  cbind(1,rbind(1,m))
  
}
