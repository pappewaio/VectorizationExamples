#So here is one I was working on last week. There is a function, acepack::ace(), that returns a result that can then 
#be used to calculate a non-linear correlation like:
# z <- acepack::ace(x, y); correlation <- cor(z$tx, z$ty
#The problem is that the ace function doesn't appear to be vectorized and therefore can only take 2 vectors at a time as input.
#I was trying to calculate the correlation matrix for the columns of a 22662x672 matrix (genes vs samples). 
#Intuitivley this was a nested for loop solution which I managed to vectorize in the end but not without having to reduce the number
#of genes used so that the matrix could be held in memory. This is potentially more a "big data" problem than a vectorization problem but
#also provides an opportunity to illustrate how to work in situations where memory (as you said, vectorizations biggest weakness) is an issue.

#My solution looked like this although I'll cut down the size of the matrix to 50x672:

#######make data
mat <- matrix(
  sample(seq(0, 10^5, 1), size=2500, replace=TRUE), 
  ncol=50,
  dimnames=list(
    paste(letters, 1:50, sep=""), 
    paste(LETTERS, 1:50, sep="")
  )
)

######nested for loop style
out1 <- matrix(NA, ncol=50, nrow=50)
for(i in 1:ncol(mat)) {
  for(y in 1:ncol(mat)) {
    t <- acepack::ace(mat[, i], mat[ ,y])
    out1[i, y] <- cor(t$tx, t$ty)
  }
}
colnames(out1) <- colnames(mat)
rownames(out1) <- colnames(mat)

#######my solution

#create a matrix to hold all comparisons (this transposes the matrix)
m1 <- matrix(
  rep(
    as.integer(mat), #integers take less memory
    ncol(mat)
  ),
  ncol=nrow(mat),
  byrow=TRUE
)
rownames(m1) <- rep(colnames(mat), ncol(mat))
    
o <- order(rownames(m1)) #prevents me from creating another large variable to hold the 2nd argument to ace().

#a unit ot do the correlation
.inner <- function(x, y) {
  out <- acepack::ace(x, y)
  cor(out$tx, out$ty)
}

#not sure if sapply strictly constitutes a vectorized solution (unclear to me as well/pappewaio)
res <- sapply(1:nrow(m1), function(j)
  .inner(m1[j,], m1[o[j],])
)

#order and rename everything so it is as expected in the final output
out2 <- matrix(res, ncol=ncol(mat))
rownames(out2) <- unique(rownames(m1)[1:nrow(m1)])
colnames(out2) <- unique(rownames(m1)[o[1:nrow(m1)]])
out2 <- out2[,match(rownames(out2), colnames(out2))]

########solutions identical?
identical(out1, out2)

#so there are 2 things with the example above:
#a) even though my solution is vectorized (or at least not a nested for loop), the code becomes very difficult to follow. This is a 
#typical reason that I sometimes choose to use for loops, i.e. readability. 
#b) tracking which correlation corresponds to which sample is not so straight forward


#########Pappewaio solution
# first aim: try to use apply to gain speed. (as my belief is that apply is faster than sapply, this needs to be tested)

#usually things are ordered in the same way as they were sent in. That is another unsaid rule when it comes to vectorization.
#if we want to minimize our memory foot print
storage.mode(mat)
mat2 <- apply(mat, c (1, 2), as.integer)
storage.mode(mat2)
# create an index of comparisons to makeo 
inx <- t(combn(1:nrow(mat2),2))
#if interested check that we have no redundancy(not necessary)
#tf <- apply(test,2,function(x){x[1]==49 & x[2]==50})
#sum(tf)
#tf <- apply(test,2,function(x){x[1]==50 & x[2]==49})
#sum(tf)

.inner.pap <- function(x, y){
  out <- acepack::ace(x, y)
  cor(out$tx, out$ty)
}

#mapply version (returns vector)
res <- mapply(inx[,1], inx[,2], FUN=function(r1, r2, m){
      .inner.pap(m[,r1], m[,r2])
  }, MoreArgs=list(m=mat2)
)
out3 <- matrix(NA, ncol=ncol(mat),nrow=nrow(mat))
out3[inx] <- res

#apply version (returns vector)
res <- apply(inx, 1, function(i, m){
      .inner.pap(m[i[1],], m[i[2],])
  }, m=mat2)

out4 <- matrix(NA, ncol=ncol(mat),nrow=nrow(mat))
out4[inx] <- res

#To get a more sofisticated version, you could rewrite the ace function to accept the matrix and index and to only calculate exactly what you are after, but in this case the ace function actually accepts a matrix (disclaimer: the matrix input seems not to yield the same result, and it seems like there is fortran code in there so not so easy to rewrite it to accept an index.).

y <- mat[,1]
x <- mat
#this one could be
out <- acepack::ace(x, y)
c1a <- cor(out$tx, out$ty)
out <- acepack::ace(t(x), y)
c1b <- cor(out$tx, out$ty)

#the same as this one
out <- acepack::ace(x[,2], y)
c2 <- cor(out$tx, out$ty)

identical(c1a[2,], c2)
identical(c1b[2,], c2)

any(as.vector(c1a) == as.vector(c2))
any(as.vector(c1b) == as.vector(c2))

#Fail.. it seems not that the matrix input yields the same result




##################################################
#So here is another nasty one.
#Below there is a matrix where each row contains some values.
#Dependant on the values we consider some of the cases (columns) related.
#We specify those values that indicate a relationship with a cutoff.
#The next step is to assemble all of the related cases and count the number of times we find them to be related. 

#set up some example data
set.seed(234)
mat <- matrix(
    sample(c(seq(0, 0.2, 10^-4), seq(0.8, 1, 10^-4)), size=15, replace=FALSE),
   ncol=5,
   dimnames=list(letters[1:3], LETTERS[1:5])
)

#specify the cutoff
cutoff <- 0.5

#by looking at the mat variable, and considering the cutoff, we can already here see what the expected result is:
#row 1 related cases: AB, AC, AE, BC, BE, CE
#row 2 related cases: AC, AD, AE, CD, CE, DE
#row 3 related cases: BC, BD, CD
#so the expected result is:
# AB: 1
# AC: 2
# AD: 1
# AE: 2
# BC: 2
# BD: 1
# BE: 1
# CD: 2
# CE: 2
# DE: 1

#We also want the counts to be the same for the reverse, i.e. if AB = 1 then BA should be 1
#All the rest should be 0.

#begin calculations
tMat <- t(mat)

#first we just find the values over the cutoff and this is vectorized so... no problem here.
sobj.bool <- tMat > cutoff

#set up a dataframe to store results
results <- data.frame(
  x=rep(cases, each=ncases),
  y=rep(cases, ncases),
  times=rep(0, ncases^2)
)

#set up variables necessary for calculations
cases <- rownames(tMat)
ncases <- length(cases)

#find related cases
for(i in 1:(dim(sobj.bool)[2])) {
  o <- which(sobj.bool[,i])
  if(length(o) > 1) {
    for(j in 1:(length(o)-1)) {
      for(k in (j+1):length(o)) {
        ind1 <- (o[j]-1)*ncases+o[k]
        results[ind1,3] <- results[ind1,3]+1
        ind2 <- (o[k]-1)*ncases+o[j]
        results[ind2,3] <- results[ind2,3]+1
      }
    }
  } 
}

##################################################
#So here is an easier one (maybe)
#I have a list where each element (a, b, c) has several other elements (A, B) like this:
l <- list(
  a=list(A=letters[1:10], B=LETTERS[11:20]), 
  b=list(A=letters[1:10], B=LETTERS[11:20]), 
  c=list(A=letters[1:10], B=LETTERS[11:20])
)

# so I want to extract each sub-element and construct a dataframe like this:
As <- data.frame()
Bs <- data.frame()
for(i in 1:length(l)) {
  As <- rbind(As, data.frame(t(l[[i]][['A']])))
  Bs <- rbind(Bs, data.frame(t(l[[i]][['B']])))      
}

#vectorized solution?
