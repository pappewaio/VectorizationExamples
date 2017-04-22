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
#storage.mode(mat)
#mat2 <- apply(mat, c (1, 2), as.integer)
#storage.mode(mat2)
# create an index of comparisons to makeo 
#inx <- t(combn(1:nrow(mat2),2))
#if interested check that we have no redundancy(not necessary)
#tf <- apply(test,2,function(x){x[1]==49 & x[2]==50})
#sum(tf)
#tf <- apply(test,2,function(x){x[1]==50 & x[2]==49})
#sum(tf)

sol3 <- function(mat){
	mat2 <- apply(mat, c (1, 2), as.integer)
	inx <- t(combn(1:nrow(mat2),2))

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
	out3
}

res3 <- sol3(mat)

sol4 <- function(mat){
	mat2 <- apply(mat, c (1, 2), as.integer)
	inx <- t(combn(1:nrow(mat2),2))

	.inner.pap <- function(x, y){
	  out <- acepack::ace(x, y)
	  cor(out$tx, out$ty)
	}
	#apply version (returns vector)
	res <- apply(inx, 1, function(i, m){
	      .inner.pap(m[i[1],], m[i[2],])
	  }, m=mat2)
	
	out4 <- matrix(NA, ncol=ncol(mat),nrow=nrow(mat))
	out4[inx] <- res
	out4
}

res4 <- sol4(mat)

#To get a more sofisticated version, you could rewrite the ace function to accept the matrix and index and to only calculate exactly what you are after, but in this case the ace function actually accepts a matrix (disclaimer: the matrix input seems not to yield the same result, and it seems like there is fortran code in there so not so easy to rewrite it to accept an index.).

#y <- mat[,1]
#x <- mat
##this one could be
#out <- acepack::ace(x, y)
#c1a <- cor(out$tx, out$ty)
#out <- acepack::ace(t(x), y)
#c1b <- cor(out$tx, out$ty)
#
##the same as this one
#out <- acepack::ace(x[,2], y)
#c2 <- cor(out$tx, out$ty)
#
#identical(c1a[2,], c2)
#identical(c1b[2,], c2)
#
#any(as.vector(c1a) == as.vector(c2))
#any(as.vector(c1b) == as.vector(c2))

#Fail.. it seems not that the matrix input yields the same result

## Another version where we calculate the cor as a second step.
sol5 <- function(mat){
	mat2 <- apply(mat, c (1, 2), as.integer)
	inx <- t(combn(1:nrow(mat2),2))

	.inner.pap <- function(x, y){
	  out <- acepack::ace(x, y)
	  cor(out$tx, out$ty)
	}

	res <- apply(inx, 1, function(i, m){
	  		out <- acepack::ace(m[i[1],], m[i[2],])
			c(out$tx, out$ty)
	  }, m=mat2)
	
	res2 <- apply(res,2, function(m2, g1){
				cor(m2[1:g1], m2[(g1+1):length(m2)])
	  }, g1=50)
	
	out5 <- matrix(NA, ncol=ncol(mat),nrow=nrow(mat))
	out5[inx] <- res2
	out5
}

res5 <- sol5(mat)

#as the cor function does not allow indexes to be used, we cant vectorize furthe
# without rewriting that function as well. 

#benchmark the solutions
library(microbenchmark)
microbenchmark(sol3(mat), sol4(mat), sol5(mat))

#Unit: milliseconds
#      expr      min       lq     mean   median       uq      max neval
# sol3(mat) 459.2650 494.2500 515.0101 504.0770 519.7576 793.7025   100
# sol4(mat) 470.1018 492.8400 514.1665 506.5867 520.3040 788.7378   100
# sol5(mat) 499.9738 529.1344 551.3592 543.3093 560.5120 705.9766   100

#The solutions don't differ that much with this example although these differences could be larger with more complex problems.
#On average sol4 is the fastest although sol3 is very close to being the same. sol5 have the worst performance on average.

#Ran it on my machine as well, same conclusion /pappewaio
#> microbenchmark(sol3(mat), sol4(mat), sol5(mat))
#Unit: milliseconds
#      expr      min       lq     mean   median       uq      max neval
# sol3(mat) 430.7390 465.7790 488.4251 480.3787 506.7110 651.4949   100
# sol4(mat) 424.9015 472.9759 494.1501 492.7945 516.0457 584.3677   100
# sol5(mat) 433.3937 482.0668 512.0837 509.7396 533.9748 684.4709   100


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

#set up variables necessary for calculations
cases <- rownames(tMat)
ncases <- length(cases)

#set up a dataframe to store results
results <- data.frame(
  x=rep(cases, each=ncases),
  y=rep(cases, ncases),
  times=rep(0, ncases^2)
)


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

##pappewaio vectorized solution
sol1 <- function(mat, cutoff){
	logic <- mat > cutoff
	totcomb <- apply(t(combn(colnames(mat),2)),1,paste,collapse="")
	funx <- function(row, totcomb){
			totcomb %in% apply(t(combn(names(row)[which(row)],2)),1,paste,collapse="")
	}
	res <- apply(apply(logic, 1, funx, totcomb),1, sum)

	#make data comparable to other solutions
	xy <- t(combn(colnames(mat),2))
	data.frame(x=xy[,1],y=xy[,2],times=res)
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
rownames(As) <- names(l)
rownames(Bs) <- names(l)

#vectorized solution?
#here is one
As2 <- data.frame(t(sapply(l, "[[", 1)))
Bs2 <- data.frame(t(sapply(l, "[[", 2)))

identical(As, As2) #fails because each factor in the columns are named in As2 and not in As; As[,1]; As2[,1]
identical(Bs, Bs2) #fails because each factor in the columns are named in Bs2 and not in Bs; Bs[,1]; Bs2[,1]
all.equal(As, As2)
all.equal(Bs, BS2)

#####
#pappewaio solution
#####

sol2 <- function(l){
	#As all elements have the same length we can just unlist everything
	ar <- array(unlist(l, use.names=FALSE),c(10,2,3))
	As3 <- as.data.frame(t(ar[,1,]))
	Bs3 <- as.data.frame(t(ar[,2,]))
	list(As3, Bs3)
}

##################################
#Challange 4
#The overall goal of this challange is pretty straight forward. We have a group of different elements that interact with each other at
#different frequencies. Constraints for the frequencies of interaction for each individual element are generally known. The challange is
#to derive a dataset of X elements that interact within the known frequency constraints.
#The challange could be divided into 2 stages. 
#The overall goal is to a) set up an interaction table that dictates the frequency of interactions between a number of elements.
#Once this is accomplished we want to b) derive a dataset comprised of x interaction pairs where the whole dataset reflects the interaction
#frequencies portrayed in the table.

#Here is an example:

#say we have 5 elements:
elem <- LETTERS[1:5]

#the constraints for the elements interactions are as follows:
#a) each element has a relativley large probability of interacting with itself.
#b) each element has a relativley large probability of interacting with a subset (1-3) of the other elements but
#c) a relativley small probability of interacting equally with all other elements.
#A probability distribution of the interactions of each element with another element might look something like this:

x <- c(rnorm(1000, 20, 10), rnorm(200, 70, 10))
hist(x[x>1 & x<100], breaks=100, main="Interaction probability", xlab="Probability x 100")

#Despite this, drawing from such a distribution is most likely to give a high proportion of elements that only react weakly 
#with all other elements and that is not desired. Thus we need a somewhat more structured approach to get the desired results.
#Such an attempt might generate the following:

a <- c(47, 34, 1, 17,  1)
b <- c(5, 56, 2, 36, 1)
c <- c(1, 3, 50, 2, 44)
d <- c(45, 1, 1, 52, 1)
e <- c(1, 1, 5, 49, 44)
intMat <- matrix(c(a,b,c,d,e), ncol=length(elem), nrow=length(elem), dimnames=list(elem, elem))
intMat

#The matrix describes the probability of reaction with each element with each other element in each column where each column sums to 100%
#thus representing all interactions for that element.
#For example we can see that element A interacts strongly with itself (47) and with B and D (34, 17) but weakly 
#with other elements.
#Although intMat describes the interactions in a manner that follows the specificed constraints. We would desire a way to do this
#dynamically and with up to 20 elements.

#In stage 2 we want to compose a dataset of interacting element sets (with up to 5 elements per set) that mirrors the frequencies of
#interactions in the interaction table.
#We can first generate all the possible desired interaction types using the following function:

allInts <- function(n, elem) {
	switch(n-1,
	       comb <- expand.grid(elem, elem),
	       comb <- expand.grid(elem, elem, elem),
	       comb <- expand.grid(elem, elem, elem, elem),
	       comb <- expand.grid(elem, elem, elem, elem, elem)
	)
	dat.sort <- t(apply(comb, 1, sort))
	dat <- comb[!duplicated(dat.sort),]
	return(dat)
}

two <- apply(allInts(2, elem), 1, paste, collapse="")
three <- apply(allInts(3, elem), 1, paste, collapse="")
four <- apply(allInts(4, elem), 1, paste, collapse="")
five <- apply(allInts(5, elem), 1, paste, collapse="")
combos <- c(two, three, four, five)

#Now the variable called combos includes all desired combinations of the elments and the goal is to pick from these elements so that
#the interaction frequency is mirrored by that in the interaction table. In addition, we would like interaction elements from each
#interaction type (two, three, four, five) in the final result.
#The following function may be useful which calculates the interactions in the input variable:

intFreq <- function(inp) {
	s2 <- strsplit(gsub("(.{1})", "\\1 ", inp), " ")
	com <- lapply(s2, function(x) combn(x, 2))
	nam <- lapply(com, function(j) sapply(1:ncol(j), function(u) paste(j[1,u], j[2,u], sep="")))
	return(table(unlist(nam)))
}

intFreq(combos)

#We can see that the dataset currently has 120 AA interactions out of a total of 600 A interactions (20%) but, 
#according to the interaction matrix we want 47% of A interactions to be AA interactions.
#Thus, the challange here is to adjust the interactions in the combos variable so that the combos variable, in the end, reflects the
#interaction frequency in the interactions table.
#Note that acheiving this goal with a small (<5) number of elements is potentially impossible. 
#Ideally we would be able to solve the problem for element numbers between 5-20.

#As an example, although with a very limited number of elements, with the following interaction table:
matrix(
	c((20/44/1)*100, (24/44/1)*100, (24/44/1)*100, (20/44/1)*100), 
	ncol=2, 
	nrow=2, 
	dimnames=list(LETTERS[1:2], LETTERS[1:2])
)

#A sucessful result would be:
res <- combos[c(1,2,6,16,17,21,56,66,126,136,121,191)]
intFreq(res)

#The results show a total AA interactions to be 20 which is 20/(20+24)*100 of all A interactions and 
#corresponds with the interaction table. The same is true for BB.
#Finally AB and BA interactions (here only shown as BA) are 24/(20+24)*100 and also correspond with the
#desired number in the interaction table.
