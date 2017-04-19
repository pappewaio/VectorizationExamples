
################################################################################
#Example 1 - simplest case
#add 1 to each element in a vector
################################################################################
#create vector 1 to 10
x0 <- 1:10

#for loop version
x1 <- x0
for (i in 1:length(x0)){
	x1[i] <- x0[i]+1
}

#vectorized version
x2 <- x0+1

identical(x1,x2)

################################################################################
#Example 2 - sequential standard operations
#multiply with 3 and then add 2 to each element in a vector
################################################################################
#create vector 1 to 10
x0 <- 1:10

#for loop version
x1 <- x0
for (i in 1:length(x0)){
	x1[i] <- (x0[i]*3)+2
}

#vectorized version
x2 <- (x0*3)+2

identical(x1,x2)

################################################################################
#Example 3 - if cases
#multiply with 3 if the element is larger than 5 for each element in a vector
################################################################################
#create vector 1 to 10
x0 <- 1:10

#for loop version
x1 <- x0
for (i in 1:length(x0)){
	if(x0[i] > 5) x1[i] <- (x0[i]*3)
}

#vectorized version
x2 <- x0
x2[x0>5] <- (x0[x0>5]*3)

identical(x1,x2)

#Note: it can be achieved in one line if we replace the elements in x0 instead
#      creating the extra vector x2. 

################################################################################
#Example 4 - nested for loop
# Add 2 to each element
################################################################################
#create vector 1 to 10
x0 <- matrix(1:100, ncol=10)

#for loop version
x1 <- x0
for (i in 1:nrow(x0)){
	for (j in 1:ncol(x0)){
		x1[i,j] <- x0[i,j]+2
	}
}

#vectorized version
x2 <- x0+2

identical(x1,x2)


################################################################################
#Example 5 - multiple manipulations
# Add 2 to each element if lower than 5
# Add 5 to each element if higher than 10
################################################################################
#create vector 1 to 10
x0 <- 1:10

#for loop version
x1 <- x0
for (i in 1:length(x0)){
	if(x0[i] < 5) x1[i] <- (x0[i]+2)
	if(x0[i] > 10) x1[i] <- (x0[i]+5)
}

#vectorized version
x2 <- x0
x2[x0<5] <- (x0[x0<5]+2)
x2[x0>10] <- (x0[x0>10]+5)

identical(x1,x2)

#Note: We had to pay the price of having two lines of code to cover the two
#      if statements, but observe there is no need to use more variables. 
#      The price we pay for vectorization is primary memory. But usually
#      it is easy to design the code so if it becomes expensive in respect
#      to RAM, then we can split up in pieces that fit into the RAM and run
#      them sequentially like a for loop, or even better to distribute to
#      different nodes on a clustes for a simulatneous evaluation. 





