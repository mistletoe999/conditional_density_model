#install.packages("restorepoint")  
#install.packages("mvtnorm")  
#install.packages("RandomFields")  

library(restorepoint)           #clone.environment
library(mvtnorm)                #dmvnorm
library(RandomFields)           #RMexp, RFsimulate

source("Integration.R")


#Return the joint pdf via the copula
calDen <- function(x, inten, sigma) {
	p <- qnorm(pexp(x, rate = inten))
	comp1 <- dmvnorm(p, sigma = sigma) 
	comp2 <- dexp(x, rate = inten)
	comp3 <- dnorm(p)

	comp1 * prod(comp2) / prod(comp3)
}




#Create initial density
createInitialDen <- function(initialDen) {
	xlen <- length(initialDen$x)
	ylen <- length(initialDen$y)
	zlen <- length(initialDen$z)
	initialDen$den <- array(0, c(xlen, ylen, zlen))
	for (i in 1:xlen)
		for (j in 1:ylen)
			for (k in 1:zlen) {
				xyz <- c(initialDen$x[i], initialDen$y[j], initialDen$z[k])
				initialDen$den[i, j, k] <- calDen(xyz, initialDen$marginInten, 
					initialDen$sigma)
			}
	initialDen$den <- initialDen$den / integration3D(initialDen)
}



#Randomize the seed of the ramdom field 
RFoptions(seed=NA)

#create a 3D random field
createRandomField <- function(randomField) {
	xlen <- length(randomField$x)
	ylen <- length(randomField$y)
	zlen <- length(randomField$z)

	xseq <- outer(outer(randomField$x, rep(1, ylen)), rep(1, zlen))
	yseq <- outer(outer(rep(1, xlen), randomField$y), rep(1, zlen))
	zseq <- outer(outer(rep(1, xlen), rep(1, ylen)), randomField$z)

	model <- RMexp(scale = 1 / randomField$rho)  
	method <- RPspectral(model, sp_lines = 100)
	sim <- RFsimulate(method, x = xseq, y = yseq, z = zseq)
	randomField$val <- array(sim$variable1, c(xlen, ylen, zlen))
}



createRandomFieldList <- function(randomFieldList) {
	randomField  <- new.env()
	randomField$x <- randomFieldList$x 
	randomField$y <- randomFieldList$y
	randomField$z <- randomFieldList$z
	randomField$rho <- randomFieldList$rho
	for (i in 1:randomFieldList$N) {
		randomFieldList$den[[i]] <- createRandomField(randomField)
	} 
	rm(randomField)
}



#create the dynamic of the density
createDynamicDen2 <- function(dynamicDen, randomFieldList) {

	dynamicDen$x <- dynamicDen$initial$x
	dynamicDen$y <- dynamicDen$initial$y
	dynamicDen$z <- dynamicDen$initial$z
	dynamicDen$den <- list()

	xlen <- length(dynamicDen$x)
	ylen <- length(dynamicDen$y)
	zlen <- length(dynamicDen$z)
	dt <- dynamicDen$T / dynamicDen$N
 
	randomField  <- new.env()
	randomField$x <- dynamicDen$x
	randomField$y <- dynamicDen$y
	randomField$z <- dynamicDen$z
	randomField$rho <- dynamicDen$rho

	den <- clone.environment(dynamicDen$initial)
	for (tt in 1:dynamicDen$N) {
		
		randomField$val <- randomFieldList$den[[tt]]	
	
		x_seq <- (dynamicDen$x - tt * dt) ^ 2
		y_seq <- (dynamicDen$y - tt * dt) ^ 2    
		z_seq <- (dynamicDen$z - tt * dt) ^ 2
 		h <- sqrt(outer(outer(x_seq, y_seq, FUN = "+"), z_seq, FUN = "+"))
		sigma <- dynamicDen$a * exp(- dynamicDen$b * h)
		change <- exp(sigma * sqrt(dt) * randomField$val - sigma ^ 2 * dt / 2)
		den$den <- den$den * change

		dynamicDen$den[[tt]] = den$den / integration3D(den)

	}

	rm(den)
}


#get density at time t
getDen <- function(t, dynamicDen) {

	den  <- new.env()
	den$x <- dynamicDen$x
	den$y <- dynamicDen$y
	den$z <- dynamicDen$z

	indicator <- which(abs(t - dynamicDen$t) < 0.0001)
	den$den <- dynamicDen$den[[indicator]]
	den
}



getMarginDen <- function(den, x = NULL, y = NULL, z = NULL) {

	marginDen  <- new.env()

	if (is.null(x) & is.null(z)) {
		marginDen$x <- den$x
		marginDen$y <- den$z
		indicator <- which(abs(y - den$y) < 0.0001)
		marginDen$den <- den$den[,indicator,]		
	}

	if (is.null(x) & is.null(y)) {
		marginDen$x <- den$x
		marginDen$y <- den$y
		indicator <- which(abs(z - den$z) < 0.0001)
		marginDen$den <- den$den[,,indicator]
	}

	marginDen$den <- marginDen$den / integration2D(marginDen)
	marginDen
}
