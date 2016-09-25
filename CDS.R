source("Integration.R")
source("Density.R")


#calculte the CDS rate at time 0
priCDS <- function(CDS, initialDen, r, recovery) {
	fun <-  function(x, y, z) exp(-r * x) 
	range <- list()
	range$x <- c(0, CDS$maturity)
	term1 <- integration3D(initialDen, fun = fun, range = range)

	N <- round(CDS$maturity / CDS$freq)

	term2 <- 0
	for(Ti in CDS$Ti) {
		range <- list()
		range$x <- c(Ti, 1000)		
		prob <- integration3D(initialDen, range = range)
		term2  <- term2 + CDS$freq * exp(- r * Ti) * prob
		#print(CDS$freq * exp(- r * Ti) * prob)
		fun <-  function(x, y, z) exp(-r * x) * (x - Ti + CDS$freq)
		range$x <- c(Ti - CDS$freq, Ti)	
		
		term2  <- term2 + integration3D(initialDen, fun = fun, range = range)	
	}

	term1 / term2 * (1- recovery)
}





evalCDS2 <- function(t, CDS, marginDen, r, recovery1) {

	term1 <- 0
	fun <-  function(x, y) exp(-r * (x - t)) 
	range <- list()
	if (t < (CDS$maturity - 0.0001)) {
		range$x <- c(t, CDS$maturity)
		range$y <- c(t, 1000)
		term1 <- term1 + integration2D(marginDen, fun = fun, range = range)
	}

	term2 <- 0
	for(Ti in CDS$Ti) {
		if (Ti >= t) {
			range <- list()
			range$x <- c(Ti, 1000)
			range$y <- c(t, 1000)	
			prob <- integration2D(marginDen, range = range)
			term2  <- term2 + CDS$freq * exp(- r * (Ti  - t)) * prob

			fun <-  function(x, y) exp(-r * (x - t)) * (x - Ti + CDS$freq)
			#print(Ti)
			#print(t)
			if (t < (Ti - 0.0001)) {
				range$x <- c(max(Ti - CDS$freq, t), Ti)	
				range$y <- c(t, 1000)	
				term2  <- term2 + integration2D(marginDen, fun = fun, range = range)
			}		

		}
	}

	range <- list()
	range$x <- c(t, 1000)
	range$y <- c(t, 1000)
	term3 <- integration2D(marginDen, range = range)

	(term1 * (1- recovery1) - term2 * CDS$rate) / term3

}



evalCVA0 <- function(CDS, dynamicDen, r, recovery) {

	CVA <- 0
	for (t in dynamicDen$t) {
		den <- getDen(t, dynamicDen)
		marginDen <- getMarginDen(den, y = t)
		CDSvalue <- evalCDS2(t, CDS, marginDen, r, recovery[1])
		#print(CDSvalue)


		k <- which(abs(t - dynamicDen$t) < 0.0001)
		Tk <- which(abs(t - dynamicDen$y) < 0.0001)

		dt <- dynamicDen$T / dynamicDen$N
		for (j in 1:k) {
			x_seq <- (dynamicDen$x - j * dt) ^ 2
			y_seq <- (dynamicDen$y - j * dt) ^ 2  
			z_seq <- (dynamicDen$z - j * dt) ^ 2
			h <- sqrt(outer(outer(x_seq, y_seq, FUN = "+"), z_seq, FUN = "+"))
			sigma <- dynamicDen$a * exp(- dynamicDen$b * h)


			if (j == 1) {
				pre <- sigma[,Tk,] * randomFieldList$den[[j]][,Tk,]
				pre <- pre - sigma[,Tk,] ^ 2 / 2 * dt
			} else {
				pre <- pre + sigma[,Tk,] * (randomFieldList$den[[j]][,Tk,] 
					- randomFieldList$den[[j - 1]][,Tk,])
				pre <- pre - sigma[,Tk,] ^ 2 / 2 * dt
			}
		

 			integrand <- new.env()
 			integrand$x <- dynamicDen$initial$x
	 		integrand$y <- dynamicDen$initial$y
			integrand$den <- dynamicDen$initial$den[,Tk,] * exp(pre)

			range <- list()
 			range$x <- c(t, 1000)	
 			range$y <- c(t, 1000)	
			integrate <- integration2D(integrand, range = range) 
			CVA <- CVA + max(CDSvalue, 0) * integrate * exp (-r * t) * dt * (1 - recovery[2]) 

		}

	}

	CVA 
}



evalDVA0 <- function(CDS, dynamicDen, r, recovery) {

	DVA <- 0
	for (t in dynamicDen$t) {
		den <- getDen(t, dynamicDen)
		marginDen <- getMarginDen(den, z = t)
		CDSvalue <- evalCDS2(t, CDS, marginDen, r, recovery[1])
		#print(CDSvalue)


		k <- which(abs(t - dynamicDen$t) < 0.0001)
		Tk <- which(abs(t - dynamicDen$z) < 0.0001)

		dt <- dynamicDen$T / dynamicDen$N
		for (j in 1:k) {
			x_seq <- (dynamicDen$x - j * dt) ^ 2
			y_seq <- (dynamicDen$y - j * dt) ^ 2  
			z_seq <- (dynamicDen$z - j * dt) ^ 2
			h <- sqrt(outer(outer(x_seq, y_seq, FUN = "+"), z_seq, FUN = "+"))
			sigma <- dynamicDen$a * exp(- dynamicDen$b * h)


			if (j == 1) {
				pre <- sigma[,,Tk] * randomFieldList$den[[j]][,,Tk]
				pre <- pre - sigma[,,Tk] ^ 2 / 2 * dt
			} else {
				pre <- pre + sigma[,,Tk] * (randomFieldList$den[[j]][,,Tk] 
					- randomFieldList$den[[j - 1]][,,Tk])
				pre <- pre - sigma[,,Tk] ^ 2 / 2 * dt
			}
		

 			integrand <- new.env()
 			integrand$x <- dynamicDen$initial$x
	 		integrand$y <- dynamicDen$initial$y
			integrand$den <- dynamicDen$initial$den[,,Tk] * exp(pre)

			range <- list()
 			range$x <- c(t, 1000)	
 			range$y <- c(t, 1000)	
			integrate <- integration2D(integrand, range = range) 
			DVA <- DVA + max(-CDSvalue, 0) * integrate * exp (-r * t) * dt * (1 - recovery[3]) 

		}

	}

	DVA 
}
