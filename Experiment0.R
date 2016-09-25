source("Integration.R")
source("Density.R")
source("CDS.R")


CVADVA <- function(initialDen, randomFieldList, paras) {
	dynamicDen  <- new.env()
	dynamicDen$T <- 5
	dynamicDen$N <- 60
	dynamicDen$a <- paras$a 
	dynamicDen$b <- paras$b 
	dynamicDen$initial <- initialDen 
	dynamicDen$rho <- randomFieldList$rho
	dynamicDen$t <- seq(dynamicDen$T / dynamicDen$N, dynamicDen$T,
		dynamicDen$T / dynamicDen$N)


	createDynamicDen2(dynamicDen, randomFieldList)



	recovery <- c(0.4, 0.4, 0.4)
	r <- 0.03

	CDS <- new.env()
	CDS$maturity <- 5
	CDS$freq <- 0.25
	CDS$Ti <- seq(CDS$freq, CDS$maturity, CDS$freq)
	CDS$rate <- priCDS(CDS, initialDen, 0.03, 0.4)


	CVA <- evalCVA0(CDS, dynamicDen, r, recovery) 
	DVA <- evalDVA0(CDS, dynamicDen, r, recovery) 
	BCVA <- CVA - DVA

	c(CVA, DVA, BCVA)

}


