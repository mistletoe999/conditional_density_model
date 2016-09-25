#setwd("C:/Users/WC/Dropbox/credit_risk/codes/")
#setwd("~/Dropbox/Credit_risk/codes/")
#setwd("~/Dropbox/QuantShare/codes/")

source("Density.R")
source("Integration.R")
source("CDS.R")


#Step 1: create initial density

coordinate <- c(c(0.001, 1/ 96, 1 / 48, 1/ 24, 3 / 48, 1 / 12),
	seq(3 / 24, 1 / 2, 1 / 24), seq(7 / 12, 5 - 1 / 12, 1 / 12),
	5 ^ (seq(1, 3.1, 2.1 / 60)))

initialDen <- new.env()
initialDen$x <- coordinate
initialDen$y <- coordinate
initialDen$z <- coordinate
initialDen$marginInten <- c(0.2, 0.06, 0.06)

for (rho in seq(0.0, 0.9, 0.1)) {
	initialDen$sigma <- matrix(c(1, rho, rho, rho, 1, rho, rho, rho, 1), nrow = 3) 
	createInitialDen(initialDen)
	filename <- paste("initialDen_", formatC(rho, 2, format = 'f'), ".RData", sep = "")
	save(initialDen, file = filename)
}


#load("initialDen_0_4.RData")

#test: the sum should be 1
#integration3D(initialDen)



#Step 2: simulate Gaussian random fields
randomField  <- new.env()
randomField$x <- coordinate
randomField$y <- coordinate
randomField$z <- coordinate
randomField$rho <- 0.0050
createRandomField(randomField)



randomFieldList <- new.env()
randomFieldList$x <- coordinate
randomFieldList$y <- coordinate
randomFieldList$z <- coordinate
randomFieldList$rho <- 0.0050
randomFieldList$N <- 60


#setwd("~/Documents/")
#setwd("D:/data/")
for (i in 5:100) {
	createRandomFieldList(randomFieldList)
	filename <- paste("randomFieldList_", 
		formatC(randomFieldList$rho, 4, format = 'f'), "_", 
		formatC(i, 3, format = 'd',  flag = "0"), ".RData", sep = "")
	save(randomFieldList, file = filename)
}



#Step 2: generate dynamic density
dynamicDen  <- new.env()
dynamicDen$T <- 5
dynamicDen$N <- 60
dynamicDen$a <- 0.2
dynamicDen$b <- 0.02
dynamicDen$initial <- initialDen 
dynamicDen$rho <- 0.005
dynamicDen$t <- seq(dynamicDen$T / dynamicDen$N, dynamicDen$T,
	dynamicDen$T / dynamicDen$N)

t <- Sys.time()
createDynamicDen(dynamicDen, randomFieldList)
Sys.time() - t



#Step 3: evaluate CDS
CDS <- new.env()
CDS$maturity <- 5
CDS$freq <- 0.25
CDS$Ti <- seq(CDS$freq, CDS$maturity, CDS$freq)
CDS$rate <- priCDS(CDS, initialDen, 0.03, 0.4)

evalCDS2(0.1, CDS, dynamicDen, 0.03, 0.4)




#Step 4: evaluate CVA on CDS



