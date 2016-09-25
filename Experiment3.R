setwd("~/data/")


source("Experiment0.R")




dataFile <- "Experiment3Data.txt"

heading_names <- c("rho0", "rho", "a", "b", "i", "CVA", "DVA", "BCVA")
write.table(t(heading_names), dataFile, append = T, col.names = F, row.names = F)



#Experiment 3: rho


rho0 <- 0.4
rhoList <- c(0.0001, 0.0005, 0.0010, 0.0050, 0.01, 0.05, 0.1, 0.5, 1)

filename <- paste("initialDen_", formatC(rho0, 2, format = 'f'), ".RData", sep = "")
load(filename)


coordinate <- c(c(0.001, 1/ 96, 1 / 48, 1/ 24, 3 / 48, 1 / 12),
	seq(3 / 24, 1 / 2, 1 / 24), seq(7 / 12, 5 - 1 / 12, 1 / 12),
	5 ^ (seq(1, 3.1, 2.1 / 60)))


for (rho in rhoList) {

	randomFieldList <- new.env()
	randomFieldList$x <- coordinate
	randomFieldList$y <- coordinate
	randomFieldList$z <- coordinate
	randomFieldList$rho <- rho
	randomFieldList$N <- 60


	createRandomFieldList(randomFieldList)

	paras <- list()
	paras$rho0 <- rho0
	paras$rho <- rho
	paras$a <- 0.4
	paras$b <- 0.02
	paras$i <- 1

	
	BCVA <- CVADVA(initialDen, randomFieldList, paras)

	data <- c(as.vector(paras), BCVA)
	write.table(t(data), dataFile, append = T, col.names = F, row.names = F)

	rm(randomFieldList)	

}

rm(initialDen)
















