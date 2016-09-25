setwd("~/data/")


source("Experiment0.R")




dataFile <- "Experiment1Data.txt"

rho0List <- c(seq(0.1, 0.9, 0.1), 0.99)


heading_names <- c("rho0", "rho", "a", "b", "i", "CVA", "DVA", "BCVA")
write.table(t(heading_names), dataFile, append = T, col.names = F, row.names = F)


#Experiment 1: rho0
for (i in 1:100) {
	filename <- paste("randomFieldList_", 
		formatC(0.0050, 4, format = 'f'), "_", 
		formatC(i, 3, format = 'd',  flag = "0"), ".RData", sep = "")
	load(filename)


	for (rho0 in rho0List) {
		filename <- paste("initialDen_", formatC(rho0, 2, format = 'f'), ".RData", sep = "")
		load(filename)

		paras <- list()
		paras$rho0 <- rho0
		paras$rho <- 0.0050
		paras$a <- 0.4
		paras$b <- 0.02
		paras$i <- i

	
		BCVA <- CVADVA(initialDen, randomFieldList, paras)


		data <- c(as.vector(paras), BCVA)
		write.table(t(data), dataFile, append = T, col.names = F, row.names = F)	

		rm(initialDen)	
		
	}

	
	rm(randomFieldList)
}













