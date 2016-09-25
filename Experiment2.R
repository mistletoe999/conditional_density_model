setwd("~/data/")


source("Experiment0.R")




dataFile <- "Experiment2Data.txt"

aList <- seq(0.1, 0.8, 0.1)


heading_names <- c("rho0", "rho", "a", "b", "i", "CVA", "DVA", "BCVA")
write.table(t(heading_names), dataFile, append = T, col.names = F, row.names = F)






#Experiment 2: a


rho0 <- 0.4

filename <- paste("initialDen_", formatC(rho0, 2, format = 'f'), ".RData", sep = "")
load(filename)


for (i in 1:100) {
	filename <- paste("randomFieldList_", 
		formatC(0.0050, 4, format = 'f'), "_", 
		formatC(i, 3, format = 'd',  flag = "0"), ".RData", sep = "")
	load(filename)

	for (a in aList) {
	

		paras <- list()
		paras$rho0 <- rho0
		paras$rho <- 0.0050
		paras$a <- a
		paras$b <- 0.02
		paras$i <- i

	
		BCVA <- CVADVA(initialDen, randomFieldList, paras)


		data <- c(as.vector(paras), BCVA)
		write.table(t(data), dataFile, append = T, col.names = F, row.names = F)	
					
	}
	
	rm(randomFieldList)
}


rm(initialDen)
























