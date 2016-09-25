setwd("~/data/")


data <- read.table("Experiment1Data.txt", header = T, colClasses = "numeric")

CVA <- aggregate(CVA ~ rho0, data, mean)
DVA <- aggregate(DVA ~ rho0, data, mean)

plot(DVA$rho0, DVA$DVA, type="o", col="blue", xlab=expression(rho[0]), ylab="CVA/DVA")
lines(CVA$rho0, CVA$CVA, type="o", pch=22, lty=2, col="red")
legend(0.65, 0.28, c("CVA on CDS","DVA on CDS"), cex=0.8, col=c("red", "blue"), pch=c(22,21), lty=c(2,1))








data <- read.table("Experiment2Data_0.9.txt", header = T, colClasses = "numeric")
data <- read.table("Experiment2Data_0.3.txt", header = T, colClasses = "numeric")
data <- read.table("Experiment2Data_0.4.txt", header = T, colClasses = "numeric")
data <- read.table("Experiment2Data.txt", header = T, colClasses = "numeric")

CVA <- aggregate(CVA ~ a, data, mean)
DVA <- aggregate(DVA ~ a, data, mean)




plot(DVA$a, DVA$DVA, type="o", col="blue", xlab="a", ylab="CVA/DVA", ylim=c(0,0.205))
lines(CVA$a, CVA$CVA, type="o", pch=22, lty=2, col="red")
legend(0.62, 0.125, c("CVA on CDS","DVA on CDS"), cex=0.8, col=c("red", "blue"), pch=c(22,21), lty=c(2,1))









