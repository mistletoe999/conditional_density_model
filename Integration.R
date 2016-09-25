library(restorepoint)           #clone.environment



cutDen2D <- function(den, index, range1) {

	xlen <- length(den$x)
	ylen <- length(den$y)

	if (index == "x") {
		newIndex <- den$x >= (range1[1] - 0.0001)  & den$x <= (range1[2] + 0.0001)
		den$x <- den$x[newIndex]
		den$den <- den$den[newIndex,]
	}


	if (index == "y") {
		newIndex <- den$y >= (range1[1] - 0.0001) & den$y <= (range1[2] + 0.0001)
		den$y <- den$y[newIndex]
		den$den <- den$den[,newIndex]
	}

}

#integrate over a surface
integration2D <- function(marginDen, fun = NULL, range = NULL) {

	newDen <- clone.environment(marginDen)

	if (!is.null(range)) {
		if(!is.null(range$x)) 
			cutDen2D(newDen, "x", range$x)
		if(!is.null(range$y)) 
			cutDen2D(newDen, "y", range$y)
	}
	
	xlen <- length(newDen$x)
	ylen <- length(newDen$y)

	if (is.null(fun)) {
			val <- newDen$den
	} else {
		xseq <- outer(newDen$x, rep(1, ylen))
		yseq <- outer(rep(1, xlen), newDen$y)
		val <- fun(xseq, yseq) * newDen$den
	}


	avg <- (val[-xlen, -ylen] + val[-1, -ylen] +
		val[-xlen, -1] + val[-1, -1] ) / 4
	xDiff <- newDen$x[-1] - newDen$x[-xlen]
	yDiff <- newDen$y[-1] - newDen$y[-xlen]

	rm(newDen)
	
	sum(avg * outer(xDiff, yDiff))

}




#cut the 3D cube
cutDen3D <- function(den, index, range) {

	xlen <- length(den$x)
	ylen <- length(den$y)
	zlen <- length(den$z)

	newIndex <- den$x > range[1]  & den$x < range[2]
	bdy <- which(newIndex)
	lower <- max(bdy[1] - 1, 1)
	upper <- min(bdy[length(bdy)] + 1, length(den$x))
	newIndex[lower] <- newIndex[upper] <- TRUE

	if (den$x[lower] < range[1]) {
		w1 <- (den$x[lower + 1] - range[1]) / (den$x[lower + 1] - den$x[lower])
		den$den[lower,,] <- w1 * den$den[lower,,] + (1 - w1) * den$den[lower + 1,,]
		den$x[lower] <- range[1]
	}
	

	if (den$x[upper] > range[2]) {	
		w1 <- (den$x[upper] - range[2]) / (den$x[upper] - den$x[upper - 1])
		den$den[upper,,] <- w1 * den$den[upper - 1,,] + (1 - w1) * den$den[upper,,]
		den$x[upper] <- range[2]
	}

	den$x <- den$x[newIndex]
	den$den <- den$den[newIndex,,]

}




#Integrate on a 3D cube
integration3D <- function(den, fun = NULL, range = NULL) {

	newDen <- clone.environment(den)

	if (!is.null(range)) {
		if(!is.null(range$x)) 
			cutDen3D(newDen, "x", range$x)
		if(!is.null(range$y)) 
			cutDen3D(newDen, "y", range$y)
		if(!is.null(range$z)) 
			cutDen3D(newDen, "z", range$z)
	}
	
	xlen <- length(newDen$x)
	ylen <- length(newDen$y)
	zlen <- length(newDen$z)

	if (is.null(fun)) {
			val <- newDen$den
	} else {
		xseq <- outer(outer(newDen$x, rep(1, ylen)), rep(1, zlen))
		yseq <- outer(outer(rep(1, xlen), newDen$y), rep(1, zlen))
		zseq <- outer(outer(rep(1, xlen), rep(1, ylen)), newDen$z)
		val <- fun(xseq, yseq, zseq) * newDen$den
	}


	avg <- (val[-xlen, -ylen, -zlen] + val[-1, -ylen, -zlen] +
		val[-xlen, -1, -zlen] + val[-xlen, -ylen, -1] +val[-1, -1, -zlen] 
		+ val[-1, -ylen, -1] + val[-xlen, -1, -1] + val[-1, -1, -1] ) / 8
	xDiff <- newDen$x[-1] - newDen$x[-xlen]
	yDiff <- newDen$y[-1] - newDen$y[-ylen]
	zDiff <- newDen$z[-1] - newDen$z[-zlen]


	rm(newDen)
	sum(avg * outer(outer(xDiff, yDiff), zDiff))
}
