library(chutils)
library(dplyr)
library(smartGridSearch)
library(nlFit)

nlFitRun <- function(targets, estimates, upperBound, lowerBound, visibleReferencePoints, conceptualReferencePoints = NULL, freeMemoryLength = FALSE, freeRangeLength = FALSE, freeNumberSensitivity = TRUE, freepIncludeConceptualPoints = FALSE, loops = 500, targetOrder = "random", plotFileName = "nlFitplot.pdf", sinkFileName = "nlFitout.txt", verbose = FALSE, yLim = NULL, xLim = NULL) {
	
	data.in <- data.frame(targets, estimates)
	
	finalLoops <- loops

	test <- FALSE
	useMultiCore <- TRUE
	multicorePackages <- "nlFit"


	#number of intervals to split each bounded parameter into in order to get the random values between the hi and low values
	numGridIntervals <- 200
	#the number of times to loop with grid search before looking for the top 10 best fits.
	numGridLoops <- ifelse(test, 50, 250) #250
	#paramtable confirmation loops
	optBoundLoops <- 20
	#The number of best fit runs from which to resize the parameter bounds.
	#Here, I set it to the sqrt(total number of runs per grid loop) or 10 at minimum.
	  #optParamListN <- ifelse ( trunc(sqrt(numGridLoops)) < 10, 10, trunc(sqrt(numGridLoops)) )
	optParamListN <- 20
	#the number of simulations to average using the final parameters.
	numSimsToAverage <- 20

	numLoops <- numGridLoops
	numIntervals <- numGridIntervals
	optParamListN <- optParamListN
	optBoundLoops <- optBoundLoops

	dataTargetCol <- "targets"
	dataEstimateCol <- "estimates"
	defaultRangeLength <- 1
	defaultNumberSensitivity <- 0

		df.out <- NULL
		df.dat <- NULL

		numTargets <- length(unique(data.in$target))
		df.tmp <- data.in %>% group_by(targets) %>% summarize(estimates = mean(estimates, na.rm = T))
				
		### nlFit
		pUpper <- list(firstEstimate = 1, accuracyPercent = 1)
		pLower <- list(firstEstimate = 0, accuracyPercent = 0)
		pInt <- list(firstEstimate = 0.001, accuracyPercent = 0.001)
				
		otherParameterList <- list(data = df.tmp, upperBound = upperBound, lowerBound = lowerBound,visibleReferencePoints = visibleReferencePoints, conceptualReferencePoints = conceptualReferencePoints, loops = loops, targetOrder = targetOrder, dataTargetCol = "targets", dataEstimateCol = "estimates", minimizeStat = 'BIC', verbose = verbose)
				
		if(freeNumberSensitivity) {
			pUpper[["numberSensitivity"]] <- 1
			pLower[["numberSensitivity"]] <- 0
			pInt[["numberSensitivity"]] <- 0.001
		} else {
			otherParameterList[["numberSensitivity"]] <- defaultNumberSensitivity
		}
		if(freeRangeLength) {
			pUpper[["rangeLength"]] <- numTargets
			pLower[["rangeLength"]] <- 1
			pInt[["rangeLength"]] <- 1
		} else {
			otherParameterList[["rangeLength"]] <- defaultRangeLength
		}
		if(freeMemoryLength) {
			pUpper[["memoryLength"]] <- numTargets
			pLower[["memoryLength"]] <- 1
			pInt[["memoryLength"]] <- 1
		} else {
			otherParameterList[["memoryLength"]] <- numTargets
		}
		if(freepIncludeConceptualPoints) {
			pUpper[["pIncludeConceptualPoints"]] <- 1
			pLower[["pIncludeConceptualPoints"]] <- 0
			pInt[["pIncludeConceptualPoints"]] <- 0.001
		} else {
			otherParameterList[["pIncludeConceptualPoints"]] <- 0
		}
		pars.n <- length(pUpper)
		otherParameterList[["pars.n"]] <- pars.n
				
		x.grid <- smartGridSearch(getOrdinalNumberlineFit, pUpper, pLower, pInt, otherParameterList, numLoops = numLoops, numIntervals = numIntervals, optParamListN = optParamListN, optBoundLoops = optBoundLoops, multicore = useMultiCore, multicorePackages = multicorePackages)
	
		statList <- x.grid$final
		statList$upperBound <- upperBound
		statList$lowerBound <- lowerBound
		if(!freeMemoryLength) statList$memoryLength <- numTargets
		if(!freeRangeLength) statList$rangeLength <- defaultRangeLength
		if(!freepIncludeConceptualPoints) statList$pIncludeConceptualPoints <- 0
		if(!freeNumberSensitivity) statList$numberSensitivity <- defaultNumberSensitivity
		statList$targetOrder <- targetOrder
		statList$visibleReferencePoints <- visibleReferencePoints
		statList$conceptualReferencePoints <- conceptualReferencePoints
		statList$loops <- finalLoops
				
		outList <- outputNLfitStats(df.tmp, statList, dataTargetCol = "targets", dataEstimateCol = "estimates", pars.n = pars.n, loops = loops, sinkFilename =  sinkFileName, appendSinkFile = TRUE, plotFileName = plotFileName, xlim = xLim, ylim = yLim)
					

		df.out <- data.frame(firstEstimate = outList$runStats$parameters$firstEstimate, accuracyPercent = outList$runStats$parameters$accuracyPercent, numberSensitivity = outList$runStats$parameters$numberSensitivity, pIncludeConceptualPoints = outList$runStats$parameters$pIncludeConceptualPoints, rangeLength = outList$runStats$parameters$rangeLength, memoryLength = outList$runStats$parameters$memoryLength, r2 = outList$runStats$fitStats$r2, BIC = outList$runStats$fitStats$BIC)
		
		returnList = list(data = outList$df.fitted, fit = df.out)
		return(returnList)
}
