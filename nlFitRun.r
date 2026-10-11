library(chutils)
library(dplyr)
library(smartGridSearch2)
library(nlFit)

nlFitRun <- function(targets, estimates, upperBound, lowerBound, visibleReferencePoints, conceptualReferencePoints = NULL, freeMemoryLength = FALSE, freeNumberSensitivity = TRUE, freepConceptual = FALSE, loops = 500, targetOrder = "random", plotFileName = "nlFitplot.pdf", sinkFileName = "nlFitout.txt", verbose = FALSE, yLim = NULL, xLim = NULL) {
	
	data.in <- data.frame(targets, estimates)
	
	finalLoops <- loops

	test <- FALSE
	useMultiCore <- TRUE
	multicorePackages <- "nlFit"


	#intervals each bounded parameter is split into for the random draws
	numGridIntervals <- 200
	#grid search loops before the best fits are kept
	numGridLoops <- ifelse(test, 50, 250)
	#confirmation loops on the parameter bounds
	optBoundLoops <- 20
	#best fit runs used to resize the parameter bounds
	optParamListN <- 20

	numLoops <- numGridLoops
	numIntervals <- numGridIntervals
	optParamListN <- optParamListN
	optBoundLoops <- optBoundLoops

	dataTargetCol <- "targets"
	dataEstimateCol <- "estimates"
	defaultNumberSensitivity <- 1

		df.out <- NULL
		df.dat <- NULL

		numTargets <- length(unique(data.in$target))
		df.tmp <- data.in %>% group_by(targets) %>% summarize(estimates = mean(estimates, na.rm = T))
				
		#free parameters and their ranges
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
		if(freeMemoryLength) {
			pUpper[["memoryLength"]] <- numTargets
			pLower[["memoryLength"]] <- 1
			pInt[["memoryLength"]] <- 1
		} else {
			otherParameterList[["memoryLength"]] <- numTargets
		}
		if(freepConceptual) {
			pUpper[["pConceptual"]] <- 1
			pLower[["pConceptual"]] <- 0
			pInt[["pConceptual"]] <- 0.001
		} else {
			otherParameterList[["pConceptual"]] <- 0
		}
		pars.n <- length(pUpper)
		otherParameterList[["pars.n"]] <- pars.n
				
		x.grid <- smartGridSearch(getOrdinalNumberlineFit, pUpper, pLower, pInt, otherParameterList, numLoops = numLoops, numIntervals = numIntervals, optParamListN = optParamListN, optBoundLoops = optBoundLoops, multicore = useMultiCore, multicorePackages = multicorePackages)
	
		statList <- x.grid$final
		statList$upperBound <- upperBound
		statList$lowerBound <- lowerBound
		if(!freeMemoryLength) statList$memoryLength <- numTargets
		if(!freepConceptual) statList$pConceptual <- 0
		if(!freeNumberSensitivity) statList$numberSensitivity <- defaultNumberSensitivity
		statList$targetOrder <- targetOrder
		statList$visibleReferencePoints <- visibleReferencePoints
		statList$conceptualReferencePoints <- conceptualReferencePoints
		statList$loops <- finalLoops
				
		outList <- outputNLfitStats(df.tmp, statList, dataTargetCol = "targets", dataEstimateCol = "estimates", pars.n = pars.n, loops = loops, sinkFilename =  sinkFileName, appendSinkFile = TRUE, plotFileName = plotFileName, xlim = xLim, ylim = yLim)
					

		df.out <- data.frame(firstEstimate = outList$runStats$parameters$firstEstimate, accuracyPercent = outList$runStats$parameters$accuracyPercent, numberSensitivity = outList$runStats$parameters$numberSensitivity, pConceptual = outList$runStats$parameters$pConceptual, memoryLength = outList$runStats$parameters$memoryLength, r2 = outList$runStats$fitStats$r2, BIC = outList$runStats$fitStats$BIC)
		
		returnList = list(data = outList$df.fitted, fit = df.out)
		return(returnList)
}
