library(nlFit)
library(chutils)
library(tidyr)
source("../nlFitRun.r")

nlFitLoops <- 250
targetOrder <- "random"
#Only run new models that have not been run before?
skipModelsWithExistingDirs <- TRUE
freepIncludeConceptualPoints <- FALSE
conceptualReferencePoints <- 0


filenamexl <- paste("samClarissaProcessedData_", fileTag, ".txt", sep="")
df.long <- read.table(filenamexl, header = T, sep="\t")

mainDirPath <- getwd()

##### Run analysis on each subject and as a group
subjects <- unique(df.long$SN)
types <- unique(df.long$nlType)
types <- c(types, "unboundedBetweenLabels")
subjects <- c(subjects, "all")

yLim <- c(min(df.long$s_minBound), max(df.long$s_maxBound))
xLim <- yLim

keepCols <- c("nlType", "s_targ", "s_est", "s_ub", "s_lb", "s_maxBound", "s_minBound", "isUniversalBounded")

	all.tmp <- NULL
	data.tmp <- NULL
	
	df.long.1 <- df.long
	subDirExists <- ch.newDir (mainDirPath, "nlFit")
	ffMilPath <- getwd()
	outFilename <- "nlGenFitData.txt"
	dataOutFilename <- "nlGenPlotData.txt"

	if(subDirExists == FALSE | skipModelsWithExistingDirs == FALSE) {
		ch.newDir (ffMilPath, "nlFitOut")
		nlFitPath <- getwd()
		setwd(mainDirPath)
		ch.newDir (ffMilPath, "nlFitPlots")
		nlFitPlotPath <- getwd()
		setwd(ffMilPath)

		for(sub in subjects) {
			if(sub == "all") {
				df.tmp.1 <- df.long.1 %>% group_by(nlType, s_targ) %>% summarize(s_est = mean(s_est, na.rm = T), s_maxBound= mean(s_maxBound, na.rm = T), s_minBound= mean(s_minBound, na.rm = T), s_lb= mean(s_lb, na.rm = T),s_ub= mean(s_ub, na.rm = T), isUniversalBounded= mean(isUniversalBounded, na.rm = T))
			} else {
				df.tmp.1 <- df.long.1[df.long.1$SN == sub, keepCols]
			}
	
			for(typ in types) {
				if(typ == "unboundedBetweenLabels") {
					df.tmp <- df.tmp.1[df.tmp.1$isUniversalBounded == TRUE, ]
				} else {
					df.tmp <- df.tmp.1[df.tmp.1$nlType == typ, ]		
				}
		 
				plotFileName.1 <- paste(sub, typ, "nlFitplot.pdf", sep="-")
				plotFileName.2 <- file.path(nlFitPlotPath, plotFileName.1)
				sinkFileName.1 <- paste(sub, typ, "nlFitout.txt", sep="-")
				sinkFileName.2 <- file.path(nlFitPath, sinkFileName.1)
				outlist.r <- nlFitRun(df.tmp$s_targ, df.tmp$s_est, upperBound = df.tmp$s_maxBound[1], lowerBound = df.tmp$s_minBound[1], visibleReferencePoints = c(df.tmp$s_lb[1], df.tmp$s_ub[1]), conceptualReferencePoints = conceptualReferencePoints, freepIncludeConceptualPoints = freepIncludeConceptualPoints, targetOrder = targetOrder, plotFileName = plotFileName.2, sinkFileName = sinkFileName.2, loops = nlFitLoops, yLim = yLim, xLim = xLim)

				df.tmp.fit <- data.frame(sn = sub, nlType = typ, nl_firstEstimate_Random =  outlist.r$fit$firstEstimate, nl_accuracyPercent_Random =  outlist.r$fit$accuracyPercent, nl_numberSensitivity_Random =  outlist.r$fit$numberSensitivity, nl_pIncludeConceptualPoints_Random =  outlist.r$fit$pIncludeConceptualPoints, nl_rangeLength_Random =  outlist.r$fit$rangeLength, nl_memoryLength_Random =  outlist.r$fit$memoryLength, nl_r2_Random =  outlist.r$fit$r2,  nl_BIC_Random =  outlist.r$fit$BIC)

				all.tmp <- ch.rbind(all.tmp, df.tmp.fit)
				
				tmpNlfitDatafit <- outlist.r$data
				tmpNlfitDatafit$sn <- sub
				tmpNlfitDatafit$nlType <- typ
				data.tmp <- ch.rbind(data.tmp, tmpNlfitDatafit)
			}
		}
		write.table(all.tmp, file=outFilename, quote=F, sep = "\t", row.names=F)
		write.table(data.tmp, file=dataOutFilename, quote=F, sep = "\t", row.names=F)
		setwd(mainDirPath)
	}

