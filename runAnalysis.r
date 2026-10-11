library(chutils)

datasets <- c("kent", "prolific", "mielicki")
#datasets <- c("kent")
#run each dataset with and without the outlier filters
outlierSets <- c("withoutOutliers", "withOutliers")
runDataPrep <- TRUE
runNlFit <- TRUE
runTradAnalysis <- TRUE
runParamAnalysis <- TRUE

mainDir.A <- getwd()
nlFitRunScript <- file.path(mainDir.A, "nlFitRun.r")
#fits are stored in folders named for the nlFit version that produced them
nlFitDirName <- paste0("nlFit_v", packageVersion("nlFit"))
nlFitGroupDirName <- paste0("nlFitGroupData_v", packageVersion("nlFit"))
#the nine AMAS items that sum to AMAS_Comp
amasItems <- paste("AMAS_", 1:9, sep="")
for(dt in datasets) {
	if(dt == "kent") {
		ch.newDir (mainDir.A, "Exp 1 Kent")
		rawDataDir <- getwd()
		filenameRawDataExl <- "SSZpilotdata_for_dale_Sep17.xlsx"
		demoCols <- c("AMAS_Comp", "Bio_Sex", "Gender", "Age", "Race", amasItems)
		demoFilename <- "SSZkentDemographics.txt"
		fileTag <- "Kent"
		#kent and prolific share the scripts in the root folder
		dataPrepScript <- file.path(mainDir.A, "samClarissaDataPrep.r")
		nlFitScript <- file.path(mainDir.A, "samClarissaNlFitAnalysis.r")
		tradScript <- file.path(mainDir.A, "samClarissaTraditionalAnalysis.r")
		paramScript <- file.path(mainDir.A, "samClarissaParameterAnalysis.r")
	} else if(dt == "prolific") {
		ch.newDir (mainDir.A, "Exp 2 Prolific")
		rawDataDir <- getwd()
		filenameRawDataExl <- "SSZ_SM1data_fordale_Sep17.xlsx"
		demoCols <- c("AMAS_Comp", "Bio_Sex", "Gender", "Race", "Age", "EdLev", amasItems)
		demoFilename <- "SSZprolificDemographics.txt"
		fileTag <- "Prolific"
		dataPrepScript <- file.path(mainDir.A, "samClarissaDataPrep.r")
		nlFitScript <- file.path(mainDir.A, "samClarissaNlFitAnalysis.r")
		tradScript <- file.path(mainDir.A, "samClarissaTraditionalAnalysis.r")
		paramScript <- file.path(mainDir.A, "samClarissaParameterAnalysis.r")
	} else if(dt == "mielicki") {
		ch.newDir (mainDir.A, "mielicki 2022")
		rawDataDir <- getwd()
		filenameRawDataExl <- "OSF/VM data FINAL.xlsx"
		#rate_anx is the mean of these five math anxiety items
		anxItems <- c(ma_gen = "Q155_1", ma_wn = "Q157_1", ma_per = "Q161_1", ma_freq = "Q361_1", ma_frac = "Q159_1")
		#fraction equivalence key: 1 = true, 2 = false
		equivKey <- list(Q386 = 1, Q387 = 1, Q388 = 1, Q389 = 2, Q390 = 2, Q391 = 2)
		#relational reasoning key: 1-4 = a-d
		torrKey <- list(Q466 = 1, Q467 = 4, Q468 = 4, Q469 = 3)
		equivMinAnswered <- 4
		torrMinAnswered <- 3
		#the mielicki raw frame is subset by demoCols alone, so SN is listed
		demoCols <- c("SN", "Age", "Gender", "MathCourse", "Race", names(anxItems), "rate_anx")
		demoFilename <- "mielickiDemo.txt"
		fileTag <- "Mielicki"
		#the mielicki scripts live in their own folder
		dataPrepScript <- file.path(rawDataDir, "mielickiDataPrep.r")
		nlFitScript <- file.path(rawDataDir, "mielickiNlFitAnalysis.r")
		tradScript <- file.path(rawDataDir, "mielickiTraditionalAnalysis.r")
		paramScript <- file.path(rawDataDir, "mielickiParameterAnalysis.r")
	} else {
		stop(paste("unknown dataset:", dt))
	}
	for(os in outlierSets) {
		removeOutliers <- (os == "withoutOutliers")
		ch.newDir (rawDataDir, os)
		dataDir <- getwd()
		if(runDataPrep) {
			source(dataPrepScript)
			setwd(dataDir)
		}
		if(runNlFit) {
			source(nlFitScript)
			setwd(dataDir)
		}
		if(runTradAnalysis) {
			source(tradScript)
			setwd(dataDir)
		}
		if(runParamAnalysis) {
			source(paramScript)
			setwd(dataDir)
		}
		setwd(rawDataDir)
	}
	setwd(mainDir.A)
}
