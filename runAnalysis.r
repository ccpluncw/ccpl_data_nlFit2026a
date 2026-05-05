library(chutils)

datasets <- c("kent", "prolific")
#datasets <- c("kent")
runDataPrep <- FALSE
runNlFit <- FALSE
runTradAnalysis <- FALSE
runParamAnalysis <- TRUE

mainDir.A <- getwd()
for(dt in datasets) {
	if(dt == "kent") {
		ch.newDir (mainDir.A, "Exp 1 Kent")
		dataDir <- getwd()
		filenameRawDataExl <- "SSZpilotdata_for_dale_Jan27.xlsx"
		demoCols <- c("AMAS_Comp", "Bio_Sex", "Gender", "Age", "Race")
		demoFilename <- "SSZkentDemographics.txt"
		fileTag <- "Kent"
	} else {
		ch.newDir (mainDir.A, "Exp 2 Polific")
		dataDir <- getwd()
		filenameRawDataExl <- "SSZ_SM1data_forD&K_Oct9.xlsx"
		demoCols <- c("AMAS_Comp", "Bio_Sex", "Gender", "Race", "Race_13_TEXT", "Age", "EdLev", "HouseInc",	"Diag_Stat", "Diag_Age",	"pLastCasto", "HCs1Casto", "HCs1Casto_3_TEXT", "HCs2Casto", "HCs2Casto_7_TEXT", "HCb1Casto", "HCb2Casto")
		demoFilename <- "SSZprolificDemographics.txt"
		fileTag <- "Prolific"
	
	}
	if(runDataPrep) {
		source("../samClarissaDataPrep.r")
		setwd(dataDir)
	}
	if(runNlFit) {
		source("../samClarissaNlFitAnalysis.r")
		setwd(dataDir)
	}
	if(runTradAnalysis) {
		source("../samClarissaTraditionalAnalysis.r")
		setwd(dataDir)
	}
	if(runParamAnalysis) {
		source("../samClarissaParameterAnalysis.r")
		setwd(mainDir.A)
	}
}