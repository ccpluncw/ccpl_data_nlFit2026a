library(chutils)
library(tidyr)
library(dplyr)

skipModelsWithExistingDirs <- FALSE

filenamexl <- paste("samClarissaProcessedData_", fileTag, ".txt", sep="")
df.long <- read.table(filenamexl, header = T, sep="\t")

mainDirPath <- getwd()

##### score each subject and the group
subjects <- unique(df.long$SN)
types <- unique(df.long$nlType)
subjects <- c(subjects, "all")

keepCols <- c("nlType", "s_targ", "s_est", "s_ub", "s_lb", "s_maxBound", "s_minBound", "isUniversalBounded")
all.tmp <- NULL

df.long.1 <- df.long
subDirExists <- ch.newDir (mainDirPath, "traditional")
tradPath <- getwd()
outFilename <- "traditionalFitData.txt"

if(subDirExists == FALSE | skipModelsWithExistingDirs == FALSE) {
	for(sub in subjects) {
		if(sub == "all") {
			df.tmp.1 <- df.long.1 %>% group_by(nlType, s_targ) %>% summarize(s_est = mean(s_est, na.rm = T), s_maxBound= mean(s_maxBound, na.rm = T), s_minBound= mean(s_minBound, na.rm = T), s_lb= mean(s_lb, na.rm = T),s_ub= mean(s_ub, na.rm = T), isUniversalBounded= mean(isUniversalBounded, na.rm = T))
		} else {
			df.tmp.1 <- df.long.1[df.long.1$SN == sub, keepCols]
		}

		for(typ in types) {
			df.tmp <- df.tmp.1[df.tmp.1$nlType == typ, ]

			#### calculate PAE
			 df.tmp$error <- df.tmp$s_est - df.tmp$s_targ
			 #scale the error by the width of the visible number line
			 PAEdivisor <- df.tmp$s_ub[1] - df.tmp$s_lb[1]
			 PAE <-  mean(abs(df.tmp$error), na.rm = T)/PAEdivisor

			 df.tmp.fit <- data.frame(sn = sub, nlType = typ, PAE = PAE)

			all.tmp <- ch.rbind(all.tmp, df.tmp.fit)
		}
	}
	write.table(all.tmp, file=outFilename, quote=F, sep = "\t", row.names=F)
	setwd(mainDirPath)
}
