library(chutils)
library(tidyr)
library(dplyr)

skipModelsWithExistingDirs <- FALSE

filenamexl <- paste("mielickiProcessedData_", fileTag, ".txt", sep="")
df.long <- read.table(filenamexl, header = T, sep="\t")

mainDirPath <- getwd()

df.long.1 <- df.long[complete.cases(df.long), ]

##### score each subject and the group
#a participant with no usable response has nothing to score
subjects <- unique(df.long.1$SN)
types <- unique(df.long.1$nlType)
subjects <- c(subjects, "all")

keepCols <- c("nlType", "target", "estimate", "s_ub", "s_lb", "s_maxBound", "s_minBound")
all.tmp <- NULL

subDirExists <- ch.newDir (mainDirPath, "traditional")
tradPath <- getwd()
outFilename <- "traditionalFitData.txt"

if(subDirExists == FALSE | skipModelsWithExistingDirs == FALSE) {
	for(sub in subjects) {
		if(sub == "all") {
			df.tmp.1 <- df.long.1 %>% group_by(nlType, target) %>% summarize(estimate = mean(estimate, na.rm = T), s_maxBound= mean(s_maxBound, na.rm = T), s_minBound= mean(s_minBound, na.rm = T), s_lb= mean(s_lb, na.rm = T),s_ub= mean(s_ub, na.rm = T))
		} else {
			df.tmp.1 <- df.long.1[df.long.1$SN == sub, keepCols]
		}

		for(typ in types) {
			df.tmp <- df.tmp.1[df.tmp.1$nlType == typ, ]

			#### calculate PAE
			 df.tmp$error <- df.tmp$estimate - df.tmp$target
			 #scale the error by the width of the number line
			 PAEdivisor <- df.tmp$s_maxBound[1] - df.tmp$s_minBound[1]
			 PAE <-  mean(abs(df.tmp$error), na.rm = T)/PAEdivisor

			 df.tmp.fit <- data.frame(sn = sub, nlType = typ, PAE = PAE)

			all.tmp <- ch.rbind(all.tmp, df.tmp.fit)
		}
	}
	write.table(all.tmp, file=outFilename, quote=F, sep = "\t", row.names=F)
	setwd(mainDirPath)
}
