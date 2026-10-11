library(tidyr)
library(readxl)
library(dplyr)

###### functions
#count the responses a threshold flags, in one nlType or in all of them
countFlagged <- function(df.in, df.flagged, estimateCol, step, nlType, removed) {
	if(nlType != "all") {
		df.in <- df.in[df.in$nlType == nlType, ]
		df.flagged <- df.flagged[df.flagged$nlType == nlType, ]
	}
	nTotal <- sum(!is.na(df.in[[estimateCol]]))
	nFlagged <- sum(!is.na(df.flagged[[estimateCol]]))
	df.out <- data.frame(step = step, nlType = nlType, nFlagged = nFlagged, nTotal = nTotal, pctFlagged = round(100*nFlagged/nTotal, 2), removed = removed)
	return(df.out)
}

############# end

boundPxSize <- 350
bufferPxSize <- 75
madOutlierThresh <- 7.5
probTargetThresh <- 1.25

xLim <- c(-3.5, 3.5)
yLim <- xLim

mainDirPath <- getwd()
df.in <- read_excel(file.path(rawDataDir, filenameRawDataExl), na = c("", "NA"))

### fix inconsistent column naming between datasets
if ("AMAS_comp" %in% colnames(df.in)) {
  colnames(df.in)[colnames(df.in) == "AMAS_comp"] <- "AMAS_Comp"
}
if ("Sex" %in% colnames(df.in)) {
  colnames(df.in)[colnames(df.in) == "Sex"] <- "Bio_Sex"
}
if ("bUNL7_n" %in% colnames(df.in)) {
  colnames(df.in)[which(names(df.in) == "bUNL7_n")] <- "bUNL7_est"
}

#the composite must be the sum of the items it is built from
if(!isTRUE(all.equal(as.numeric(df.in$AMAS_Comp), as.numeric(rowSums(df.in[, amasItems]))))) {
	stop("AMAS_Comp does not equal the sum of the AMAS items")
}

df.demo <- df.in[, c("SN", demoCols)]
#the demographics are the same in both outlier sets, so save them once
write.table(df.demo, file=file.path(rawDataDir, demoFilename), quote=F, sep = "\t", row.names=F)

df.in <- df.in[, !names(df.in) %in% demoCols]

#reshape from wide to long: condition, trial and DV come from the column names
df.long.1 <- df.in %>%
  pivot_longer(
    cols = -c(SN),
    names_to = c("condition", "trial", "DV"),
    names_pattern = "([A-Za-z]+)(\\d+)_(.*)",
    values_to = "value"
  )

#spread the DVs into their own columns (ub, lb, targ, est)
df.long <- df.long.1 %>%
  pivot_wider(
    names_from = DV,
    values_from = value
  )

	### recode conditions so it is bounded vs unbounded
	df.long$nlType <- ifelse(df.long$condition == "TNL", "bounded", "unbounded")
	df.long$isUniversalBounded <- ifelse(df.long$condition == "bUNL", TRUE, FALSE)

#### put all number lines on the same scale
df.long$halfBoundRange <- (df.long$ub - df.long$lb)/2
df.long$maxBound <- df.long$ub + 2*df.long$halfBoundRange
df.long$minBound <- df.long$lb - 2*df.long$halfBoundRange
df.long$unitsPerPixel <- (2*df.long$halfBoundRange)/boundPxSize
df.long$bufferUnits <- df.long$unitsPerPixel * bufferPxSize

df.long$s_ub <- (df.long$ub - df.long$halfBoundRange - df.long$lb)/df.long$halfBoundRange
df.long$s_lb <- (df.long$lb - df.long$halfBoundRange - df.long$lb)/df.long$halfBoundRange
df.long$s_targ <- (df.long$targ - df.long$halfBoundRange - df.long$lb)/df.long$halfBoundRange
df.long$s_est <- (df.long$est - df.long$halfBoundRange - df.long$lb)/df.long$halfBoundRange
#the buffer adds the 75 pixels on either side of the target region
df.long$s_maxBound <- (df.long$maxBound - df.long$halfBoundRange - df.long$lb + df.long$bufferUnits)/df.long$halfBoundRange
df.long$s_minBound <- (df.long$minBound - df.long$halfBoundRange - df.long$lb - df.long$bufferUnits)/df.long$halfBoundRange
###### end same scale

#the bounded number line allows responses only between its labels
df.long$s_maxBound <- ifelse(df.long$nlType == "bounded", df.long$s_ub, df.long$s_maxBound)
df.long$s_minBound <- ifelse(df.long$nlType == "bounded", df.long$s_lb, df.long$s_minBound)

#### filter outliers
#first identify and remove problem targets
df.long <- df.long %>% group_by(nlType, SN) %>% mutate(z_diff = scale(s_targ - s_est))
df.zDiff <- df.long %>% group_by(targ) %>% summarize(meanZdiff = mean(z_diff, na.rm = T))
problemTargets <- df.zDiff[abs(df.zDiff$meanZdiff) > probTargetThresh, ]
df.long.problemTargets <- df.long[df.long$targ %in% problemTargets$targ, ]
if(removeOutliers) {
	df.long.rmProblemTargets <- df.long[!(df.long$targ %in% problemTargets$targ), ]
} else {
	df.long.rmProblemTargets <- df.long
}
df.long.rmProblemTargets$z_diff <- NULL

#now identify and remove individual outliers
df.long.rmProblemTargets$diff <- df.long.rmProblemTargets$s_targ - df.long.rmProblemTargets$s_est
#scale each difference by the participant's median absolute difference
df.long.rmProblemTargets <- df.long.rmProblemTargets %>% group_by(nlType, SN) %>% mutate(median_diff = median(abs(diff), na.rm = T), mad_diff = diff/median_diff)
#filter by a threshold
df.long.IndividualOutliers <- df.long.rmProblemTargets[abs(df.long.rmProblemTargets$mad_diff) >= madOutlierThresh, ]
if(removeOutliers) {
	df.long.rmProblemTargets.filtered <- df.long.rmProblemTargets[abs(df.long.rmProblemTargets$mad_diff) < madOutlierThresh, ]
} else {
	df.long.rmProblemTargets.filtered <- df.long.rmProblemTargets
}

#tally what each threshold flags, whether or not it is removed
df.exclusions <- NULL
for(ty in c(unique(df.long$nlType), "all")) {
	df.exclusions <- rbind(df.exclusions, countFlagged(df.long, df.long.problemTargets, "s_est", "problemTargets", ty, removeOutliers))
	df.exclusions <- rbind(df.exclusions, countFlagged(df.long.rmProblemTargets, df.long.IndividualOutliers, "s_est", "individualOutliers", ty, removeOutliers))
}

pdf("plotsDescribingData.pdf")
pTargs <- unique(df.long$targ)
for(pt in pTargs) {
	df.tmp.pt <- df.long[df.long$targ == pt,]
	types <- unique(df.tmp.pt$nlType)
	for(ty in types) {
		df.tmp.pt.1 <- df.tmp.pt[df.tmp.pt$nlType == ty, ]
		sTargs <- unlist(unique(df.tmp.pt.1[, "s_targ"]))
		for(st in sTargs) {
			df.tmp.pt.2 <- df.tmp.pt.1[df.tmp.pt.1$s_targ == st, ]
			with(df.tmp.pt.2, hist(s_est, breaks = 30, main = paste(ty, format(pt, scientific = FALSE, big.mark = ",")), , xlim = c(-3.5, 3.5), cex.axis = 1.5,lwd = 2))
			abline(v=st)
		}
	}
}

subjects <- unique(df.long.rmProblemTargets.filtered$SN)
types <- unique(df.long.rmProblemTargets.filtered$nlType)
for(sn in subjects) {
	for(ty in types) {
		df.tmp.1 <- df.long.rmProblemTargets.filtered[df.long.rmProblemTargets.filtered$SN == sn & df.long.rmProblemTargets.filtered$nlType == ty, ]
		with(df.tmp.1, plot(s_est ~ s_targ, main = sn, xlim = xLim, ylim = yLim))
		abline(a=0, b=1)
	}
}
dev.off()

filename <- paste("samClarissaProcessedData_", fileTag, ".txt", sep = "")
write.table(df.long.rmProblemTargets.filtered, file=filename, quote=F, sep = "\t", row.names=F)
filename <- paste("samClarissaProblemTargets_", fileTag, ".txt", sep = "")
write.table(problemTargets, file=filename, quote=F, sep = "\t", row.names=F)
filename <- paste("samClarissaIndividualOutliers_", fileTag, ".txt", sep = "")
write.table(df.long.IndividualOutliers, file=filename, quote=F, sep = "\t", row.names=F)
filename <- paste("samClarissaExclusions_", fileTag, ".txt", sep = "")
write.table(df.exclusions, file=filename, quote=F, sep = "\t", row.names=F)
