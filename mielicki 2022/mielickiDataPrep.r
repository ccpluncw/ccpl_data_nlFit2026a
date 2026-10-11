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

#score forced-choice items against a key, a missing response stays missing
scoreKeyed <- function(df.in, key, prefix) {
	df.out <- data.frame(SN = as.numeric(df.in$SN))
	for(nm in names(key)) {
		resp <- suppressWarnings(as.numeric(df.in[[nm]]))
		df.out[[paste(prefix, nm, sep="")]] <- ifelse(is.na(resp), NA, as.numeric(resp == key[[nm]]))
	}
	return(df.out)
}

#proportion correct, NA when too few items were answered
propCorrect <- function(df.items, minAnswered) {
	nAnswered <- rowSums(!is.na(df.items))
	pc <- rowMeans(df.items, na.rm = TRUE)
	pc[nAnswered < minAnswered] <- NA
	return(as.numeric(pc))
}

#MathCourse lists course codes 1-10, code 11 is other and is not a course
countMathCourses <- function(x) {
	#the raw column is text, so a blank answer arrives as the string "NA"
	if(is.na(x) || x == "NA") return(NA_integer_)
	codes <- suppressWarnings(as.numeric(unlist(strsplit(as.character(x), ","))))
	codes <- codes[!is.na(codes) & codes >= 1 & codes <= 10]
	return(length(unique(codes)))
}

############# end

madOutlierThresh <- 7.5
probTargetThresh <- 1.25

nlCols <- c("Lg_WNL_1", "Lg_WNL_2", "Lg_WNL_3", "Lg_WNL_4", "Lg_WNL_5", "Lg_WNL_6", "Lg_WNL_7", "Lg_WNL_8", "Lg_WNL_9", "Lg_WNL_10", "Lg_WNL_11", "Lg_WNL_12", "Lg_WNL_13", "Lg_WNL_14", "Lg_WNL_15",  "Lg_WNL_16", "Lg_WNL_17", "Lg_WNL_18", "Lg_WNL_19", "Lg_WNL_20", "Lg_WNL_21", "Lg_WNL_22")

mainDirPath <- getwd()
df.in <- read_excel(file.path(rawDataDir, filenameRawDataExl))

### fix inconsistent column naming between datasets
for(nm in names(anxItems)) {
  colnames(df.in)[colnames(df.in) == anxItems[nm]] <- nm
  df.in[[nm]] <- as.numeric(ifelse(df.in[[nm]] == "NA", NA, df.in[[nm]]))
}
df.in$rate_anx <- rowMeans(df.in[, names(anxItems)])
if ("AgeSurvey" %in% colnames(df.in)) {
  colnames(df.in)[colnames(df.in) == "AgeSurvey"] <- "Age"
}
if ("id" %in% colnames(df.in)) {
  colnames(df.in)[colnames(df.in) == "id"] <- "SN"
}
df.demo <- df.in[, demoCols]

#### math covariates
#scored and standardized before the exclusions, so the scale does not move
df.equiv.items <- scoreKeyed(df.in, equivKey, "eq_")
df.torr.items <- scoreKeyed(df.in, torrKey, "tr_")
equivCols <- names(df.equiv.items)[names(df.equiv.items) != "SN"]
torrCols <- names(df.torr.items)[names(df.torr.items) != "SN"]
#the scored frames keep the row order of df.in, so they line up with df.demo
df.demo <- data.frame(df.demo, df.equiv.items[, equivCols], df.torr.items[, torrCols])
df.demo$equiv <- propCorrect(df.demo[, equivCols], equivMinAnswered)
df.demo$torr <- propCorrect(df.demo[, torrCols], torrMinAnswered)
df.demo$nMath <- as.numeric(sapply(df.demo$MathCourse, countMathCourses))
df.demo$equivZ <- as.numeric(scale(df.demo$equiv))
df.demo$torrZ <- as.numeric(scale(df.demo$torr))
df.demo$nMathZ <- as.numeric(scale(df.demo$nMath))
#anxiety is standardized before the exclusions too
df.demo$rate_anxZ <- as.numeric(scale(df.demo$rate_anx))

#the demographics are the same in both outlier sets, so save them once
write.table(df.demo, file=file.path(rawDataDir, demoFilename), quote=F, sep = "\t", row.names=F)

df.in <- df.in[, names(df.in) %in% c("SN", nlCols)]

df.in <- df.in %>%
  rename(
		target_150000000 = Lg_WNL_1,
		target_722000000 = Lg_WNL_2,
		target_2000000 = Lg_WNL_3,
		target_56000000 = Lg_WNL_4,
		target_606000000 = Lg_WNL_5,
		target_738000000 = Lg_WNL_6,
		target_147000000 = Lg_WNL_7,
		target_163000000 = Lg_WNL_8,
		target_818000000 = Lg_WNL_9,
		target_100000000 = Lg_WNL_10,
		target_725000000 = Lg_WNL_11,
		target_486000000 = Lg_WNL_12,
		target_122000000 = Lg_WNL_13,
		target_18000000 = Lg_WNL_14,
		target_366000000 = Lg_WNL_15,
		target_34000000 = Lg_WNL_16,
		target_5000000 = Lg_WNL_17,
		target_246000000 = Lg_WNL_18,
		target_179000000 = Lg_WNL_19,
		target_754000000 = Lg_WNL_20,
		target_938000000 = Lg_WNL_21,
		target_78000000 = Lg_WNL_22
	)

#reshape from wide to long, the column name carries the target
df.long <- df.in %>%
  pivot_longer(
		cols = starts_with("target_"),
		names_to = "target",
		names_prefix = "target_",
		values_to = "estimate"
	)

	### every mielicki number line is bounded from 0 to 1 billion
	df.long$nlType <- "bounded"
	df.long$isUniversalBounded <- FALSE
	df.long$s_ub <- 1000000000
	df.long$s_lb <- 0
	df.long$s_maxBound <- 1000000000
	df.long$s_minBound <- 0
	
	df.long$target <- ifelse(df.long$target == "NA", NA, as.numeric(as.character(df.long$target)))
	df.long$estimate <- ifelse(df.long$estimate == "NA", NA, as.numeric(as.character(df.long$estimate)))


#### filter outliers
#first identify and remove problem targets
df.long <- df.long %>% group_by(nlType, SN) %>% mutate(z_diff = scale(target - estimate))
df.zDiff <- df.long %>% group_by(target) %>% summarize(meanZdiff = mean(z_diff, na.rm = T))
problemTargets <- df.zDiff[abs(df.zDiff$meanZdiff) > probTargetThresh, ]
df.long.problemTargets <- df.long[df.long$target %in% problemTargets$target, ]
if(removeOutliers) {
	df.long.rmProblemTargets <- df.long[!(df.long$target %in% problemTargets$target), ]
} else {
	df.long.rmProblemTargets <- df.long
}
df.long.rmProblemTargets$z_diff <- NULL

#now identify and remove individual outliers
df.long.rmProblemTargets$diff <- df.long.rmProblemTargets$target - df.long.rmProblemTargets$estimate
#scale each difference by the participant's median absolute difference
df.long.rmProblemTargets <- df.long.rmProblemTargets %>% group_by(nlType, SN) %>% mutate(median_diff = median(abs(diff), na.rm = T), mad_diff = diff/median_diff)
#filter by a threshold
df.long.IndividualOutliers <- df.long.rmProblemTargets[abs(df.long.rmProblemTargets$mad_diff) >= madOutlierThresh, ]
df.long.IndividualOutliers <- df.long.IndividualOutliers[!is.na(df.long.IndividualOutliers$SN), ]
if(removeOutliers) {
	df.long.rmProblemTargets.filtered <- df.long.rmProblemTargets[abs(df.long.rmProblemTargets$mad_diff) < madOutlierThresh, ]
} else {
	df.long.rmProblemTargets.filtered <- df.long.rmProblemTargets
}

#tally what each threshold flags, whether or not it is removed
df.exclusions <- NULL
for(ty in c(unique(df.long$nlType), "all")) {
	df.exclusions <- rbind(df.exclusions, countFlagged(df.long, df.long.problemTargets, "estimate", "problemTargets", ty, removeOutliers))
	df.exclusions <- rbind(df.exclusions, countFlagged(df.long.rmProblemTargets, df.long.IndividualOutliers, "estimate", "individualOutliers", ty, removeOutliers))
}

xLim = c(0, 1000000000)
yLim = c(0, 1000000000)

pdf("plotsDescribingData.pdf")
pTargs <- unique(df.long$target)
for(pt in pTargs) {
	df.tmp.pt <- df.long[df.long$target == pt,]
	types <- unique(df.tmp.pt$nlType)
	for(ty in types) {
		df.tmp.pt.1 <- df.tmp.pt[df.tmp.pt$nlType == ty, ]
		sTargs <- unlist(unique(df.tmp.pt.1[, "target"]))
		for(st in sTargs) {
			df.tmp.pt.2 <- df.tmp.pt.1[df.tmp.pt.1$target == st, ]
			with(df.tmp.pt.2, hist(estimate, breaks = 30, main = paste(ty, format(pt, scientific = FALSE, big.mark = ",")), , xlim = xLim, cex.axis = 1.5,lwd = 2))
			abline(v=st)
		}
	}
}

subjects <- unique(df.long.rmProblemTargets.filtered$SN)
types <- unique(df.long.rmProblemTargets.filtered$nlType)
types <- types[!is.na(types)]
for(sn in subjects) {
	for(ty in types) {
		df.tmp.1 <- df.long.rmProblemTargets.filtered[df.long.rmProblemTargets.filtered$SN == sn & df.long.rmProblemTargets.filtered$nlType == ty, ]
		with(df.tmp.1, plot(estimate ~ target, main = sn, xlim = xLim, ylim = yLim))
		abline(a=0, b=1)
	}
}
dev.off()

filename <- paste("mielickiProcessedData_", fileTag, ".txt", sep = "")
write.table(df.long.rmProblemTargets.filtered, file=filename, quote=F, sep = "\t", row.names=F)
filename <- paste("mielickiProblemTargets_", fileTag, ".txt", sep = "")
write.table(problemTargets, file=filename, quote=F, sep = "\t", row.names=F)
filename <- paste("mielickiIndividualOutliers_", fileTag, ".txt", sep = "")
write.table(df.long.IndividualOutliers, file=filename, quote=F, sep = "\t", row.names=F)
filename <- paste("mielickiExclusions_", fileTag, ".txt", sep = "")
write.table(df.exclusions, file=filename, quote=F, sep = "\t", row.names=F)

