library(tidyr)
library(readxl)
library(dplyr)

zOutlierThresh <- 3
madOutlierThresh <- 7.5
probTargetThresh <- 1.25
fileTag <- "Kent"

demoFilename <- "mielickiDemo.txt"
filenameRawDataExl <- "./OSF/VM data FINAL.xlsx"
nlCols <- c("Lg_WNL_1", "Lg_WNL_2", "Lg_WNL_3", "Lg_WNL_4", "Lg_WNL_5", "Lg_WNL_6", "Lg_WNL_7", "Lg_WNL_8", "Lg_WNL_9", "Lg_WNL_10", "Lg_WNL_11", "Lg_WNL_12", "Lg_WNL_13", "Lg_WNL_14", "Lg_WNL_15",  "Lg_WNL_16", "Lg_WNL_17", "Lg_WNL_18", "Lg_WNL_19", "Lg_WNL_20", "Lg_WNL_21", "Lg_WNL_22")

mainDirPath <- getwd()
df.in <- read_excel(filenameRawDataExl)

### fix inconsistent column naming between datasets
if ("Q155_1" %in% colnames(df.in)) {
  colnames(df.in)[colnames(df.in) == "Q155_1"] <- "rate_anx"
}
if ("AgeSurvey" %in% colnames(df.in)) {
  colnames(df.in)[colnames(df.in) == "AgeSurvey"] <- "Age"
}
if ("id" %in% colnames(df.in)) {
  colnames(df.in)[colnames(df.in) == "id"] <- "SN"
}
demoCols <- c("SN", "Age", "Gender", "MathCourse", "Race", "rate_anx")

df.demo <- df.in[, demoCols]
write.table(df.demo, file=demoFilename, quote=F, sep = "\t", row.names=F)

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

# Reshape the data from wide to long format
df.long <- df.in %>%
  pivot_longer(
		cols = starts_with("target_"),  # Selects all columns that start with "target_"
		names_to = "target",            # Name of the new column storing the column names
		names_prefix = "target_",       # Strips "target_" away, leaving just the number
		values_to = "estimate"          # Name of the new column storing the 
	)

	### recode conditions so it is bounded vs unbounded
	df.long$nlType <- "bounded"
	df.long$isUniversalBounded <- FALSE
	df.long$s_ub <- 1000000000
	df.long$s_lb <- 0
	df.long$s_maxBound <- 1000000000
	df.long$s_minBound <- 0
	
	df.long$target <- ifelse(df.long$target == "NA", NA, as.numeric(as.character(df.long$target)))
	df.long$estimate <- ifelse(df.long$estimate == "NA", NA, as.numeric(as.character(df.long$estimate)))


#### Filter outliers - 
#first identify and remove problem targets 
df.long <- df.long %>% group_by(nlType, SN) %>% mutate(z_diff = scale(target - estimate))
df.zDiff <- df.long %>% group_by(target) %>% summarize(meanZdiff = mean(z_diff, na.rm = T))
problemTargets <- df.zDiff[abs(df.zDiff$meanZdiff) > probTargetThresh, ]
df.long.rmProblemTargets <- df.long[!(df.long$target %in% problemTargets$target), ]
df.long.problemTargets <- df.long[df.long$target %in% problemTargets$target, ]
df.long.rmProblemTargets$z_diff <- NULL

#now identify and remove individual outliers
df.long.rmProblemTargets$diff <- df.long.rmProblemTargets$target - df.long.rmProblemTargets$estimate
#get the median absolute difference and, for each point, calculate the number of MADs each difference it. 
df.long.rmProblemTargets <- df.long.rmProblemTargets %>% group_by(nlType, SN) %>% mutate(median_diff = median(abs(diff), na.rm = T), mad_diff = diff/median_diff)
#filter by a threshold 
df.long.rmProblemTargets.filtered <- df.long.rmProblemTargets[abs(df.long.rmProblemTargets$mad_diff) < madOutlierThresh, ]
df.long.IndividualOutliers <- df.long.rmProblemTargets[abs(df.long.rmProblemTargets$mad_diff) >= madOutlierThresh, ]
df.long.IndividualOutliers <- df.long.IndividualOutliers[!is.na(df.long.IndividualOutliers$SN), ]

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

