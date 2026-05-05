library(tidyr)
library(readxl)
library(dplyr)

boundPxSize <- 350
bufferPxSize <- 75
zOutlierThresh <- 3
madOutlierThresh <- 7.5
probTargetThresh <- 1.25

makeAllDataPositive <- FALSE
scaleValue <- 10

mainDirPath <- getwd()
df.in <- read_excel(filenameRawDataExl)

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

df.demo <- df.in[, c("SN", demoCols)]
write.table(df.demo, file=demoFilename, quote=F, sep = "\t", row.names=F)

df.in <- df.in[, !names(df.in) %in% demoCols]

# Reshape the data from wide to long format
df.long.1 <- df.in %>%
  pivot_longer(
    cols = -c(SN), # Exclude SN columns from pivoting
    names_to = c("condition", "trial", "DV"), # New columns to create from parts of the old column names
    names_pattern = "([A-Za-z]+)(\\d+)_(.*)", # Regular expression to capture parts
    values_to = "value" # Name of the column to store the values
  )

# spread the DVs into their own columns (ub, lb, targ, est)
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

df.long$halfBoundRange <- (df.long$ub - df.long$lb)/2
df.long$s_ub <- (df.long$ub - df.long$halfBoundRange - df.long$lb)/df.long$halfBoundRange
df.long$s_lb <- (df.long$lb - df.long$halfBoundRange - df.long$lb)/df.long$halfBoundRange
df.long$s_targ <- (df.long$targ - df.long$halfBoundRange - df.long$lb)/df.long$halfBoundRange
df.long$s_est <- (df.long$est - df.long$halfBoundRange - df.long$lb)/df.long$halfBoundRange
#### we add the bufferUnits to the target region to capture the 75 pixels on either side
df.long$s_maxBound <- (df.long$maxBound - df.long$halfBoundRange - df.long$lb + df.long$bufferUnits)/df.long$halfBoundRange
df.long$s_minBound <- (df.long$minBound - df.long$halfBoundRange - df.long$lb - df.long$bufferUnits)/df.long$halfBoundRange
###### end same scale

##### now change the s_maxBound and s_minBound to be correct for the bounded number line
df.long$s_maxBound <- ifelse(df.long$nlType == "bounded", df.long$s_ub, df.long$s_maxBound)
df.long$s_minBound <- ifelse(df.long$nlType == "bounded", df.long$s_lb, df.long$s_minBound)

if(makeAllDataPositive) {
		shift <- abs(min(df.long$s_minBound))
		df.long$s_targ <- (df.long$s_targ + shift) * scaleValue
		df.long$s_est <- (df.long$s_est + shift) * scaleValue
		df.long$s_ub <- (df.long$s_ub + shift) * scaleValue
		df.long$s_lb <- (df.long$s_lb + shift) * scaleValue
		df.long$s_minBound <- (df.long$s_minBound + shift) * scaleValue
		df.long$s_maxBound <- (df.long$s_maxBound + shift) * scaleValue
		xLim <- c(min(df.long$s_minBound), max(df.long$s_maxBound))
		yLim <- xLim
} else {
	xLim <- c(-3.5, 3.5)
	yLim <- xLim
}

#### Filter outliers - 
#first identify and remove problem targets 
df.long <- df.long %>% group_by(nlType, SN) %>% mutate(z_diff = scale(s_targ - s_est))
df.zDiff <- df.long %>% group_by(targ) %>% summarize(meanZdiff = mean(z_diff, na.rm = T))
problemTargets <- df.zDiff[abs(df.zDiff$meanZdiff) > probTargetThresh, ]
df.long.rmProblemTargets <- df.long[!(df.long$targ %in% problemTargets$targ), ]
df.long.problemTargets <- df.long[df.long$targ %in% problemTargets$targ, ]
df.long.rmProblemTargets$z_diff <- NULL

#now identify and remove individual outliers
df.long.rmProblemTargets$diff <- df.long.rmProblemTargets$s_targ - df.long.rmProblemTargets$s_est
#get the median absolute difference and, for each point, calculate the number of MADs each difference it. 
df.long.rmProblemTargets <- df.long.rmProblemTargets %>% group_by(nlType, SN) %>% mutate(median_diff = median(abs(diff), na.rm = T), mad_diff = diff/median_diff)
#filter by a threshold 
df.long.rmProblemTargets.filtered <- df.long.rmProblemTargets[abs(df.long.rmProblemTargets$mad_diff) < madOutlierThresh, ]
df.long.IndividualOutliers <- df.long.rmProblemTargets[abs(df.long.rmProblemTargets$mad_diff) >= madOutlierThresh, ]

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

