library(nlFit)
library(chutils)
library(dplyr)
library(readxl)
library(tidyr)
library(lmerTest)
library(emmeans)

options(tibble.width = Inf)

source(nlFitRunScript)

###### functions
getAlpha <- function(df.items) {
	df.items <- df.items[complete.cases(df.items), ]
	k <- ncol(df.items)
	alpha <- (k/(k-1)) * (1 - sum(apply(df.items, 2, var))/var(rowSums(df.items)))
	df.out <- data.frame(kItems = k, nComplete = nrow(df.items), alpha = alpha)
	return(df.out)
}

#anova table in long form
aovLong <- function(aov.in, param, n) {
	df.out <- data.frame(param = param, term = rownames(aov.in), df = aov.in$Df, SS = aov.in$"Sum Sq", F = aov.in$"F value", p = aov.in$"Pr(>F)", n = n)
	return(df.out)
}

#the sex by anxiety model with the simple slopes and the sex gap it implies
anxietyModel <- function(df.in, yNm, param) {
	fr.anx <- as.formula(paste(yNm, baseFormula))
	lm.anx <- lm(fr.anx, data = df.in)
	aov.anx <- anova(lm.anx)
	#male is the first factor level, so the slopes come out male then female
	sl <- as.data.frame(test(emtrends(lm.anx, ~ Bio_Sex, var = "AMAS_CompZ")))
	#anxiety is standardized, so this is the sex gap at mean anxiety
	gp <- as.data.frame(pairs(emmeans(lm.anx, ~ Bio_Sex, at = list(AMAS_CompZ = 0))))
	df.slope <- data.frame(param = param, Bio_Sex = sl$Bio_Sex, slope = sl$AMAS_CompZ.trend, se = sl$SE, df = sl$df, t = sl$t.ratio, p = sl$p.value)
	df.gap <- data.frame(param = param, contrast = gp$contrast, estimate = gp$estimate, se = gp$SE, df = gp$df, t = gp$t.ratio, p = gp$p.value)
	#print the formula that was fit, not the expression that built it
	lm.anx$call$formula <- fr.anx
	return(list(lm = lm.anx, aov = aov.anx,
				reg = aovLong(aov.anx, param, nobs(lm.anx)),
				slopes = df.slope,
				gap = df.gap))
}

############# end

nlFitLoops <- 250
targetOrder <- "random"
#skip the group fit if its folder already exists
skipModelsWithExistingDirs <- TRUE
freepConceptual <- FALSE
conceptualReferencePoints <- 0

runNLfit <- TRUE
belowR2threshold <- 0
statsOutFilename <- paste(fileTag, "StatsOut.txt", sep="")
mainDirPathPA <- getwd()

#parameters in the number line type by sex mixed models
aovParams <- c("FE", "AP", "NS", "R2")
#parameters in the sex by anxiety regressions
anxParams <- c("AP", "FE", "NS", "BIC")
#anxiety is standardized in the regressions, the parameters are not
baseFormula <- " ~ Bio_Sex + AMAS_CompZ + Bio_Sex:AMAS_CompZ"

nlFitFilename <- file.path(".", nlFitDirName, "nlGenFitData.txt")
tradFilename <- "./traditional/traditionalFitData.txt"
exclusionsFilename <- paste("samClarissaExclusions_", fileTag, ".txt", sep="")
nlDataFilename <- paste("samClarissaProcessedData_", fileTag, ".txt", sep="")

df.raw <- read_excel(file.path(rawDataDir, filenameRawDataExl), na = c("", "NA"))
### fix inconsistent column naming between datasets
if ("AMAS_comp" %in% colnames(df.raw)) {
  colnames(df.raw)[colnames(df.raw) == "AMAS_comp"] <- "AMAS_Comp"
}
if ("Sex" %in% colnames(df.raw)) {
  colnames(df.raw)[colnames(df.raw) == "Sex"] <- "Bio_Sex"
}
df.demo <- df.raw[, c("SN", demoCols)]
#standardized before the exclusions, so the scale does not move with them
df.demo$AMAS_CompZ <- as.numeric(scale(df.demo$AMAS_Comp))

df.param.nlFit.in <- read.table(nlFitFilename, header = T, sep="\t")
df.param.nlFit.in <- df.param.nlFit.in %>% rename(
											FE = nl_firstEstimate_Random,
											AP = nl_accuracyPercent_Random,
											NS = nl_numberSensitivity_Random,
											R2 = nl_r2_Random,
											BIC = nl_BIC_Random
											)
df.param.nlFit.in <- df.param.nlFit.in[, c("sn", "nlType", "FE", "AP", "NS", "R2", "BIC")]
nlFitParams <- names(df.param.nlFit.in)
nlFitParams <- nlFitParams[!(nlFitParams %in% c("sn", "nlType"))]

df.param.trad.in <- read.table(tradFilename, header = T, sep="\t")
tradParams <- names(df.param.trad.in)
tradParams <- tradParams[!(tradParams %in% c("sn", "nlType"))]
allParams <- c(nlFitParams, tradParams)

df.param.long <- merge(df.param.nlFit.in, df.param.trad.in, by = c("sn", "nlType"))
nlTypes <- unique(df.param.long$nlType)


############## get demographics
sink(statsOutFilename, append = F)
	cat("\n\n\n\n **************** Demographics before Filtering ****************\n\n")
	cat("\n Total N: ", nrow(df.demo), "\n")
	if("Bio_Sex" %in% names(df.demo)) print(with(df.demo, table(list(Bio_Sex, Gender))))
	if("Race" %in% names(df.demo)) print(with(df.demo, table(Race)))
	if("Age" %in% names(df.demo)) cat("\n\n Mean Age: ",mean(df.demo$Age))
	if("Age" %in% names(df.demo)) cat("\n SD Age: ",sd(df.demo$Age), "\n\n")
sink(NULL)
##### remove bad subjects
nPreFilter <- nrow(df.demo)
badSubs <- unique(df.param.long[df.param.long$R2 < belowR2threshold, "sn"])
df.param.long <- df.param.long[!(df.param.long$sn %in% badSubs), ]
df.demo <- df.demo[!(df.demo$SN %in% badSubs), ]

##### remove nonbinary subjects
nonBinary <- unique(unlist(df.demo[df.demo$Bio_Sex != df.demo$Gender, "SN"]))
df.demo <- df.demo[!(df.demo$SN %in% nonBinary), ]
df.param.long <- df.param.long[!(df.param.long$sn %in% nonBinary), ]

############## get demographics after filtering
sink(statsOutFilename, append = T)
	cat("\n\n\n\n **************** Demographics after Filtering ****************\n\n")
	cat("\n Final N: ", nrow(df.demo), "\n")
	cat("\n Removed because R2 < ", belowR2threshold, ": ", length(badSubs), " (", round(100*length(badSubs)/nPreFilter, 2), "% of ", nPreFilter, ")\n")
	print(badSubs)
	cat("\n\n Removed because Bio_Sex != Gender: ", length(nonBinary), " (", round(100*length(nonBinary)/nPreFilter, 2), "% of ", nPreFilter, ")\n\n")
	print(nonBinary)
	cat("\n\n")
	if("Bio_Sex" %in% names(df.demo)) print(with(df.demo, table(list(Bio_Sex, Gender))))
	if("Race" %in% names(df.demo)) print(with(df.demo, table(Race)))
	if("Age" %in% names(df.demo)) cat("\n\n Mean Age: ",mean(df.demo$Age))
	if("Age" %in% names(df.demo)) cat("\n SD Age: ",sd(df.demo$Age), "\n\n")

	cat("\n\n\n\n **************** Response Exclusions ****************\n\n")
	print(read.table(exclusionsFilename, header = T, sep="\t"))

	cat("\n\n\n\n **************** Gender by AMAS_Comp ****************\n\n")
	print(with(df.demo, t.test(AMAS_Comp ~ Bio_Sex)))
sink(NULL)

df.param <- merge(df.param.long, df.demo, by.x = c("sn"), by.y = c("SN"))
goodSubs <- as.numeric(unique(df.param$sn))
#male first, so every sex effect reads as female minus male
df.param$Bio_Sex <- factor(df.param$Bio_Sex, levels = c("Male", "Female"))

#### fit the group data for the published plot
if(runNLfit) {
	df.nlData <- read.table(nlDataFilename, header = T, sep="\t")
	df.nlData <- df.nlData[df.nlData$SN %in% goodSubs, ]
	yLim <- c(min(df.nlData$s_minBound), max(df.nlData$s_maxBound))
	xLim <- yLim

	subDirExists <- ch.newDir (mainDirPathPA, nlFitGroupDirName)
	nlFitPath <- getwd()

	if(subDirExists == FALSE | skipModelsWithExistingDirs == FALSE) {
		df.tmp.1 <- df.nlData %>% group_by(nlType, s_targ) %>% summarize(s_est = mean(s_est, na.rm = T), s_maxBound= mean(s_maxBound, na.rm = T), s_minBound= mean(s_minBound, na.rm = T), s_lb= mean(s_lb, na.rm = T),s_ub= mean(s_ub, na.rm = T), isUniversalBounded= mean(isUniversalBounded, na.rm = T))

		types <- unique(df.nlData$nlType)
		for(typ in types) {
			df.tmp <- df.tmp.1[df.tmp.1$nlType == typ, ]

			plotFileName.1 <- paste("all", typ, "nlFitplot.pdf", sep="-")
			sinkFileName.1 <- paste("all", typ, "nlFitout.txt", sep="-")
			outlist.r <- nlFitRun(df.tmp$s_targ, df.tmp$s_est, upperBound = df.tmp$s_maxBound[1], lowerBound = df.tmp$s_minBound[1], visibleReferencePoints = c(df.tmp$s_lb[1], df.tmp$s_ub[1]), conceptualReferencePoints = conceptualReferencePoints, freepConceptual = freepConceptual, targetOrder = targetOrder, plotFileName = plotFileName.1, sinkFileName = sinkFileName.1, loops = nlFitLoops, yLim = yLim, xLim = xLim)
		}
	}

	setwd(mainDirPathPA)
}

df.param.wide <- df.param %>% pivot_wider(id_cols = c("sn", "AMAS_Comp", "Bio_Sex"), names_from = c(nlType), values_from = all_of(allParams), names_glue = "{nlType}_{.value}")

#### summary statistics
#mean and sd of every parameter, named mAP, sdAP and so on
df.sum <- df.param %>%
  group_by(Bio_Sex, nlType) %>%
  summarize(
    across(
      all_of(allParams),
      list(m = ~mean(., na.rm = TRUE),
           sd = ~sd(., na.rm = TRUE)),
      .names = "{.fn}{.col}"
    ),
    .groups = "drop"
  )

df.sum.sex <- df.param %>%
  group_by(Bio_Sex) %>%
  summarize(
    across(
      all_of(allParams),
      list(m = ~mean(., na.rm = TRUE),
           sd = ~sd(., na.rm = TRUE)),
      .names = "{.fn}{.col}"
    ),
    .groups = "drop"
  )

df.sum.typ <- df.param %>%
  group_by(nlType) %>%
  summarize(
    across(
      all_of(allParams),
      list(m = ~mean(., na.rm = TRUE),
           sd = ~sd(., na.rm = TRUE)),
      .names = "{.fn}{.col}"
    ),
    .groups = "drop"
  )

#### anxiety descriptives, one row per participant
df.anx <- unique(df.param[, c("sn", "Bio_Sex", "AMAS_Comp")])
df.anx.out <- NULL
#c() on a factor drops the labels, so loop over the labels themselves
for(sx in c(as.character(unique(df.anx$Bio_Sex)), "all")) {
	if(sx == "all") {
		df.tmp <- df.anx
	} else {
		df.tmp <- df.anx[df.anx$Bio_Sex == sx, ]
	}
	df.anx.out <- ch.rbind(df.anx.out, data.frame(anxiety = "AMAS_Comp", Bio_Sex = sx, n = nrow(df.tmp), mean = mean(df.tmp$AMAS_Comp, na.rm = T), sd = sd(df.tmp$AMAS_Comp, na.rm = T), min = min(df.tmp$AMAS_Comp, na.rm = T), max = max(df.tmp$AMAS_Comp, na.rm = T)))
}

#### reliability of the absolute error across the number line items
df.items <- read.table(nlDataFilename, header = T, sep="\t")
df.items <- df.items[df.items$SN %in% goodSubs, ]
df.items$absError <- abs(df.items$s_est - df.items$s_targ)
df.items$item <- paste(df.items$condition, df.items$trial, sep="")
df.rel.out <- NULL
for(nt in nlTypes) {
	df.tmp <- df.items[df.items$nlType == nt, ]
	df.tmp.wide <- df.tmp[, c("SN", "item", "absError")] %>% pivot_wider(id_cols = "SN", names_from = "item", values_from = "absError")
	df.tmp.wide <- as.data.frame(df.tmp.wide[, names(df.tmp.wide) != "SN"])
	df.rel.out <- ch.rbind(df.rel.out, data.frame(nlType = nt, measure = "absError", getAlpha(df.tmp.wide)))
}

#### reliability of the nine AMAS items in AMAS_Comp
df.amasItems <- df.demo[df.demo$SN %in% goodSubs, amasItems]
df.rel.out <- ch.rbind(df.rel.out, data.frame(nlType = "all", measure = "AMAS_Comp", getAlpha(df.amasItems)))


sink(statsOutFilename, append = T)
	cat("\n\n\n\n **************** Summary ****************\n\n")
	print(df.sum.sex)
	print(df.sum.typ)
	print(df.sum)
	write.table(df.sum.sex, file=paste(fileTag, "SummaryStatsBySex.txt", sep=""), quote=F, sep = "\t", row.names=F)
	write.table(df.sum.typ, file=paste(fileTag, "SummaryStatsByType.txt", sep=""), quote=F, sep = "\t", row.names=F)
	write.table(df.sum, file=paste(fileTag, "SummaryStatsBySexByType.txt", sep=""), quote=F, sep = "\t", row.names=F)

	cat("\n\n\n\n **************** Anxiety Descriptives ****************\n\n")
	print(df.anx.out)
	write.table(df.anx.out, file=paste(fileTag, "AnxietyDescriptives.txt", sep=""), quote=F, sep = "\t", row.names=F)

	cat("\n\n\n\n **************** Reliability ****************\n\n")
	print(df.rel.out)
	write.table(df.rel.out, file=paste(fileTag, "Reliability.txt", sep=""), quote=F, sep = "\t", row.names=F)

	#### each parameter correlated across the two number line types
	df.cor.out <- NULL
	pdf("nlTypeParams.pdf")
	for(pr in nlFitParams) {
		for(idxNt in 2:length(nlTypes)) {
			nt1 <- nlTypes[idxNt - 1]
			nt2 <- nlTypes[idxNt]
			col1 <- paste(nt1, pr, sep = "_")
			col2 <- paste(nt2, pr, sep = "_")
			lm.tmp <- lm(df.param.wide[[col1]] ~ df.param.wide[[col2]])
			plot(df.param.wide[[col1]] ~ df.param.wide[[col2]], main = paste(nt1, "v", nt2, pr))
			abline(lm.tmp)
			cor.res.1.e <- round(cor.test(df.param.wide[[col1]] , df.param.wide[[col2]])$estimate, 2)
			cor.res.1.p <- round(cor.test(df.param.wide[[col1]] , df.param.wide[[col2]])$p.value, 4)
			cat("\n\n\n\n **************** ", pr, nt1, " v ", nt2, " ****************\n\n")
			print(summary(lm.tmp))
			print(cor.res.1.e)
			print(cor.res.1.p)

			df.cor.out <- ch.rbind(df.cor.out, data.frame(param = pr, nlType1 = nt1, nlType2 = nt2, cor = cor.res.1.e, pVal = cor.res.1.p))

		}
	}
	cat("\n\n\n\n **************** Correlations Across Tasks ****************\n\n")
	print(df.cor.out)
	write.table(df.cor.out, file=paste(fileTag, "CorAcrossTasks.txt", sep=""), quote=F, sep = "\t", row.names=F)

#### number line type by sex on each parameter
df.anova.out <- NULL
for (pr in aovParams) {
	tmp.fr <- as.formula(paste(pr, "~ Bio_Sex*nlType + (1|sn)"))
	lme.out <- lmerTest::lmer(tmp.fr, data = df.param)

	cat("\n\n\n\n **************** ", pr, " Mixed Effect Anova ****************\n\n")
	aov.tmp <- anova(lme.out)
	print(aov.tmp)
	df.anova.out <- ch.rbind(df.anova.out, data.frame(param = pr, term = rownames(aov.tmp), numDF = aov.tmp$NumDF, denDF = aov.tmp$DenDF, F = aov.tmp$"F value", p = aov.tmp$"Pr(>F)"))
	if(anova(lme.out)["nlType", "Pr(>F)"] <= 0.05) {
		post_hoc_tukey <- emmeans(lme.out, specs = pairwise ~ nlType, adjust = "tukey")
		cat("\n\n\n******* nlType: Post Hoc Tukey \n\n")
		print(post_hoc_tukey)
	}
	if(anova(lme.out)["Bio_Sex", "Pr(>F)"] <= 0.05) {
		post_hoc_tukey <- emmeans(lme.out, specs = pairwise ~ Bio_Sex, adjust = "tukey")
		cat("\n\n\n******* Bio_Sex: Post Hoc Tukey \n\n")
		print(post_hoc_tukey)
	}
	if(anova(lme.out)["Bio_Sex:nlType", "Pr(>F)"] <= 0.05) {
		emm_interaction <- emmeans(lme.out, pairwise ~ nlType | Bio_Sex, adjust = "Tukey")
		cat("\n\n\n******* Interaction: Post Hoc Tukey \n\n")
		print(emm_interaction)
	}

	box.fr <- as.formula(paste(pr, "~ nlType"))
	boxplot(box.fr, data = df.param, frame = F, cex.axis = 1.25, main = pr, ylim = c(0,1))

	box.fr <- as.formula(paste(pr, "~ nlType*Bio_Sex"))
	boxplot(box.fr, data = df.param, frame = F, cex.axis = 1.25, main = pr, ylim = c(0,1))
}
cat("\n\n\n\n **************** Sex By nlType Anovas ****************\n\n")
print(df.anova.out)
write.table(df.anova.out, file=paste(fileTag, "AnovaSexByType.txt", sep=""), quote=F, sep = "\t", row.names=F)
dev.off()

#### sex, anxiety and their interaction on each parameter
df.reg.out <- NULL
cat("\n\n\n\n **************** Sex by Anxiety ****************\n\n")
for(nt in nlTypes) {
	df.tmp <- df.param[df.param$nlType == nt, ]
	for(pr in anxParams) {
		am <- anxietyModel(df.tmp, pr, pr)
		cat("\n\n\n ******** nlType: ", nt, " Parameter ", pr, " By Anxiety\n")
		print(am$aov)
		print(summary(am$lm))
		cat("\n simple slope of anxiety within each sex\n")
		print(am$slopes)
		cat("\n sex gap at mean anxiety\n")
		print(am$gap)
		df.reg.out <- ch.rbind(df.reg.out, data.frame(nlType = nt, am$reg))
	}
}

cat("\n\n\n ******** Sex By Anxiety Summary ********\n\n")
print(df.reg.out)
write.table(df.reg.out, file=paste(fileTag, "RegSexByAnxiety.txt", sep=""), quote=F, sep = "\t", row.names=F)

sink(NULL)
