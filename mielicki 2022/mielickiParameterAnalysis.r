library(chutils)
library(dplyr)
library(tidyr)
library(emmeans)

options(tibble.width = Inf)

###### functions
getAlpha <- function(df.items) {
	df.items <- df.items[complete.cases(df.items), ]
	k <- ncol(df.items)
	alpha <- (k/(k-1)) * (1 - sum(apply(df.items, 2, var))/var(rowSums(df.items)))
	df.out <- data.frame(kItems = k, nComplete = nrow(df.items), alpha = alpha)
	return(df.out)
}

#descriptives for one variable
descriptives <- function(x, varName) {
	x <- x[!is.na(x)]
	data.frame(variable = varName, n = length(x), mean = mean(x), sd = sd(x), min = min(x), max = max(x))
}

#name the covariate model when a parameter is fit more than one way
insertModel <- function(df.in, model) {
	if(is.null(model)) return(df.in)
	return(data.frame(df.in[, "param", drop = F], model = model, df.in[, setdiff(names(df.in), "param")]))
}

#anova table in long form
aovLong <- function(aov.in, param, n, model = NULL) {
	df.out <- data.frame(param = param, term = rownames(aov.in), df = aov.in$Df, SS = aov.in$"Sum Sq", F = aov.in$"F value", p = aov.in$"Pr(>F)", n = n)
	return(insertModel(df.out, model))
}

#gender by anxiety model, with simple slopes and the gender gap
anxietyModel <- function(df.in, yNm, param, model = NULL) {
	fr.anx <- as.formula(paste(yNm, baseFormula))
	lm.anx <- lm(fr.anx, data = df.in)
	aov.anx <- anova(lm.anx)
	#male is the first factor level, so the slopes come out male then female
	sl <- as.data.frame(test(emtrends(lm.anx, ~ Gender, var = "rate_anxZ")))
	#anxiety is standardized, so this is the gender gap at mean anxiety
	gp <- as.data.frame(pairs(emmeans(lm.anx, ~ Gender, at = list(rate_anxZ = 0))))
	df.slope <- data.frame(param = param, Gender = sl$Gender, slope = sl$rate_anxZ.trend, se = sl$SE, df = sl$df, t = sl$t.ratio, p = sl$p.value)
	df.gap <- data.frame(param = param, contrast = gp$contrast, estimate = gp$estimate, se = gp$SE, df = gp$df, t = gp$t.ratio, p = gp$p.value)
	#print the formula that was fit, not the expression that built it
	lm.anx$call$formula <- fr.anx
	return(list(lm = lm.anx, aov = aov.anx,
				reg = aovLong(aov.anx, param, nobs(lm.anx), model),
				slopes = insertModel(df.slope, model),
				gap = insertModel(df.gap, model)))
}

############# end

belowR2threshold <- 0
statsOutFilename <- paste(fileTag, "StatsOut.txt", sep="")

anxParams <- c("AP", "FE", "NS", "BIC")
covariates <- c("equiv", "torr", "nMath")
#the matrix is on the raw columns, standardizing cannot change a correlation
corVars <- c("AP", "FE", "NS", "BIC", "rate_anx", covariates)
#entry order: equivalence, then reasoning, then course count
hierOrder <- c("equivZ", "torrZ", "nMathZ")
#the reduced model keeps the predictors significant in the ordered entry
reducedAlpha <- 0.05
#every predictor is standardized, the parameters are not
baseFormula <- " ~ Gender + rate_anxZ + Gender:rate_anxZ"

nlFitFilename <- file.path(".", nlFitDirName, "nlGenFitData.txt")
tradFilename <- "./traditional/traditionalFitData.txt"
exclusionsFilename <- paste("mielickiExclusions_", fileTag, ".txt", sep="")
nlDataFilename <- paste("mielickiProcessedData_", fileTag, ".txt", sep="")

df.demo <- read.table(file.path(rawDataDir, demoFilename), header = T, sep = "\t")

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
	if("Gender" %in% names(df.demo)) print(with(df.demo, table(Gender, useNA = "ifany")))
	if("Race" %in% names(df.demo)) print(with(df.demo, table(Race), useNA = "ifany"))
	if("Age" %in% names(df.demo)) cat("\n\n Mean Age: ",mean(df.demo$Age, na.rm=T))
	if("Age" %in% names(df.demo)) cat("\n SD Age: ",sd(df.demo$Age, na.rm=T), "\n\n")
sink(NULL)
##### remove bad subjects
nPreFilter <- nrow(df.demo)
badSubs <- unique(df.param.long[df.param.long$R2 < belowR2threshold, "sn"])
df.param.long <- df.param.long[!(df.param.long$sn %in% badSubs), ]
df.demo <- df.demo[!(df.demo$SN %in% badSubs), ]

##### remove subjects who made no estimates
#a participant with no row in the parameter file never made an estimate
noEstimates <- df.demo$SN[!(df.demo$SN %in% df.param.long$sn)]
df.demo <- df.demo[!(df.demo$SN %in% noEstimates), ]

##### remove nonbinary subjects
#gender: 1 male, 2 female, 3 nonbinary, 4 different, 9 rather not report
nonBinary <- unique(unlist(df.demo[which(df.demo$Gender > 2), "SN"]))
Nas <- unique(unlist(df.demo[is.na(df.demo$Gender), "SN"]))
nonBinary <- c(nonBinary, Nas)
df.demo <- df.demo[!(df.demo$SN %in% nonBinary), ]
df.param.long <- df.param.long[!(df.param.long$sn %in% nonBinary), ]

############## get demographics after filtering
sink(statsOutFilename, append = T)
	cat("\n\n\n\n **************** Demographics after Filtering ****************\n\n")
	cat("\n Final N: ", nrow(df.demo), "\n")
	cat("\n Removed because R2 < ", belowR2threshold, ": ", length(badSubs), " (", round(100*length(badSubs)/nPreFilter, 2), "% of ", nPreFilter, ")\n")
	print(badSubs)
	cat("\n\n Removed because they made no estimates: ", length(noEstimates), " (", round(100*length(noEstimates)/nPreFilter, 2), "% of ", nPreFilter, ")\n\n")
	print(noEstimates)
	cat("\n\n Removed because gender was not male or female: ", length(nonBinary), " (", round(100*length(nonBinary)/nPreFilter, 2), "% of ", nPreFilter, ")\n\n")
	print(nonBinary)
	cat("\n\n")
	if("Gender" %in% names(df.demo)) print(with(df.demo, table(Gender)))
	if("Race" %in% names(df.demo)) print(with(df.demo, table(Race)))
	if("Age" %in% names(df.demo)) cat("\n\n Mean Age: ",mean(df.demo$Age))
	if("Age" %in% names(df.demo)) cat("\n SD Age: ",sd(df.demo$Age), "\n\n")

	cat("\n\n\n\n **************** Response Exclusions ****************\n\n")
	print(read.table(exclusionsFilename, header = T, sep="\t"))

	cat("\n\n\n\n **************** Gender by rate_anx ****************\n\n")
	print(with(df.demo, t.test(rate_anx ~ Gender)))
sink(NULL)

df.param <- merge(df.param.long, df.demo, by.x = c("sn"), by.y = c("SN"))
goodSubs <- as.numeric(unique(df.param$sn))
#male first, so every gender effect reads as female minus male
df.param$Gender <- factor(df.param$Gender, levels = c(1, 2), labels = c("male", "female"))

#### summary statistics
#n, mean and sd of every parameter, named mAP, sdAP and so on
df.sum.all <- df.param %>%
  summarize(
		n = n(),
    across(
      all_of(allParams),
      list(m = ~mean(., na.rm = TRUE),
           sd = ~sd(., na.rm = TRUE)),
      .names = "{.fn}{.col}"
    ),
    .groups = "drop"
  )

df.sum.sex <- df.param %>%
  group_by(Gender) %>%
  summarize(
		n = n(),
    across(
      all_of(allParams),
      list(m = ~mean(., na.rm = TRUE),
           sd = ~sd(., na.rm = TRUE)),
      .names = "{.fn}{.col}"
    ),
    .groups = "drop"
  )

#### anxiety descriptives, one row per participant
df.anx <- unique(df.param[, c("sn", "Gender", "rate_anx")])
df.anx.out <- NULL
#c() on a factor drops the labels, so loop over the labels themselves
for(gd in c(as.character(unique(df.anx$Gender)), "all")) {
	if(gd == "all") {
		df.tmp <- df.anx
	} else {
		df.tmp <- df.anx[df.anx$Gender == gd, ]
	}
	df.anx.out <- ch.rbind(df.anx.out, data.frame(anxiety = "rate_anx", Gender = gd, n = nrow(df.tmp), mean = mean(df.tmp$rate_anx, na.rm = T), sd = sd(df.tmp$rate_anx, na.rm = T), min = min(df.tmp$rate_anx, na.rm = T), max = max(df.tmp$rate_anx, na.rm = T)))
}

#### reliability of the absolute error across the number line items
df.items <- read.table(nlDataFilename, header = T, sep="\t")
df.items <- df.items[df.items$SN %in% goodSubs, ]
df.items$absError <- abs(df.items$estimate - df.items$target)
df.items$item <- df.items$target
df.rel.out <- NULL
for(nt in nlTypes) {
	df.tmp <- df.items[df.items$nlType == nt, ]
	df.tmp.wide <- df.tmp[, c("SN", "item", "absError")] %>% pivot_wider(id_cols = "SN", names_from = "item", values_from = "absError")
	df.tmp.wide <- as.data.frame(df.tmp.wide[, names(df.tmp.wide) != "SN"])
	df.rel.out <- ch.rbind(df.rel.out, data.frame(nlType = nt, measure = "absError", getAlpha(df.tmp.wide)))
}

#### reliability of the five math anxiety items in rate_anx
df.anxItems <- df.demo[df.demo$SN %in% goodSubs, names(anxItems)]
df.rel.out <- ch.rbind(df.rel.out, data.frame(nlType = "all", measure = "rate_anx", getAlpha(df.anxItems)))

#### reliability of the scored covariate items, on the same participants
df.equivItems <- df.demo[df.demo$SN %in% goodSubs, paste("eq_", names(equivKey), sep="")]
df.rel.out <- ch.rbind(df.rel.out, data.frame(nlType = "all", measure = "equiv", getAlpha(df.equivItems)))
df.torrItems <- df.demo[df.demo$SN %in% goodSubs, paste("tr_", names(torrKey), sep="")]
df.rel.out <- ch.rbind(df.rel.out, data.frame(nlType = "all", measure = "torr", getAlpha(df.torrItems)))


sink(statsOutFilename, append = T)
	cat("\n\n\n\n **************** Summary ****************\n\n")
	print(df.sum.all)
	print(df.sum.sex)
	write.table(df.sum.all, file=paste(fileTag, "SummaryStatsAll.txt", sep=""), quote=F, sep = "\t", row.names=F)
	write.table(df.sum.sex, file=paste(fileTag, "SummaryStatsByGender.txt", sep=""), quote=F, sep = "\t", row.names=F)

	cat("\n\n\n\n **************** Anxiety Descriptives ****************\n\n")
	print(df.anx.out)
	write.table(df.anx.out, file=paste(fileTag, "AnxietyDescriptives.txt", sep=""), quote=F, sep = "\t", row.names=F)

	cat("\n\n\n\n **************** Reliability ****************\n\n")
	print(df.rel.out)
	write.table(df.rel.out, file=paste(fileTag, "Reliability.txt", sep=""), quote=F, sep = "\t", row.names=F)

#the parameter frame is long by nlType, so collapse to one row per participant
df.cov.p <- df.param[!duplicated(df.param$sn), ]

#### correlations among the parameters, anxiety and the covariates
df.desc.out <- NULL
for(cv in corVars) {
	df.desc.out <- ch.rbind(df.desc.out, descriptives(df.cov.p[[cv]], cv))
}
#pairwise complete, so a missing score does not cost the other pairs their cases
#Hmisc is not attached because it would mask the dplyr summarize
cor.out <- Hmisc::rcorr(as.matrix(df.cov.p[, corVars]))

cat("\n\n\n\n **************** Correlation Matrix ****************\n\n")
print(df.desc.out)
cat("\n")
print(round(cor.out$r, 2))
cat("\n")
print(round(cor.out$P, 4))
cat("\n n complete on all: ", sum(complete.cases(df.cov.p[, corVars])), "\n")
write.table(data.frame(variable = corVars, cor.out$r), file=paste(fileTag, "CorMatrix.txt", sep=""), quote=F, sep = "\t", row.names=F)
write.table(data.frame(variable = corVars, cor.out$P), file=paste(fileTag, "CorMatrixP.txt", sep=""), quote=F, sep = "\t", row.names=F)

#### gender, anxiety and their interaction on each parameter
df.anxReg.out <- NULL
cat("\n\n\n\n **************** Gender by Anxiety ****************\n\n")
for(pr in anxParams) {
	am <- anxietyModel(df.param, pr, pr)
	cat("\n\n\n ******** Parameter ", pr, " By Anxiety\n")
	print(am$aov)
	print(summary(am$lm))
	cat("\n simple slope of anxiety within each gender\n")
	print(am$slopes)
	cat("\n gender gap at mean anxiety\n")
	print(am$gap)
	df.anxReg.out <- ch.rbind(df.anxReg.out, am$reg)
}

cat("\n\n\n ******** Gender By Anxiety Summary ********\n\n")
#lead with nlType, as the kent and prolific tables do
df.anxReg.out <- data.frame(nlType = nlTypes[1], df.anxReg.out)
print(df.anxReg.out)
write.table(df.anxReg.out, file=paste(fileTag, "RegSexByAnxiety.txt", sep=""), quote=F, sep = "\t", row.names=F)

#### hierarchical regression of each parameter on the math covariates
df.hier.out <- NULL
df.red.out <- NULL
lm.reduced <- list()
reducedRhs <- list()
cat("\n\n\n\n **************** Hierarchical Regression on Math Covariates ****************\n\n")
for(pr in anxParams) {
	df.cc.p <- df.cov.p[complete.cases(df.cov.p[, c(pr, covariates)]), ]
	fr.hier <- as.formula(paste(pr, "~", paste(hierOrder, collapse=" + ")))
	lm.hier <- lm(fr.hier, data = df.cc.p)
	lm.hier$call$formula <- fr.hier
	aov.hier <- anova(lm.hier)
	#each step's sequential SS as a proportion of the total SS
	st <- rownames(aov.hier) != "Residuals"
	r2ch <- aov.hier$"Sum Sq"[st]/sum(aov.hier$"Sum Sq")
	df.step <- data.frame(param = pr, step = 1:sum(st), term = rownames(aov.hier)[st], df = aov.hier$Df[st], SS = aov.hier$"Sum Sq"[st], F = aov.hier$"F value"[st], p = aov.hier$"Pr(>F)"[st], R2change = r2ch, R2cum = cumsum(r2ch), n = nobs(lm.hier))
	df.hier.out <- ch.rbind(df.hier.out, df.step)

	cat("\n\n\n ******** Parameter ", pr, " By Math Covariates\n")
	print(df.step)

	#the reduced model keeps the steps that cleared alpha, in order
	keep <- hierOrder[aov.hier$"Pr(>F)"[st] < reducedAlpha]
	#a parameter that keeps nothing gets the intercept only model
	rhs <- if(length(keep) > 0) paste(keep, collapse=" + ") else "1"
	fr.red <- as.formula(paste(pr, "~", rhs))
	lm.red <- lm(fr.red, data = df.cc.p)
	lm.red$call$formula <- fr.red
	aov.red <- anova(lm.red)
	#the residualized models use what these models leave over
	lm.reduced[[pr]] <- lm.red
	reducedRhs[[pr]] <- rhs
	df.red.out <- ch.rbind(df.red.out, data.frame(param = pr, model = rhs, term = rownames(aov.red), df = aov.red$Df, SS = aov.red$"Sum Sq", F = aov.red$"F value", p = aov.red$"Pr(>F)", R2 = summary(lm.red)$r.squared, n = nobs(lm.red)))

	cat("\n reduced model: ", pr, " ~ ", rhs, "\n")
	print(aov.red)
	print(summary(lm.red))
}

cat("\n\n\n ******** Hierarchical Steps Summary ********\n\n")
print(df.hier.out)
write.table(df.hier.out, file=paste(fileTag, "HierCovariates.txt", sep=""), quote=F, sep = "\t", row.names=F)

cat("\n\n\n ******** Reduced Covariate Models Summary ********\n\n")
print(df.red.out)
write.table(df.red.out, file=paste(fileTag, "ReducedCovariates.txt", sep=""), quote=F, sep = "\t", row.names=F)

#### gender by anxiety on the math-residualized parameters
df.res.out <- NULL
df.resCor.out <- NULL
cat("\n\n\n\n **************** Gender by Anxiety on Math-Residualized Parameters ****************\n\n")
for(pr in anxParams) {
	rNm <- paste(pr, "resid", sep="")
	df.cc.p <- df.cov.p[complete.cases(df.cov.p[, c(pr, covariates)]), ]
	#the reduced model was fit on these same rows, so the residuals line up
	df.cc.p[[rNm]] <- resid(lm.reduced[[pr]])
	am <- anxietyModel(df.cc.p, rNm, pr, reducedRhs[[pr]])
	cat("\n\n\n ******** Parameter ", rNm, " By Anxiety, from ", pr, " ~ ", reducedRhs[[pr]], "\n")
	print(am$aov)
	print(summary(am$lm))
	cat("\n simple slope of anxiety within each gender\n")
	print(am$slopes)
	cat("\n gender gap at mean anxiety\n")
	print(am$gap)
	df.res.out <- ch.rbind(df.res.out, am$reg)

	ct <- cor.test(df.cc.p[[rNm]], df.cc.p$rate_anxZ)
	df.resCor.out <- ch.rbind(df.resCor.out, data.frame(param = pr, r = as.numeric(ct$estimate), p = ct$p.value, n = as.numeric(ct$parameter) + 2))
}

cat("\n\n\n ******** Residual By Anxiety Summary ********\n\n")
print(df.res.out)
write.table(df.res.out, file=paste(fileTag, "ResidAnxietyReg.txt", sep=""), quote=F, sep = "\t", row.names=F)

cat("\n residual with anxiety\n\n")
print(df.resCor.out)

sink(NULL)
