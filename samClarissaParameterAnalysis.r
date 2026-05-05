library(nlFit)
library(chutils)
library(dplyr)
library(readxl)
library(tidyr)
library(lmerTest)
library(emmeans)

options(tibble.width = Inf)

source("../nlFitRun.r")

nlFitLoops <- 250
targetOrder <- "random"
#Only run new models that have not been run before?
freepIncludeConceptualPoints <- FALSE
conceptualReferencePoints <- 0

runNLfit <- FALSE
includeUBR <- FALSE
belowR2threshold <- 0
statsOutFilename <- paste(fileTag, "StatsOut.txt", sep="")
mainDirPathPA <- getwd()

nlFitFilename <- "./nlFit/nlGenFitData.txt"
tradFilename <- "./traditional/traditionalFitData.txt"

#demographicCols <- c("SN", "AMAS_Comp", "Bio_Sex", "Gender", "Race", "Race_13_TEXT", "Age", "EdLev", "HouseInc", "Diag_Stat", "Diag_Age", "pLastCasto", "HCs1Casto", "HCs1Casto_3_TEXT", "HCs2Casto", "HCs2Casto_7_TEXT", "HCb1Casto", "HCb2Casto")
df.raw <- read_excel(filenameRawDataExl)
### fix inconsistent column naming between datasets
if ("AMAS_comp" %in% colnames(df.raw)) {
  colnames(df.raw)[colnames(df.raw) == "AMAS_comp"] <- "AMAS_Comp"
}
if ("Sex" %in% colnames(df.raw)) {
  colnames(df.raw)[colnames(df.raw) == "Sex"] <- "Bio_Sex"
}
df.demo <- df.raw[, c("SN", demoCols)]

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

df.param.long.1 <- merge(df.param.nlFit.in, df.param.trad.in, by = c("sn", "nlType"))

if(includeUBR) {
	df.param.long <- df.param.long.1
} else {
	df.param.long <- df.param.long.1[df.param.long.1$nlType != "unboundedBetweenLabels", ]
}
nlTypes <- unique(df.param.long$nlType)


############## get demographics
sink(statsOutFilename, append = F)
	cat("\n\n **************** Demographics before Filtering ****************\n\n")
	cat("\n Total N: ", nrow(df.demo), "\n")
	if("Bio_Sex" %in% names(df.demo)) print(with(df.demo, table(list(Bio_Sex, Gender))))
	if("Race" %in% names(df.demo)) print(with(df.demo, table(Race)))
	if("Age" %in% names(df.demo)) cat("\n\n Mean Age: ",mean(df.demo$Age))
	if("Age" %in% names(df.demo)) cat("\n SD Age: ",sd(df.demo$Age), "\n\n")
sink(NULL)
##### remove bad subjects
badSubs <- unique(df.param.long[df.param.long$R2 < belowR2threshold, "sn"])
df.param.long <- df.param.long[!(df.param.long$sn %in% badSubs), ]
df.demo <- df.demo[!(df.demo$SN %in% badSubs), ]

##### remove nlnbinay subjects
nonBinary <- unique(unlist(df.demo[df.demo$Bio_Sex != df.demo$Gender, "SN"]))
df.demo <- df.demo[!(df.demo$SN %in% nonBinary), ]
df.param.long <- df.param.long[!(df.param.long$sn %in% nonBinary), ]

############## get demographics after filtering
sink(statsOutFilename, append = T)
	cat("\n\n **************** Demographics after Filtering ****************\n\n")
	cat("\n Final N: ", nrow(df.demo), "\n")
	cat("\n Removed because R2 < ", belowR2threshold, ": ", length(badSubs), "\n")
	print(badSubs)
	cat("\n\n Removed because Bio_Sex != Gender: ", length(nonBinary), "\n\n")
	print(nonBinary)
	cat("\n\n")
	if("Bio_Sex" %in% names(df.demo)) print(with(df.demo, table(list(Bio_Sex, Gender))))
	if("Race" %in% names(df.demo)) print(with(df.demo, table(Race)))
	if("Age" %in% names(df.demo)) cat("\n\n Mean Age: ",mean(df.demo$Age))
	if("Age" %in% names(df.demo)) cat("\n SD Age: ",sd(df.demo$Age), "\n\n")
		
		
	cat("\n\n **************** Gender by AMAS_Comp ****************\n\n")
	print(with(df.demo, t.test(AMAS_Comp ~ Bio_Sex)))
sink(NULL)

df.param <- merge(df.param.long, df.demo, by.x = c("sn"), by.y = c("SN"))
goodSubs <- as.numeric(unique(df.param$sn))

	#### run analysis on group data to generate plot for publication
if(runNLfit) {
	filenameNLdata <- paste("samClarissaProcessedData_", fileTag, ".txt", sep="")
	df.nlData <- read.table(filenameNLdata, header = T, sep="\t")
	df.nlData <- df.nlData[df.nlData$SN %in% goodSubs, ]
	yLim <- c(min(df.nlData$s_minBound), max(df.nlData$s_maxBound))
	xLim <- yLim

	ch.newDir (mainDirPathPA, "nlFitGroupData")
	nlFitPath <- getwd()

	df.tmp.1 <- df.nlData %>% group_by(nlType, s_targ) %>% summarize(s_est = mean(s_est, na.rm = T), s_maxBound= mean(s_maxBound, na.rm = T), s_minBound= mean(s_minBound, na.rm = T), s_lb= mean(s_lb, na.rm = T),s_ub= mean(s_ub, na.rm = T), isUniversalBounded= mean(isUniversalBounded, na.rm = T))

	types <- unique(df.nlData$nlType)
	for(typ in types) {
		df.tmp <- df.tmp.1[df.tmp.1$nlType == typ, ]		

		plotFileName.1 <- paste("all", typ, "nlFitplot.pdf", sep="-")
		sinkFileName.1 <- paste("all", typ, "nlFitout.txt", sep="-")
		outlist.r <- nlFitRun(df.tmp$s_targ, df.tmp$s_est, upperBound = df.tmp$s_maxBound[1], lowerBound = df.tmp$s_minBound[1], visibleReferencePoints = c(df.tmp$s_lb[1], df.tmp$s_ub[1]), conceptualReferencePoints = conceptualReferencePoints, freepIncludeConceptualPoints = freepIncludeConceptualPoints, targetOrder = targetOrder, plotFileName = plotFileName.1, sinkFileName = sinkFileName.1, loops = nlFitLoops, yLim = yLim, xLim = xLim)
	}

	setwd(mainDirPathPA)
}
	#############

df.param.wide <- df.param %>% pivot_wider(id_cols = c("sn", "AMAS_Comp", "Bio_Sex"), names_from = c(nlType), values_from = all_of(allParams), names_glue = "{nlType}_{.value}")

#### Get Summary statistics
df.sum <- df.param %>%
  group_by(Bio_Sex, nlType) %>%
  summarize(
    # Apply functions to every column in your vector
    across(
      all_of(allParams),
      
      # List of functions to apply
      list(m = ~mean(., na.rm = TRUE), 
           sd = ~sd(., na.rm = TRUE)),
           
      # This creates the custom names like "mAP", "sdAP"
      .names = "{.fn}{.col}"
    ),
    
    # Good practice to ungroup after summarizing
    .groups = "drop" 
  )

df.sum.sex <- df.param %>%
  group_by(Bio_Sex) %>%
  summarize(
    # Apply functions to every column in your vector
    across(
      all_of(allParams),
  
      # List of functions to apply
      list(m = ~mean(., na.rm = TRUE), 
           sd = ~sd(., na.rm = TRUE)),
       
      # This creates the custom names like "mAP", "sdAP"
      .names = "{.fn}{.col}"
    ),

    # Good practice to ungroup after summarizing
    .groups = "drop" 
  )

df.sum.typ <- df.param %>%
  group_by(nlType) %>%
  summarize(
    # Apply functions to every column in your vector
    across(
      all_of(allParams),
    
      # List of functions to apply
      list(m = ~mean(., na.rm = TRUE), 
           sd = ~sd(., na.rm = TRUE)),
         
      # This creates the custom names like "mAP", "sdAP"
      .names = "{.fn}{.col}"
    ),
  
    # Good practice to ungroup after summarizing
    .groups = "drop" 
  )


sink(statsOutFilename, append = T)
	cat("\n\n number of participants who responded randomly: ", length(unique(badSubs)), "\n")
	print(badSubs)
	cat("\n\n **************** Summary ****************\n\n")
	print(df.sum.sex)
	print(df.sum.typ)
	print(df.sum)

	df.cor.out <- NULL
	pdf("nlTypeParams.pdf")
	for(pr in allParams) {
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
			cat("\n\n **************** ", pr, nt1, " v ", nt2, " ****************\n\n")
			print(summary(lm.tmp))
			print(cor.res.1.e)
			print(cor.res.1.p)
			
			df.cor.out <- ch.rbind(df.cor.out, data.frame(param = pr, nlType1 = nt1, nlType2 = nt2, cor = cor.res.1.e, pVal = cor.res.1.p))
			
		}
	}

#### relation between PAE and AP
df.cor2.out <- NULL
nlTypes <- unique(df.param$nlType)
for(pr in nlFitParams) {
	for(nt in nlTypes) {
		for(tp in tradParams) {
			cat("\n\n **************** ", nt, tp, " ~ ", pr, " ****************\n\n")
			form.tmp <- paste(tp, " ~ ", pr)
			lm.tmp <- lm(form.tmp, data = df.param[df.param$nlType == nt,])
			
			cor.res.e <- round(cor.test(df.param[df.param$nlType == nt, tp] , df.param[df.param$nlType == nt, pr])$estimate, 2)
			cor.res.p <- round(cor.test(df.param[df.param$nlType == nt, tp] , df.param[df.param$nlType == nt, pr])$p.value, 4)
			df.cor2.out <- ch.rbind(df.cor2.out, data.frame(nlType = nt, nlParam = pr, tradParam = tp, cor = cor.res.e, pValue = cor.res.p))
			
			print(summary(lm.tmp))
			print(cor.res.e)
			print(cor.res.p)
			with(df.param[df.param$nlType == nt,], plot(as.formula(form.tmp), main = paste(nt, tp, " ~ ", pr)))
			abline(lm.tmp)
			if(tp == "dPAE") abline(h = 0, col = "grey")
			
		}
	}
}


for (pr in allParams) {
	tmp.fr <- as.formula(paste(pr, "~ Bio_Sex*nlType + (1|sn)"))
	lme.out <- lmerTest::lmer(tmp.fr, data = df.param)

	cat("\n\n **************** ", pr, " Mixed Effect Anova ****************\n\n")
	print(anova(lme.out))
	if(anova(lme.out)["nlType", "Pr(>F)"] <= 0.05) {
		post_hoc_tukey <- emmeans(lme.out, specs = pairwise ~ nlType, adjust = "tukey")		
		cat("\n\n******* nlType: Post Hoc Tukey \n\n")		
		print(post_hoc_tukey)
	}
	if(anova(lme.out)["Bio_Sex", "Pr(>F)"] <= 0.05) {
		post_hoc_tukey <- emmeans(lme.out, specs = pairwise ~ Bio_Sex, adjust = "tukey")		
		cat("\n\n******* Bio_Sex: Post Hoc Tukey \n\n")		
		print(post_hoc_tukey)
	}
	if(anova(lme.out)["Bio_Sex:nlType", "Pr(>F)"] <= 0.05) {
		emm_interaction <- emmeans(lme.out, pairwise ~ nlType | Bio_Sex, adjust = "Tukey")
		cat("\n\n******* Interaction: Post Hoc Tukey \n\n")		
		print(emm_interaction)
	}
	
	if(pr %in% c("FE", "AP", "NS", "R2")) {
		box.fr <- as.formula(paste(pr, "~ nlType"))
		boxplot(box.fr, data = df.param, frame = F, cex.axis = 1.25, main = pr, ylim = c(0,1))

		box.fr <- as.formula(paste(pr, "~ nlType*Bio_Sex"))
		boxplot(box.fr, data = df.param, frame = F, cex.axis = 1.25, main = pr, ylim = c(0,1))
	}
}

if(dt == "prolific"){
	df.param.hc <- df.param[df.param$Bio_Sex == "Female",]
	df.param.hc$HCuse <- ifelse(df.param.hc$HCs1Casto == "I AM currently using a prescription contraceptive", TRUE, 				ifelse(df.param.hc$HCs1Casto == "I am NOT currently using any prescription contraceptives", FALSE,NA ))
	for (pr in allParams) {
		tmp.fr <- as.formula(paste(pr, "~ HCuse*nlType + (1|sn)"))
		lme.out <- lmerTest::lmer(tmp.fr, data = df.param.hc)

		cat("\n\n **************** ", pr, " HC use Mixed Effect Anova ****************\n\n")
		print(anova(lme.out))
	}
	
}
dev.off()


nlTypes <- unique(df.param$nlType)
nlTypes <- nlTypes[nlTypes != "unboundedBetweenLabels"]
sigParams <- NULL
df.cor3.out <- NULL
pdf("anxParamsRelation.pdf")
for(nt in nlTypes) {
	df.tmp <- df.param[df.param$nlType == nt, ]
	for(pr in allParams) {
		tmp.fr <- paste(pr, " ~ AMAS_Comp")
		lm.tmp <- lm(tmp.fr, data = df.tmp)
		if (summary(lm.tmp)$coefficients[2,4] < 0.05) {
			sigParams <- rbind(sigParams, data.frame(nlType = nt, param = pr))
		}
		
		cor.res.e <- round(cor.test(df.param[df.param$nlType == nt, "AMAS_Comp"] , df.param[df.param$nlType == nt, pr])$estimate, 2)
		cor.res.p <- round(cor.test(df.param[df.param$nlType == nt, "AMAS_Comp"] , df.param[df.param$nlType == nt, pr])$p.value, 4)
		df.cor3.out <- ch.rbind(df.cor3.out, data.frame(nlType = nt, anx = "AMAS_Comp", param = pr, cor = cor.res.e, pValue = cor.res.p))
		
		with(df.tmp, plot(as.formula(tmp.fr), main = paste(nt, pr)))
		abline(lm.tmp)
		cat("\n ******** nlType: ", nt, " Parameter ", pr, " By Anxiety\n")
		print(summary(lm.tmp))
		print(cor.res.e)
		print(cor.res.p)
	}
}
dev.off()


nlTypes <- unique(df.param$nlType)
pdf("PAE_BIC_Relation.pdf")
	lm.tmp <- lm(df.param$PAE ~ df.param$BIC)
	plot(df.param$PAE ~ df.param$BIC, main = paste("PAE BIC"))
	abline(lm.tmp)
	cat("\n ******** PAE By BIC\n")
	print(summary(lm.tmp))
dev.off()

testParams <- sigParams[!(sigParams$param %in% grep('BIC|R2', sigParams$param, value=TRUE)),]
nlTypes <- unique(testParams$nlType)
for(nt in nlTypes) {
	df.tmp <- df.param[df.param$nlType == nt, ]
	tmpParams <- testParams[testParams$nlType == nt, "param"]
	#remove BIC 
	lm.BIC <- lm("AMAS_Comp ~ BIC", data = df.tmp)
	df.tmp$resAnx <- resid(lm.BIC)
	for(pr in tmpParams) {
		tmp.fr <- paste(pr, " ~ resAnx")
		lm.tmp <- lm(tmp.fr, data = df.tmp)
		cat("\n ******** nlType: ", nt, " Parameter ", pr, " by BIC residual Anxiety\n")
		print(summary(lm.tmp))
	}
}

anxParams <- c("AP","BIC")
for(nt in nlTypes) {
	df.tmp <- df.param[df.param$nlType == nt, ]
	for(pr in anxParams) {
		tmp.fr <- paste(pr, " ~ Bio_Sex + AMAS_Comp + Bio_Sex:AMAS_Comp")
		lm.tmp <- lm(tmp.fr, data = df.tmp)
		
		cat("\n ******** nlType: ", nt, " Parameter ", pr, " By Anxiety\n")
		print(anova(lm.tmp))
		print(summary(lm.tmp))
	}
}



sink(NULL)
