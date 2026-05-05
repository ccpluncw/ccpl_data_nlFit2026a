library(nlFit)
library(chutils)
library(dplyr)
library(readxl)
library(tidyr)
library(lmerTest)
library(emmeans)

options(tibble.width = Inf)

fileTag <- "Kent"
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
demoFilename <- "mielickiDemo.txt"

df.demo <- read.table(demoFilename, header = T, sep = "\t")

df.param.nlFit.in <- read.table(nlFitFilename, header = T, sep="\t")
df.param.nlFit.in <- df.param.nlFit.in %>% rename(
											FE = nl_firstEstimate_Random,
											AP = nl_accuracyPercent_Random,
											NS = nl_numberSensitivity_Random,
											R2 = nl_r2_Random,
											BIC = nl_BIC_Random
											)
df.param.long <- df.param.nlFit.in[, c("sn", "nlType", "FE", "AP", "NS", "R2", "BIC")]
nlFitParams <- names(df.param.long)
nlFitParams <- nlFitParams[!(nlFitParams %in% c("sn", "nlType"))]

nlTypes <- unique(df.param.long$nlType)

############## get demographics
sink(statsOutFilename, append = F)
	cat("\n\n **************** Demographics before Filtering ****************\n\n")
	cat("\n Total N: ", nrow(df.demo), "\n")
	if("Gender" %in% names(df.demo)) print(with(df.demo, table(Gender, useNA = "ifany")))
	if("Race" %in% names(df.demo)) print(with(df.demo, table(Race), useNA = "ifany"))
	if("Age" %in% names(df.demo)) cat("\n\n Mean Age: ",mean(df.demo$Age, na.rm=T))
	if("Age" %in% names(df.demo)) cat("\n SD Age: ",sd(df.demo$Age, na.rm=T), "\n\n")
sink(NULL)
##### remove bad subjects
badSubs <- unique(df.param.long[df.param.long$R2 < belowR2threshold, "sn"])
df.param.long <- df.param.long[!(df.param.long$sn %in% badSubs), ]
df.demo <- df.demo[!(df.demo$SN %in% badSubs), ]

##### remove nlnbinay subjects
nonBinary <- unique(unlist(df.demo[df.demo$Gender > 2, "SN"]))
Nas <- unique(unlist(df.demo[is.na(df.demo$Gender), "SN"]))
nonBinary <- c(nonBinary, Nas)
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
	if("Gender" %in% names(df.demo)) print(with(df.demo, table(Gender)))
	if("Race" %in% names(df.demo)) print(with(df.demo, table(Race)))
	if("Age" %in% names(df.demo)) cat("\n\n Mean Age: ",mean(df.demo$Age))
	if("Age" %in% names(df.demo)) cat("\n SD Age: ",sd(df.demo$Age), "\n\n")
		
		
	cat("\n\n **************** Gender by AMAS_Comp ****************\n\n")
	print(with(df.demo, t.test(rate_anx ~ Gender)))
sink(NULL)

df.param <- merge(df.param.long, df.demo, by.x = c("sn"), by.y = c("SN"))
goodSubs <- as.numeric(unique(df.param$sn))

df.param.wide <- df.param %>% pivot_wider(id_cols = c("sn", "rate_anx", "Gender"), names_from = c(nlType), values_from = all_of(nlFitParams), names_glue = "{nlType}_{.value}")

#### Get Summary statistics
df.sum.all <- df.param %>%
  summarize(
		n = n(),
    # Apply functions to every column in your vector
    across(
      all_of(nlFitParams),
  
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
  group_by(Gender) %>%
  summarize(
		n = n(),
    # Apply functions to every column in your vector
    across(
      all_of(nlFitParams),

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
	print(df.sum.all)
	print(df.sum.sex)


pdf("GenderBoxplots.pdf")
for (pr in nlFitParams) {
	tmp.fr <- as.formula(paste(pr, "~ Gender"))
	t.out <- t.test(tmp.fr, data = df.param)
	boxplot(tmp.fr, data = df.param, frame = F, main = pr)

	cat("\n\n **************** ", pr, " T.test by gender ****************\n\n")
	print((t.out))
}
dev.off()

sigParams <- NULL
df.cor3.out <- NULL
pdf("anxParamsRelation.pdf")
	for(pr in nlFitParams) {
		tmp.fr <- paste(pr, " ~ rate_anx")
		lm.tmp <- lm(tmp.fr, data = df.param)
		if (summary(lm.tmp)$coefficients[2,4] < 0.05) {
			sigParams <- rbind(sigParams, data.frame(param = pr))
		}
		
		cor.res.e <- round(cor.test(df.param[, "rate_anx"] , df.param[, pr])$estimate, 2)
		cor.res.p <- round(cor.test(df.param[, "rate_anx"] , df.param[, pr])$p.value, 4)
		df.cor3.out <- ch.rbind(df.cor3.out, data.frame(anx = "rate_anx", param = pr, cor = cor.res.e, pValue = cor.res.p))
		
		with(df.param, plot(as.formula(tmp.fr), main = paste(pr)))
		abline(lm.tmp)
		cat("\n ******** Parameter ", pr, " By Anxiety\n")
		print(summary(lm.tmp))
		print(cor.res.e)
		print(cor.res.p)
	}

dev.off()

anxParams <- c("AP","FE", "NS", "BIC")
	for(pr in anxParams) {
		tmp.fr <- paste(pr, " ~ Gender + rate_anx + Gender:rate_anx")
		lm.tmp <- lm(tmp.fr, data = df.param)

		cat("\n ******** Parameter ", pr, " By Anxiety\n")
		print(anova(lm.tmp))
		print(summary(lm.tmp))
	}



sink(NULL)
