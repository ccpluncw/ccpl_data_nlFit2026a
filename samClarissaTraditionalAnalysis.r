library(chutils)
library(tidyr)

###### functions
getMLLM <- function(targetVals, estimateVals, upperBound, constrained = TRUE) {
	df.tmp <- data.frame(target = targetVals, estimate = estimateVals)
	#add constrained LL model
	df.tmp$logTarget <- (upperBound/log(upperBound)) * log(df.tmp$target) 
	
	if(constrained) {
		####### LL without the a
		fit.ll <- tryCatch ({
			 nls(estimate ~ ( ((1-lambda)*target) + (lambda*logTarget)), data=df.tmp, start=list(lambda = 1), algorithm="port", lower=c(0), upper=c(1), control = nls.control(minFactor=1/10000000, maxiter=10000, warnOnly = FALSE))
		}, error = function(x) {
		            print(paste("Constrained LL Model Failed"))								
		            # Choose a return value in case of error
		            NULL
		})

		if(!is.null(fit.ll)) {
			df.tmp$fit.ll <- predict(fit.ll, newdata = df.tmp)
			coef.lam.ll <- round(coef(fit.ll)[[1]], 2)
			r2.ll <- round(ch.R2(df.tmp$estimate, df.tmp$fit.ll, standardize = FALSE),2)
			BIC.ll <- round(ch.BIC(df.tmp$estimate, df.tmp$fit.ll, 1, standardize = FALSE), 0)
			coef.int.ll <- 1
		} else {
			coef.int.ll <- NA
			coef.lam.ll <- NA
			r2.ll <- NA
			BIC.ll <- NA
		}
	} else {

		### LL model with the a
		####### LL without the a
		fit.ll <- tryCatch ({
			 nls(estimate ~ a*( ((1-lambda)*target) + (lambda*logTarget)), data=df.tmp, start=list(a=1, lambda = 1), 					algorithm="port", lower=c(-Inf, 0), upper=c(Inf, 1), control = nls.control(minFactor=1/10000000, maxiter=10000, warnOnly = FALSE))
		}, error = function(x) {
		            print(paste ("Full LL Model Failed on"))								
		            # Choose a return value in case of error
		            NULL
		})

		if(!is.null(fit.ll)) {
			df.tmp$fit.ll <- predict(fit.ll, newdata = df.tmp)
			coef.int.ll <- round(coef(fit.ll)[[1]],2)
			coef.lam.ll <- round(coef(fit.ll)[[2]], 2)
			r2.ll <- round(ch.R2(df.tmp$estimate, df.tmp$fit.ll, standardize = FALSE),2)
			BIC.ll <- round(ch.BIC(df.tmp$estimate, df.tmp$fit.ll, 2, standardize = FALSE), 0)
		} else {
			coef.int.ll <- NA
			coef.lam.ll <- NA
			r2.ll <- NA
			BIC.ll <- NA
	
		}
	}
	
	df.out <- data.frame(intercept = coef.int.ll, lambda = coef.lam.ll, r2 = r2.ll, BIC = BIC.ll)
	return(df.out)
	
}

############# end

skipModelsWithExistingDirs <- FALSE
runMLLM <- FALSE
scaleValue <- 10
PAEdivisor <- 6

filenamexl <- paste("samClarissaProcessedData_", fileTag, ".txt", sep="")
df.long <- read.table(filenamexl, header = T, sep="\t")

mainDirPath <- getwd()

##### Run analysis on each subject and as a group
subjects <- unique(df.long$SN)
types <- unique(df.long$nlType)
types <- c(types, "unboundedBetweenLabels")
subjects <- c(subjects, "all")
keep55mil <- c(TRUE, FALSE)

keepCols <- c("nlType", "s_targ", "s_est", "s_ub", "s_lb", "s_maxBound", "s_minBound", "isUniversalBounded")
all.tmp <- NULL

df.long.1 <- df.long
subDirExists <- ch.newDir (mainDirPath, "traditional")
tradPath <- getwd()
outFilename <- "traditionalFitData.txt"
outPlotFileName <- "traditionalFitData.pdf"

xHigh <- (abs(min(df.long.1$s_minBound)) + max(df.long.1$s_maxBound)) * scaleValue
xLim <- c(0, xHigh)
yLim <- xLim
if(subDirExists == FALSE | skipModelsWithExistingDirs == FALSE) {
	if(runMLLM) pdf(outPlotFileName)

	for(sub in subjects) {
		if(sub == "all") {
			df.tmp.1 <- df.long.1 %>% group_by(nlType, s_targ) %>% summarize(s_est = mean(s_est, na.rm = T), s_maxBound= mean(s_maxBound, na.rm = T), s_minBound= mean(s_minBound, na.rm = T), s_lb= mean(s_lb, na.rm = T),s_ub= mean(s_ub, na.rm = T), isUniversalBounded= mean(isUniversalBounded, na.rm = T))
		} else {
			df.tmp.1 <- df.long.1[df.long.1$SN == sub, keepCols]
		}

		for(typ in types) {

			if(typ == "unboundedBetweenLabels") {
				df.tmp <- df.tmp.1[df.tmp.1$isUniversalBounded == TRUE, ]
				shift <- 1
				upperBound <- (df.tmp$s_ub[1] + shift) * scaleValue
			} else {
				df.tmp <- df.tmp.1[df.tmp.1$nlType == typ, ]		
				shift <- abs(df.tmp$s_minBound[1])
				upperBound <- (df.tmp$s_maxBound[1] + shift) * scaleValue
			}

			#### calculate PAE
			 df.tmp$error <- df.tmp$s_est - df.tmp$s_targ
			 ##use constant PAEdivisor to make bounded and unbounded comparable
			 PAE <-  mean(abs(df.tmp$error), na.rm = T)/PAEdivisor
			 dPAE <-  mean(df.tmp$error, na.rm = T)/PAEdivisor
	
			 df.tmp.fit <- data.frame(sn = sub, nlType = typ, PAE = PAE, dPAE = dPAE)

			if(runMLLM) {
				#### add constant to estimate and target to make all targets positive
				df.tmp$s_targ <- (df.tmp$s_targ + shift) * scaleValue
				df.tmp$s_est <- (df.tmp$s_est + shift) * scaleValue
				df.tmp$s_ub <- (df.tmp$s_ub + shift) * scaleValue
				df.tmp$s_lb <- (df.tmp$s_lb + shift) * scaleValue
 
				#### calculate MLLM Constrained
				df.constrained <- getMLLM(df.tmp$s_targ, df.tmp$s_est, upperBound, constrained = TRUE) 
				df.tmp$logTarget <- (upperBound/log(upperBound)) * log(df.tmp$s_targ) 
				df.tmp$fit.c <- df.constrained$intercept*( ((1-df.constrained$lambda)*df.tmp$s_targ) + (df.constrained$lambda*df.tmp$logTarget))
				#### calculate MLLM Unconstrained
				df.unconstrained <- getMLLM(df.tmp$s_targ, df.tmp$s_est, upperBound, constrained = FALSE) 
				df.tmp$fit.uc <- df.unconstrained$intercept*( ((1-df.unconstrained$lambda)*df.tmp$s_targ) + (df.unconstrained$lambda*df.tmp$logTarget))

				df.tmp <- df.tmp[order(df.tmp$s_targ), ]
				with(df.tmp, plot(s_est ~ s_targ, ylim = yLim, xlim = xLim, main = paste(sub, typ)))
				abline(a = 0, b = 1, col = "grey75", lty=2)
				with(df.tmp, lines(fit.c ~ s_targ, col = "black"))
				with(df.tmp, lines(fit.uc ~ s_targ, col = "grey50"))

				df.tmp.fit <- cbind(df.tmp.fit, data.frame(Unc_int =  df.unconstrained$intercept, Unc_Lam =  df.unconstrained$lambda, Unc_R2 =  df.unconstrained$r2, Unc_BIC =  df.unconstrained$BIC, Con_Lam =  df.constrained$lambda, Con_R2 =  df.constrained$r2, Con_BIC =  df.constrained$BIC))
			}

			all.tmp <- ch.rbind(all.tmp, df.tmp.fit)
		}
	}
	write.table(all.tmp, file=outFilename, quote=F, sep = "\t", row.names=F)
	setwd(mainDirPath)
	if(runMLLM) dev.off()
}

