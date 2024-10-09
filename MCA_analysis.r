library("FactoMineR")
library("factoextra")
#library("missMDA")

rm(list=ls())

# open data
data <- read.csv2(file = "../CIRCOM060324.txt", sep="\t")

############################################################
# Make variables as categories (by quartiles)
make_cat <- function(x, n) {
	x <- as.numeric(x)
	qt1 <- quantile(x , na.rm=TRUE)
	
	res <- integer(0)
	for (i in which(!is.na(x))){
		if (x[i] >= qt1[1] & x[i] <= qt1[2]) {res[i] <- paste(n, "(q1)", sep="")}
		if (x[i] > qt1[2] & x[i] <= qt1[3]) {res[i] <- paste(n, "(q2)", sep="")}
		if (x[i] > qt1[3] & x[i] <= qt1[4]) {res[i] <- paste(n, "(q3)", sep="")}
		if (x[i] > qt1[4] & x[i] <= qt1[5]) {res[i] <- paste(n, "(q4)", sep="")}
		}
	res[is.na(x)] <- NA	
	return (res)
	}

data$AntibioticTreatment6Months[data$AntibioticTreatment6Months == 1 & !is.na(data$AntibioticTreatment6Months)] <- "antibiotic intake"
data$AntibioticTreatment6Months[data$AntibioticTreatment6Months == 0 & !is.na(data$AntibioticTreatment6Months)] <- "no antibiotic intake"

data$ESBLPresence[data$ESBLPresence == 1 & !is.na(data$ESBLPresence)] <- "ESBL presence"
data$ESBLPresence[data$ESBLPresence == 0 & !is.na(data$ESBLPresence)] <- "ESBL absence"

## MCA
# Clinical data
data2.1 <- data[!is.na(data$AntibioticTreatment6Months) & !is.na(data$AssociatedSeripousnessGravity) & !is.na(data$AssociatedRiskFactorESBL) & !is.na(data$ESBLPresence),]

Gravity2 <- make_cat(data2.1$AssociatedSeripousnessGravity, "Gravity")
data2.1["Gravity2"] <- Gravity2

ESBLRisk2 <- make_cat(data2.1$AssociatedRiskFactorESBL, "ESBLRisk")
data2.1["ESBLRisk2"] <- ESBLRisk2

data2.2 <- data2.1[c("CasTem", "ESBLRisk2", "AntibioticTreatment6Months", "ESBLPresence")]
data2.2 <- apply(data2.2, 2, as.factor)
table(data2.2[,1])

res.mca = MCA(data2.2)
round(res.mca$eig, 3)
dimdesc(res.mca)
summary(res.mca)
round(res.mca$var$eta, 3)

### logistic regression to assess the effects
data2.1["y"] <- as.numeric(as.factor(data2.1$CasTem))
data2.1$y <- ifelse(data2.1$y == 1, 1, 0) # Control = 0 / Case = 1

binom1 <- glm(y ~ ESBLPresence, data=data2.1, family=binomial(link=logit))
summary(binom1) # ns
binom2 <- glm(y ~ AntibioticTreatment6Months, data=data2.1, family=binomial(link=logit))
summary(binom2) # 0.00497
binom3 <- glm(y ~ AssociatedRiskFactorESBL, data=data2.1, family=binomial(link=logit))
summary(binom3) # 0.0337


# Lifestyle data
data2.1 <- data[!is.na(data$SportScore) & !is.na(data$PrecariousnessScore) & !is.na(data$Shannon) & !is.na(data$NutScore),] # 24 lignes

NutScore2 <- make_cat(data2.1 $NutScore, "NutScore")
data2.1["NutScore2"] <- c(NutScore2)

ScorePreca2 <- make_cat(data2.1 $PrecariousnessScore, "PrecariousnessScore")
data2.1["PrecariousnessScore2"] <- ScorePreca2

DivAlpha2 <- make_cat(data2.1 $Shannon, "AlphaDiversity")
data2.1["Shannon2"] <- DivAlpha2

ScoreSport2 <- make_cat(data2.1 $SportScore, "SportScore")
data2.1["SportScore2"] <- ScoreSport2

data2.2 <- data2.1[c("CasTem", "SportScore2", "Shannon2", "NutScore2","PrecariousnessScore2")]
data2.2 <- apply(data2.2, 2, as.factor)
table(data2.2[,1])

res.mca = MCA(data2.2)
round(res.mca$eig, 3)
dimdesc(res.mca)
summary(res.mca)
round(res.mca$var$eta, 3)

### logistic regression to assess the effects
data2.1["y"] <- as.numeric(as.factor(data2.1$CasTem))
data2.1$y <- ifelse(data2.1$y == 1, 1, 0) # Control = 0 / Case = 1
data2.1[,c(7:9, 13)] <- apply(data2.1[,c(7:9, 13)], 2, as.numeric) 

binom1 <- glm(y ~ PrecariousnessScore, data=data2.1, family=binomial(link=logit))
summary(binom1) # ns
binom2 <- glm(y ~ SportScore, data=data2.1, family=binomial(link=logit))
summary(binom2) # ns
binom3 <- glm(y ~ NutScore, data=data2.1, family=binomial(link=logit))
summary(binom3) # ns
binom4 <- glm(y ~ Shannon, data=data2.1, family=binomial(link=logit))
summary(binom4) # 0.0124

