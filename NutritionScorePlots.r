library("FactoMineR")
library("factoextra")
library("ggpubr")
library("pwr")
library("rstatix")

data <- read.csv2(file = "scores.txt", sep="\t")
div <- read.csv2(file = "../Diversity/recapDiversity.csv", sep=",")

#### Cas Tem
tab <- as.data.frame(data[,c(2,4,7)])
tab$NutScore <- as.numeric(tab$NutScore)

# NutScore
w_test <- compare_means(NutScore ~ CasTem,  data = tab, method = "wilcox.test", p.adjust.method="hochberg")
w_test <- w_test %>% mutate(y.position = c(12, 14, 16))
ggboxplot(tab, x = "CasTem", y = "NutScore", color = "CasTem", add = "dotplot")+ 
         stat_pvalue_manual(w_test, label="p.adj")

## Power analysis
# compute effect size
model <- aov(NutScore ~ CasTem, data = tab)
f <- eta_squared(model)

# anova (k = nb of groups, n = nb of obs per group, f = effect size, power = 1 minus Type II error probability)
p1 <- pwr.anova.test(k=3, n=table(tab$CasTem), f=f)
mean(p1$power) # 0.1795279
pwr.anova.test(k=3, n=rep(100,3), f=f)

#### Type cirrhose
tab <- as.data.frame(data[,c(10,4,7)])
tab$CirrhosisType2[data$CasTem == "Healthy"] <- "Healthy"
tab$NutScore <- as.numeric(tab$NutScore)

# NutScore
w_test <- compare_means(NutScore ~ CirrhosisType2,  data = tab, method = "wilcox.test", p.adjust.method="hochberg")
w_test <- w_test %>% mutate(y.position = seq(12, 30, 2))
ggboxplot(tab, x = "CirrhosisType2", y = "NutScore", color = "CirrhosisType2", add = "dotplot")+ 
         stat_pvalue_manual(w_test, label="p.adj")

## Power analysis
# compute effect size
model <- aov(BodyMassIndex ~ CirrhosisType2, data = tab)
f <- eta_squared(model)

# anova (k = nb of groups, n = nb of obs per group, f = effect size, power = 1 minus Type II error probability)
p2 <- pwr.anova.test(k=5, n=table(tab$CirrhosisType2), f=f)
mean(p2$power) # 0.6188778
pwr.anova.test(k=5, n=rep(100,5), f=f)





