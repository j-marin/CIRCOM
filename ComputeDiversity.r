library("phyloseq")
library("ggplot2")
library("vegan")
library("DESeq2")
library("stringr")
library("dplyr")

# Open mothur results
a <- import_biom("../../Mothur/mothur/data/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.1.biom")
df_b <- as.data.frame(otu_table(a))

# Open data
data <- read.csv("../recap.csv", sep = ";", header = TRUE, na.strings = "", fileEncoding = "UTF-8-BOM")
attach(data)
dim(data)

CorrTab <- read.csv2(file = "../../../Stages/Younes/Analyse diversité/CorrTable2.csv")

# Filter data
data_filtre <- filter(data, Remarque != "MELANGE ESWAB")
data_supp <- filter(data, Remarque == "MELANGE ESWAB")
df_b1 = df_b %>%
  select(colnames=-data_supp$Sample)
  
# rarefaction curve
rarecurve(t(df_b1), step=50, cex=0.5) #### HERE
    
# sample size
plot_bar(a, fill="Rank6") + theme(legend.position = "none")
ss <- colSums(otu_table(a))

#harmonize abondance
plot_bar(a, fill="Rank6") + theme(legend.position = "none")
a.rarefied = rarefy_even_depth(a, rngseed=1, sample.size=1000, replace=F)
plot_bar(a.rarefied, fill="Rank6") + theme(legend.position = "none")

# compute alpha diversity (Observed and Shannon)
alphaD <- cbind(Sample=str_replace(rownames(alphaD), "X", ""),estimate_richness(a.rarefied, measures = "Observed"),estimate_richness(a.rarefied, measures = "Shannon"))

