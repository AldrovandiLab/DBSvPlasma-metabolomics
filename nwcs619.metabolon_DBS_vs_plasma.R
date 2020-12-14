#!/usr/bin/Rscript

library(ggplot2)
library(ape)
library(plyr)
library(reshape2)
library(cluster)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(gplots)
library(vegan)
library(irr)
library(useful)
library(pscl)
library(parallel)
library(igraph)
library(randomForest)
library(ROCR)
library(stringi)
library(mixOmics)
library(ggfortify)
library(Rtsne)
library(ggforce)
library(emmeans)
library(tableone)
library(abind)
library(limma)
library(psych)
library(glmnet)

source("utils.R")
source("mcc.R")

cols.sig <- c("black", "red", "grey"); names(cols.sig) <- c("ns", "sig", "NA")
cols.cohort <- c("#808080", "#cccccc", "#0571b0", "#92c5de", "#ca0020", "#f4a582", "#ffa500", "#ffdb99"); names(cols.cohort) <- c("Control.untreated", "Case.untreated", "Control.zdv", "Case.zdv", "Control.PI-ART", "Case.PI-ART", "Control.other", "Case.other")
siglevel <- 0.05
dircolors <- c("blue", "red", "grey"); names(dircolors) <- c("down", "up", "NS")

mapping_fn <- "data/DBSvPlasma-metabolomics_Mapping.121420.txt"
mapping <- read.table(mapping_fn, header=T, as.is=T, sep="\t")
metadata_variables <- read.table("data/metadata_variables.121420.txt", header=T, as.is=T, sep="\t", row.names=1)
for (mvar in rownames(metadata_variables)) {
	if (metadata_variables[mvar, "type"] == "factor") {
		mapping[,mvar] <- factor(mapping[,mvar])
		if (metadata_variables[mvar, "baseline"] != "") {
			mapping[,mvar] <- relevel(mapping[,mvar], metadata_variables[mvar, "baseline"])
		}
	} else if (metadata_variables[mvar, "type"] == "ordered") {
		if (metadata_variables[mvar, "baseline"] == "") {
			lvls <- unique(mapping[,mvar])
		} else {
			lvls <- unlist(strsplit(metadata_variables[mvar, "baseline"], ","))
		}
		mapping[,mvar] <- ordered(mapping[,mvar], levels=lvls)
	} else if (metadata_variables[mvar, "type"] == "numeric") {
		mapping[,mvar] <- as.numeric(as.character(mapping[,mvar]))
	}
}

#########################################################################################################
### read in Metabolon data
metabolite_levels <- c("BIOCHEMICAL", "SUB.PATHWAY", "SUPER.PATHWAY")
## parse metabolon data and Z-transform
df.metabolon <- list()
metabolon_map <- data.frame()
metabolon_sortorder <- list()
for (st in c("DBS", "Plasma")) {
	df.metabolon[[st]] <- list()
	metabolon <- read.table(sprintf("data/ScaledImpData.%s.121420.txt", st), header=T, as.is=T, sep="\t", quote="", comment.char="")
	colnames(metabolon) <- gsub("^X", "", colnames(metabolon))
	metabolon <- subset(metabolon, SUPER.PATHWAY!="") # remove uncharacterized molecules
	metabolon_map <- rbind(metabolon_map, metabolon[,c("BIOCHEMICAL", "SUB.PATHWAY", "SUPER.PATHWAY")])
	sel <- setdiff(colnames(metabolon), metabolite_levels)
	sel <- intersect(sel, c(as.character(mapping$patid)))
	for (mlevel in metabolite_levels) {
		tmp <- metabolon[,c(mlevel, sel)]
		labels <- tmp[, mlevel]; tmp <- tmp[, setdiff(colnames(tmp), mlevel)]; ids <- colnames(tmp)
		tmp <- as.data.frame(t(apply(tmp, 1, function(x) {
			x <- as.numeric(gsub(",", "", x))
			to_impute <- which(is.na(x))
			x[to_impute] <- min(x, na.rm=T)
			x
		})))
		colnames(tmp) <- ids; tmp[, mlevel] <- labels
		agg <- aggregate(as.formula(sprintf(". ~ %s", mlevel)), tmp, sum); rownames(agg) <- agg[,mlevel]; agg <- agg[,-1]
		agg <- agg[, as.character(intersect(colnames(agg), mapping$patid))] # filter to just the samples in the mapping file
		agg <- t(log(agg))
		agg <- agg[, setdiff(1:ncol(agg), which(is.na(apply(agg, 2, sd))))] # remove entries with zero variation
		df.metabolon[[st]][[length(df.metabolon[[st]])+1]] <- agg
	}
	names(df.metabolon[[st]]) <- metabolite_levels
}
names(df.metabolon) <- c("DBS", "Plasma")
metabolon_map <- unique(metabolon_map); rownames(metabolon_map) <- metabolon_map$BIOCHEMICAL
cols.superpathway <- c(brewer.pal(length(unique(metabolon_map$SUPER.PATHWAY)), "Set1"), "#bbbbbb"); names(cols.superpathway) <- c(unique(metabolon_map$SUPER.PATHWAY), "NOT_METABOLITE")

#######################
## parse raw metabolon data for detection statistics
df.metabolon.raw <- list()
for (st in c("DBS", "Plasma")) {
	df.metabolon.raw[[st]] <- list()
	metabolon <- read.table(sprintf("data/OrigScale.%s.121420.txt", st), header=T, as.is=T, sep="\t", quote="", comment.char="")
	colnames(metabolon) <- gsub("^X", "", colnames(metabolon))
	metabolon <- subset(metabolon, SUPER.PATHWAY!="") # remove uncharacterized molecules
#	metabolon_map <- rbind(metabolon_map, metabolon[,c("BIOCHEMICAL", "SUB.PATHWAY", "SUPER.PATHWAY")])
	sel <- setdiff(colnames(metabolon), metabolite_levels)
	sel <- intersect(sel, c(as.character(mapping$patid)))
	for (mlevel in c("BIOCHEMICAL")) {
		tmp <- metabolon[,c(mlevel, sel)]
		labels <- tmp[, mlevel]; tmp <- tmp[, setdiff(colnames(tmp), mlevel)]; ids <- colnames(tmp)
		tmp <- as.data.frame(t(apply(tmp, 1, function(x) {
			x <- as.numeric(gsub(",", "", x))
			x
		})))
		to_remove <- which(apply(tmp, 1, function(x) all(is.na(x)))) # remove rows that are all NA
		tmp <- tmp[setdiff(rownames(tmp), to_remove),,drop=F]; labels <- labels[setdiff(1:length(labels), to_remove)]
		colnames(tmp) <- ids; tmp[, mlevel] <- labels
		# fix metabolite names as necessary
#		colnames(agg) <- gsub("\\*", "", colnames(agg))
		tmp <- tmp[, as.character(intersect(colnames(tmp), c(mapping$patid)))]
		rownames(tmp) <- labels
		df.metabolon.raw[[st]][[length(df.metabolon.raw[[st]])+1]] <- tmp
	}
	names(df.metabolon.raw[[st]]) <- c("BIOCHEMICAL")
}
names(df.metabolon.raw) <- c("DBS", "Plasma")

out_pdf <- sprintf("output/metabolon_analysis_DBS_vs_plasma.%s.pdf", format(Sys.Date(), "%m%d%y"))
pdf(out_pdf, width=12)

#########################################################################################################
### Analysis
p <- ggplot(mapping) + geom_blank() + theme_classic() + annotate("text", x=1, y=1, label=sprintf("Analysis of DBS vs Plasma metabolomics"), size=16)
print(p)

### QC data about metabolomics
## number of metabolites detected in each sample type, Venn diagrams
mlevel <- "BIOCHEMICAL"
qc <- {}; merged <- data.frame(BIOCHEMICAL=rownames(metabolon_map), detected.maternal_DBS=NA, detected.maternal_plasma=NA, median.maternal_DBS=NA, median.maternal_plasma=NA); rownames(merged) <- merged$BIOCHEMICAL
for (st in c("DBS", "Plasma")) {
	tmp <- df.metabolon.raw[[st]][[mlevel]]
	# count detectable as any non-NA value
	counts.maternal <- apply(tmp[, as.character(mapping$patid)], 1, function(x) length(which(!is.na(x))))
	median.maternal <- apply(tmp[, as.character(mapping$patid)], 1, function(x) median(x,na.rm=T))
	out <- data.frame(BIOCHEMICAL=rownames(tmp), subtype=st, detected.maternal=counts.maternal, median.maternal=median.maternal)
	if (st == "DBS") {
		merged[as.character(out$BIOCHEMICAL), c("BIOCHEMICAL", "detected.maternal_DBS", "median.maternal_DBS")] <- out[, c("BIOCHEMICAL", "detected.maternal", "median.maternal")]
	} else {
		merged[as.character(out$BIOCHEMICAL), c("BIOCHEMICAL", "detected.maternal_plasma", "median.maternal_plasma")] <- out[, c("BIOCHEMICAL", "detected.maternal", "median.maternal")]
	}
	write.table(out, file=sprintf("output/detected_counts.%s.txt", st), quote=F, sep="\t", row.names=F, col.names=T)
}
merged$FLAG.maternal_DBS <- ifelse(is.na(merged$detected.maternal_DBS), FALSE, merged$detected.maternal_DBS > 0)
merged$FLAG.maternal_plasma <- ifelse(is.na(merged$detected.maternal_plasma), FALSE, merged$detected.maternal_plasma > 0)
merged$Group <- ifelse(merged$FLAG.maternal_plasma & merged$FLAG.maternal_DBS, "both", ifelse(merged$FLAG.maternal_DBS, "DBS", ifelse(merged$FLAG.maternal_plasma, "Plasma", "other")))
merged <- subset(merged, Group != "other") # remove BIOCHEMICALS with all NA values
merged$SUB.PATHWAY <- metabolon_map[rownames(merged), "SUB.PATHWAY"]
merged$SUPER.PATHWAY <- metabolon_map[rownames(merged), "SUPER.PATHWAY"]
vc <- vennCounts(merged[, c("FLAG.maternal_DBS", "FLAG.maternal_plasma")])
vennDiagram(vc, cex=c(1,1))
write.table(merged, file=sprintf("output/detected_counts.%s.txt", "merged"), quote=F, sep="\t", row.names=F, col.names=T)
# pie charts of BIOCHEMICAL classes in each group
for (gr in unique(merged$Group)) {
	p <- plot_metabolite_breakdown(rownames(subset(merged, Group==gr)), metabolon_map) + ggtitle(sprintf("Detected metabolites (%s)", gr))
	print(p)
}
# chisq test for enrichment of SUPER.PATHWAY classes
res <- {}
counts <- table(merged[,c("SUPER.PATHWAY", "Group")])
for (cl in unique(merged$SUPER.PATHWAY)) {
	for (st in c("DBS", "Plasma")) {
		tmp <- counts[cl, c(st, "both")]
		sums <- colSums(counts[, c(st, "both")])
		tab <- rbind(tmp, sums-tmp); rownames(tab) <- c(cl, sprintf("not %s", cl))
		test <- chisq.test(tab)
		tmp <- c(as.vector(test$observed), as.vector(test$expected))
		names(tmp) <- c(sprintf("%s - %s (observed)", rep(colnames(test$observed), each=2), rep(rownames(test$observed), 2)), sprintf("%s - %s (expected)", rep(colnames(test$expected), each=2), rep(rownames(test$expected), 2)))
		res <- rbind(res, c(cl, st, test$p.value, tmp))
	}
}
res <- as.data.frame(res)
colnames(res)[1:3] <- c("SUPER.PATHWAY", "SampleType", "pval")
write.table(res, file=sprintf("output/detected_metabolites_by_class.%s.txt", "merged"), quote=F, sep="\t", row.names=F, col.names=T)
# heatmap of SUPER.PATHWAY class counts by group, points scaled to absolute count, colored by enrichment relative to counts in the 'both' group
counts.rel <- normalizeByCols(counts); counts.rel <- log2(counts.rel / counts.rel[,1])
df <- melt(counts); tmp <- melt(counts.rel); df$rel <- tmp[,3]
colnames(df) <- c("SUPER.PATHWAY", "Group", "count", "rel")
df2 <- merge(df, res, by.x=c("SUPER.PATHWAY", "Group"), by.y=c("SUPER.PATHWAY", "SampleType"), all.x=T)[, 1:5]
df2$pval <- as.numeric(as.character(df2$pval)); df2$padj <- p.adjust(df2$pval, method="fdr")
df2$dir <- ifelse(is.na(df2$padj), "NS", ifelse(df2$padj < 0.05, "SIG", "NS"))
p <- ggplot(df2, aes(x=Group, y=SUPER.PATHWAY, fill=rel)) + geom_point(aes(size=count, color=dir), shape=21) + geom_text(aes(label=count), size=2, hjust=0.5, vjust=0.5) + theme_classic() + ggtitle(sprintf("Enrichment of lipid classes")) + scale_fill_gradient2(low="blue", mid="white", high="red") + scale_color_manual(values=c("black", "purple"))
print(p)

## summary statistics and dynamic range for metabolites by DBS/plasma
mapping.sel <- mapping; rownames(mapping.sel) <- mapping$patid
res <- {}
for (st in c("DBS", "Plasma")) {
	sel.metabolites <- rownames(subset(merged, Group %in% c(st, "both")))
	sel.ids <- intersect(c(as.character(mapping$patid), as.character(mapping$cpatid)), rownames(df.metabolon[[st]][[mlevel]]))
	data <- df.metabolon[[st]][[mlevel]][sel.ids, sel.metabolites]
	tmp <- t(apply(data, 2, function(x) {
		c(mean(x, na.rm=T), median(x, na.rm=T), min(x, na.rm=T), max(x, na.rm=T), IQR(x, na.rm=T))
	}))
	colnames(tmp) <- c("mean", "median", "min", "max", "IQR")
	tmp <- as.data.frame(tmp)
	tmp$Assay <- st
	tmp$Group <- merged[rownames(tmp), "Group"]
	tmp$Group2 <- sprintf("%s-%s", tmp$Group, tmp$Assay)
	tmp$metabolite <- rownames(tmp)
	res <- rbind(res, tmp)
}
for (metric in c("mean", "median", "min", "max", "IQR")) {
	test.both <- wilcox.test(as.formula(sprintf("%s ~ Group2", metric)), subset(res, Group2 %in% c("both-DBS", "both-Plasma")))
	test.separate <- wilcox.test(as.formula(sprintf("%s ~ Group2", metric)), subset(res, Group2 %in% c("DBS-DBS", "Plasma-Plasma")))
	p <- ggplot(res, aes_string(x="Group2", y=metric)) + geom_boxplot() + theme_classic() + ggtitle(sprintf("%s by Group+Assay (Wilcox p=%.4g and p=%.4g for both/separate)", metric, test.both$p.value, test.separate$p.value))
	print(p)
}


## ICC and Spearman correlation in paired maternal plasma/DBS samples
subtype <- "maternal"
mapping.sel <- mapping; rownames(mapping.sel) <- mapping$patid
sel.metabolites <- rownames(subset(merged, FLAG.maternal_plasma & FLAG.maternal_DBS))
data <- matrix()
data.plasma <- df.metabolon[["Plasma"]][[mlevel]][as.character(mapping$patid), sel.metabolites]
data.dbs <- df.metabolon[["DBS"]][[mlevel]][as.character(mapping$patid), sel.metabolites]
data <- abind(data.plasma, data.dbs, along=0); dimnames(data)[[1]] <- c("Plasma", "DBS")
## ICC by metabolite
final <- apply(data, 3, function(x) {
	test <- cor.test(x[1,], x[2,], method="spearman")
	unlist(c(mean(x[1,]), mean(x[2,]), sd(x[1,]), sd(x[2,]), mean(abs(x[1,]-x[2,])), ICC(t(x))$results["Single_raters_absolute", c("ICC", "p", "lower bound", "upper bound")], test$estimate, test$p.value))
#	icc(t(x))$value
}); final <- as.data.frame(t(final)); colnames(final) <- c("mean.plasma", "mean.DBS", "sd.plasma", "sd.DBS", "mean_difference", "ICC", "p", "lower_bound", "upper_bound", "rho", "spearman_p")
final$padj <- p.adjust(final$p, method="fdr")
final$spearman_padj <- p.adjust(final$spearman_p, method="fdr")
final <- final[order(final$p), ]
write.table(final, file="output/ICC.metabolite.txt", quote=F, sep="\t", row.names=T, col.names=T)
final <- subset(final, mean.plasma > -5 & mean.DBS > -5) # remove a few really low abundance outliers
test <- cor.test(~ mean.plasma + ICC, final, method="spearman")
p <- ggplot(final, aes(x=mean.plasma, y=ICC)) + geom_point() + stat_smooth(method="lm") + theme_classic() + ggtitle(sprintf("ICC vs mean.plasma (by metabolite, Spearman rho=%.4g, p=%.4g)", test$estimate, test$p.value))
print(p)
test <- cor.test(~ mean.DBS + ICC, final, method="spearman")
p <- ggplot(final, aes(x=mean.DBS, y=ICC)) + geom_point() + stat_smooth(method="lm") + theme_classic() + ggtitle(sprintf("ICC vs mean.DBS (by metabolite, Spearman rho=%.4g, p=%.4g)", test$estimate, test$p.value))
print(p)
p <- ggplot(final, aes(x=ICC)) + geom_histogram(bins=50) + theme_classic() + ggtitle(sprintf("Distribution of ICC (by metabolite, mean=%.2g, median=%.2g)", mean(final$ICC), median(final$ICC))) + geom_vline(xintercept=mean(final$ICC), color="red")
print(p)
test <- cor.test(~ mean.plasma + rho, final, method="spearman")
p <- ggplot(final, aes(x=mean.plasma, y=rho)) + geom_point() + stat_smooth(method="lm") + theme_classic() + ggtitle(sprintf("Spearman rho vs mean.plasma (by metabolite, Spearman rho=%.4g, p=%.4g)", test$estimate, test$p.value))
print(p)
test <- cor.test(~ mean.DBS + rho, final, method="spearman")
p <- ggplot(final, aes(x=mean.DBS, y=rho)) + geom_point() + stat_smooth(method="lm") + theme_classic() + ggtitle(sprintf("Spearman rho vs mean.DBS (by metabolite, Spearman rho=%.4g, p=%.4g)", test$estimate, test$p.value))
print(p)
p <- ggplot(final, aes(x=ICC)) + geom_histogram(bins=50) + theme_classic() + ggtitle(sprintf("Distribution of ICC (by metabolite, mean=%.2g, median=%.2g)", mean(final$ICC), median(final$ICC))) + geom_vline(xintercept=mean(final$ICC), color="red")
print(p)
p <- ggplot(final, aes(x=rho)) + geom_histogram(bins=50) + theme_classic() + ggtitle(sprintf("Distribution of Spearman rho (by metabolite, mean=%.2g, median=%.2g)", mean(final$rho, na.rm=T), median(final$rho, na.rm=T))) + geom_vline(xintercept=mean(final$rho, na.rm=T), color="blue")
print(p)
df <- melt(final[,c("ICC", "rho")])
p <- ggplot(df, aes(x=value, fill=variable)) + geom_histogram(bins=50, position=position_dodge()) + theme_classic() + ggtitle(sprintf("Distribution of ICC+Spearman rho (by metabolite)")) + geom_vline(xintercept=mean(final$ICC, na.rm=T), color="red") + geom_vline(xintercept=mean(final$rho, na.rm=T), color="blue") + scale_fill_brewer(palette="Set1")
print(p)
thresholds <- seq(from=0.05, to=0.95, by=0.05)
df <- sapply(thresholds, function(thresh) {
	length(which(final$ICC >= thresh))
}); df <- data.frame(threshold=thresholds, n=df); df$pct <- 100*(df$n / nrow(final))
p <- ggplot(df, aes(x=threshold, y=pct)) + geom_point() + geom_line() + geom_text(label=sprintf("%d (%.2g%%)", df$n, df$pct)) + theme_classic() + ggtitle(sprintf("Number of metabolites at ICC threshold (n=%d padj<0.05)", length(which(final$padj<0.05))))
print(p)
p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("Number of metabolites at ICC threshold")) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(df), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
print(p)
df <- as.matrix(summary(final$rho))
p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("ICC by metabolite summary statistics")) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(df), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
print(p)
df <- sapply(thresholds, function(thresh) {
	length(which(final$rho >= thresh))
}); df <- data.frame(threshold=thresholds, n=df); df$pct <- 100*(df$n / nrow(final))
p <- ggplot(df, aes(x=threshold, y=pct)) + geom_point() + geom_line() + geom_text(label=sprintf("%d (%.2g%%)", df$n, df$pct)) + theme_classic() + ggtitle(sprintf("Number of metabolites at rho threshold (n=%d padj<0.05)", length(which(final$spearman_padj<0.05))))
print(p)
p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("Number of metabolites at rho threshold")) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(df), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
print(p)
df <- as.matrix(summary(final$ICC))
p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("rho by metabolite summary statistics")) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(df), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
print(p)
df <- sapply(thresholds, function(thresh) {
	c(length(which(final$ICC >= thresh)), length(which(final$rho >= thresh)))
}); df <- as.data.frame(t(df)); colnames(df) <- c("icc", "spearman"); df$threshold <- thresholds
df <- melt(df, id.vars=c("threshold")); colnames(df) <- c("threshold", "method", "n")
df$pct <- 100*(df$n / nrow(final))
p <- ggplot(df, aes(x=threshold, y=pct, color=method)) + geom_point() + geom_line() + geom_text(label=sprintf("%d (%.2g%%)", df$n, df$pct)) + theme_classic() + ggtitle(sprintf("Number of metabolites at ICC/rho threshold")) + scale_color_brewer(palette="Set1")
print(p)
# tables and distribution plots of highly reproducible metabolites
out <- subset(final, ICC>=0.9)
out <- merge(out, metabolon_map, by="row.names")
write.table(out, file="output/high_ICC_metabolites.txt", row.names=F, col.names=T, sep="\t")
for (metabolite in out$Row.names) {
	tmp <- rbind(data.plasma[,metabolite], data.dbs[,metabolite]); rownames(tmp) <- c("plasma", "DBS")
	df <- melt(tmp); colnames(df) <- c("SampleType", "PID", "value"); df$PID <- factor(df$PID)
	p <- ggplot(df, aes(x=SampleType, y=value, group=PID)) + geom_point() + geom_line() + theme_classic() + ggtitle(sprintf("%s in paired samples", metabolite))
	print(p)
}

## plot of mean metabolite values (axes are values in DBS and plasma), colored by ICC/Spearman rho
# maybe also draw ellipse for each point as SD of each metabolite?
p <- ggplot(final, aes(x=mean.plasma, y=mean.DBS, color=ICC)) + geom_point() + theme_classic() + ggtitle(sprintf("mean metabolite values and ICC")) + scale_color_gradient(low="black", high="green") + geom_vline(xintercept=0) + geom_hline(yintercept=0) + geom_abline(slope=1)
print(p)
p <- ggplot(final, aes(x=mean.plasma, y=mean.DBS, color=ICC)) + geom_point() + geom_errorbarh(aes(xmin=mean.plasma-sd.plasma, xmax=mean.plasma+sd.plasma), height=0.4) + geom_errorbar(aes(ymin=mean.DBS-sd.DBS, ymax=mean.DBS+sd.DBS), width=0.4) + theme_classic() + ggtitle(sprintf("mean metabolite values and ICC")) + scale_color_gradient(low="black", high="green")
print(p)
p <- ggplot(final, aes(x=mean.plasma, y=mean.DBS, color=rho)) + geom_point() + theme_classic() + ggtitle(sprintf("mean metabolite values and rho")) + scale_color_gradient(low="black", high="green") + geom_vline(xintercept=0) + geom_hline(yintercept=0) + geom_abline(slope=1)
print(p)
p <- ggplot(final, aes(x=mean.plasma, y=mean.DBS, color=rho)) + geom_point() + geom_errorbarh(aes(xmin=mean.plasma-sd.plasma, xmax=mean.plasma+sd.plasma), height=0.4) + geom_errorbar(aes(ymin=mean.DBS-sd.DBS, ymax=mean.DBS+sd.DBS), width=0.4) + theme_classic() + ggtitle(sprintf("mean metabolite values and rho")) + scale_color_gradient(low="black", high="green")
print(p)

## metabolite class distribution of highly reproducible metabolites
for (thresh in thresholds) {
	metlist <- rownames(subset(final, ICC>=thresh))
	if (length(metlist) > 0) {
		p <- plot_metabolite_breakdown(metlist, metabolon_map) + ggtitle(sprintf("Reproducible metabolites (ICC >= %.2g)", thresh))
		print(p)
	}
}
for (thresh in thresholds) {
	metlist <- rownames(subset(final, rho>=thresh))
	if (length(metlist) > 0) {
		p <- plot_metabolite_breakdown(metlist, metabolon_map) + ggtitle(sprintf("Reproducible metabolites (rho >= %.2g)", thresh))
		print(p)
	}
}

## metabolite class and ICC distribution of random forests selected features
featurelist <- {}
subtype <- "maternal"
mlevel <- "BIOCHEMICAL"
for (st in c("DBS", "Plasma")) {
	for (regi in c("untreated", "zdv", "PI-ART")) {
		tmp <- read.table(sprintf("output/randomForest_METABOLITE.%s.%s.%s.%s.features.txt", regi, subtype, mlevel, st), header=F, as.is=T, sep="\t", quote=""); colnames(tmp) <- c("feature", "importance")
		featurelist <- c(featurelist, tmp$feature)
	}
}
featurelist <- unique(featurelist)
df <- matrix(0, ncol=6, nrow=length(featurelist)); rownames(df) <- featurelist; colnames(df) <- c("Plasma - untreated", "DBS - untreated", "Plasma - zdv", "DBS - zdv", "Plasma - PI-ART", "DBS - PI-ART")
for (st in c("DBS", "Plasma")) {
	for (regi in c("untreated", "zdv", "PI-ART")) {
		tmp <- read.table(sprintf("output/randomForest_METABOLITE.%s.%s.%s.%s.features.txt", regi, subtype, mlevel, st), header=F, as.is=T, sep="\t", quote=""); colnames(tmp) <- c("feature", "importance")
		df[tmp$feature, sprintf("%s - %s", st, regi)] <- tmp$importance
	}
}
heatmap.2(df, Rowv=T, Colv=F, dendrogram="row", trace="none", col=colorRampPalette(c("white","blue"))(150), margins=c(6,18), cexCol=0.8, cexRow=0.6, main=sprintf("All RF features (n=%d)", nrow(df)))
heatmap.2(df, Rowv=T, Colv=T, dendrogram="both", trace="none", col=colorRampPalette(c("white","blue"))(150), margins=c(6,18), cexCol=0.8, cexRow=0.6, main=sprintf("All RF features (n=%d)", nrow(df)))
sel <- names(which(rowSums(df>0)>1)) # metabolites found in >1 model
df <- df>0; df.nonunique <- df[sel,] # convert to binary flag
sel <- c(names(which(df.nonunique[,"DBS - untreated"] & df.nonunique[,"Plasma - untreated"])), names(which(df.nonunique[,"DBS - zdv"] & df.nonunique[,"Plasma - zdv"])), names(which(df.nonunique[,"DBS - PI-ART"] & df.nonunique[,"Plasma - PI-ART"])))
heatmap.2(df.nonunique[sel,]+1-1, Rowv=T, Colv=F, dendrogram="row", trace="none", col=colorRampPalette(c("white","blue"))(150), margins=c(6,18), cexCol=0.8, cexRow=0.6, main=sprintf("RF features found consistently (n=%d)", length(sel)))
# metabolite class distributions
metlist <- rownames(df)
p <- plot_metabolite_breakdown(metlist, metabolon_map) + ggtitle(sprintf("All %d RF features", length(metlist)))
print(p)
metlist <- rownames(df.nonunique)
p <- plot_metabolite_breakdown(metlist, metabolon_map) + ggtitle(sprintf("%d RF features in >1 model", length(metlist)))
print(p)
metlist <- sel
p <- plot_metabolite_breakdown(metlist, metabolon_map) + ggtitle(sprintf("%d RF features found consistently", length(metlist)))
print(p)
out <- merge(df, metabolon_map, by="row.names"); out$Flag <- ifelse(out$Row.names %in% sel, "Consistent", "Not consistent")
rownames(out) <- out$Row.names; out <- out[, -1]; out <- merge(out, final, by="row.names", all.x=T)
test <- wilcox.test(ICC ~ Flag, out)
write.table(out, file="output/RF_features_in_multiple_models.txt", row.names=F, col.names=T, sep="\t", quote=F)
p <- ggplot(out, aes(x=Flag, y=ICC)) + geom_boxplot() + theme_classic() + ggtitle(sprintf("ICC by RF consistency (Wilcoxon p=%.4g)", test$p.value))
print(p)

## ICC and Spearman correlation in paired maternal plasma/DBS samples (stratified by Country)
subtype <- "maternal"
mapping.sel <- mapping; rownames(mapping.sel) <- mapping$patid
sel.metabolites <- rownames(subset(merged, FLAG.maternal_plasma & FLAG.maternal_DBS))
data <- matrix()
data.plasma <- df.metabolon[["Plasma"]][[mlevel]][as.character(mapping$patid), sel.metabolites]
data.dbs <- df.metabolon[["DBS"]][[mlevel]][as.character(mapping$patid), sel.metabolites]
data <- abind(data.plasma, data.dbs, along=0); dimnames(data)[[1]] <- c("Plasma", "DBS")
final <- {}
for (country in levels(mapping.sel$Country)) {
	mapping.country <- subset(mapping.sel, Country==country)
	data.country <- data[,rownames(mapping.country),]
	tmp <- t(apply(data.country, 3, function(x) {
		test <- try(cor.test(x[1,], x[2,], method="spearman"), silent=T)
		test.icc <- try(ICC(t(x)), silent=T)
		if (class(test) != "try-error" && class(test.icc) != "try-error") {
			unlist(c(country, mean(x[1,]), mean(x[2,]), sd(x[1,]), sd(x[2,]), mean(abs(x[1,]-x[2,])), test.icc$results["Single_raters_absolute", c("ICC", "p", "lower bound", "upper bound")], test$estimate, test$p.value))
		} else {
			unlist(c(country, mean(x[1,]), mean(x[2,]), sd(x[1,]), sd(x[2,]), mean(abs(x[1,]-x[2,])), NA, NA, NA, NA, NA, NA))
		}
	}))
	tmp <- as.data.frame(tmp)
	tmp$metabolite <- rownames(tmp)
	final <- rbind(final, tmp)
}
final <- as.data.frame(final); colnames(final) <- c("country", "mean.plasma", "mean.DBS", "sd.plasma", "sd.DBS", "mean_difference", "ICC", "p", "lower_bound", "upper_bound", "rho", "spearman_p", "metabolite")
final$p <- as.numeric(as.character(final$p)); final$spearman_p <- as.numeric(as.character(final$spearman_p)); final$mean.plasma <- as.numeric(as.character(final$mean.plasma)); final$mean.DBS <- as.numeric(as.character(final$mean.DBS)); final$ICC <- as.numeric(as.character(final$ICC)); final$rho <- as.numeric(as.character(final$rho))
final$padj <- p.adjust(final$p, method="fdr")
final$spearman_padj <- p.adjust(final$spearman_p, method="fdr")
final <- final[order(final$p), ]
write.table(final, file="output/ICC_by_Country.metabolite.txt", quote=F, sep="\t", row.names=T, col.names=T)
final <- subset(final, mean.plasma > -5 & mean.DBS > -5) # remove a few really low abundance outliers
p <- ggplot(final, aes(x=mean.plasma, y=ICC, color=country)) + geom_point() + stat_smooth(method="lm") + theme_classic() + ggtitle(sprintf("ICC vs mean.plasma (by metabolite, stratified by Country)")) + scale_color_brewer(palette="Set1")
print(p)
p <- ggplot(final, aes(x=mean.DBS, y=ICC, color=country)) + geom_point() + stat_smooth(method="lm") + theme_classic() + ggtitle(sprintf("ICC vs mean.DBS (by metabolite, stratified by Country)")) + scale_color_brewer(palette="Set1")
print(p)
p <- ggplot(final, aes(x=ICC, color=country)) + geom_density() + theme_classic() + ggtitle(sprintf("Distribution of ICC (by metabolite, stratified by Country)")) + scale_color_brewer(palette="Set1")
print(p)
summ <- aggregate(ICC~country, final, summary)
p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("ICC by metabolite+Country summary statistics")) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(t(summ)), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
print(p)

## ICC and Spearman correlation in paired maternal plasma/DBS samples (stratified by Group)
subtype <- "maternal"
mapping.sel <- mapping; rownames(mapping.sel) <- mapping$patid
sel.metabolites <- rownames(subset(merged, FLAG.maternal_plasma & FLAG.maternal_DBS))
data <- matrix()
data.plasma <- df.metabolon[["Plasma"]][[mlevel]][as.character(mapping$patid), sel.metabolites]
data.dbs <- df.metabolon[["DBS"]][[mlevel]][as.character(mapping$patid), sel.metabolites]
data <- abind(data.plasma, data.dbs, along=0); dimnames(data)[[1]] <- c("Plasma", "DBS")
final <- {}
for (gr in setdiff(levels(mapping.sel$Group), c("Case.other", "Control.other"))) {
	mapping.gr <- subset(mapping.sel, Group==gr)
	data.gr <- data[,rownames(mapping.gr),]
	tmp <- t(apply(data.gr, 3, function(x) {
		test <- try(cor.test(x[1,], x[2,], method="spearman"), silent=T)
		test.icc <- try(ICC(t(x)), silent=T)
		if (class(test) != "try-error" && class(test.icc) != "try-error") {
			unlist(c(gr, mean(x[1,]), mean(x[2,]), sd(x[1,]), sd(x[2,]), mean(abs(x[1,]-x[2,])), test.icc$results["Single_raters_absolute", c("ICC", "p", "lower bound", "upper bound")], test$estimate, test$p.value))
		} else {
			unlist(c(gr, mean(x[1,]), mean(x[2,]), sd(x[1,]), sd(x[2,]), mean(abs(x[1,]-x[2,])), NA, NA, NA, NA, NA, NA))
		}
	}))
	tmp <- as.data.frame(tmp)
	tmp$metabolite <- rownames(tmp)
	final <- rbind(final, tmp)
}
final <- as.data.frame(final); colnames(final) <- c("Group", "mean.plasma", "mean.DBS", "sd.plasma", "sd.DBS", "mean_difference", "ICC", "p", "lower_bound", "upper_bound", "rho", "spearman_p", "metabolite")
final$p <- as.numeric(as.character(final$p)); final$spearman_p <- as.numeric(as.character(final$spearman_p)); final$mean.plasma <- as.numeric(as.character(final$mean.plasma)); final$mean.DBS <- as.numeric(as.character(final$mean.DBS)); final$ICC <- as.numeric(as.character(final$ICC)); final$rho <- as.numeric(as.character(final$rho))
final$padj <- p.adjust(final$p, method="fdr")
final$spearman_padj <- p.adjust(final$spearman_p, method="fdr")
final <- final[order(final$p), ]
write.table(final, file="output/ICC_by_Group.metabolite.txt", quote=F, sep="\t", row.names=T, col.names=T)
final <- subset(final, mean.plasma > -5 & mean.DBS > -5) # remove a few really low abundance outliers
p <- ggplot(final, aes(x=mean.plasma, y=ICC, color=Group)) + geom_point() + stat_smooth(method="lm") + theme_classic() + ggtitle(sprintf("ICC vs mean.plasma (by metabolite, stratified by Group)")) + scale_color_brewer(palette="Set1")
print(p)
p <- ggplot(final, aes(x=mean.DBS, y=ICC, color=Group)) + geom_point() + stat_smooth(method="lm") + theme_classic() + ggtitle(sprintf("ICC vs mean.DBS (by metabolite, stratified by Group)")) + scale_color_brewer(palette="Set1")
print(p)
p <- ggplot(final, aes(x=ICC, color=Group)) + geom_density() + theme_classic() + ggtitle(sprintf("Distribution of ICC (by metabolite, stratified by Group)")) + scale_color_brewer(palette="Set1")
print(p)
p <- ggplot(final, aes(x=Group, y=ICC)) + geom_boxplot() + theme_classic() + ggtitle(sprintf("Distribution of ICC (by metabolite, stratified by Group)"))
print(p)
summ <- aggregate(ICC~Group, final, summary)
p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("ICC by metabolite+Group summary statistics")) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(t(summ)), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
print(p)

## ICC by subject
final <- apply(data, 2, function(x) {
	unlist(c(mean(x[1,]), mean(x[2,]), mean(abs(x[1,]-x[2,])), ICC(t(x))$results["Single_raters_absolute", c("ICC", "p", "lower bound", "upper bound")]))
#	icc(t(x))$value
}); final <- as.data.frame(t(final)); colnames(final) <- c("mean.plasma", "mean.DBS", "mean_difference", "ICC", "p", "lower_bound", "upper_bound")
final$padj <- p.adjust(final$p, method="fdr")
final <- final[order(final$p), ]
write.table(final, file="output/ICC.subject.txt", quote=F, sep="\t", row.names=T, col.names=T)
test <- cor.test(~ mean.plasma + ICC, final, method="spearman")
p <- ggplot(final, aes(x=mean.plasma, y=ICC)) + geom_point() + stat_smooth(method="lm") + theme_classic() + ggtitle(sprintf("ICC vs mean.plasma (by subject, Spearman rho=%.4g, p=%.4g)", test$estimate, test$p.value))
print(p)
test <- cor.test(~ mean.plasma + ICC, final, method="spearman")
p <- ggplot(final, aes(x=mean.DBS, y=ICC)) + geom_point() + stat_smooth(method="lm") + theme_classic() + ggtitle(sprintf("ICC vs mean.DBS (by subject, Spearman rho=%.4g, p=%.4g)", test$estimate, test$p.value))
print(p)
p <- ggplot(final, aes(x=ICC)) + geom_histogram(bins=50) + theme_classic() + ggtitle(sprintf("Distribution of ICC (by subject, mean=%.2g, median=%.2g)", mean(final$ICC), median(final$ICC))) + geom_vline(xintercept=mean(final$ICC), color="red")
print(p)
thresholds <- seq(from=0.05, to=0.95, by=0.05)
df <- sapply(thresholds, function(thresh) {
	length(which(final$ICC >= thresh))
}); df <- data.frame(threshold=thresholds, n=df); df$pct <- 100*(df$n / nrow(final))
p <- ggplot(df, aes(x=threshold, y=pct)) + geom_point() + geom_line() + geom_text(label=sprintf("%d (%.2g%%)", df$n, df$pct)) + theme_classic() + ggtitle(sprintf("Number of metabolites at ICC threshold (by subject, n=%d padj<0.05)", length(which(final$padj<0.05))))
print(p)
p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("Number of metabolites at ICC threshold (by subject)")) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(df), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
print(p)
df <- as.matrix(summary(final$ICC))
p <- qplot(1:10, 1:10, geom = "blank") + theme_bw() + ggtitle(sprintf("ICC by metabolite summary statistics (by subject)")) + theme(line = element_blank()) + annotation_custom(grob = tableGrob(df), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
print(p)
## tSNE by subject with both plasma and DBS
rownames(data.plasma) <- sprintf("%s.Plasma", rownames(data.plasma)); rownames(data.dbs) <- sprintf("%s.DBS", rownames(data.dbs))
data2 <- rbind(data.plasma, data.dbs)
tsne.out <- Rtsne(data2, perplexity=20)
df <- data.frame(PC1 = tsne.out$Y[,1], PC2 = tsne.out$Y[,2], SampleID=rownames(data2)); rownames(df) <- df$SampleID
df$SubjectID <- unlist(lapply(rownames(df), function(x) unlist(strsplit(x, "\\."))[1]))
df$SampleType <- unlist(lapply(rownames(df), function(x) unlist(strsplit(x, "\\."))[2]))
p <- ggplot(df, aes_string(x="PC1", y="PC2", colour="SampleType")) + geom_point() + geom_text(aes(label=SampleID), vjust="inward", hjust="inward") + theme_classic() + ggtitle(sprintf("tSNE by subject (%s)", mvar)) + scale_color_brewer(palette="Set1")
print(p)
for (mvar in c("Group", "Delivery")) {
	df[, mvar] <- mapping.sel[df$SubjectID, mvar]
	p <- ggplot(df, aes_string(x="PC1", y="PC2", colour=mvar)) + geom_point() + geom_text(aes(label=SampleID), vjust="inward", hjust="inward") + theme_classic() + ggtitle(sprintf("tSNE by subject (%s)", mvar)) + scale_color_brewer(palette="Set1")
	print(p)
}
## PCA
pca <- prcomp(data2, center=F, scale=F)
eigs <- pca$sdev^2
pvar <- 100*(eigs / sum(eigs))
df.pca <- df # from tSNE above
df.pca$PC1 <- pca$x[,1]; df.pca$PC2 <- pca$x[,2]
p <- ggplot(df.pca, aes_string(x="PC1", y="PC2", colour="SampleType")) + geom_point() + theme_classic() + ggtitle(sprintf("PCoA by SampleType (Euclidean distance)")) + xlab(sprintf("PC1 [%.1f%%]", pvar[1])) + ylab(sprintf("PC2 [%.1f%%]", pvar[2])) + scale_color_brewer(palette="Set1") + stat_ellipse(type="t")
print(p)

## PERMANOVA
res <- adonis2(data2 ~ SampleType + Group + SubjectID, data=df, permutations=999, method='euclidean')
sink(sprintf("output/metabolon_PERMANOVA_by_SampleType.txt"))
print(res)
sink()
## distance boxplots
distance_metric <- "euclidean"
dm <- as.matrix(dist(data2, method=distance_metric))
pair_df <- as.data.frame(t(combn(rownames(df), 2))); colnames(pair_df) <- c("s1", "s2")
pair_df$s1 <- as.character(pair_df$s1); pair_df$s2 <- as.character(pair_df$s2)
pair_df$id1 <- df[pair_df$s1, "SubjectID"]; pair_df$id2 <- df[pair_df$s2, "SubjectID"]
pair_df$Group <- ifelse(pair_df$id1 == pair_df$id2, "Within-Subject", "Between-Subject")
pair_df$dist <- dm[cbind(pair_df$s1, pair_df$s2)]
test <- wilcox.test(dist ~ Group, pair_df)
p <- ggplot(pair_df, aes(x=Group, y=dist)) + geom_boxplot() + theme_classic() + ggtitle(sprintf("%s distance comparison (Wilcoxon p=%.4g)", distance_metric, test$p.value))
print(p)

## differential abundance
# t-test/Wilcoxon
subtype <- "maternal"
mapping.sel <- mapping; rownames(mapping.sel) <- mapping$patid
sel.metabolites <- rownames(subset(merged, FLAG.maternal_plasma & FLAG.maternal_DBS))
data <- matrix()
data.plasma <- df.metabolon[["Plasma"]][[mlevel]][as.character(mapping$patid), sel.metabolites]
data.dbs <- df.metabolon[["DBS"]][[mlevel]][as.character(mapping$patid), sel.metabolites]
data <- abind(data.plasma, data.dbs, along=0); dimnames(data)[[1]] <- c("Plasma", "DBS")
final <- apply(data, 3, function(x) {
	test <- t.test(x[1,], x[2,])
	test2 <- wilcox.test(x[1,], x[2,])
	unlist(c(mean(x[1,]), mean(x[2,]), sd(x[1,]), sd(x[2,]), mean(abs(x[1,]-x[2,])), test$statistic, test$p.value, test2$statistic, test2$p.value))
}); final <- as.data.frame(t(final)); colnames(final) <- c("mean.plasma", "mean.DBS", "sd.plasma", "sd.DBS", "mean_difference", "t", "p", "W", "p.wilcox")
final$padj <- p.adjust(final$p, method="fdr")
final$padj.wilcox <- p.adjust(final$p.wilcox, method="fdr")
final <- final[order(final$p), ]
write.table(final, file="output/DA.metabolite.txt", quote=F, sep="\t", row.names=T, col.names=T)
metlist <- rownames(subset(final, padj<siglevel))
p <- plot_metabolite_breakdown(metlist, metabolon_map) + ggtitle(sprintf("DA metabolites (t-test padj<%.2g)", siglevel))
print(p)
metlist <- rownames(subset(final, padj.wilcox<siglevel))
p <- plot_metabolite_breakdown(metlist, metabolon_map) + ggtitle(sprintf("DA metabolites (Wilcoxon padj<%.2g)", siglevel))
print(p)


dev.off()





