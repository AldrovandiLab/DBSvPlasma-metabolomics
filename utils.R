multiplot <- function(..., plotlist=NULL, file, cols=1, rows=1) {
  require(grid)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
	
	i = 1
	while (i < numPlots) {
		numToPlot <- min(numPlots-i+1, cols*rows)
		# Make the panel
		# ncol: Number of columns of plots
		# nrow: Number of rows needed, calculated from # of cols
		layout <- matrix(seq(i, i+cols*rows-1), ncol = cols, nrow = rows, byrow=T)
		if (numToPlot==1) {
		  print(plots[[i]])
		} else {
		  # Set up the page
		  grid.newpage()
		  pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
		  # Make each plot, in the correct location
		  for (j in i:(i+numToPlot-1)) {
		    # Get the i,j matrix positions of the regions that contain this subplot
		    matchidx <- as.data.frame(which(layout == j, arr.ind = TRUE))
		    print(plots[[j]], vp = viewport(layout.pos.row = matchidx$row,
		                                    layout.pos.col = matchidx$col))
		  }
		}
		i <- i+numToPlot
  }
}

normalizeByRows <- function (df, rsum=1)
{
	while (any(abs((rowSums(df)-rsum))>1e-13)) {
		df <- rsum*(df / rowSums(df))
	}
	return(df)
}
normalizeByCols <- function (df, csum=1, level=NULL, delim="\\|")
{
	if (is.null(level)) {
		while (any(abs((colSums(df)-csum))>1e-13 & colSums(df)!=0, na.rm=T)) {
			missing <- which(colSums(df)==0)
			df <- sweep(df, 2, colSums(df)/csum, "/")
			df[,missing] <- 0
		}
	} else {
	 tmp <- df
	 tmp$taxa <- rownames(tmp)
	 tmp$splitter <- factor(unlist(lapply(rownames(tmp), function(x) unlist(strsplit(x, delim))[level])))
	 names <- rownames(tmp)[order(tmp$splitter)]
	 tmp <- ddply(tmp, .(splitter), function(x) {
	 		x <- x[, setdiff(colnames(x), c("taxa", "splitter"))]
			while (any(abs((colSums(x)-csum))>1e-13 & colSums(df)!=0, na.rm=T)) {
				x <- sweep(x, 2, colSums(x)/csum, "/")
			}
			x
		})
		rownames(tmp) <- names
		df <- tmp[, setdiff(colnames(tmp), "splitter")]
	}
	return(df)
}

renameLevelsWithCounts <- function(fvec, originalLevelsAsNames=FALSE) {
	tab <- table(fvec)
	retval <- sprintf("%s (n=%d)", fvec, tab[unlist(lapply(fvec, function(x) match(x, names(tab))))])
#	newlevels <- sprintf("%s (n=%d)", levels(fvec), tab[levels(fvec)])
	newlevels <- sprintf("%s (n=%d)", levels(fvec), tab[unlist(lapply(names(tab), function(x) which(levels(fvec)==x)))])
	retval <- factor(retval, levels=newlevels)
	if (originalLevelsAsNames) {
		names(retval) <- fvec
	}
	return(retval)
}

# get selected taxonomy
# DEFAULT: do not include entries that are not classified to the requested level (return NA instead)
getTaxonomy <- function(otus, tax_tab, level, na_str = c("unclassified", "unidentified", "NA", ""), includeUnclassified = FALSE) {
	# coerce taxonomyTable to data.frame (see https://github.com/joey711/phyloseq/issues/983), otherwise hangs as of phyloseq v1.26
	if (class(tax_tab) == "taxonomyTable") {
		tax_tab <- as.data.frame(tax_tab@.Data)
	}
	ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
	sel <- ranks[1:match(level, ranks)]
	inds <- apply(tax_tab[otus,sel], 1, function(x) max(which(!(x %in% na_str | is.na(x)))))
	if (includeUnclassified) {
		retval <- as.data.frame(tax_tab)[cbind(otus, ranks[inds])]
		inds.nmatch <- inds!=match(level, ranks)
		retval[inds.nmatch] <- paste(na_str[1], retval[inds.nmatch], sep=" ")
		if (level == "Species") {
			tmp <- as.data.frame(tax_tab)[cbind(otus, "Genus")]
			retval[!inds.nmatch] <- sprintf("%s %s", tmp[!inds.nmatch], retval[!inds.nmatch])
		}
	} else {
		retval <- as.data.frame(tax_tab)[cbind(otus, level)]
		if (level == "Species") {
			tmp <- as.data.frame(tax_tab)[cbind(otus, "Genus")]
			retval[!is.na(retval)] <- sprintf("%s %s", tmp[!is.na(retval)], retval[!is.na(retval)])
		}
	}
	retval <- gsub("\\[|\\]", "", retval)
	return(retval)
}


# select best available model from list of glmmTMB models, or return NA if none are valid
# returns a one-element list with names as provided
selectBestModel <- function(models) {
	
	valid <- lapply(models, class) != "try-error"
	AICs <- lapply(models, function(x) {
		ifelse(class(x)=="try-error", NA, AIC(x))
	})
	sel <- which.min(AICs)
	# if there is a valid minimum AIC model, return it; otherwise return NA
	if (length(sel)==1 & valid[sel]) {
		retval <- models[sel]
	} else {
		retval <- NA
	}
	
	return(retval)
}

# stratified sampling from a factor vector
sample.stratified <- function(groups, pct=0.75) {
	retval <- sapply(levels(groups), function(strata) {
		x <- which(groups == strata)
		sample(x, size=ceiling(pct*length(x)))
	})
	unlist(retval)
}

# plot pie chart of BIOCHEMICAL classes
plot_metabolite_breakdown <- function(feature_list, metabolon_map) {
	cols.superpathway <- brewer.pal(9, "Set1"); names(cols.superpathway) <- c("Amino Acid", "Carbohydrate", "Cofactors and Vitamins", "Energy", "Lipid", "Nucleotide", "Partially Characterized Molecules", "Peptide", "Xenobiotics")
	df <- melt(table(metabolon_map[feature_list, "SUPER.PATHWAY"]))
	colnames(df) <- c("Class", "count")
	p <- ggplot(df, aes(x=factor(1), y=count, fill=Class)) + geom_bar(stat="identity", color="black", width=1) + geom_text(aes(label=count), position=position_stack(vjust=0.5)) + theme_classic() + coord_polar(theta="y") + theme(axis.ticks=element_blank(), axis.text.x=element_blank(), axis.title=element_blank(), axis.line=element_blank()) + scale_fill_manual(values=cols.superpathway)
	p
}


scale_by_group <- function(vec, groups) {
	for (gr in unique(groups)) {
		inds <- which(groups==gr)
		vec[inds] <- vec[inds] / median(vec[inds])
	}
	vec
}
	
