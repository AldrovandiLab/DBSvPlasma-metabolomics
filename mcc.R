# X is indicator matrix of predictions, Y is indicator matrix of truth
# columns are classes, rows are samples
mcc <- function(preds=NULL, actuals=NULL, x=NULL, y=NULL) {
	# if preds and actuals are provided, x and y will be ignored
	if (!is.null(preds)) {
		nclasses <- length(union(preds, actuals))
		x <- matrix(0, nrow=length(preds), ncol=nclasses)
		y <- matrix(0, nrow=length(actuals), ncol=nclasses)
		x[cbind(1:nrow(x), preds+1)] <- 1
		y[cbind(1:nrow(y), actuals+1)] <- 1
	}
	if (!all(dim(x) == dim(y))) {
		stop("X and Y must have the same dimensions")
	}
	
	cov_biased <- function(x, y) {
		sum(sapply(1:ncol(x), function(k) {
			cov(x[,k], y[,k]) # unbiased estimate with (n-1) denominator as opposed to (n), but cancels out anyways so identical result
		}))
	}
	numerator <- cov_biased(x,y)
	denominator <- sqrt(cov_biased(x,x) * cov_biased(y,y))
	numerator / denominator	
}
