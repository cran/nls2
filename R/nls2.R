
nls2 <- function(formula, data = parent.frame(), start, control = nls.control(),
	algorithm = c("default", "plinear", "port", "brute-force"), ...,
	all = FALSE) { 

	L <- (list(formula = formula, data = data, control = control))
	if (!missing(start)) L$start <- start
	L <- append(L, list(...))
	algorithm <- match.arg(algorithm)
	if (algorithm == "brute-force") {
	   nls <- function(formula, data, start, ...) {
	      nlsModel <- stats:::nlsModel
	      environment(nlsModel) <- environment()
	      #  disable nlsModel gradient error
	      stop <- function(...) {
	        msg <- "singular gradient matrix at initial parameter estimates"
	        if (list(...)[[1]] == msg) return()
	        stop(...)
	      }
	      structure(list(m = nlsModel(formula, data, start), 
	         call = list(algorithm = algorithm), 
	         convInfo = list(isConv = TRUE, finIter = 0, finTol = 0)), 
	         class = "nls")
	   }
	} else L$algorithm <- algorithm

	if (missing(start)) return(do.call(nls, L))
	else L$start <- as.data.frame(as.list(start))

	if (NROW(L$start) == 1) return(do.call(nls, L))

	result <- apply(L$start, 1, function(start) {
		L$start <- start
		do.call(nls, L)
	})
	if (all) result else {
		ss <- lapply(result, function(x) sum(resid(x)^2))
		result[[which.min(ss)]]
	}
}

