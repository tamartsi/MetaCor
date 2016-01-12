check.column.order.input.MetaCor <-
function(beta, var, cov){
	if (!is.matrix(cov)) stop("cov should be a matrix. Try cov[,, drop = F]")
	S <- ncol(beta)
	strata.names <- sapply(colnames(beta), function(x){
		strsplit(x, split = "Beta.", fixed = TRUE)[[1]][2]
	})
	
	if (!all(colnames(var) == paste0("var.", strata.names))) return(list(good.input = FALSE, message = "column names of var matrix are not var.strat_name"))
	
	cov.colnames.should.be <- rep("", (S-1)*S/2)
	cur.ind <- 1
	for (i in 1:(S-1)){
		cov.colnames.should.be[cur.ind:(cur.ind + S - i -1)] <- paste0("cov.", strata.names[i], ":", strata.names[(i+1):S])
		cur.ind <- cur.ind + S - i
	}
	
	if (!all(colnames(cov) == cov.colnames.should.be)) return(list(good.input = FALSE, message = "column names of cov matrix are not cov.strat1_name:strat2_name"))
	
	return(list(good.input = TRUE))
}
