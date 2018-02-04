# function takes two arguments
# abs is sample x wavelength matrix of absorbances
# conc is sample x analyte matrix of concentration
# function returns one result
# eb is analyte x wavelength matrix of predicted eb values
# as.matrix assures dataframe is treated as matrix
findeb = function(abs, conc){
	abs.m = as.matrix(abs)
	conc.m = as.matrix(conc)
	ct = t(conc.m)
	ctc = ct %*% conc.m
	invctc = solve(ctc)
	eb = invctc %*% ct %*% abs.m
	output = eb
	invisible(output)
}

# function takes two arguments
# abs is sample x wavelength matrix of absorbances
# eb is analyte x wavelength matrix of eb values
# function returns one result
# pred.conc is sample x analyte matrix of predicted concentrations
# as.matrix assures dataframe is treated as matrix
# round trucates all values to four decimal place
findconc = function(abs, eb){
	abs.m = as.matrix(abs)
	eb.m = as.matrix(eb)
	ebt = t(eb.m)
	ebebt = eb %*% ebt
	invebebt = solve(ebebt)
	pred_conc = round(abs.m %*% ebt %*% invebebt, digits = 5)
	output = pred_conc
	invisible(output)
}
