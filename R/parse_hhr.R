
parse_hhr <- function (db_name, hhrF, tabF, aliF) {

	# read *.hhr file into vector (one line - one element of the vector)
	con <- file(hhrF, open='r')
	lin <- readLines(con)
	close(con)


	# TABLE
	# table with hits is preceded and followed by the first two empty lines in *.hhr file
	idxT <- which( grepl('^\\s*$', lin) )[1:2]
	linT <- lin[ (idxT[1]+1):(idxT[2]-1) ]


	# extracted table requires complex parsing, because its second and last columns ('Hit' and 'Template HMM') may or may not contain spaces
	L <- lapply(linT[2:length(linT)], function (x) {

						templ_coo <- x
						templ_coo <- sub('^.+\\s(\\d+-\\d+\\s*\\(\\d+\\))\\s*$', '\\1', templ_coo)
						x <- sub(templ_coo, '', x, fixed=TRUE)
						templ_coo <- gsub(' ', '', templ_coo)

						v <- strsplit(x,'\\s+')[[1]]
						if (v[1]=='') { v <- v[-1] }
						if (v[length(v)]=='') { v <- v[-length(v)] }
						l <- length(v)

						hit_name  <- v[2]

						r <- c(db_name, hit_name, v[(l-6):l], templ_coo, aliF, v[1])
						return(r)
	})

	df <- data.frame( matrix(unlist(L), nrow=length(L), byrow=TRUE), stringsAsFactors=FALSE)

	linT[1] <- sub('^\\s*', '', linT[1])
	colnames(df) <- c('DB', strsplit(linT[1],'\\s+')[[1]][2:8], 'QueryHMM', 'TemplateHMM', 'aliF', 'aliNo')
	colnames(df)[colnames(df)=='E-value'] <- 'E_value'


	# write table extracted from *.hhr file
	write.table(df, file=tabF, sep='\t', quote=FALSE, row.names=FALSE)


	# ALIGNMENTS
	idxA <- which( grepl('^No\\s\\d+\\s*$', lin) )

	A <- list()

	for (j in seq_along(idxA)) {

		no <- sub('^No\\s(\\d+)\\s*$', '\\1', lin[idxA[j]])

		from <- idxA[j] + 1
		to <- ifelse(j==length(idxA), length(lin), idxA[j+1]) - 1

		A[[no]] <- lin[from:to]
	}

	saveRDS(A, file=aliF)

}
