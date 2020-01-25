
make_cumulative_hits_table <- function(q_id, from, to, tabFs, cmF) {

	from <- as.numeric(from)
	to   <- as.numeric(to)


	# combine tables extracted from *.hhr files
	DF <- data.frame(NULL, stringsAsFactors=FALSE)
	for (f in tabFs) {
		df <- read.table(f, header=TRUE, stringsAsFactors=FALSE)
		DF <- rbind(DF, df)
	}


	# query fragment coordinates in the query sequence
	DF$q_id   <- q_id
	DF$q_from <- from
	DF$q_to   <- to
	DF$q_len  <- to - from + 1


	# hit coordinates in the query sequence
	L <- strsplit(DF$QueryHMM,'-')
	L <- lapply(L, as.numeric)

	DF$h_from <- sapply(L, function (v) v[1]) + from - 1
	DF$h_to   <- sapply(L, function (v) v[2]) + from - 1
	DF$h_len  <- DF$h_to - DF$h_from + 1


	# order hits in decreasing order of their Probability, increasing order of their E-value
	DF <- DF[ order(-DF$Prob, DF$E_value), ]


	# write cumulative table
	DF <- DF[,c('q_id', 'q_from', 'q_to', 'q_len', 'DB', 'Hit', 'Prob', 'E_value', 'Score', 'h_from', 'h_to', 'h_len', 'TemplateHMM', 'aliF', 'aliNo')]
	write.table(DF, file=cmF, sep='\t', quote=FALSE, row.names=FALSE)
}
