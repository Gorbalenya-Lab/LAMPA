
predict_tm_regions <- function(seqF, tmhmmpath, tm_gap, tmF_raw, tmF_coo, tmF_sum) {

	tm_gap <- as.numeric(tm_gap)

	# predict
	system(paste(tmhmmpath, seqF, ">", tmF_raw))

	# TM helices coordinates
	tab <- read.table(tmF_raw, stringsAsFactors=FALSE)
	idx <- which(tab[,3]=='TMhelix')

	if (length(idx) > 0) {

		tab <- tab[idx,4:5]

		# cluster TM helices, separated by distance < tm_gap
		tm_ir <- IRanges::IRanges(start=tab[,1], end=tab[,2])
		tm_cl <- IRanges::reduce(tm_ir, min.gapwidth=tm_gap)

		# write output tables
		coo <- as.data.frame( tm_cl )[, c('start','end')]
		colnames(coo) <- c('tm_cl_from', 'tm_cl_to')
		coo$tm_cl_index <- 1:nrow(coo)
		write.table(coo, file=tmF_coo, sep='\t', quote=FALSE, row.names=FALSE)

		colnames(tab) <- c('tm_helix_from', 'tm_helix_to')
		tab$tm_cl_index <- as.data.frame( IRanges::findOverlaps(tm_ir, tm_cl) )$subjectHits
		write.table(tab, file=tmF_sum, sep='\t', quote=FALSE, row.names=FALSE)

	}
}
