
split_seq_by_hits <- function(clustF, qp_len, cooF) {

	DF <- data.frame(NULL, stringsAsFactors=FALSE)


	if (file.exists(clustF)) {

		CL <- read.table(clustF, header=TRUE, stringsAsFactors=FALSE)

		for (x in unique(CL$q_from)) { # for each region, to which there were hits...

			idx <- which(CL$q_from==x)
			from <- CL$q_from[idx][1]
			to <- CL$q_to[idx][1]

			# clusters of overlapping hits, belonging to the region
			ir <- IRanges::IRanges(start=CL$cl_from[idx], end=CL$cl_to[idx])

			# gaps between clusters of hits
			gaps <- IRanges::gaps(ir, start=from, end=to)
			gaps <- as.data.frame(gaps)
			gaps <- gaps[gaps$width >= qp_len, c('start','end')]

			# store regions between clusters of hits
			DF <- rbind(DF, gaps)
		}
	}


	if (nrow(DF)>0) {
		write.table(DF, file=cooF, sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)
	}
}
