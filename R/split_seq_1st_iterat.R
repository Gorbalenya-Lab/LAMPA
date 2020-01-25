
split_seq_1st_iterat <- function(whole_pr_cooF, noF, tmF_coo, clustF, qp_len, cooF) {

	# whole protein coordinates
	from <- 1
	to <- read.table(whole_pr_cooF, stringsAsFactors=FALSE)[1,2]


	# clusters of TM helices
	tm_bool <- ((!is.null(tmF_coo)) && file.exists(tmF_coo))
	if (tm_bool) {
		tab <- read.table(tmF_coo, header=TRUE, stringsAsFactors=FALSE)
		tm_ir <- IRanges::IRanges(start=tab$tm_cl_from, end=tab$tm_cl_to)
	}


	# clusters of hits
	hit_bool <- file.exists(clustF)
	if (hit_bool) {
		CL <- read.table(clustF, header=TRUE, stringsAsFactors=FALSE)
		hit_ir <- IRanges::IRanges(start=CL$cl_from, end=CL$cl_to)
	}


	# gaps between clusters of hits and/or TM helices
	if (tm_bool || hit_bool) {

		if (tm_bool && hit_bool) { ir <- IRanges::union(tm_ir, hit_ir) }

		if (tm_bool && (!hit_bool)) { ir <- tm_ir; unlink(noF) }

		if ((!tm_bool) && hit_bool) { ir <- hit_ir }


		gaps <- IRanges::gaps(ir, start=from, end=to)
		gaps <- as.data.frame(gaps)
		gaps <- gaps[gaps$width >= qp_len, c('start','end')]


		if (nrow(gaps)>0) {
			write.table(gaps, file=cooF, sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)
		}	
	}
}
