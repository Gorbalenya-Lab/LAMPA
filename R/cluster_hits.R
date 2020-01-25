
cluster_hits <- function(d, P, E, L, noF, clustF, prefix, iterat_type, iterat_num) {

	P <- as.numeric(P)
	E <- as.numeric(E)
	L <- as.numeric(L)


	# join all cumulative tables from the iteration into one
	t <- data.frame(NULL, stringsAsFactors=FALSE)


	cmFs <- list.files(path=d, pattern='_cumulative_hits.tsv', full.names=TRUE)
	for (f in cmFs) {

		# read cumulative table
		cm <- read.table(f, header=TRUE, colClasses=c('aliNo'='character'), stringsAsFactors=FALSE)

		# filter hits
		sele <- which( (cm$Prob > P)   &   (cm$E_value < E)   &   (cm$h_len > L) )

		# no hits on region under consideration, satisfying the criteria? print region's coordinates to noF!
		if (length(sele)==0) {
			cat(cm$q_from[1], '\t', cm$q_to[1], '\n', sep='', file=noF, append=TRUE)
		} else {
			t <- rbind(t, cm[sele,])
		}
	}


	# cluster overlapping hits
	CL <- data.frame(NULL, stringsAsFactors=FALSE)

	ir <- IRanges::IRanges(start=t$h_from, end=t$h_to)
	clust <- IRanges::reduce(ir, min.gapwidth=0)

	# to which cluster each hit belongs?
	hit2clust <- as.data.frame( IRanges::findOverlaps(ir, clust) )$subjectHits

	# save hits constituting each cluster as a table
	for (n in unique(hit2clust)) {

		cl_index <- paste0(iterat_num,'.',n)
		tab <- data.frame(iterat_num, iterat_type, cl_index, stringsAsFactors=FALSE)

		tab <- cbind(tab, clust[n,])
		colnames(tab)[4:6] <- c('cl_from', 'cl_to', 'cl_len')

		idx <- which(hit2clust==n)

		# alignments
		lin <- c()
		for (i in idx) {
			f <- paste0(d, '/', t$aliF[i])
			A <- readRDS(f)
			a <- A[[ t$aliNo[i] ]]
			lin <- c(lin, a)
		}
		con <- file(paste0(prefix, cl_index, '_alignments.txt'), 'w')
		writeLines(lin, con=con)
		close(con)

		# table
		tab <- cbind(tab, t[idx,])
		q_cols <- c('q_id', 'q_from', 'q_to', 'q_len')
		tab <- tab[,c(q_cols, colnames(tab)[!colnames(tab)%in%q_cols])]
		tab$aliF  <- NULL
		tab$aliNo <- NULL
		write.table(tab, file=paste0(prefix, cl_index, '.tsv'), sep='\t', quote=FALSE, row.names=FALSE)

		# top hit in cluster
		CL <- rbind(CL, tab[1,])
	}


	# save top hits in clusters as a table
	if (nrow(CL) > 0) {
		CL <- CL[order(CL$cl_from),]
		write.table(CL, file=clustF, sep='\t', quote=FALSE, row.names=FALSE)
	}
}
