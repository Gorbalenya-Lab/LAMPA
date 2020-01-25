
split_seq_by_len <- function(noF, ap_len, f1, f2) {

	ap_len <- as.numeric(ap_len)


	# get regions w/o hits
	tab <- read.table(noF, stringsAsFactors=FALSE)

	tab[,3] <- tab[,2] - tab[,1] + 1
	idx <- which(tab[,3] > ap_len)
	tab <- tab[idx,]


	# split regions w/o hits
	if (nrow(tab) == 0) {

		# exit if no suitable regions
		print("No regions without hits exceeding ap_len length threshold, hence no Average-protein-size-specific iterations will be conducted.")

	} else {

		# split starting from N-terminus with and without inset
		for (mode in c('N-term', 'Inset')) {

			DF <- data.frame(NULL, stringsAsFactors=FALSE)

			for (i in 1:nrow(tab)) {

				from <- tab[i,1] + ifelse(mode=='N-term', 0, (ap_len%/%2))
				to   <- tab[i,2]
				brks <- seq(from=from-1, to=to, by=ap_len)

				if (length(brks) > 1) {

					# fragments of size ap_len
					df <- data.frame(start=head(brks,-1)+1, end=tail(brks,-1), stringsAsFactors=FALSE)

					# extend C-terminal fragment to include the remaining part of the region under consideration
					if ( ((to - df$end[nrow(df)]) < (ap_len%/%2)) && (df$start[nrow(df)] != from) ) { df$end[nrow(df)] <- to }

					# store info about the deliniated fragments
					DF <- rbind(DF, df)

				}

			}

			if (nrow(DF) > 0) {
				write.table(DF, file=ifelse(mode=='N-term', f1, f2), sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)
			}
		}
	}
}
