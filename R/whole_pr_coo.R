
whole_pr_coo <- function(seqF, cooF) {

	# input seq len
	seq_len <- length(seqinr::read.fasta(seqF, seqtype='AA')[[1]])

	# data frame with coords
	df <- data.frame(start=1, end=seq_len)

	# output
	write.table(df, file=cooF, sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)
}
