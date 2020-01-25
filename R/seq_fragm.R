
seq_fragm <- function(seqF, from, to, fragmF) {

	from <- as.numeric(from)
	to   <- as.numeric(to)

	whole_seq <- seqinr::read.fasta(seqF, seqtype='AA')
	fragm_seq <- lapply(whole_seq, function (v) v[from:to])

	nm <- names(fragm_seq)
	nm[length(nm)] <- paste0(nm[length(nm)], '_', from, '-', to)

	seqinr::write.fasta(sequences=fragm_seq, names=nm, file.out=fragmF)

}
