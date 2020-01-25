
run_hhsearch <- function(d, q_id, seqF, cooF, DBs, hhmakepath, hhsearchpath, cpu) {

	coo <- read.table(cooF, header=FALSE, stringsAsFactors=FALSE)

	# for each query sequence fragment ...
	for (i in 1:nrow(coo)) {

		from <- coo[i,1]
		to  <- coo[i,2]


		# file names
		baseF <- paste0(q_id, "_", from, "-", to)
		fragmF <- paste0(baseF, ifelse(grepl('\\.a3m$', seqF), '.a3m', '.fasta'))
		hhmF <- paste0(baseF, ".hhm")


		# prepare file with query sequence fragment
		seq_fragm(seqF, from, to, fragmF)


		# build query profile
		system(paste(hhmakepath, "-M first -name", baseF, "-i", fragmF, "-o", hhmF, "-v 0"))


		# for each target database ...
		for (db_name in names(DBs)) {

			# run HHsearch
			hhrF <- paste0(baseF, "_VS_", db_name, ".hhr")
			system(paste(hhsearchpath, "-i", hhmF, "-d", DBs[[db_name]], "-o", hhrF, "-cpu", cpu, "-v 0"))

			# extract hits table & alignments
			tabF <- paste0(baseF, "_VS_", db_name, "_hits.tsv")
			aliF <- paste0(baseF, "_VS_", db_name, "_alignments.rds")
			parse_hhr(db_name, hhrF, tabF, aliF)
		}


		# make cumulative hits table (one query sequence fragment, multiple databases)
		cmF <- paste0(baseF, "_cumulative_hits.tsv")

		pattern <- paste0(baseF, "_VS_[^\\s]+_hits.tsv")
		tabFs <- list.files(path=d, pattern=pattern)

		make_cumulative_hits_table(q_id, from, to, tabFs, cmF)
	}
}
