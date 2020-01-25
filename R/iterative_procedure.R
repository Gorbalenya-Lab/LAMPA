
iterative_procedure <- function(seqF, DBs, hhmakepath, hhsearchpath, addsspath, tmhmmpath, tm_gap, qp_len, ap_len, P, E, L, cpu) {

	# define qp_len
	if (qp_len < (L+1)) { qp_len <- L+1 }


	# save path to working directory
	work_dir <- getwd()


	# extract query sequence ID from FASTA file
	q_id <- attr(seqinr::read.fasta(seqF, seqtype="AA")[[1]], which="name")


	# predict transmembrane regions
	tmF_coo <- NULL
	if (!is.null(tmhmmpath)) {
		tmF_raw <- paste0(work_dir, "/", q_id, "_TMHMM.txt")		# file for raw TMHMM output
		tmF_coo <- paste0(work_dir, "/", q_id, "_TM_cl_coo.tsv")	# file for coordinates of TM helices clusters
		tmF_sum <- paste0(work_dir, "/", q_id, "_TM.tsv")		# file for TM prediction summary
		predict_tm_regions(seqF, tmhmmpath, tm_gap, tmF_raw, tmF_coo, tmF_sum)
	}


	# predict secondary structure
	if (!is.null(addsspath)) {
		a3mF <- paste0(work_dir, "/", q_id, ".a3m")
		system(paste(addsspath, seqF, a3mF, "-fas"))
		seqF <- a3mF
	}


	# prepare file with query coordinates for the first iteration, query = whole protein
	whole_pr_cooF <- paste0(work_dir, "/", q_id, "_query_regions_coo_iteration_1.tsv")
	whole_pr_coo(seqF, whole_pr_cooF)
	cooF <- whole_pr_cooF


	# file for coordinates of regions w/o hits
	noF <- paste0(work_dir, "/", q_id, "_regions_without_hits.tsv")


	# main loop
	I = 1
	STATUS = "RUN"
	while (TRUE) {

		# if iteration w/o significant hits is reached, prepare for Average-protein-size-specific iterations
		if (!file.exists(cooF) && (STATUS == "RUN")) {

			if (file.exists(noF)) {
				f1 <- paste0(work_dir, "/", q_id, "_query_regions_coo_iteration_", I, ".tsv")
				f2 <- paste0(work_dir, "/", q_id, "_query_regions_coo_iteration_", (I + 1), ".tsv")
				split_seq_by_len(noF, ap_len, f1, f2)
			}
			else {
				print("No regions without hits left, hence no Average-protein-size-specific iterations will be conducted.")
			}
			STATUS = "END"
		}


		# if final iteration of the procedure is reached, summarize results, clean up and exit
		if (!file.exists(cooF) && (STATUS == "END")) {

			plotF <- paste0(work_dir, "/", q_id, "_annotation_plot.pdf")
			tabF <- paste0(work_dir, "/", q_id, "_annotation_table.tsv")
			plot_and_table((I-1), q_id, whole_pr_cooF, tmF_coo, plotF, tabF)

			toMove <- list.files(work_dir)
			toMove <- grep("_TM.tsv$|_hits_cluster_|_annotation_table.tsv$|_annotation_plot.pdf$", toMove, value=TRUE, invert=TRUE)
			toDir <- paste0(work_dir, "/utility_data/")
			dir.create(toDir)
			file.copy(from=toMove, to=toDir, recursive=TRUE)
			unlink(toMove, recursive=TRUE)

			break
		}


		# iteration folder
		d <- paste0(work_dir, "/iteration_", I)
		dir.create(d)


		# run HHsearch
		setwd(d)
		run_hhsearch(d, q_id, seqF, cooF, DBs, hhmakepath, hhsearchpath, cpu)
		setwd(work_dir)


		# cluster overlapping hits
		clustF <- paste0(work_dir, "/", q_id, "_clusters_of_hits_iteration_", I, ".tsv")
		prefix <- paste0(work_dir, "/", q_id, "_hits_cluster_")	# prefix for files with info about each individual cluster

		if (I == 1) { iterat_type <- 1 }
		if ((I > 1) && (STATUS == "RUN")) { iterat_type <- 2 }
		if ((I > 1) && (STATUS == "END")) { iterat_type <- 3 }
		cluster_hits(d, P, E, L, noF, clustF, prefix, iterat_type, I)


		# prepare coords file for the next iteration
		cooF <- paste0(work_dir, "/", q_id, "_query_regions_coo_iteration_", (I+1), ".tsv")

		if (STATUS == "RUN") {
			if (I == 1) {
				split_seq_1st_iterat(whole_pr_cooF, noF, tmF_coo, clustF, qp_len, cooF)
			}
			else {
				split_seq_by_hits(clustF, qp_len, cooF)
			}
		}


		# i++
		I = I + 1
	}
}
