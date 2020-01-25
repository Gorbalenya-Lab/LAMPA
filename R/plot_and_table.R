
plot_and_table <- function(I, q_id, whole_pr_cooF, tmF_coo, plotF, tabF) {

	# preparation
	I <- as.numeric(I)
	seq_len <- read.table(whole_pr_cooF, stringsAsFactors=FALSE)[1,2]
	tm_bool <- (!is.null(tmF_coo))


	# initialize plot
	pdf(plotF, height=pdf.options()$height*0.15*(I+tm_bool), width=pdf.options()$width*2)
	par(mar=c(1,0,5,0))
	plot(NULL, xlim=seq_len*c(-0.1, 1.05), ylim=c(I,-1*tm_bool), ann=FALSE, xaxt='n', yaxt='n', bty='n')
	title(main=q_id, cex=1.25)


	# plot axes
	axis(side=3, at=seq(from=0,to=seq_len,by=1000), pos=-1*tm_bool)
	text(x=-seq_len*0.05, y=1:I-0.4, labels=1:I, cex=1.25)
	if (tm_bool) { text(x=-seq_len*0.05, y=-0.4, labels='TM', cex=1.25) }


	# plot TM regions
	if (tm_bool) {
		rect(xleft=0, ybottom=0, xright=seq_len, ytop=-0.8, col='grey85', border='black')
		if (file.exists(tmF_coo)) {
			tab <- read.table(tmF_coo, header=TRUE, stringsAsFactors=FALSE)
			segments(x0=tab$tm_cl_from-1, x1=tab$tm_cl_to, y0=-0.3, lwd=5, col='darkred', lend='butt')
			text(x=(tab$tm_cl_from+tab$tm_cl_to)/2, y=-0.5, labels=tab$tm_cl_index, col='darkred')
		}
	} 


	# plot query regions & build table
	TAB <- data.frame(NULL, stringsAsFactors=FALSE)
	for (i in 1:I) {
		qF <- paste0(q_id, '_query_regions_coo_iteration_',i,'.tsv')
		qt <- read.table(qF, stringsAsFactors=FALSE)
		rect(xleft=qt[,1]-1, ybottom=i, xright=qt[,2], ytop=i-0.8, col='grey85', border='black')

		clustF <- paste0(q_id, '_clusters_of_hits_iteration_',i,'.tsv')
		if (!file.exists(clustF)) { next }
		CL <- read.table(clustF, header=TRUE, colClasses=c('cl_index'='character'), stringsAsFactors=FALSE)
		TAB <- rbind(TAB, CL)
	}


	# plot hits & write table
	if (nrow(TAB) > 0) {
		TAB <- TAB[order(TAB$cl_from),]

		segments(x0=TAB$cl_from-1, x1=TAB$cl_to, y0=TAB$iterat_num-0.3, lwd=5, col='royalblue', lend='butt')
		text(x=(TAB$cl_from+TAB$cl_to)/2, y=TAB$iterat_num-0.5, labels=TAB$cl_index, col='royalblue')

		write.table(TAB, file=tabF, sep='\t', quote=FALSE, row.names=FALSE)
	}
	else { print('No hits found!') }


	dev.off()
}
