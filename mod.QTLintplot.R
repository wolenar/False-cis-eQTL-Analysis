# Script to make a QTL plot-like graph using the
# marker -log10 p-values from the affyGG probeset QTL

# plot.int: logical, if true the probe marker interaction is also plotted
# plot.sig: logical, if true the pvalue for the most significant qtl is plotted

mod.QTLintplot<-function(probeset, probesetQtlProfile, markerPos, chrOffsets, gene.gmb, plot.int, plot.sig) {
	
	if(missing(plot.int)){
		plot.int<-FALSE
	}
	if(missing(plot.sig)){
		plot.sig<-FALSE
	}

	markerPos.gmb <- markerPos$gmb

	text.size<-1
	
	# Data
	qtl.pvalues<-probesetQtlProfile$minlog10pmarker

	# Get y-axis limits
	###################
	ylow<-0
	# Determine qtl.pvalues or int.pvalues has largest value
	# then use that number as the max of y-axis
	yhigh<-max(qtl.pvalues)
	
	# Main title
	title <-paste(probeset, "Probe-set QTL plot")
	
	
	# Add parameters necessary for plotting interaction data
	if(!missing(plot.int) & plot.int==TRUE){
		int.pvalues<-probesetQtlProfile$minlog10pinteraction
		yhigh<-ifelse(max(int.pvalues) > max(qtl.pvalues), max(int.pvalues),yhigh)
		title<-sub("plot","interaction analysis",title)
	}
	

	yhigh<-ceiling(yhigh)+1
	# Determine the yaxis interval
	ystep<-ifelse(yhigh>=10, 2, 1)

	# Main Plot
	###############
	plot(markerPos.gmb, qtl.pvalues, type="n", axes=FALSE, main=title, 
	ylab="-log10(p-value)", xlab="Megabases", ylim=c(ylow,yhigh), lwd=1.75, 
	font.lab=1, cex.lab=1.2, cex.main=1.5, cex.axis=1.2)
	
	# Axes
	##########
	# X-axis	
	axis(1, at=(seq(0,2600,200)),
	labels=c("0", "200", "400", "600", "800", "1,000", "1,200", "1,400", 
	"1,600", "1,800", "2,000", "2,200", "2,400", 	"2,600" ), cex.axis=1)
	# Y-axis,
	axis(2, at=(seq(ylow, yhigh, ystep)),cex.axis=1)

	# Add triangle symbol to indicate gene's location
	points(gene.gmb, -.1, pch=17, col="purple")

	# Chromosome dividers
	#########################
	nChr<-20
	text.size=1
	even<-seq(2,nChr,2)	# To find even chr's for shading
    
	for (i in 1:nChr) {	
		# Shaded chromosomes
		if (sum(i==even)==1) {
			polygon(c(chrOffsets[i], chrOffsets[i], chrOffsets[i+1], 
			chrOffsets[i+1]), c(0, yhigh, yhigh, 0), 
			border=FALSE, col = rgb(0,0,0,.05))
		}
		# Chromosome labels 1 - 19 	   
		if (i < nChr) {
			text(median(c(chrOffsets[i], chrOffsets[i+1])), 
			yhigh*.97, i, cex=text.size)
		}
		else if (i == nChr) {
			# Text for X chr
			text(median(c(chrOffsets[i], 2649)), 
			yhigh*.97, "X", cex=text.size)
			# Shading for X chr  
			polygon(c(chrOffsets[i], chrOffsets[i], 2649, 
			2649), c(0, yhigh, yhigh, 0), 
			border=FALSE, col = rgb(0,0,0,.05))
		}
	}
	
	if(plot.int==TRUE){
		lines(markerPos.gmb, int.pvalues, lwd=1.75, col=rgb(.9,0,0,.8))
		legend(x="top", legend=c("QTL p-value", "Interaction p-value"),
		col=c("blue", "red"), lty=c(1,1), lwd=2, box.lty=0, horiz=TRUE, cex=.8)
	}
	
	lines(markerPos.gmb, qtl.pvalues, lwd=1.75, col=rgb(0,0,.9,.8))
	
	if(plot.sig==TRUE){
		# Find peak marker
	peak.marker<-probesetQtlProfile$markers[probesetQtlProfile$minlog10pmarker==max(probesetQtlProfile$minlog10pmarker)][1]
		# Look-up peak markers mb location
		peak.mb<-markerPos[rownames(markerPos)==peak.marker,"gmb"]
		# Find peak pvalue
		peak.pvalue<-min(probesetQtlProfile$pmarker)
		# Round peak pvalue to 4 signficant numbers
		peak.pvalue<-signif(peak.pvalue, digits=4)
		# Add star above peak QTL
		points(peak.mb, max(qtl.pvalues)*1.03, pch=8)
		
		legend(x="bottomleft", legend=paste("peak p-value: ", peak.pvalue, sep=""), pch=8, box.lty=0, cex=.7)	
	}

}