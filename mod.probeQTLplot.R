# Custom probe-level QTL plot

mod.probeQTLplot<-function (probeset, probeQtlProfiles, markerPos, chrOffsets) {

	markerPos.gmb <- markerPos[, "gmb"]
	nChr <- length(chrOffsets)
	nProbes <- length(unique(probeQtlProfiles$probenr))
	maxheight <- max(probeQtlProfiles[, 5])
	stepsize <- maxheight
	qtlThresh <- 3 # Suggested by the paper, I didn't calculate this.
	
	# Main title
	title<-paste(probeset, " Probe-level QTL analysis", sep="")
	
	plot(markerPos.gmb, probeQtlProfiles[probeQtlProfiles$probenr == 1, 5],
		 type = "n", yaxt = "n", main = title, xlab = "Marker", ylab = "-10log p-value", 
		ylim = range(0:(maxheight * 110 * nProbes)/100), cex = 1)


	totmaxheight <- (maxheight * 110 * nProbes)/100

	for (i in 1:nChr) {
		lines(c(chrOffsets[i], chrOffsets[i]), c(0, totmaxheight), 
		lty = 1, lwd = 1, col = "grey")
	
		
		if (i == 1) {
			text((chrOffsets[i]+chrOffsets[i+1])/2, totmaxheight, "Chr 1",)
		}
			
		else if (i < nChr && i > 1) {
			text((chrOffsets[i]+chrOffsets[i+1])/2, totmaxheight, i)
		}
			
		else if (i == nChr) {
			lines(c(max(markerPos.gmb), max(markerPos.gmb)), c(0, totmaxheight), 
			lty = 1, lwd = 1, col = "grey")
			text((chrOffsets[i]+max(markerPos.gmb))/2, totmaxheight, "X")
		}
	}
	
	
	mycols <- rep(c("orange1", "orange4", "lightblue4", "blue4"), 6)

	for (i in 1:nProbes) {
		pp <- data.frame(probeQtlProfiles[probeQtlProfiles$probenr == i, 5])
		lines(markerPos.gmb, pp[, 1] + (i - 1) * stepsize, col = mycols[i])
		text(24, (i - 1) * stepsize, i, col = mycols[i], pos = 2)
		lines(c(0, max(markerPos.gmb)), c((i - 1) * stepsize + qtlThresh, 
		(i - 1) * stepsize + qtlThresh), lty = 2, col = mycols[i])
		}
		
	# Marker location ticks
	for (idxMarkers in 1:length(markerPos.gmb)) {
		lines(c(markerPos.gmb[idxMarkers], markerPos.gmb[idxMarkers]), 
		#c(0.4 - stepsize, 2.5 - stepsize), col = "darkgrey")
		c(-2, -1), col = "darkgrey")
		}

}
