# Plots probe signals for a given probe set. Adapted from
# affyGG package. 

# traits: matrix, probe-level data for a probeset. Rows=probes, cols=samples.
# genotypes: matrix of genotype data. Rows=markers, cols=samples. 
# marker: character, marker where strains are grouped by genotype
# alleles: character vector of all alleles. First two are used to color lines. 
# probesPos: numeric vector of each probes relative location
# Script will color lines based on genotype at the QTL marker
# (1 = ISS, 2 = ILS)
# or
# (1=B6, 2=D2, 0=H)
mod.probePlot<-function (traits, probeset, marker, genotypes, alleles, probes.pos) {
	
	# Get y-axis limits
	ylow=floor(range(traits)[1])-1
	yhigh=floor(range(traits)[2])+1

	# Allele color sequence
	allele.seq<-as.numeric(genotypes[marker,])
	# Convert 1 to blue, 2 to red, and 3 to green
	allelecolors<-ifelse(allele.seq==1,"blue",ifelse(allele.seq==2,"red","green"))

	# Main Title
	title <-paste(probeset, " Probe-Level RMA Values", sep="")

	# Adjust probe positions for graph by starting probes at 1
	probesPos <- probes.pos - min(probes.pos) + 1
	lengthprobe <- 25
	
	# Top plot
	###########
	
	# mai gives margin size in inches
	par(mai = c(0, 1, 0.2, 0.3))
	# fig sets coordinates for the top plot
	par(fig = c(0, 1, 0.79, 1))

	x<-c(probesPos, max(probesPos) + lengthprobe + 2)
	y<-c(1:length(probesPos), length(probesPos)+2)

	plot(x, y, xlab = paste("probeset:", probeset), 
		ylab = "", type = "n", xaxt = "n", yaxt="n", main=title)

	for (i in 1:length(probesPos)) {
		lines(c(probesPos[i], probesPos[i] + lengthprobe), c(i,i), lwd = 2, col = "blue")
		text(probesPos[i] + 1, i + 0.9, i, cex=.8)
		}

	text(probesPos[length(probesPos)] + 25, 2, "Relative probe position (bp)", pos = 2, cex=.8)
	
	# Bottom plot
	#############
	# mai gives margin size in inches
	par(mai = c(1.2, 1, 0.1, 0.3))
	# fig sets coordinates for the top plot
	par(fig = c(0, 1, 0, 0.8), new = T, xpd = T)
	myxlab <- paste("Probes from", as.character(probeset), sep=" ")
	matplot(traits, type = "l", main = "", xlab = myxlab, xaxt = "n", 
		ylab = "log2 intensity", lty = 1, ylim = c(ylow,yhigh), 
		col = allelecolors, axes=FALSE)

	axis(1, 1:nrow(traits))
	axis(2)

	#text(10,ylow, paste("Genotypes at", marker),cex=.9, font=3)
	
	# Allele legend
	legend(x="bottomleft", legend=alleles, col=c("blue","red","green"), lty=1, lwd=2,
		box.lty=0, horiz=TRUE, title=paste("Genotypes at", marker))

}
