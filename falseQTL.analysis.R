falseQTL.analysis <- function(probesets, geno.file, cel.dir, strain, label, 
						cis.buffer=5, batch=NULL) {
	######################
	# Functions          #
	######################
	# RMA Preprocessing
	mod.rma.preprocessing<-function(A){
		A2 <- bg.correct.rma(A)
		A3 <- normalize.AffyBatch.quantiles(A2)
		pmperprobe<-cbind(probeset=probeNames(A3), 
			round(log2(data.frame(pm(A2))), 3))
		return(pmperprobe)
	}
	
	# Print current time
	print.time<-function(){paste("\n",format(Sys.time(),"%D %l:%M %p"),sep="")}
	
	# Convert nucleotides to megabases
	toMb<-function(x){
		return(x/1000000)
	}
	
	######################
	# Main Program       #
	######################
	
	orig.wd<-getwd()
	
	writeLines(paste(print.time(), "Analysis started..."))

	# Load necessary libraries
	require(affy)
	require(affyGG)
	
	
	chrs<-c(1:19,"X")

	strain<-toupper(strain)

	# Load genotype data
	geno<-read.csv(geno.file, row.names=1)

	# Extract marker positions from genotype data
	if(length(grep("chr", colnames(geno), ignore.case=TRUE)>0)){
		markerPos<-data.frame(chr=geno[,grep("chr", colnames(geno),
		ignore.case=TRUE)], mb=geno[,grep("mb", colnames(geno),
		ignore.case=TRUE)], row.names=rownames(geno))
  } else {
      stop("\nGenotype file missing marker position columns: Chr and/or Mb.")
      }
	  
	# Create genotype dataframe without marker positions
	geno<-data.frame(geno[,grep(strain,colnames(geno))])
	strains<-colnames(geno)
  
	# Create numeric genotype matrix
	#################################
	# Count how many times each unique allele, ignore case
	allele.counts<-table(toupper(unlist(unlist(geno))))
	# Sort alleles by frequency.
	allele.counts<-sort(allele.counts, dec=TRUE)
	# Create a key for alleles and their numeric code in the genotype matrix
	alleles<-data.frame(allele=names(allele.counts), code=1:length(names(allele.counts)))
	# Replace most common allele with 1's, second most common with 2's, etc...
	geno.mat<-geno
	for(i in 1:nrow(alleles)){
		geno.mat<-apply(geno.mat, MARGIN=1, FUN=function(x) 
			#replace(x, list=which(x==alleles[i,"allele"]), values=as.character(i)))
			replace(x, list=grep(alleles[i,"allele"], x, ignore.case=TRUE), 
				values=as.character(i)))
	}
	geno.mat<-t(geno.mat)
	
	# Read in CEL files to determine chip type
	writeLines(paste(print.time(),"Loading CEL files..."))
	cel.names<-list.celfiles(cel.dir)
	affy.object<-ReadAffy(celfile.path=cel.dir, filenames=cel.names)
	array.type<-cdfName(affy.object)
	cdf.name<-cleancdfname(array.type)
	array.name<-sub("cdf","",cdf.name)
	
	# Reconciling genotype data and celfiles strains
	################################################
	keep.strains<-sapply(strains, 
		FUN=function(x) length(grep(paste(x,"\\D",sep=""), 
		cel.names, perl=TRUE)))
	# Only keep strains with a matching cel.file
	keep.strains<-keep.strains[keep.strains>0]
	# Check to make sure no strains matched multiple cel.files
	if(sum(keep.strains>1)>0){
		stop("At least one Strain ID from genotype data file matches more than one CEL file. Please ensure all genotype columns and CEL files are labeled correctly")
	}
	
	# Exclude strains with no available celfiles
	geno<-geno[,colnames(geno) %in% names(keep.strains)]
	geno.mat<-geno.mat[,colnames(geno.mat) %in% names(keep.strains)]
	writeLines(paste("\n", ncol(geno), 
		"strains from genotype data file have corresponding CEL files.\n"))
		
	# RMA Preprocessing
	###########################
	writeLines(paste(print.time(),"Calculating probe signals..."))
	probesignals<-mod.rma.preprocessing(affy.object)

	# Load array specific probe and annotation data from bioconductor
	#################################################################
	writeLines(paste(print.time(),"Loading sequence and annotation data..."))

	# Bioconductor packages
	probe.pkg<-paste(array.name,"probe",sep="")
	annot.pkg<-paste(array.name,".db",sep="")

	# Check for necessary bioconductor packages
	if(sum(c(probe.pkg, annot.pkg) %in% rownames(installed.packages()))<2){
		stop(paste("\nMissing required bioconductor package. Please ensure that",
			probe.pkg, "and", annot.pkg, "are installed before retrying."))
	}

	# Load array probe sequence info
	require(package=probe.pkg,character.only=TRUE, quietly=TRUE)

	# Load array annotation data
	require(package=annot.pkg, character.only=TRUE)
	
	# Calculate gigabases
	########################	
	# Load chromosome lengths
	chrLengths<-eval(as.name(paste(array.name,"CHRLENGTHS",sep="")))
	chrLengths<-toMb(chrLengths[names(chrLengths) %in% chrs])
	chrOffsets<-c(0,cumsum(chrLengths)[1:19])
	names(chrOffsets)=chrs
	
	# Add GMb to markerPos
	markerPos<-cbind(markerPos,gmb=markerPos$mb)
	for(i in 1:length(chrs)){
		markerPos$gmb[markerPos$chr==chrs[i]] <-
			markerPos$gmb[markerPos$chr==chrs[i]]+chrOffsets[i]
	}
	
	# Repeat for each provided probeset
	#######################################
	for(probeset.num in 1:length(probesets)){
		probeset<-probesets[probeset.num]
		writeLines(paste(print.time(),paste("Begining analysis for",probeset)))
		# Extract probe positions for probeset
		probes.pos<-data.frame(subset(eval(as.name(probe.pkg)),
			subset=Probe.Set.Name==probeset, select=Probe.Interrogation.Position))
		probes.pos<-as.numeric(probes.pos[,1])
		
		# Gene symbol
		gene.sym<-get(probeset,
			envir=eval(as.name(paste(array.name,"SYMBOL",sep=""))))
		# If no gene symbol is available just use probeset ID
		if(is.na(gene.sym)){
			writeLines(paste("No gene symbol available for",probeset))
			gene.sym<-probeset
		}
		# Get full gene name
		gene.name<-get(probeset,
			envir=eval(as.name(paste(array.name,"GENENAME",sep=""))))
		# Get gene's chromosomal location
		gene.chr<-get(probeset,
			envir=eval(as.name(paste(array.name,"CHR",sep=""))))
		# Get gene's chromosomal start location
		gene.start.mb<-get(probeset,
			envir=eval(as.name(paste(array.name,"CHRLOC",sep=""))))[1]
		gene.strand<-ifelse(gene.start.mb<0, "-", "+")
		gene.start.mb<-toMb(abs(gene.start.mb))
		# Get gene's chromosomal stop location
		gene.end.mb<-get(probeset,
			envir=eval(as.name(paste(array.name,"CHRLOCEND",sep=""))))[1]
		gene.end.mb<-toMb(abs(gene.end.mb))
		
		# Gene's genomic location
		gene.start.gmb <- gene.start.mb+chrOffsets[gene.chr]
		gene.end.gmb <- gene.end.mb+chrOffsets[gene.chr]
		
		# Gene's transcription start location
		gene.tx.mb<-ifelse(gene.strand=="+", gene.start.mb, gene.end.mb)	
		gene.tx.gmb<-gene.tx.mb+chrOffsets[gene.chr]

		# Prefix to add to all output files
		results.prefix<-paste(gene.sym, probeset, label, sep="_")

		# Create a new folder in wd to store results of analysis
		# then set that as the working directory
		results.folder<-paste(results.prefix, "results", sep="_")
		system(paste("mkdir ", results.folder, sep=""))
		setwd(results.folder)
	
		# Select probe-level data for one probe set
		traits <- as.matrix(probesignals[probesignals$probeset==probeset,
			2:ncol(probesignals)])
		colnames(traits)=names(keep.strains)
		rownames(traits)=paste(probeset,1:nrow(traits),sep="")
		write.csv(traits, file=paste(results.prefix, "probesignals.csv",sep="_"))
	
		# Perform QTL Analysis on probe level
		#####################################
		# Outputs a data.frame containing Fmarkerperprobe and
		# pmarkerperprobe for each probe in a probe-set at all markers
		writeLines(paste(print.time(),"Probe-level QTL analysis..."))
		qtlmap.probe<-qtlMap.xProbe(genotypes=geno.mat, 
			traits=traits, batch=batch)

		# Add -log10 pmarker per probe values for mapping
		qtlmap.probe<-cbind(qtlmap.probe,
			minlog10pmarkerperprobe=-log10(qtlmap.probe[,4]))
		# Export probe-level QTL analysis results
		write.csv(qtlmap.probe, file=paste(results.prefix,
			"probe_level_QTL_results.csv", sep="_"))

		# Perform QTL analysis on Probeset-level
		########################################
		writeLines(paste(print.time(),"ProbeSet-level QTL analysis..."))
		qtlmap.probeset<-qtlMap.xProbeSet(genotypes=geno.mat, 
			traits=traits, batch=batch)
		## 4.5 min on quad-core desktop
		qtlProbeset <- -log10(qtlmap.probeset$pmarker)
		intProbeset <- -log10(qtlmap.probeset$pinteraction)

		# Add -log10 values
		qtlmap.probeset<-cbind(qtlmap.probeset, minlog10pmarker=qtlProbeset,
			minlog10pinteraction=intProbeset)

		# Export probeset-level QTL analysis results
		write.csv(qtlmap.probeset, file=paste(results.prefix,
			"probeset_QTL_results.csv", sep="_"))

		# Identify marker where qtlProbeset is highest
		###############################################
		# If there's a tie for peak marker
		# select the marker closest to gene's
		# transcription start location	

		peak.genome.marker<-as.character(qtlmap.probeset$marker[which
			(qtlmap.probeset$pmarker==min(qtlmap.probeset$pmarker))])
		
		if(length(peak.genome.marker)>1){
			peak.genome.marker<-peak.genome.marker[
				which.min(abs(markerPos$gmb[match(peak.genome.marker, 
				rownames(markerPos))]-gene.tx.gmb))]
		}
	
		# Identify all markers within the specified cis.buffer
		cis.markers<-subset(markerPos, subset=markerPos$chr==gene.chr & 
			markerPos$mb >= (gene.start.mb - cis.buffer) &
			markerPos$mb <=(gene.end.mb + cis.buffer))
		cis.markers<-rownames(cis.markers)
	
		# Subset of qtl data containing only cis markers
		cis.probeset<-qtlmap.probeset[match(cis.markers,
			qtlmap.probeset$marker),]
	
		# Ensure there is a marker with a pvalue <.05
		if(sum(cis.probeset$pmarker<=.05)==0){
			stop(paste("No markers within the specified cis eQTL buffer ",
			"region are significantly associated with ",
			probeset,"'s expression.", sep=""), call.=FALSE)
		}
	
		peak.cis.marker<-as.character(cis.probeset$marker[which
			(cis.probeset$pmarker==min(cis.probeset$pmarker))])
	
		if(length(peak.cis.marker)>1){
			peak.cis.marker<-peak.cis.marker[
				which.min(abs(markerPos$gmb[match(peak.cis.marker, 
				rownames(markerPos))]-gene.tx.gmb))]
		}
	
		peak.cis.pval<-as.numeric(subset(qtlmap.probeset, 
			subset=markers==peak.cis.marker, select=pmarker))
		
		peak.geno.pval<-as.numeric(subset(qtlmap.probeset, 
			subset=markers==peak.genome.marker, select=pmarker))
	
		if(peak.cis.marker==peak.genome.marker){
			writeLines(paste("\nPeak marker for",probeset,
			"cis-eQTL is",peak.cis.marker,"\non",
			"Chr", markerPos[peak.cis.marker,"chr"], "at",
			 markerPos[peak.cis.marker,"mb"],
			"Mb with a p-value of", formatC(peak.cis.pval, form="e", dig=2)))
		} else {
			writeLines(paste("\nPeak marker for",probeset,
			"is a trans eQTL at",peak.genome.marker,
			"\non", "Chr", markerPos[peak.genome.marker,"chr"], "at", 
			markerPos[peak.genome.marker,"mb"], "Mb with a p-value of", 
			formatC(peak.geno.pval, form="e", dig=2),
			"\n\nThe most significant cis-eQTL is",
			peak.cis.marker,"\non","Chr", markerPos[peak.cis.marker,"chr"], "at", 
			markerPos[peak.cis.marker,"mb"], "Mb with a p-value of", 
			formatC(peak.cis.pval, form="e", dig=2)))
		}

		# Perform probe elimination assay
		####################################
		writeLines(paste(print.time(),"Probe elimination analysis..."))
		pe <- probeElimination(probesetName=probeset,
			markerName=peak.cis.marker, genotypes=geno.mat, 
			traits=traits, batch=batch)

		# Create spread sheet with probe-elimination results
		write.csv(pe, file=paste(results.prefix, 
			"probe_elimination_results.csv", sep="_"))


		# Plots
		########
		writeLines(paste(print.time(),"Generating graphs..."))
		pdf(paste(results.prefix,"probePlot.pdf", sep="_"), 
			height = 8, width = 10.5)
		mod.probePlot(traits=traits, probeset=probeset, marker=peak.cis.marker, 
			genotypes=geno.mat, alleles=alleles[,"allele"], probes.pos=probes.pos,
			probeset.loc=gene.tx.mb, marker.loc=markerPos[peak.cis.marker,"mb"])
		dev.off()
		
		# Export data needed to recreate probePlot
		probePlot.input<-as.list(traits, probeset, peak.cis.marker, geno.mat,
			alleles[,"allele"], probes.pos, gene.tx.mb, markerPos[peak.cis.marker,"mb"])

		# Create custom probe-level QTL plot
		pdf(paste(results.prefix,"probe-level_qtlPlot.pdf", sep="_"), 
			height = 8, width = 10.5)
		mod.probeQTLplot(probeset=probeset, probeQtlProfiles=qtlmap.probe, 
			markerPos=markerPos, chrOffsets=chrOffsets)
		dev.off()

		# Create custom probeset QTL interaction plot
		pdf(paste(results.prefix, "probeset_interaction_qtlPlot.pdf", sep="_"), 
			height = 8, width = 10.5)
		mod.QTLintplot(probeset=probeset, probesetQtlProfile=qtlmap.probeset,
			markerPos=markerPos, 
			chrOffsets=chrOffsets, gene.gmb=gene.tx.gmb, plot.int=TRUE)
		dev.off()
	
		# Create custom probeset QTLplot
		pdf(paste(results.prefix, "probeset_qtlPlot.pdf", sep="_"), 
			height = 8, width = 10.5)
		mod.QTLintplot(probeset=probeset, probesetQtlProfile=qtlmap.probeset,
			markerPos=markerPos, chrOffsets=chrOffsets, gene.gmb=gene.tx.gmb, 
			plot.int=FALSE, plot.sig=TRUE)
		dev.off()	
		setwd(orig.wd)
	}

	writeLines(paste(print.time(),"Analysis complete. Huzzah!"))
}
