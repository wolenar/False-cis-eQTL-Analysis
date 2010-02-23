geno.file<-"/Users/milesg5users/Dropbox/MilesLab Share/RI Datasets/BXD datasets/BXD Genotype data/BXD Genotypes.csv"
cel.dir<-"/Users/milesg5users/Documents/Aaron/CEL_files/BXD CEL Files/PFC Saline"

strain<-"bxd"
# Currently strain is only necessary to cidentify which columns in 
# genotype file contain actual data

falseQTL.analysis<- function(probeset, geno.file, cel.dir, label, cis.buffer, batch) {
	
  print.time<-function(){format(Sys.time(),"%D %l:%M %p")}
	writeLines(paste("Analysis started: ",print.time(), sep=""))
	
	# Set batch to null if none is provided
	if(missing(batch)){
		batch<-NULL
	}

	# Load necessary libraries
	require(affy)
	require(affyGG)
	chrs<-c(1:19,"X")

	# Convert strain to uppercase
	# strain<-toupper(strain) UNCECESSARY?

	# Load genotype data
	geno<-read.csv(geno.file, row.names=1)

	# Extract marker positions from genotype data
	if(length(grep("chr", colnames(geno), ignore.case=TRUE)>0)){
     markerPos<-data.frame(chr=geno[,grep("chr", colnames(geno), ignore.case=TRUE)],
      mb=geno[,grep("mb", colnames(geno), ignore.case=TRUE)],
      row.names=rownames(geno))
  } else {
      stop("Genotype file missing marker position columns: Chr and/or Mb.")
      }
	  
	# Create genotype dataframe without marker positions
	geno<-data.frame(geno[,grep(strain,colnames(geno))])
  
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

	# Read in CEL files to determine chip type
	writeLines("Reading in CEL files...")
	writeLines(print.time())
	cel.names<-list.celfiles(cel.dir)
	affy.object<-ReadAffy(celfile.path=cel.dir, filenames=cel.names)
	array.type<-cdfName(affy.object)
	cdf.name<-cleancdfname(array.name)
	array.name<-sub("cdf","",cdf.name)

	# Load array specific probe and annotation data from bioconductor
	#################################################################
	writeLines("Loading sequence and annotation data...")
	writeLines(print.time())

	# Load array probe sequence info
	require(package=paste(array.name,"probe",sep=""), character.only=TRUE)
	probes.pos<-mouse430a2probe[mouse430a2probe$Probe.Set.Name==probeset,5]
	
	probes.pos<-data.frame(subset(eval(as.name(paste(array.name,"probe",sep=""))),
		subset=Probe.Set.Name==probeset, select=Probe.Interrogation.Position))
	probes.pos<-as.numeric(probes.pos[,1])

	# Load array annotation data
	require(package=paste(array.name,".db",sep=""), character.only=TRUE)
	# Gene symbol
	gene.sym<-get(probeset, envir=eval(as.name(paste(array.name,"SYMBOL",sep=""))))
	# If no gene symbol is available just use probeset ID
	if(is.na(gene.sym)){
		writeLines(paste("No gene symbol available for",probeset))
		gene.sym<-probeset
	}
	# Get full gene name
	gene.name<-
		get(probeset, envir=eval(as.name(paste(array.name,"GENENAME",sep=""))))
	# Get gene's chromosomal location
	gene.chr<-
		get(probeset, envir=eval(as.name(paste(array.name,"CHR",sep=""))))
	# Get gene's chromosomal start location
	gene.start<-
		get(probeset, envir=eval(as.name(paste(array.name,"CHRLOC",sep=""))))
	# Get gene's chromosomal stop location
	gene.end<-
		get(probeset, envir=eval(as.name(paste(array.name,"CHRLOCEND",sep=""))))

	# Calculate gigabases
	########################	
	# Load chromosome lengths
	chrLengths<-eval(as.name(paste(array.name,"CHRLENGTHS",sep="")))
	chrLengths<-chrLengths[names(chrLengths) %in% chrs]/1000000
	chrOffsets<-c(0,cumsum(chrLengths)[1:19])
	names(chrOffsets)=chrs
		gene.gmb<-as.numeric(as.character(chip.annot[probeset, "GMb"]))
		
	# Add GMb to markerPos
	markerPos<-cbind(markerPos,gmb=markerPos$mb)
	for(i in 1:length(chrs)){
		markerPos$gmb[markerPos$chr==chrs[i]] <-
			markerPos$gmb[markerPos$chr==chrs[i]]+chrOffsets[i]
	}
		

	# Prefix to add to all output files
	results.prefix<-paste(gene.sym, probeset, label, sep="_")

	# Check to make sure annotation information was available
	# for current probeset. If it's not, prompt user to enter it.
	#if(is.na(gene.sym)==TRUE) {
	#	system("printf '\a'")
	#	gene.sym<-readline(paste("I'm not familiar with ",probeset,
	#			". Please enter its gene name here: ",sep=""))
	#	}

	#if(is.na(gene.mb)==TRUE) {
	#	system("printf '\a'")
	#	gene.chr<-readline(paste("What chromosome is ",probeset,
	#			" located on? ",sep=""))

	#	gene.mb<-readline(paste("And where is ", probeset,
	#			" located on chromosome ", chr, " (eg 74324567)? ",sep=""))
	#	}

	# Create a new folder in wd to store results of analysis
	# then set that as the working directory
	results.folder<-paste(results.prefix, "results", sep="_")
	system(paste("mkdir ", results.folder, sep=""))
	setwd(results.folder)

	# Set batch to NULL if batch.factor=FALSE
	if (!exists("batch")) {
		batch<-as.null()
		}

	# Reconciling genotype and celfile strains
	strains<-colnames(geno)
	celfiles<-list.celfiles(cel.dir)

	keep.strains<-as.logical()
	for(i in 1:length(strains)){
		keep.strains[i]<-length(grep(strains[i], celfiles))==1
		}
		
	# Exclude strains with no available celfiles
	strains<-strains[keep.strains]
	geno<-geno[,keep.strains]
	geno.mat<-geno.mat[,keep.strains]
  
	writeLines("Calculating probe signals...")
	writeLines(print.time())
	probesignals<-rma.preprocessing(cel.dir, celfiles)
	
	# Select probe-level data for one probe set
	traits <- as.matrix(probesignals[probesignals$probeset==probeset,2:ncol(probesignals)])
	colnames(traits)=strains
	rownames(traits)=paste(probeset,1:nrow(traits),sep="")
	write.csv(traits, file=paste(results.prefix, "probesignals.csv",sep="_"))
	
	# Perform QTL Analysis on probe level
	#####################################
	# Outputs a data.frame containing Fmarkerperprobe and
	# pmarkerperprobe for each probe in a probe-set at all markers
	writeLines("Probe-level QTL Analysis...")
	writeLines(print.time())
	qtlmap.probe<-qtlMap.xProbe(genotypes=geno.mat, traits=traits, batch=batch)
	# Update data frame to include -log10 pmarker per probe values for mapping
	qtlmap.probe<-cbind(qtlmap.probe,minlog10pmarkerperprobe=-log10(qtlmap.probe[,4]))

	# Create spread sheet with probe-level QTL analysis results
	write.table(qtlmap.probe, file=paste(results.prefix, "probe_level_QTL_results.csv", sep="_"),
				sep=",", col.names=TRUE, row.names=FALSE)

	# Perform QTL analysis on Probeset-level
	########################################
	writeLines("Probeset-level QTL Analysis...")
	writeLines(print.time())
	qtlmap.Probeset<-qtlMap.xProbeSet(genotypes=geno.mat, traits=traits, batch=batch)
	qtlProbeset <- -log10(qtlmap.Probeset$pmarker)
	intProbeset <- -log10(qtlmap.Probeset$pinteraction)

	# Update data frame to include the -log10 values for mapping
	qtlmap.Probeset<-cbind(qtlmap.Probeset, minlog10pmarker=qtlProbeset,
		minlog10pinteraction=intProbeset)

	# Create spread sheet with probeset-level QTL analysis results
	write.table(qtlmap.Probeset, file=paste(results.prefix, "probeset_QTL_results.csv", sep="_"),
				sep=",", col.names=TRUE, row.names=FALSE)

	# Identify marker where qtlProbeset is highest
	###############################################
	# If cis.buffer is NULL search entire genome
	# Otherwise limit search to cisrange of gene
	if (is.null("cis.buffer")) {
		peak.marker <- as.character(qtlmap.Probeset$marker[qtlmap.Probeset$pmarker == min(qtlmap.Probeset$pmarker)])
		}

	if (!is.null("cis.buffer")) {
		# Identify markers inside the cis QTL range
		# Create dataframe with only markerPos from probeset's home chr
		markerPos.chr<-subset(markerPos, subset=markerPos$chr==gene.chr)
		cis.markers<-rownames(markerPos.chr[markerPos.chr$mb<=gene.mb+cis.buffer &
		           markerPos.chr$mb>=gene.mb-cis.buffer,])
		# Create new data frame with markers as rownames
		cis.qtlmap.Probeset<-data.frame(qtlmap.Probeset[,2:7], row.names=qtlmap.Probeset[,1])
		# Keep only data for markers inside QTL range
		cis.qtlmap.Probeset<-subset(cis.qtlmap.Probeset, subset=rownames(cis.qtlmap.Probeset) %in% cis.markers)
		# Identify peak markers in this range
		peak.marker <- cis.qtlmap.Probeset[cis.qtlmap.Probeset$pmarker == min(cis.qtlmap.Probeset$pmarker),]
		peak.marker<-rownames(peak.marker)
		# If multipe markers are tied for lowest pvalue, 
		# grab marker nearest to probeset
		if(length(peak.marker)>1){
			peak.marker<-peak.marker[which.min(abs(markerPos.chr[peak.marker,"mb"]-gene.mb))]
			}
		rm(cis.qtlmap.Probeset)
		}

	# Get peak marker's location
	peak.markerPos<-markerPos[peak.marker,]
	
	writeLines(paste("\nPeak marker for",probeset,"is",peak.marker,":\n",
		peak.markerPos,"\n"))


	# Eliminiate deviating probes using the statisical  model
	#########################################################
	writeLines("Performing probe elimination assay...")
	writeLines(print.time())
	pe <- probeElimination(probesetName=probeset, markerName=peak.marker,
		genotypes=geno.mat, traits=traits, batch=batch)
	pe
	# Create spread sheet with probe-elimination results
	write.table(pe, file=paste(gene.sym,probeset,strain,label,
				"probe_elimination_results.csv", sep="_"),
				sep=",", col.names=TRUE, row.names=FALSE)

	# Create workspace .RData file
	#print("Saving workspace...", sep="")
	#workspace.name<-paste(gene.sym,"_",treatment,"_affyGG.RData", sep="")
	#save.image(file=workspace.name)

	# Plots
	########
	writeLines("Creating plots...")
	writeLines(print.time())
	
	pdf(paste(results.prefix,"probePlot.pdf", sep="_"), 
		height = 8, width = 10.5)
	mod.probePlot(traits=traits, probeset=probeset, marker=peak.marker, 
		genotypes=geno.mat, alleles=alleles, probes.pos=probes.pos)
	dev.off()

	# Create custom probe-level QTL plot
	pdf(paste(results.prefix,"probe-level_qtlPlot.pdf", sep="_"), 
		height = 8, width = 10.5)
	mod.probeQTLplot(probeset=probeset, probeQtlProfiles=qtlmap.probe, 
		markerPos=markerPos, chrOffsets=chrOffsets)
	dev.off()

	# Create custom probeset QTL interaction plot
	pdf(paste(results.prefix, "probeset_interaction_qtlPlot.pdf", sep="_"), 
		height = 8, width = 10.5)
	mod.QTLintplot(probeset=probeset, probesetQtlProfile=qtlmap.Probeset,
		markerPos=markerPos, 
		chrOffsets=chrOffsets, gene.gmb=gene.gmb, plot.int=TRUE)
	dev.off()
	
	# Create custom probeset QTLplot
	pdf(paste(results.prefix, "probeset_qtlPlot.pdf", sep="_"), 
		height = 8, width = 10.5)
	mod.QTLintplot(probeset=probeset, probesetQtlProfile=qtlmap.Probeset,
		markerPos=markerPos, chrOffsets=chrOffsets, gene.gmb=gene.gmb, 
		plot.int=FALSE, plot.sig=TRUE)
	dev.off()

	writeLines(paste("Analysis finished: ", print.time(), sep=""))
}