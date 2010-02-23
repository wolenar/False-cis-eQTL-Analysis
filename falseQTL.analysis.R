probeset<-"1428600_at"
geno.file<-"/home/aaron/Dropbox/MilesLab Share/RI Datasets/BXD datasets/BXD Genotype data/BXD Genotypes.csv"
#cel.dir<-"/Users/milesg5users/Documents/Aaron/CEL_files/BXD CEL Files/PFC Saline"
cel.dir<-'/media/F82823B028236CB6/Users/Aaron Wolen/Documents/VCU/Miles Lab/BXD CEL Files/PFC Saline'
label<-"test"
strain<-"bxd"

qtlmap.probe<-read.csv(file="")

# Currently strain is only necessary to cidentify which columns in 
# genotype file contain actual data

falseQTL.analysis <- function(probeset, geno.file, cel.dir, strain, label, 
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
	
	# Variables
	######################
	
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
	     markerPos<-data.frame(chr=geno[,grep("chr", colnames(geno), ignore.case=TRUE)], 
	     mb=geno[,grep("mb", colnames(geno), ignore.case=TRUE)], row.names=rownames(geno))
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
	# Extract probe positions for probeset
	probes.pos<-data.frame(subset(eval(as.name(probe.pkg)),
		subset=Probe.Set.Name==probeset, select=Probe.Interrogation.Position))
	probes.pos<-as.numeric(probes.pos[,1])

	# Load array annotation data
	require(package=annot.pkg, character.only=TRUE)
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
	gene.start.mb<-
		get(probeset, envir=eval(as.name(paste(array.name,"CHRLOC",sep=""))))[1]
	gene.strand<-ifelse(gene.start.mb<0, "-", "+")
	gene.start.mb<-toMb(abs(gene.start.mb))
	# Get gene's chromosomal stop location
	gene.end.mb<-
		get(probeset, envir=eval(as.name(paste(array.name,"CHRLOCEND",sep=""))))[1]
	gene.end.mb<-toMb(abs(gene.end.mb))

	# Calculate gigabases
	########################	
	# Load chromosome lengths
	chrLengths<-eval(as.name(paste(array.name,"CHRLENGTHS",sep="")))
	chrLengths<-toMb(chrLengths[names(chrLengths) %in% chrs])
	chrOffsets<-c(0,cumsum(chrLengths)[1:19])
	names(chrOffsets)=chrs
	
	# Gene's genomic location
	gene.start.gmb <- gene.start.mb+chrOffsets[gene.chr]
	gene.end.gmb <- gene.end.mb+chrOffsets[gene.chr]
		
	# Gene's transcription start location
	gene.tx.mb<-ifelse(gene.strand=="+", gene.start.mb, gene.end.mb)	
	gene.tx.gmb<-gene.tx.mb+chrOffsets[gene.chr]
		
	# Add GMb to markerPos
	markerPos<-cbind(markerPos,gmb=markerPos$mb)
	for(i in 1:length(chrs)){
		markerPos$gmb[markerPos$chr==chrs[i]] <-
			markerPos$gmb[markerPos$chr==chrs[i]]+chrOffsets[i]
	}
		

	# Prefix to add to all output files
	results.prefix<-paste(gene.sym, probeset, label, sep="_")


	# Create a new folder in wd to store results of analysis
	# then set that as the working directory
	results.folder<-paste(results.prefix, "results", sep="_")
	system(paste("mkdir ", results.folder, sep=""))
	setwd(results.folder)
	
	
	# Reconciling genotype data and celfiles strains
	################################################

	keep.strains<-sapply(strains, 
		FUN=function(x) length(grep(paste(x,"\\D",sep=""), cel.names, perl=TRUE)))
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
	
	# Select probe-level data for one probe set
	traits <- as.matrix(probesignals[probesignals$probeset==probeset,2:ncol(probesignals)])
	colnames(traits)=names(keep.strains)
	rownames(traits)=paste(probeset,1:nrow(traits),sep="")
	write.csv(traits, file=paste(results.prefix, "probesignals.csv",sep="_"))
	
	# Perform QTL Analysis on probe level
	#####################################
	# Outputs a data.frame containing Fmarkerperprobe and
	# pmarkerperprobe for each probe in a probe-set at all markers
	writeLines(paste(print.time(),"Probe-level QTL analysis..."))
	qtlmap.probe<-qtlMap.xProbe(genotypes=geno.mat, traits=traits, batch=batch)
	## 8.17 min on quad-core desktop
	# Update data frame to include -log10 pmarker per probe values for mapping
	qtlmap.probe<-cbind(qtlmap.probe,minlog10pmarkerperprobe=-log10(qtlmap.probe[,4]))
	# Create spread sheet with probe-level QTL analysis results
	write.csv(qtlmap.probe, file=paste(results.prefix, "probe_level_QTL_results.csv", sep="_"))

	# Perform QTL analysis on Probeset-level
	########################################
	writeLines(paste(print.time(),"ProbeSet-level QTL analysis..."))
	qtlmap.probeset<-qtlMap.xProbeSet(genotypes=geno.mat, traits=traits, batch=batch)
	## 4.5 min on quad-core desktop
	qtlProbeset <- -log10(qtlmap.probeset$pmarker)
	intProbeset <- -log10(qtlmap.probeset$pinteraction)

	# Update data frame to include the -log10 values for mapping
	qtlmap.probeset<-cbind(qtlmap.probeset, minlog10pmarker=qtlProbeset,
		minlog10pinteraction=intProbeset)

	# Create spread sheet with probeset-level QTL analysis results
	write.csv(qtlmap.probeset, file=paste(results.prefix, "probeset_QTL_results.csv", sep="_"))

	# Identify marker where qtlProbeset is highest
	###############################################
	# If there's a tie for peak marker, select the marker closest to the genes
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
		markerPos$mb >= (gene.tx.mb - cis.buffer) &
		markerPos$mb <=(gene.tx.mb + cis.buffer))
	cis.markers<-rownames(cis.markers)
	
	# Subset of qtl data containing only cis markers
	cis.probeset<-qtlmap.probeset[match(cis.markers, qtlmap.probeset$marker),]
	
	peak.cis.marker<-as.character(cis.probeset$marker[which
		(cis.probeset$pmarker==min(cis.probeset$pmarker))])
	
	if(length(peak.cis.marker)>1){
		peak.cis.marker<-peak.cis.marker[
			which.min(abs(markerPos$gmb[match(peak.cis.marker, 
			rownames(markerPos))]-gene.tx.gmb))]
	}


	# Get peak marker's location
	peak.markerPos<-markerPos[peak.marker,]
	
	writeLines(paste("\nPeak marker for",probeset,"is",peak.marker,":\n",
		peak.markerPos,"\n"))


	# Eliminiate deviating probes using the statisical  model
	#########################################################
	writeLines(paste(print.time(),"Probe elimination analysis..."))
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
	
	# Export data
	########################

	# Plots
	########
	writeLines(paste(print.time(),"Generating graphs..."))
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

	writeLines(paste(print.time(),"Analysis complete. Huzzah!"))
	
}
