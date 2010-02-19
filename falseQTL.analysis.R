#################
# AffyGG Script #
#################
# Only original CEL files are required
# Batch groups represent RNA synthesis groupings

# cis.buffer allows user to specify the Mb range that defines a cis QTL
# It will then look for the peak linkage markers within that range of the gene
# for creating the probe plot and executing the probe elimination analysis
# If cis.buffer is Null it will search genome wide

# probeset: character, probeset analysis will be performed on.

# genotype: csv file. Row names=marker names, Col1=Chr, Col2=Mb, remaing columns
# contain genotypes for individual strains

# strain: character, "lxs" or "bxd", used to select proper gene annotations
# cel.dir: character, contains directory location
# label: character object tacked on to results folder and other output
# batch: needs to be implented. Took out because caused errors on bach.

falseQTL.analysis<- function(probeset, geno.file, cel.dir, strain, label, cis.buffer, batch) {
  print.time<-function(){format(Sys.time(),"%D %l:%M %p")}
	writeLines(paste("Analysis started: ",print.time(), sep=""))
	
	# Set batch to null if none is provided
	if(missing(batch)){
		batch<-NULL
	}

	# Load necessary libraries
	require(affyGG)
	require(org.Mm.eg.db)
	chrs<-c(1:19,"X")

	# Convert strain to uppercase
	strain<-toupper(strain)

	# Load genotype data
	geno<-read.csv(geno.file, row.names=1)

	# Extract marker positions from genotype data
	if(length(grep("chr", colnames(geno), ignore.case=TRUE)>0)){
     markerPos<-data.frame(chr=geno[,grep("chr", colnames(geno), ignore.case=TRUE)],
      mb=geno[,grep("mb", colnames(geno), ignore.case=TRUE)],
      row.names=rownames(geno))
  } else {
      stop("Genotype file missing marker position columns.")
      }
	  
	# Add gmb to markerPos
	chrLengths<-org.Mm.egCHRLENGTHS
	chrLengths<-chrLengths[names(chrLengths) %in% chrs]/1000000
	chrOffsets<-c(0,cumsum(chrLengths)[1:19])
	names(chrOffsets)=chrs
	
	markerPos<-cbind(markerPos,gmb=markerPos$mb)
	for(i in 1:length(chrs)){
		markerPos$gmb[markerPos$chr==chrs[i]]<-markerPos$gmb[markerPos$chr==chrs[i]]+chrOffsets[i]
	}
	
	# Create genotype dataframe without marker positions
	geno<-data.frame(geno[,grep(strain,colnames(geno))])
  
	# Create genotype matrix
	#########################
	allele.counts<-table(unlist(unlist(geno)))
	# Sort alleles by frequency.
	allele.counts<-sort(allele.counts, dec=TRUE)
	alleles<-names(allele.counts)
	# Replace most common allele with 1's, second most common with 2's, etc...
	geno.mat<-apply(geno, MARGIN=2, FUN=function(x) 
			ifelse(x==alleles[1],1,	ifelse(x==alleles[2],2,3)))

	# LXS Data
	############
	writeLines("Loading data...")
	writeLines(print.time())
	if (strain=='LXS') {

    # Load webqtl probeset annotations
    chip.annot<-read.csv(file="/home/wolenar/RI_Datasets/LXS_datasets/webqtl_mouse430a2_annotations.csv",
      row.names=1)
		# Load chip database
		library(mouse430a2.db)
		# Get gene symbol
		gene.sym<-chip.annot[probeset, "Symbol"]
		# Get full gene description
		gene.description<-chip.annot[probeset, "Description"]
		# Get gene's chromosomal location
		gene.chr<-as.numeric(as.character(chip.annot[probeset, "Chr"]))
		# Get gene's chromosomal and genomic location
		gene.mb<-as.numeric(as.character(chip.annot[probeset, "Mb"]))
		gene.gmb<-as.numeric(as.character(chip.annot[probeset, "GMb"]))
		# Load chip probe info and get probe positions
		library(mouse430a2probe)
		probes.pos<-mouse430a2probe[mouse430a2probe$Probe.Set.Name==probeset,5]
		}
		
	# Prefix to add to all output files
	results.prefix<-paste(gene.sym, probeset, strain, label, sep="_")

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
	
	pdf(paste(results.prefix,"probePlot.pdf", sep="_"), height = 8, width = 10.5)
	myprobePlot(traits=traits, probeset=probeset, marker=peak.marker, 
		genotypes=geno.mat, alleles=alleles, probes.pos=probes.pos)
	dev.off()

	# Create custom probe-level QTL plot
	pdf(paste(results.prefix,"probe-level_qtlPlot.pdf", sep="_"), height = 8, width = 10.5)
	myprobeQTLplot(probeset=probeset, probeQtlProfiles=qtlmap.probe, markerPos=markerPos, chrOffsets=chrOffsets)
	dev.off()

	# Create custom probeset QTL interaction plot
	pdf(paste(results.prefix, "probeset_interaction_qtlPlot.pdf", sep="_"), height = 8, width = 10.5)
	myQTLintplot(probeset=probeset, probesetQtlProfile=qtlmap.Probeset, markerPos=markerPos, 
		chrOffsets=chrOffsets, gene.gmb=gene.gmb, plot.int=TRUE)
	dev.off()
	
	# Create custom probeset QTLplot
	pdf(paste(results.prefix, "probeset_qtlPlot.pdf", sep="_"), height = 8, width = 10.5)
	myQTLintplot(probeset=probeset, probesetQtlProfile=qtlmap.Probeset, markerPos=markerPos, 
		chrOffsets=chrOffsets, gene.gmb=gene.gmb, plot.int=FALSE, plot.sig=TRUE)
	dev.off()

	writeLines(paste("Analysis finished: ", print.time(), sep=""))
}