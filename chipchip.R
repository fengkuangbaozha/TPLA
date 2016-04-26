



##########prepare raw data############
setwd("~/sunyd/identify/tomato_rnaseq/SRRm82/chip-chip/GSE49125")
print(list.files())
source("/psc/bioinformatics/sunyd/identify/script/chipchip_source.R")
samples <- read.delim("sample_key.txt", header=T,comment.char="#", stringsAsFactors=F) ###read in sample info
signals <- read.pairs(samples$file,samples$sample.name)              #########read in the paired text files,get the probe expression info
layout.name <- get.layout.info(samples$file[1])              ######read in the ndf and pos files
layout <- read.layout(layout.name)                         ##########read in the probe id info

#########normalization############
require(vsn)
normalized.signals <- exprs(vsn2(signals,subsample=20000, verbose=F))

#############Identification of enriched probes###########
pmean <- probe.mean(normalized.signals,samples$sample.type)    
pstat <- probe.stat(normalized.signals,samples$sample.type)

###########Identification of enriched regions###########
hmm.out <- stateHMM(pstat, layout)   ######get the position and the score for regions
pstate <- hmm.out$probe.state
regions <- hmm.out$regions
write.gff.regions(regions, "bound_regions.gff")
##x <- list(pmean[state==1], pmean[state!=1])
