# load the predefined functions that are in the noe_source.R file
# in the current working folder of R (that should be your project folder)
source("noe_source.R")

###############################################################################
## READ NIMBLEGEN DATA
###############################################################################

# read the sample description file
# the file is a tab-delimited table. columns "file", "sample.type" and "sample.name" are
# required for this protocol
samples <- read.delim("sample_key.txt", header=T, comment.char="#", stringsAsFactors=F)

# display sample description file
samples

# read signal intensities from all sample files into a matrix, only experimental signals are kept
# matrix columns are signals from indiv. files. column names correspond to file names.
# matrix rows correspond to probe ids
signals <- read.pairs(samples$file, sub(".pair", "", samples$file))

# display the first rows of the signal intensity table
head(signals)

# get layout name from nimblegen raw data file
layout.name <- get.layout.info(samples$file[1])

# read layout and position information into data frame, only experimental signals are kept
# instead of passing layout.name one could enter the name of layout directly
layout <- read.layout(layout.name)

# display the first rows of the layout table 
head(layout)

# save the tables created so far in R objects that can rapidly be re-loaded into objects of the same
# name using load("whatever.rda")
save(samples, file="samples.rda")
save(signals, file="signals.rda")
save(layout, file="layout.rda")

###############################################################################
## INSPECTING THE RAW DATA, QUALITY CONTROL THE RAW DATA 
###############################################################################

# density plot of the raw signals of the 1st channel
# for the second channels use signals[,2], the third signals[,3] and so on
plot(density(signals[,1]))

# prepare a layout that places 2 plots into a single row
par(mfrow=c(1,2))
# density plot of the log2 transformed signals of the 1st channel
plot(density(log2(signals[,1])))
# qq plot of the log2 transformed signals of the 1st channel
qqnorm(log2(signals[,1]), pch=".")
qqline(log2(signals[,1]), lty=2, col=2)

par(mfrow=c(1,2))
# A is the average signal intensity of each probe of the 2 channels of the 1st array
A = (log2(signals[,2]) + log2(signals[,1]))/2
# M is the signal ratio log2(channel2/channel1) each probe of the 1st array
M = (log2(signals[,2]) - log2(signals[,1]))

require(affy)
# MA plot
ma.plot(A, M, plot.method="smoothScatter", cex=0.8, main="MA plot")

require(vsn)
# SD versus signal strength plot
meanSdPlot(log2(signals[,1:2]), ranks=TRUE, main="SD versus rank")

# all quality plots for array 145215 on one page
quality.control(signals, samples, "145215")

# array image reconstruction for the 1st channel
nimble.image(log2(signals[,1]),layout)

# pairwise correlations of all channels
correlate.batch(log2(signals))

###############################################################################
## NORMALIZATION
###############################################################################

# load the vsn package needed for normalisation
require(vsn)

# normalisation of all channels using the vsn2 function
normalized.signals <- exprs(vsn2(signals, subsample=20000, verbose=F))

# signal distribution boxplots before and after normalisation
par(mfrow=c(1,2))
boxplot(log2(signals)~col(signals), main="before normalization")
boxplot(normalized.signals~col(normalized.signals), main="after normalization")

# MA before and after normalization displayed for the 3rd array (channels 5 and 6)
require(affy)
par(mfrow=c(1,2))
A = (log2(signals[,6]) + log2(signals[,5]))/2
M = (log2(signals[,6]) - log2(signals[,5]))
ma.plot(A, M, plot.method="smoothScatter", main="before normalization", cex=0.8)
A = (normalized.signals[,6] + normalized.signals[,5])/2
M = (normalized.signals[,6] - normalized.signals[,5])
ma.plot(A, M, plot.method="smoothScatter", main="after normalization", cex=0.8)

###############################################################################
## PROBE LEVEL STATISTICS
###############################################################################

# compute the mean ratios of all 3 replicates, takes dye swaps into account
pmean <- probe.mean(normalized.signals, samples$sample.type)
# density plot of mean ratios
plot(density(pmean), main="mean")
							
# compute the sam statistic for the probes
pstat <- probe.stat(normalized.signals, samples$sample.type)
plot(density(pstat), main="sam statistic")

require(locfdr)
par(mfrow=c(1,1))
# compute fdr for each probe
# generates a plot to display the fits
lfdr <- locfdr(pstat)$fdr

# lets define an fdr cutoff of 0.02
cutoff <- lfdr < 0.02

require(geneplotter)
# generation of a stacked histogram to show the ratios of probes above and below cutoff
x <- list(pmean[cutoff], pmean[!cutoff])
histStack(x, breaks=seq(floor(min(pmean)), ceiling(max(pmean)), length=51), col=c(2,0),
	xlab="mean log2(ratio)", main="fdr threshold < 0.02")
legend("topright", c("enriched", "non-enriched"), fill=c(2,0), bty="n")

###############################################################################
## REGION SUMMARY
###############################################################################

# compute HMM states for each probe based on the sam statistcs (pstat)
hmm.out <- stateHMM(pstat, layout)
# retrieve the state for each probe
pstate <- hmm.out$probe.state
# get a table of bound regions
regions <- hmm.out$regions

# generation of a stacked histogram to show the ratios of enriched and non-enriched probes
x <- list(pmean[pstate==1], pmean[pstate!=1])
histStack(x, breaks=seq(floor(min(pmean)), ceiling(max(pmean)), length=51), col=c(2,0),
	xlab="mean log2(ratio)", main="hmm enriched")
legend("topright", c("enriched", "non-enriched"), fill=c(2,0), bty="n")


###############################################################################
## EXPORT DATA
###############################################################################

# create a file in wiggle format that contains the mean signal ratios of all replicates
write.wiggle(pmean, layout, "mean.wig", descr="mean")
# create a file in wiggle format that contains the sam statistic
write.wiggle(pstat, layout, "stat.wig", descr="stat")
# create a file in wiggle format that contains the HMM states (1 - enriched, 0 - non-enriched)
write.wiggle(pstate, layout, "state.wig", descr="hmmstate")
# create a gff file of the bound regions identified by tileHMM
write.gff.regions(regions, "bound_regions.gff")

