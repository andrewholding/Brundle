pkgname <- "Brundle"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "Brundle-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('Brundle')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("Brundle")
### * Brundle

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Brundle
### Title: Brundle
### Aliases: Brundle
### Keywords: DESeq2 Diffbind

### ** Examples

data(dbaControl,package="Brundle")
data(dbaExperiment,package="Brundle")
Brundle(dbaExperiment,dbaControl,"Fulvestrant","none")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Brundle", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("jg.applyNormalisation")
### * jg.applyNormalisation

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: jg.applyNormalisation
### Title: jg.applyNormalisation
### Aliases: jg.applyNormalisation
### Keywords: correction normalisation peakset

### ** Examples

data(jg.experimentPeakset, package="Brundle")
jg.experimentPeaksetNormalised<-jg.applyNormalisation(jg.experimentPeakset,
                                                       1.267618,
                                                       0.6616886,
                                                       c("1b", "2b", "3b"))




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("jg.applyNormalisation", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("jg.convertPeakset")
### * jg.convertPeakset

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: jg.convertPeakset
### Title: jg.convertPeakset
### Aliases: jg.convertPeakset
### Keywords: Convert DESeq2 DiffBind

### ** Examples

jg.convertPeakset(jg.controlPeakset)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("jg.convertPeakset", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("jg.correctDBASizeFactors")
### * jg.correctDBASizeFactors

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: jg.correctDBASizeFactors
### Title: jg.correctDBASizeFactors
### Aliases: jg.correctDBASizeFactors
### Keywords: DESeq2 Diffbind

### ** Examples

data(jg.controlPeaksetDeSeq,package="Brundle")
data(dbaExperiment,package="Brundle")
jg.controlSizeFactors = estimateSizeFactorsForMatrix(jg.controlPeaksetDeSeq)
jg.correctDBASizeFactors(dbaExperiment,jg.controlSizeFactors)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("jg.correctDBASizeFactors", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("jg.dbaGetPeakset")
### * jg.dbaGetPeakset

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: jg.dbaGetPeakset
### Title: dbaGetPeakset
### Aliases: jg.dbaGetPeakset
### Keywords: DiffBind counts peakset

### ** Examples

data(dbaExperiment, package="Brundle")
jg.experimentPeakset <- jg.dbaGetPeakset(dbaExperiment)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("jg.dbaGetPeakset", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("jg.getControlCounts")
### * jg.getControlCounts

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: jg.getControlCounts
### Title: jg.getControlCounts
### Aliases: jg.getControlCounts
### Keywords: bam bamFiles reads

### ** Examples

data(jg.controlPeakset, package="Brundle")
fpath <- system.file("extdata", "samplesheet_SLX14438_hs_CTCF_DBA.csv",package="Brundle")
jg.controlSampleSheet<-fpath
jg.controlCountsTreated<-jg.getControlCounts(jg.controlPeakset, jg.controlSampleSheet,"Fulvestrant")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("jg.getControlCounts", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("jg.getCorrectionFactor")
### * jg.getCorrectionFactor

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: jg.getCorrectionFactor
### Title: jg.getCorrectionFactor
### Aliases: jg.getCorrectionFactor
### Keywords: DiffBind correction normalisation

### ** Examples

data(jg.controlCountsTreated, package="Brundle")
data(jg.controlCountsUntreated, package="Brundle")
jg.coefficient<-jg.getNormalizationCoefficient(jg.controlCountsTreated, jg.controlCountsUntreated)





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("jg.getCorrectionFactor", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("jg.getNormalizationCoefficient")
### * jg.getNormalizationCoefficient

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: jg.getNormalizationCoefficient
### Title: jg.getNormalizationCoefficient
### Aliases: jg.getNormalizationCoefficient
### Keywords: normalization

### ** Examples

data(jg.controlCountsTreated, package="Brundle")
data(jg.controlCountsUntreated, package="Brundle")
jg.coefficient<-jg.getNormalizationCoefficient(jg.controlCountsTreated,
                                                  jg.controlCountsUntreated)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("jg.getNormalizationCoefficient", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("jg.plotDeSeq")
### * jg.plotDeSeq

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: jg.plotDeSeq
### Title: jg.plotDeSeq
### Aliases: jg.plotDeSeq
### Keywords: DESeq2 data plot

### ** Examples

 data(jg.experimentResultsDeseq,package="Brundle")
 jg.plotDeSeq(jg.experimentResultsDeseq,
  p=0.01,
  title.main="Fold-change in ER binding",
  flip=TRUE
 )




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("jg.plotDeSeq", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("jg.plotDeSeqCombined")
### * jg.plotDeSeqCombined

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: jg.plotDeSeqCombined
### Title: jg.plotDeSeqCombined
### Aliases: jg.plotDeSeqCombined
### Keywords: DESeq2 data plot

### ** Examples

data(jg.controlResultsDeseq,package="Brundle")
data(jg.experimentResultsDeseq,package="Brundle")

jg.plotDeSeqCombined(jg.controlResultsDeseq,
                    jg.experimentResultsDeseq,
                    title.main="ER and CTCF Binding Folding changes on ER treatment",
                    p=0.01,flip=TRUE)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("jg.plotDeSeqCombined", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("jg.plotNormalization")
### * jg.plotNormalization

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: jg.plotNormalization
### Title: jg.plotNormalization
### Aliases: jg.plotNormalization
### Keywords: normalization plot

### ** Examples

data(jg.controlCountsTreated, package="Brundle")
data(jg.controlCountsUntreated, package="Brundle")
jg.plotNormalization(jg.controlCountsTreated, jg.controlCountsUntreated)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("jg.plotNormalization", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("jg.runDeSeq")
### * jg.runDeSeq

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: jg.runDeSeq
### Title: jg.runDeSeq
### Aliases: jg.runDeSeq
### Keywords: DESeq2

### ** Examples

data(jg.controlPeaksetDeSeq,package="Brundle")
data(jg.conditions,package="Brundle")
jg.controlSizeFactors = estimateSizeFactorsForMatrix(jg.controlPeaksetDeSeq)
jg.runDeSeq(jg.controlPeaksetDeSeq,jg.conditions, jg.SizeFactors = NULL)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("jg.runDeSeq", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
