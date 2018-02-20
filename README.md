# Brundle
Brundle is an R package that provides a series of functions for the normalisation of ChIP-Seq data
to internal or external controls. It can be installed from the [Brundle package on CRAN](https://CRAN.R-project.org/package=Brundle) using the install.packages("Brundle") command in R.

It is supported by the data package [BrudleData](https://github.com/andrewholding/BrundleData).

[Brundle_Example](https://github.com/andrewholding/Brundle_Example) provides worked examples and preprocessing scripts for people who want to use Brundle in their own research. The examples are also packaged as a [Docker Container](http://dockerhub.com/andrewholding/brundle) pre-installed with all the tools to run the examples and pipeline. You can find [instructions on running the container](https://github.com/andrewholding/Brundle_Example/blob/master/README.md#using-docker-container) in the [Brundle_Example](https://github.com/andrewholding/Brundle_Example/blob/master/README.md) readme. The Dockerfile and other relavent code for generating the container can be found in the [BrundleDocker](https://github.com/andrewholding/BrundleDocker) repository. 

## Quick Start

Below is a quick example that will generate a normalised MA plot.

```R
# Install Package
library(devtools)
install_github("andrewholding/Brundle")
# You may also need to install packages from BioConductor e.g. DiffBind

# Load package
library(Brundle)

# Load Example Diffbind Object
data(dbaExperiment,package="Brundle")
data(dbaControl,package="Brundle")

# Load Example Samplesheets
exptCSV <- system.file("extdata", "samplesheet_SLX14438_hs_ER_DBA.csv",   package="Brundle")
ctrlCSV <- system.file("extdata", "samplesheet_SLX14438_hs_CTCF_DBA.csv", package="Brundle")
jg.ExperimentSampleSheet<-exptCSV
jg.ControlSampleSheet<-ctrlCSV


# Normalise with Brundle
jg.experimentPeaksetNormalised<-Brundle(dbaExperiment,
                                        dbaControl,
                                        "Fulvestrant",
                                        "none",
                                        jg.ExperimentSampleSheet,
                                        jg.ControlSampleSheet,
                                        jg.noBAMs=TRUE
                                        )

# Insert data back into DiffBind object
dba <- DiffBind:::pv.resetCounts(dbaExperiment, jg.experimentPeaksetNormalised)

# Process with DiffBind as normal
dba<-dba.analyze(dba)
dba.plotMA(dba, bFlip=TRUE, bSmooth=FALSE)

```
## Workflow

![Workflow](https://cdn.rawgit.com/andrewholding/Brundle_Example/master/images/workflow.svg)

## Example Data on UCSC Gene Browser
[CTCF Spike-in data](http://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=andrewholding&hgS_otherUserSessionName=Brundle)

[H2Av Spike-in data - Human](https://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=andrewholding&hgS_otherUserSessionName=ER%2FH2av)

[H2Av Spike-in data -Drosophilia](https://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=andrewholding&hgS_otherUserSessionName=ER%2FH2av%20dm3)

[hsER/mmER Spike-in data - Human](https://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=andrewholding&hgS_otherUserSessionName=hsER%2FmmER)

[hsER/mmER Spike-in data - Mouse](https://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=andrewholding&hgS_otherUserSessionName=HsER%2FmmER%20mm9)

[CTCF Spike-in (+/-E2) ER data](https://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=andrewholding&hgS_otherUserSessionName=ER%2FCTCF)

[CTCF Spike-in (+/-E2) H4k12ac data](https://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=andrewholding&hgS_otherUserSessionName=H4K12ac%2FCTCF)

[ER+ Breast Cancer PDX Data]
(https://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=andrewholding&hgS_otherUserSessionName=PDX)
