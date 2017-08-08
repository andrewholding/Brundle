

######################
#
# Brundle - ChIP-seq Normalisation Convienice Functions
# Andrew Holding 2017
#
######################

######################
#
# Functions
#
######################

#' jg.countAlignedMReads
#'
#' This function counts the number of aligned reads in millions from a list of
#' bam files. It returns these in as a list of numbers in the same order.
#' @param jg.bamFiles is a list of bam files to count.
#' @keywords bamFiles bam reads
#' @import Rsamtools

jg.countAlignedMReads <- function(jg.bamFiles) {
    jg.counts <- numeric()
    for (jg.bam in jg.bamFiles) {
        jg.counts <-
            append(jg.counts, colSums(idxstatsBam(jg.bam)["mapped"]) / 1E6)
    }
    return(jg.counts)
}

#' jg.getControlCounts
#'
#' This function counts the number of aligned reads in millions from a list of
#' bam files. It returns these in as a list of numbers in the same order.
#' @param jg.control is a peakset extracted by jg.dbaGetPeakset.
#' @param jg.controlSampleSheet is the samplesheet supplied to DiffBind to
#' generate jg.control.
#' @param jg.Condition is the condition we the counts for as specficied in the
#'  samplesheet.
#' @keywords bamFiles bam reads
#' @export
#' @examples
#' data(jg.controlPeakset, package="Brundle")
#' fpath <- system.file("extdata", "samplesheet_SLX14438_hs_CTCF_DBA.csv",package="Brundle")
#' jg.controlSampleSheet<-fpath
#' jg.controlCountsTreated<-jg.getControlCounts(jg.controlPeakset, jg.controlSampleSheet,"Fulvestrant")


jg.getControlCounts <-
    function(jg.control,
             jg.controlSampleSheet,
             jg.Condition)
    {
        jg.controlCounts <- jg.control[, -c(1:3)]
        temp <-
            read.csv(file = jg.controlSampleSheet,
                     header = TRUE,
                     sep = ",")['Condition'] == jg.Condition
        return(jg.controlCounts[, temp])
    }

#' jg.plotNormalization
#'
#' This function allows the user to visualize the normalisation. It is not
#' needed for the pipeline but provides a helpful illustration of the process.
#' @param jg.controlCountsTreated Control counts extracted from the Diffbind object for the treated condition using jg.getControlCounts
#' @param jg.controlCountsUntreated Control ounts extracted from the Diffbind object for the untreated condition using jg.getControlCounts
#' @keywords plot normalization
#' @export
#' @examples
#' data(jg.controlCountsTreated, package="Brundle")
#' data(jg.controlCountsUntreated, package="Brundle")
#' jg.plotNormalization(jg.controlCountsTreated, jg.controlCountsUntreated)


jg.plotNormalization <-
    function(jg.controlCountsTreated,
             jg.controlCountsUntreated)
    {
        plot(
            rowMeans(jg.controlCountsTreated),
            rowMeans(jg.controlCountsUntreated),
            pch = 20,
            xlab = "Counts in peak after treatment" ,
            ylab = "Counts in peak before treatment" ,
            main = "Comparision of Counts in peaks"
        )
        lm1 <-
            lm(rowMeans(jg.controlCountsUntreated) ~ 0 + rowMeans(jg.controlCountsTreated))

        abline(c(0, lm1$coef), col = "red3")
        print(lm1$coefficients)
        angularcoeff <- lm1$coef[1]

        points(
            rowMeans(jg.controlCountsTreated) * angularcoeff,
            rowMeans(jg.controlCountsUntreated),
            pch = 20,
            col = "royalblue3"
        )
        treatment_fit <- rowMeans(jg.controlCountsTreated) * angularcoeff
        lm1 <- lm(treatment_fit ~ 0 + rowMeans(jg.controlCountsUntreated))
        abline(c(0, lm1$coef), col = "purple")
        legend(
            "topleft",
            legend = c("Raw", "Normalised"),
            pch = 20,
            col = c("black", "royalblue3")
        )
    }

#' jg.getNormalizationCoefficient
#'
#' This function allows the user to caryy out the normalisation and returns a
#' coefficient by using a linear fit to the control data.
#' @param jg.controlCountsTreated Control counts extracted from the Diffbind object for the treated condition using jg.getControlCounts
#' @param jg.controlCountsUntreated Control ounts extracted from the Diffbind object for the untreated condition using jg.getControlCounts
#' @keywords normalization
#' @export
#' @examples
#' data(jg.controlCountsTreated, package="Brundle")
#' data(jg.controlCountsUntreated, package="Brundle")
#' jg.coefficient<-jg.getNormalizationCoefficient(jg.controlCountsTreated,
#'                                                   jg.controlCountsUntreated)

jg.getNormalizationCoefficient <-
    function(jg.controlCountsTreated,
             jg.controlCountsUntreated)
    {
        lm1 <-
            lm(rowMeans(jg.controlCountsUntreated) ~ 0 + rowMeans(jg.controlCountsTreated))
        return(lm1$coef[1])
    }


#' jg.plotMA
#'
#' This function plots both the control and experimental data on an MA plot.
#' It also allows for the user to provide a normalisation coefficent for the data.
#' @param jg.experimentPeakset is the peakset of the experimental data extracted from a DiffBind ojbect with jg.dbaGetPeakset
#' @param jg.controlPeakset is the peakset of the control data extracted from a DiffBind ojbect with jg.dbaGetPeakset
#' @param jg.untreatedNames is a list of sample names for the control or untreated conditions
#' @param jg.treatedNames is a list of sample samples for the treated conditions
#' @param jg.coefficient is a normalisation coefficent for the data that can be generated via the pipeline. Can be set to 1 to view before normalisation.
#' @param
#' @keywords plot normalization

jg.plotMA <-
    function(jg.experimentPeakset,
             jg.controlPeakset,
             jg.untreatedNames,
             jg.treatedNames,
             jg.coefficient)
    {
        M_corrected <- apply(jg.experimentPeakset[-c(1:3)], 1, function(x) {
            untreated <- mean(x[jg.untreatedNames])
            treated <- jg.coefficient * mean(x[jg.treatedNames])
            fc <- mean(treated) / mean(untreated)
            log2fc <- log2(fc)
            return(log2fc)
        })
        A_corrected <- apply(jg.experimentPeakset[-c(1:3)], 1, function(x) {
            untreated <- mean(x[jg.untreatedNames])
            treated <- jg.coefficient * mean(x[jg.treatedNames])
            return(log10(sum(treated + untreated)))
        })


        M_dm_corrected <- apply(jg.controlPeakset[-c(1:3)], 1, function(x) {
            untreated <- mean(x[jg.untreatedNames])
            treated <- jg.coefficient * mean(x[jg.treatedNames])
            fc <- mean(treated) / mean(untreated)
            log2fc <- log2(fc)
            return(log2fc)
        })


        A_dm_corrected <- apply(jg.controlPeakset[-c(1:3)], 1, function(x) {
            untreated <- mean(x[jg.untreatedNames])
            treated <- jg.coefficient * mean(x[jg.treatedNames])
            return(log10(sum(treated + untreated)))
        })

        plot(
            A_corrected,
            M_corrected,
            pch = 20,
            xlab = "A, log10(counts)",
            ylab = "M, log2FC(treatment)",
            main = "Normalised aligned reads"
        )
        points(A_dm_corrected,
               M_dm_corrected,
               pch = 20,
               col = "cornflowerblue")
        lm1 <- lm(M_dm_corrected ~ A_dm_corrected)
        abline(lm1$coef, col = "red4")
        abline(h = 0)

    }


#' jg.getDba
#'
#' Generates a DiffBind object from a valid SampleSheet with the required
#' data for normalisation. No examples are provided as BAM files are not
#' included in this package.
#'
#' @param jg.experimentSampleSheet is the filename of samplesheet to be loaded
#' @param dbaSummits is the peak width in bp from summits (optional)
#' @keywords DiffBind counts load samplesheet
#' @export
#' @import DiffBind


jg.getDba <- function (jg.experimentSampleSheet,
                       dbaSummits,
                       ...)
{
    dba <- dba(sampleSheet = jg.experimentSampleSheet)
    if (exists("dbaSummits"))
    {
        dba <- dba.count(dba, summits = dbaSummits, ...)
    } else {
        dba <- dba.count(dba)
    }
    dba <- dba.count(dba, peaks = NULL, score = DBA_SCORE_READS)
    return(dba)
}

#' dbaGetPeakset
#'
#' Extracts a peakset from a dba object.
#'
#' @param dba is the name of the DiffBind object
#' @keywords DiffBind counts peakset
#' @export
#' @examples
#' data(dbaExperiment, package="Brundle")
#' jg.experimentPeakset <- jg.dbaGetPeakset(dbaExperiment)
#' @import DiffBind

jg.dbaGetPeakset <- function(dba)
{
    jg.peakset <-
        dba.peakset(dba, bRetrieve = TRUE, DataType = DBA_DATA_FRAME)
    #Correct sample names back to that in sample sheet,
    #as DiffBind changes them on export.

    jg.sampleIds <- dba$samples[, 'SampleID']
    names(jg.peakset) <- c("CHR", "START", "END", jg.sampleIds)
    return(jg.peakset)
}

#' jg.getSampleIds
#'
#' Extracts the sample Id from DiffBind formatted SampleSheet in csv format.
#'
#' @param jg.controlSampleSheet is the filename of the samplesheet
#' @keywords DiffBind samplesheet sample

jg.getSampleIds <- function(jg.controlSampleSheet)
{
    jg.sampleIds <-
        as.character(read.csv(
            file = jg.controlSampleSheet,
            header = TRUE,
            sep = ","
        )[, 'SampleID'])
    return(jg.sampleIds)
}

#' jg.getCorrectionFactor
#'
#' Generates a correction factor that is applied before reinserting the data
#' into the DiffBind object for analysis.
#'
#' @param jg.experimentSampleSheet is the csv samplesheet used to load the data into DiffBind
#' @param jg.treatedNames is a list of the names of samples that are treated
#' @param jg.untreatedNames is a list of the names of samples that are untreated
#' @keywords DiffBind correction normalisation
#' @export
#' @examples
#' data(jg.controlCountsTreated, package="Brundle")
#' data(jg.controlCountsUntreated, package="Brundle")
#' jg.coefficient<-jg.getNormalizationCoefficient(jg.controlCountsTreated, jg.controlCountsUntreated)
#'
#'


jg.getCorrectionFactor <-
    function (jg.experimentSampleSheet,
              jg.treatedNames,
              jg.untreatedNames)
    {
        #Load aligned reads for experiment Bams.

        #Get list of experimentBams
        jg.experimentBams <-
            as.character(read.csv(
                file = jg.experimentSampleSheet,
                header = TRUE,
                sep = ","
            )[, 'bamReads'])
        jg.sampleIds <-
            as.character(read.csv(
                file = jg.experimentSampleSheet,
                header = TRUE,
                sep = ","
            )[, 'SampleID'])
        jg.experimentAligned <- jg.countAlignedMReads(jg.experimentBams)
        names(jg.experimentAligned) <- jg.sampleIds
        #Take ratio of treated:untreated aligned reads and create correction factor.
        jg.correctionFactor <-
            sum(jg.experimentAligned[jg.treatedNames]) / sum(jg.experimentAligned[jg.untreatedNames])
        return(jg.correctionFactor)
    }


#' jg.applyNormalisation
#'
#' Takes the experimental peakset and applies the calculated coefficient and
#' correction factor.
#'
#' @param jg.experimentPeakset is the peakset extracted from the Diffbind object
#' @param jg.coefficient is the coefficient calculated by jg.getNormalizationCoefficient
#' @param jg.correctionFactor is the correction factor calculated by jg.getCorrectionFactor
#' @param jg.treatedNames is the names of the treated samples
#' @keywords peakset correction normalisation
#' @export
#' @examples
#' data(jg.experimentPeakset, package="Brundle")
#' jg.experimentPeaksetNormalised<-jg.applyNormalisation(jg.experimentPeakset, 1.267618, 0.6616886, c("1b", "2b", "3b"))
#'

jg.applyNormalisation <-
    function(jg.experimentPeakset,
             jg.coefficient,
             jg.correctionFactor,
             jg.treatedNames)
    {
        jg.experimentPeaksetNormalised <- jg.experimentPeakset
        jg.experimentPeaksetNormalised[jg.treatedNames] <-
            (jg.coefficient * jg.correctionFactor * jg.experimentPeakset[jg.treatedNames])
        return(jg.experimentPeaksetNormalised)
    }


#' jg.plotDeSeq
#'
#' Plots the output from DESeq2 for the Brundle pipeline
#'
#' @param ma.df is the result Dataframe from DESeq2
#' @param p is the minium FDR to highlight as significant
#' @param title is the plot title
#' @param log2fold is the minimum log2 fold change for highlighted points
#' @param flip when set to TRUE flips the data
#' @keywords DESeq2 data plot
#' @export
#' @examples
#'  data(jg.experimentResultsDeseq,package="Brundle")
#'  jg.plotDeSeq(jg.experimentResultsDeseq,
#'   p=0.01,
#'   title.main="Fold-change in ER binding",
#'   flip=TRUE
#'  )
#'
#' @import lattice
#'
jg.plotDeSeq <-
    function(ma.df,
             p = 0.01,
             title.main = "Differential ChIP",
             log2fold = 0.5,
             flip = FALSE)
    {



        if (flip == TRUE)
        {
            ma.df$log2FoldChange <- -ma.df$log2FoldChange
        }

        xyplot(
            ma.df$log2FoldChange ~ log(ma.df$baseMean, base = 10),
            groups = (
                ma.df$padj < p &
                    abs(ma.df$log2FoldChange) > log2fold & !is.na(ma.df$padj)
            ),
            col = c("black", "red"),
            main = title.main,
            scales = "free",
            aspect = 1,
            pch = 20,
            cex = 0.5,
            ylab = expression("log"[2] ~ "ChIP fold change"),
            xlab = expression("log"[10] ~ "Mean of Normalized Counts"),
            par.settings = list(
                par.xlab.text = list(cex = 1.1, font = 2),
                par.ylab.text = list(cex = 1.1, font = 2)
            )
        )


    }

#' jg.plotDeSeqCombined
#'
#' Overlays the plots from the output from DESeq2 for the Brundle pipeline
#'
#' @param jg.controlResultsDeseq is the result Dataframe from DESeq2 for the control conditions
#' @param jg.experimentResultsDeseqis the result Dataframe from DESeq2 for the experimental conditions
#' @param title.main is the plot title
#' @param padjX is the minium FDR to highlight as significant
#' @param flip when set to TRUE flips the data
#' @keywords DESeq2 data plot
#' @export
#' @examples
#' data(jg.controlResultsDeseq,package="Brundle")
#' data(jg.experimentResultsDeseq,package="Brundle")
#'
#' jg.plotDeSeqCombined(jg.controlResultsDeseq,
#'                     jg.experimentResultsDeseq,
#'                     title.main="ER and CTCF Binding Folding changes on ER treatment",
#'                     p=0.01,flip=TRUE)
#'
#' @import lattice

jg.plotDeSeqCombined <-
    function(jg.controlResultsDeseq,
             jg.experimentResultsDeseq,
             title.main,
             padjX,
             flip = FALSE)
    {
        jg.controlResultsDeseq$group = 'a'
        jg.experimentResultsDeseq$group = 'b'

        if (flip == TRUE)
        {
            jg.controlResultsDeseq$log2FoldChange <-
                -jg.controlResultsDeseq$log2FoldChange
            jg.experimentResultsDeseq$log2FoldChange <-
                -jg.experimentResultsDeseq$log2FoldChange
        }

        for (i in 1:length(jg.experimentResultsDeseq$group)) {
            if (!is.na(jg.experimentResultsDeseq$padj[i]) &
                !is.na(jg.experimentResultsDeseq$log2FoldChange[i]) &
                jg.experimentResultsDeseq$padj[i] < padjX &
                jg.experimentResultsDeseq$log2FoldChange[i] < 0) {
                jg.experimentResultsDeseq$group[i] <- 'd'
            }
            else if (!is.na(jg.experimentResultsDeseq$padj[i]) &
                     !is.na(jg.experimentResultsDeseq$log2FoldChange[i]) &
                     jg.experimentResultsDeseq$padj[i] < padjX &
                     jg.experimentResultsDeseq$log2FoldChange[i] > 0) {
                jg.experimentResultsDeseq$group[i] <- 'c'
            }
        }

        for (i in 1:length(jg.controlResultsDeseq$group)) {
            if (!is.na(jg.controlResultsDeseq$padj[i]) &
                !is.na(jg.controlResultsDeseq$log2FoldChange[i]) &
                jg.controlResultsDeseq$padj[i] < padjX &
                jg.controlResultsDeseq$log2FoldChange[i] < 0) {
                jg.controlResultsDeseq$group[i] <- 'f'
            }
            else if (!is.na(jg.controlResultsDeseq$padj[i]) &
                     !is.na(jg.controlResultsDeseq$log2FoldChange[i]) &
                     jg.controlResultsDeseq$padj[i] < padjX &
                     jg.controlResultsDeseq$log2FoldChange[i] > 0) {
                jg.controlResultsDeseq$group[i] <- 'e'
            }
        }

        full.res = rbind(jg.controlResultsDeseq, jg.experimentResultsDeseq)

        colours<-
            c(
                "grey40",
                "grey80",
                "#5480ff",
                "#ff5454",
                "#08298a",
                "#750505"
            )
        names(colours)<-c("a","b","c","d","e","f")
        colours<-colours[sort(unique(full.res$group))]

        xyplot(
            full.res$log2FoldChange ~ log(full.res$baseMean, base = 10),
            data = full.res,
            groups = full.res$group,
            col = colours,
            ylab = expression('log'[2] * ' Differential ChIP'),
            xlab = expression("log"[10] ~ "Mean of Normalized Counts"),
            aspect = 1.0,
            pch = 16,
            cex = 0.5,
            main = title.main,
            scales = list(
                x = list(cex = 0.8, relation = "free"),
                y = list(cex = 0.8, relation = "free")
            ),
            between = list(y = 0.5, x = 0.5),
            auto.key = TRUE,
            key = list(
                corner = c(0, 1),
                cex = 0.75,
                points = list(
                    col = c(
                        "white",
                        "gray80",
                        "gray40",
                        "#ff5454",
                        "#5480ff",
                        "#750505",
                        "#08298a"
                    ),
                    pch = 20
                ),
                text = list(
                    c(
                        " ",
                        "Target Binding",
                        "Control Binding",
                        "Target Decreased",
                        "Target Increased",
                        "Control Decreased",
                        "Control Increased"
                    )
                )
            )
        )
    }

#' jg.convertPeakset
#'
#' Converts a DiffBind object into a DESeq2 compatible form for the workflow.
#'
#' @param jg.controlPeakset is the name of the DiffBind object to convert
#' @keywords DESeq2 DiffBind Convert
#' @export
#' @examples jg.convertPeakset(jg.controlPeakset)

jg.convertPeakset <- function(jg.controlPeakset)
{
    jg.controlPeaksetDeSeq <- jg.controlPeakset[-c(1:3)]
    row.names(jg.controlPeaksetDeSeq) <-
        paste(jg.controlPeakset[, 1],
              ':',
              jg.controlPeakset[, 2],
              '-',
              jg.controlPeakset[, 3],
              sep = '')

    return(jg.controlPeaksetDeSeq)

}

#' jg.runDeSeq
#'
#' Runs DESeq2 on our peakset after we have obtained the normalised size factors.
#'
#' @param jg.PeaksetDeSeq is the experimental peakset formatted for DESeq2
#' @param jg.conditions is the list of conditions to be compared
#' @param jg.SizeFactors is the size factors generated from the control samples
#' @keywords DESeq2
#' @export
#' @examples jg.runDeSeq(jg.PeaksetDeSeq,jg.conditions, jg.SizeFactors = NULL)
#' @import DESeq2
jg.runDeSeq <-
    function(jg.PeaksetDeSeq,
             jg.conditions,
             jg.SizeFactors = NULL)
    {
        jg.DeSeq = DESeqDataSetFromMatrix(jg.PeaksetDeSeq, jg.conditions, ~ Condition)
        if (is.null(jg.SizeFactors))
        {
            jg.DeSeq = estimateSizeFactors(jg.DeSeq)
        } else {
            sizeFactors(jg.DeSeq) <- jg.SizeFactors
        }
        jg.DeSeq = estimateDispersions(jg.DeSeq)
        jg.DeSeq = nbinomWaldTest(jg.DeSeq)
        return(jg.DeSeq)
    }

#' jg.correctDBASizeFactors
#'
#' Correct the size factors in a DiffBind object using our DESeq2 pipeline for
#' normalisation.
#'
#' @param dba
#' @param jg.controlSizeFactors
#' @keywords DESeq2 Diffbind
#' @export
#' @examples jg.correctDBASizeFactors(dba,jg.controlSizeFactors)


jg.correctDBASizeFactors <- function(dba, jg.controlSizeFactors)
{
    jg.libsizes <- as.numeric(dba$class["Reads", ])
    names(jg.libsizes) <- names(dba$class["Reads", ])

    dba.correctedReads <- jg.controlSizeFactors

    dba$class["Reads", ] <- dba.correctedReads
    return(dba)
}
