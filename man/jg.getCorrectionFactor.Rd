% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/brundle.R
\name{jg.getCorrectionFactor}
\alias{jg.getCorrectionFactor}
\title{jg.getCorrectionFactor}
\usage{
jg.getCorrectionFactor(jg.experimentSampleSheet, jg.treatedNames,
  jg.untreatedNames)
}
\arguments{
\item{jg.experimentSampleSheet}{is the csv samplesheet used to load the data into DiffBind}

\item{jg.treatedNames}{is a list of the names of samples that are treated}

\item{jg.untreatedNames}{is a list of the names of samples that are untreated}
}
\description{
Generates a correction factor that is applied before reinserting the data
into the DiffBind object for analysis.
}
\examples{
data(jg.controlCountsTreated, package="Brundle")
data(jg.controlCountsUntreated, package="Brundle")
jg.coefficient<-jg.getNormalizationCoefficient(jg.controlCountsTreated, jg.controlCountsUntreated)


}
\keyword{DiffBind}
\keyword{correction}
\keyword{normalisation}
