% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read.process.R
\name{read.process}
\alias{read.process}
\title{Read and Process Methylation Data into Bins}
\usage{
read.process(start.pos = 1, end.pos = max(pos), pos, norm.m, norm.um, los.m,
  los.um, bin.size)
}
\arguments{
\item{start.pos}{Starting position for data analysis. Default to 1.}

\item{end.pos}{End position for data analysis. Default to the maximum position in the parameter 'pos'.}

\item{pos}{CpG positions}

\item{norm.m}{A matrix contains number of methylated reads of the normal group, each column is a biological sample}

\item{norm.um}{A matrix contains number of un-methylated reads of the normal group}

\item{los.m}{A matrix contains number of methylated reads of the abnormal group}

\item{los.um}{A matrix contains number of un-methylated reads of the abnormal group}

\item{bin.size}{User-defined bin size. Default to 40}
}
\description{
This function reads in CpG methylation data and transforms observations into bins.
Users define starting position, end position, and bin size.
}