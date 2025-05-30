% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/funfea.R
\name{load_kog_annotation}
\alias{load_kog_annotation}
\title{Load COG/KOG Annotation}
\usage{
load_kog_annotation(path)
}
\arguments{
\item{path}{Path to COG/KOG annotation}
}
\value{
A dataframe with COG/KOG annotation information
}
\description{
This function loads a COG/KOG (Clusters of Orthologous Genes) annotation file into a dataframe.
}
\examples{
kog_annotation_df <- load_kog_annotation(path/to/kog/annotation)
}
