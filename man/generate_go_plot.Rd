% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/funfea.R
\name{generate_go_plot}
\alias{generate_go_plot}
\title{Generate GO Terms Enrichment Plot}
\usage{
generate_go_plot(
  enrichment_df,
  gotermType = "biological_process",
  n = 10,
  significant = "#E58601",
  plot_title = NA,
  plot_type = "bar"
)
}
\arguments{
\item{enrichment_df}{A list of dataframes with GO terms enrichment statistics per GO term type}

\item{gotermType}{A GO term type ("biological_process", "molecular_function" or "cellular_component"; default = "biological_process")}

\item{n}{A number of top-ranking GO terms to include in the plot (default = 10)}

\item{significant}{A color assignment for statistically significant GO terms (default = "#E58601")}

\item{plot_title}{A plot title (default = NA)}

\item{plot_type}{A plot type option ("bar" or "lollipop") (default = "bar")}
}
\value{
A ggplot2 bar plot object
}
\description{
This function generates a bar plot for visualization of GO terms enrichment statistics.
}
\examples{
generate_go_plot(go_enrichment_df)
generate_go_plot(go_enrichment_df, gotermType = "biological_process", n = 15, plot_title = "GO Enrichment: Biological Process")
generate_go_plot(go_enrichment_df, gotermType = "molecular_function")
generate_go_plot(go_enrichment_df, gotermType = "cellular_component")
generate_go_plot(go_enrichment_df, gotermType = "biological_process", n = 10, significant = "blue")
generate_go_plot(go_enrichment_df, gotermType = "biological_process", plot_type = "lollipop")
}
