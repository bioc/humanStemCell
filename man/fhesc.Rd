\name{fhesc}
\Rdversion{1.1}
\alias{fhesc}
\docType{data}
\title{
Data from a simple experiment on Human stem cells.
}
\description{
Human stem cells were assayed using Affymetrix 133plus 2 arrays. There
were six arrays, three were biological replicates for undifferentiated
cells, the other three were biological replicates for differentiated cells.
}
\usage{data(fhesc)}
\format{
The data are in the form of an \code{ExpressionSet} instance.
}
\details{
Human Embryonic Stem Cells, H1 Line were cultured under feeder-free conditions.
Undifferentiated samples were taken from this pool. The differentiated samples
were obtained by maintaining the cells in culture for 10 - 14 days in the 
absence of basic fibroblast growth factor and conditioned medium.
}
\source{
The data were obtained from Dr. M. Tewari.
}
\references{
These data were used to prepare the book chapter, 
R and Bioconductor packages in bioinformatics: towards
systems biology, by Nolwenn LeMeur, Michael Lawrence, Merav Bar, 
Muneesh Tewari and Robert Gentleman
}
\examples{
data(fhesc)
}
\keyword{datasets}
