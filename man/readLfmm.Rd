\name{read.lfmm}
\alias{read.lfmm}

\title{Read files in the \code{\link{lfmm}} format}
\description{
Read a file in the \code{\link{lfmm}} format.
}
\usage{
read.lfmm(input.file)
}

\arguments{
\item{input.file}{
A character string containing a path to the input file,
a genotypic matrix in the \code{lfmm} format.
}
}
\value{
\item{R}{
A matrix containing the genotypes with one line per individual and one column 
per SNP.
}
}

\details{
Numeric records are read in blocks after checking that their column counts
agree. Spaces, tabs, blank lines, LF and CRLF line endings, and a final record
without a newline are supported. Missing-value codes 9 and -9 remain numeric
codes; \code{NA} is read as an R missing value. Numeric dosage values are
supported. Storage type and matrix names follow the historical reader.

Short files and nonstandard input, such as quoted fields, comments, and
header or row-name conventions, use the historical \code{read.table} path.
Malformed records are never repaired by repartitioning a flat token stream.
The input file must remain unchanged during reading.
}

\author{
    Eric Frichot
}
\seealso{
\code{\link{write.lfmm}}
\code{\link{lfmm.data}}
\code{\link{lfmm}}
\code{\link{geno2lfmm}}
\code{\link{lfmm2geno}}
\code{\link{ancestrymap2lfmm}}
\code{\link{ped2lfmm}}
}
\examples{
# tutorial contains a matrix of genotypes R with 1000 SNPs for 165 individuals.
# and a matrix with an environmental variable C.
data("tutorial")

# write R in a file called "genotypes.lfmm"
# Create file:    "genotypes.lfmm".
write.lfmm(tutorial.R,"genotypes.lfmm")

# read the file "genotypes.lfmm".
R = read.lfmm("genotypes.lfmm")
}


\keyword{read/write}

