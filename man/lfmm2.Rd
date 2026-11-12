\name{lfmm2}
\alias{lfmm2}
\alias{lfmm2.test,lfmm2Class-method}




\title{Fitting Latent Factor Mixed Models (Least squares algorithm)}
\description{Latent Factor Mixed Models (LFMMs) are factor regression models in which the response variable is a genotypic matrix, and the explanatory variables are environmental measures of ecological interest or trait values. The \code{\link{lfmm2}} function estimates latent factors based on an exact least-squares approach. The resulting object can be used by the function \code{\link{lfmm2.test}} to identify genetic polymorphisms exhibiting association with ecological gradients or phenotypes, while correcting for unobserved confounders. An MCMC estimation algorithm is implemented in the function \code{\link{lfmm}}, but this version should be prefered.
}

\usage{
lfmm2 (input, env, K, lambda = 1e-5, effect.sizes = FALSE, formula = NULL)
}

\arguments{
\item{input}{A genotypic matrix or a character string containing a path to the input file. The genotypic  matrix must be in the \code{\link{lfmm}{lfmm_format}} format without missing values (9 or \code{NA}). See \code{\link{impute}} for completion based on nonnegative matrix factorization and consider R packages for reading large matrices.
}
\item{env}{A numeric matrix, a data frame containing only ordinary numeric columns, or a character string containing a path to the environmental file. The environment matrix must be in the \code{\link{env}} format without missing values. Response variables must be encoded as \code{numeric}.
}
\item{K}{
    An integer corresponding to the number of latent factors.  The number of latent factors  could be estimated from the elbow point in the PCA screeplot for the genotype matrix.
}
\item{lambda}{
    A positive numeric value for a ridge regularization parameter. The default value is set to 1e-5.
}
\item{effect.sizes}{ 
    A logical value that indicates whether the matrix of effect sizes should be returned or not. The default value is set to \code{FALSE} for saving memory space. 
}
\item{formula}{An optional one-sided formula, allowed only with a numeric data frame in \code{env}. The default \code{NULL} uses every column in its original order. Supported terms are numeric variables, \code{I(x^k)} for positive literal integer powers, \code{log(x)} and \code{sqrt(x)}, joined by \code{+}. Interactions, factors, arbitrary calls and \code{.} are not supported. Non-numeric data-frame columns cause an error even if not selected.}
}


\details{
Formula transformations are applied before the existing column centering.
No intercept column is added to the environmental matrix. This does not
remove the intercept from the association regressions in \code{lfmm2.test}.
The prepared matrix is not standardized. Explicit \code{~ 0 + x} and
\code{~ x - 1} are equivalent to \code{~ x} for environmental preparation.
Formula predictors must be finite and nonconstant. Rank is diagnosed on a
separate centered, unit-norm representation, with QR tolerance equal to
machine epsilon times the larger environmental dimension. An anchor is
subtracted before diagnostic centering to retain small variation around a
large absolute level. No diagnostic transformation changes the matrix
returned for fitting or association. Dependencies at floating-point
precision are rejected, and no predictor is silently dropped.
Arguments in transformations must be unnamed. Each term produces one
column, in formula order.

For objects fitted with a formula, association testing also checks that
the environmental coefficients are estimable conditional on the intercept
and latent factors. Checks of fitted models use their own QR tolerances,
which need not equal the preparatory rank tolerance. The global test also
checks absorption by the retained nuisance design at machine precision on
a centered, normalized diagnostic copy, before checking the fitted
residual regression. Unestimable predictors cause an
explicit error rather than substitution of a latent-factor statistic.

The formula is normalized and stored without its original R environment
in a small internal \code{env.spec} attribute. The five statistical S4 slots
are unchanged. No environmental matrix or numerical fingerprint is stored.
The fitted formula is reconstructed automatically by \code{lfmm2.test};
that function has no formula argument. With a formula, the same predictors,
values and individual order must be supplied when testing. Dimensions and
predictor structure are checked, but numerical identity and individual
identity cannot be verified. They remain the user's responsibility.
Old objects and objects fitted without a formula retain the historical
association path. No restriction on the identity of tested loci is added.
}

\value{
\code{lfmm2} returns an object of class \code{lfmm2Class} that contains $K$ estimated latent factors \code{@U} and their loadings  \code{@V}.

The following method can be applied to an object of class \code{lfmm2Class}:

\item{\code{\link{lfmm2.test}}}{
     P-values adjusted for the $K$ latent factors computed by lfmm2.
}
}
\author{
    Olivier Francois
}
\seealso{
\code{\link{lfmm.data}}
\code{\link{impute}}
\code{\link{lfmm2.test}}
\code{\link{pca}}
\code{\link{lfmm}}
\link{tutorial}
}

\references{
Caye K, Jumentier B, Lepeule J, Francois O. (2019). LFMM 2: fast and accurate inference of gene-environment associations in genome-wide studies. Molecular biology and evolution, 36(4), 852-860.

Gain C, Francois O. (2021). LEA 3: Factor models in population genetics and ecological genomics with R. Molecular Ecology Resources. doi: 10.1111/1755-0998.13366

}
\examples{
### Example of analysis using lfmm2 ###


# Simulation with 10 target loci, with effect sizes ranging between -10 an 10 
# n = 100 individuals and L = 1000 loci

X <- as.matrix(rnorm(100)) # causal environmental variable
B <- rep(0, 1000) 
target <- sample(1:1000, 10) # target loci
B[target] <- runif(10, -10, +10) # effect sizes

# Creating hidden factors and loadings

U <- t(tcrossprod(as.matrix(c(-1,0.5,1.5)), X))+ matrix(rnorm(300), ncol = 3)
V <- matrix(rnorm(3000), ncol = 3)

# Simulating a binarized matrix containing haploid genotypes 
# Simulation performed with the generative LFMM

Y <- tcrossprod(as.matrix(X), B) + tcrossprod(U, V) + matrix(rnorm(100000, sd = .5), nrow = 100)
Y <- matrix(as.numeric(Y > 0), ncol = 1000)

######################################
# Fitting an LFMM with K = 3 factors #
######################################

mod2 <- lfmm2(input = Y, env = X, K = 3)

# Computing P-values and plotting their minus log10 values 
# Target loci are highlighted

pv <- lfmm2.test(object = mod2, input = Y, env = X, linear = TRUE)
plot(-log10(pv$pvalues), col = "grey", cex = .4, pch = 19)
points(target, -log10(pv$pvalues[target]), col = "red")

#rm(list = ls())
}


\keyword{lfmm2}
\keyword{tutorial}
