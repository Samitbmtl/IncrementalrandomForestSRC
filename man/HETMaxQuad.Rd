\name{HETMaxQuad}
\alias{HETMaxQuad}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
HETMaxQuad
}
\description{
This function is for incremental treatment effect for continuous treatment and outcome. Quadratic maximum
}
\usage{
HETMaxQuad(formula, 
    data, 
    ...,
    incremental.CenterResponse = FALSE,
    incremental.CenterTreatment = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
formule in the form of f(Outcome, Treatment) ~ Covariates
}
  \item{data}{
Data frame used building the trees.
}
  \item{...}{
All arguments used in the rfsrc function from the randomForestSRC package version 2.9.2
}
  
  \item{incremental.CenterResponse}{
When set to TRUE, it modifies the outcome variable by centing its values.
}
  \item{incremental.CenterTreatment}{
When set to TRUE, it modifies the treatment variable by centing its values.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
predict returns maxtrix of (gmax, tau)
}
\references{

}
\author{
Sami Tabib
}
\note{

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

#example dataset
data(nmes_data, package = "causaldrf")

#dataset cleaning
dat <- nmes_data
dat <- dat[dat$TOTALEXP != 0,]
dat$uid <- seq.int(nrow(dat))
dat$beltuse <- as.numeric(unclass(dat$beltuse))
dat$POVSTALB <- as.numeric(unclass(dat$POVSTALB))
dat$Y = log(dat$TOTALEXP+1)
dat$G = dat$packyears
dat$G = dat$G - min(dat$packyears)
dat$G = dat$G / max(dat$G)

#build forest
obj <- HETMaxQuad(f(Y, G) ~ AGESMOKE+LASTAGE+MALE+RACE3+beltuse+educate+marital+POVSTALB, data = dat ,ntree = 100,samptype ="swr")

#prediction for the same dataset
pred <- predict(obj, data=dat,OOB = TRUE, obsUniqueInTree = TRUE)
pred$predicted

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
