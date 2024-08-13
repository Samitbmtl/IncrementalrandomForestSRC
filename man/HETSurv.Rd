\name{HETSurv}
\alias{HETSurv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
HETSurv
}
\description{
This function is for individual treatment effect estimation for survival data. 
}
\usage{
HETSurv(formula, 
    data, 
    ...,
    incremental.tailType = 1,
    incremental.minObsControl = 10,
    incremental.minObsTreatment = 10,
    incremental.TreatmentLowBound = 0.3,
    incremental.TreatmentUpperBound = 0.7)
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
  
  \item{incremental.tailType}{
tail type : 0 no tail or 1 Exponential tail.
}
  \item{incremental.minObsControl}{
minimum control observations in each node. Should be between 1 and nodesize
}
  \item{incremental.minObsTreatment}{
minimum Treatment observations in each node. Should be between 1 and nodesize
}
  \item{incremental.TreatmentLowBound}{
minimum proportion of treatment observations in each node. Should be between 0 and incremental.TreatmentUpperBound
}
  \item{incremental.TreatmentUpperBound}{
Maximum proportion of treatment observations in each node. Should be between TreatmentLowBound and 1
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
predict returns diference of estimation of Kaplan-Meier of treatment and control data in the terminal dataset.
}
\references{
https://academic.oup.com/bioinformatics/article/36/2/629/5542949
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
library(survival)
dat <- survival::colon

#datacleaning
dat <- subset(dat,etype==2 & rx != 'Lev+5FU')
dat <- na.omit(dat)
dat$G <- ifelse(dat$rx == "Obs",0,1)
dat$T <- dat$time 
dat$C <- dat$status

#build forest
obj <- HETSurv(f(T, G, C) ~ sex+age+obstruct+perfor+adhere+nodes+differ+extent+surg+node4, data = dat ,ntree = 100,samptype ="swr")

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