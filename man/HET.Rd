\name{HET}
\alias{HET}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
HET
}
\description{
This function is for incremental treatment effect for continuous treatment and outcome. increase the heterogeneity of the CATE estimation.
}
\usage{
HET(formula, 
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

When making prediction with this method, there is 2 arguments that could be used :

OOB = TRUE or FALSE\cr
This argument define how the Big Bag of prediction (BBOP) will be constitued. 
Let's remember that a Bag of prediction (BOP) is the formed by the training observationa of the terminal node where a predicted observation falls in the tree. For a forest, we have multiple trees, so multiple BOP, this is what we call Big Bag of prediction (BBOP).
for more information about the BBOP, please refer to https://academic.oup.com/bioinformatics/article/36/2/629/5542949\cr
If OOB = TRUE, the package will use only the tree's out of bag observations to constitue the BBOP.   \cr
If OOB = FALSE, the package will use the InBag observations to constitute the BBOP. 


obsUniqueInTree = TRUE or FALSE\cr 
This argument is used to define how the package manage repetitions in the BBOP.
For each tree, we have a BOP where the predicted observations falls in the tree.\cr
If obsUniqueInTree = TRUE, the package will only keep one repetition of a training observation in the dataset. \cr
If obsUniqueInTree = FALSE, the package will only keep all training observation in the dataset (no mater how many times they are in the BOP). 

}
\value{
predict returns estimation of beta1 vector.
}
\references{
article here
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
obj <- HET(f(Y, G) ~ AGESMOKE+LASTAGE+MALE+RACE3+beltuse+educate+marital+POVSTALB, data = dat ,ntree = 100,samptype ="swr")

#prediction for the same dataset
pred <- predict(obj, data=dat,OOB = TRUE, obsUniqueInTree = TRUE)
pred$predicted

#prediction with new data
newDataSize = 10
set.seed(123456789)
newData <- data.frame(AGESMOKE=sample(5:200, replace=TRUE, size=newDataSize),
                      LASTAGE=sample(18:100, replace=TRUE, size=newDataSize),
                      MALE=sample(c(0,1), replace=TRUE, size=newDataSize),
                      RACE3=factor(sample(c(1,2,3), replace=TRUE, size=newDataSize)),
                      beltuse=sample(c(1,2,3), replace=TRUE, size=newDataSize),
                      educate=factor(sample(1:4, replace=TRUE, size=newDataSize)),
                      marital=factor(sample(1:5, replace=TRUE, size=newDataSize)),
                      POVSTALB=sample(1:5, replace=TRUE, size=newDataSize))


pred <- predict(obj, newdata=newData,OOB = FALSE, obsUniqueInTree = TRUE)
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
