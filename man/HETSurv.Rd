\name{HETSurv}
\alias{HETSurv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
HETSurv Help
}
\description{
Random forest to estimate the CATE (conditional average treatment effect) for survival data with a binary treatment.
}
\usage{
HETSurv(formula, 
    data, 
    incremental.tailType = 1,
    incremental.minObsControl = 10,
    incremental.minObsTreatment = 10,
    incremental.TreatmentLowBound = 0.3,
    incremental.TreatmentUpperBound = 0.7,
    ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
Formula in the form f(response, treatment) ~ covariates 
}
  \item{data}{
Data frame used for building the trees.
}
  \item{incremental.tailType}{
Tail type : 0 no tail or 1 Exponential tail. 
}
  \item{incremental.minObsControl}{
Minimum number of control observations in each node. Should be between 1 and nodesize 
}
  \item{incremental.minObsTreatment}{
Minimum number of treatment observations in each node. Should be between 1 and nodesize 
}
  \item{incremental.TreatmentLowBound}{
Minimum proportion of treatment observations in each node. Should be between 0 and incremental.TreatmentUpperBound
}
  \item{incremental.TreatmentUpperBound}{
Maximum proportion of treatment observations in each node. Should be between incremental.TreatmentLowBound and 1 
}
  \item{...}{
Other arguments used in the rfsrc function from the randomForestSRC package version 2.9.2 
}
}
\details{

Refer to the article below for more details. This is the method called ITE SRF in the article. \cr
We have a binary treatment W (a value of 0 represents a control and a value of 1 represents a treatment) and a survival response Y. The CATE (or ITE for individual treatment effect) is defined by tau(x) = E[Y | X=x,W=1] -  E[Y | X=x,W=0].\cr
Trees are built using the RandomForestSRC package and leveraging its custom split rule feature. \cr
The estimated expected survival times are obtained by integrating the Kaplan-Meier (KM) estimates. Since a KM estimate can sometimes be undefined past a certain value, the argument incremental.tailType controls what to do in this case. When set to 1, an estimated exponential tail is used to complete the KM.\cr
The arguments incremental.minObsControl, incremental.minObsTreatment, incremental.TreatmentLowBound, incremental.TreatmentUpperBound control the number of control and treatment observations in each node to make sure there are enough of each group since expected survival times for the control and treatment groups are estimated separately in each node.  \cr
Predictions from a fitted forest return estimated values for tau(x). \cr 
When making prediction two arguments can be used : \cr

OOB = TRUE or FALSE\cr
This argument defines how the Bag of prediction (BOP) will be constructed. The BOP is the set of locally similar observations from the training data that we use to estimate tau(x). \cr
If OOB = TRUE, only the tree's out of bag observations are used to construct the BOP. \cr
If OOB = FALSE, the In-Bag observations are used to construct the BOP.\cr

obsUniqueInTree = TRUE or FALSE\cr
This argument is used to define how the function manages multiple copies of observations in a terminal node (which is possible is we use bootstrap samples to build the trees). \cr
If obsUniqueInTree = TRUE, the function will only keep one copy of a training observation in each tree to construct the BOP.\cr
If obsUniqueInTree = FALSE, the function will keep all copies of a training observation in each tree to construct the BOP.\cr
 

}
\value{
Returns a fitted random forest object. 
}
\references{
Tabib, S. and Larocque, D. (2020). Non-Parametric Individual Treatment Effect Estimation for Survival Data with Random Forests. Bioinformatics 36, 629–636.
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

# example dataset
library(survival)
dat <- survival::colon

# data cleaning
dat <- subset(dat,etype==2 & rx != 'Lev+5FU')
dat <- na.omit(dat)
dat$G <- ifelse(dat$rx == "Obs",0,1)
dat$T <- dat$time 
dat$C <- dat$status

# training data
datrain=dat[1:(nrow(dat)-50),]

# new data
datest=dat[(nrow(dat)-50+1):nrow(dat),]

# build forest
obj <- HETSurv(f(T, C, G) ~ sex+age+obstruct+perfor+adhere+nodes+differ+extent+surg+node4, data = datrain ,ntree = 100,samptype ="swr")

# prediction for the training data
pred <- predict(obj, data=dat,OOB = TRUE, obsUniqueInTree = TRUE)
pred$predicted

# prediction with new data
pred <- predict(obj, newdata=datest,OOB = FALSE, obsUniqueInTree = TRUE)
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