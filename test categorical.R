library(IncrementalrandomForestSRC)

getwd()
setwd("C:\\Article3a")

dgpCat=function(n)
{
  # n =sample size
  x1 = runif(n)*10
  x2 = runif(n)*10 
  x3 = runif(n)*10 
  x4 = runif(n)*10 
  x5 = runif(n)*10
  
  g = sample.int(5, n,replace = TRUE) - 1
  
  #y= runif(n)*10 + runif(n)*10*g
  y= runif(n)*10*(g+1)
  
  y= (x1 < 5) * (g==1) * 3 + (x1 >= 5) *(g==4) * 6 + rnorm(n, sd=0.1)
  
  # if (x1 < 5)
  # {
  #   y= (g==1) * 3 + rnorm(n)
  # }
  # else
  # {
  #   y= (g==4) * 6 + rnorm(n)
  # }
  
  dat <- data.frame(y , g , x5,x1,x2,x3,x4)
  
  dat$p=5
  dat$scenario=1
  dat$idobs=1:n
  dat
}

set.seed(0)
d_training <- dgpCat(1000)
  
  
v.obj1 <- irfsrc(indCat(y, g) ~ x1+x2+x3+x4+x5, data = d_training, 
                 ntree = 100, block.size = 1,nodesize = 50, do.trace = FALSE, membership = TRUE, statistics = TRUE)

v.obj1 <- irfsrc(indCat(y, g) ~ x1+x2+x3+x4+x5, data = d_training, nsplit = 0,mtry = 5 ,bootstrap = "none",
                 ntree = 1, block.size = 1,nodesize = 50, do.trace = TRUE, membership = TRUE, statistics = TRUE)



######### test #######

ds_left <- subset(d_training,x1<0.39249)
ds_right <- subset(d_training,x1>=0.39249)

mean(subset(ds_left,g==0)$y)
mean(subset(ds_left,g==1)$y)
mean(subset(ds_left,g==2)$y)
mean(subset(ds_left,g==3)$y)
mean(subset(ds_left,g==4)$y)

length(subset(ds_left,g==0)$y)
length(subset(ds_left,g==1)$y)
length(subset(ds_left,g==2)$y)
length(subset(ds_left,g==3)$y)
length(subset(ds_left,g==4)$y)

tau1 <- mean(subset(ds_left,g==1)$y) - mean(subset(ds_left,g==0)$y)
tau2 <- mean(subset(ds_left,g==2)$y) - mean(subset(ds_left,g==0)$y)
tau3 <- mean(subset(ds_left,g==3)$y) - mean(subset(ds_left,g==0)$y)
tau4 <- mean(subset(ds_left,g==4)$y) - mean(subset(ds_left,g==0)$y)

print (paste(tau1,tau2,tau3,tau4))

###

mean(subset(ds_right,g==0)$y)
mean(subset(ds_right,g==1)$y)
mean(subset(ds_right,g==2)$y)
mean(subset(ds_right,g==3)$y)
mean(subset(ds_right,g==4)$y)

length(subset(ds_right,g==0)$y)
length(subset(ds_right,g==1)$y)
length(subset(ds_right,g==2)$y)
length(subset(ds_right,g==3)$y)
length(subset(ds_right,g==4)$y)

tau1 <- mean(subset(ds_right,g==1)$y) - mean(subset(ds_right,g==0)$y)
tau2 <- mean(subset(ds_right,g==2)$y) - mean(subset(ds_right,g==0)$y)
tau3 <- mean(subset(ds_right,g==3)$y) - mean(subset(ds_right,g==0)$y)
tau4 <- mean(subset(ds_right,g==4)$y) - mean(subset(ds_right,g==0)$y)

print (paste(tau1,tau2,tau3,tau4))


sqrt(35*965) * abs(15.237132752547 - 19.7143248868903)
######################

pr1 <- predict(v.obj1, membership = TRUE)

pr1$predicted

# agg <- aggregate(pr1$predicted$y, list(pr1$predicted$TestId, pr1$predicted$g), FUN=mean) 
# names(agg) <- c("TestId","g","y")
# 
# agg0 <- subset(agg,g==0)
# agg0$y0 <- agg0$y
# agg0$y <- NULL
# agg0$g <- NULL
# 
# resAgg <- merge(agg,agg0,by=c("TestId")) 
# resAgg$tau <- resAgg$y - resAgg$y0
# resAgg <- subset(resAgg,g!=0)
# 
# resAgg$y <- NULL
# resAgg$y0 <- NULL
# 
# resAgg <- resAgg[order(resAgg$g),]
# 
# res <- reshape(resAgg, direction = "wide", idvar = c("TestId"), timevar = "g")
