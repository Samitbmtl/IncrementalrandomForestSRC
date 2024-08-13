library(IncrementalrandomForestSRC)

dgpquad=function(n,typeg,typelift,typex2)  
{
  x5 = runif(n)*10
  x1 = runif(n)*10
  x3 = runif(n)*10
  x4 = runif(n)*10
  
  
  if(typeg==1){g = runif(n)}
  if(typeg==2){
    g=-x1/50 + .6 + rnorm(n,0,.3)
    g=.01*(g<=.01) + .99*(g>=.99) + g*(g>.01)*(g<.99) }
  
  if(typelift==1){lift=1*(x5<5) + 5*(x5>=5)}
  if(typelift==2){lift=x5/2}
  if(typelift==3){lift=((x5-5)^2)/5}
  
  
  if(typex2==1){x2=runif(n)*10}
  #  if(typex2==2){x2=.5*((x1<=5)*(x3<=5)*1 + (x1<=5)*(x3>5)*2+(x1>5)*(x4<=5)*3 + (x1>5)*(x4>5)*4) + 10*g + (-lift-5) +  rnorm(n,0,.5) }
  
  b0=.5*((x1<=5)*(x3<=5)*1 + (x1<=5)*(x3>5)*2+(x1>5)*(x4<=5)*3 + (x1>5)*(x4>5)*4)
  b1=10
  b2=(-lift-5)
  
  if(typex2==2){x2= b0 + b1*g + b2*g^2 +  rnorm(n,0,.5) }
  
  maxeffect=b0-b1^2/(4*b2)   # because we always have 0 < -b1/(2 b2) < 1 and b2 < 0. 
  bestdose=-b1/(2*b2)
  
  y=b0 + b1*g + b2*g^2 +rnorm(n,0,1)
  
  dat <- data.frame(y , g , x5,x1,x2,x3,x4,b0,b1,b2,bestdose,maxeffect)
}



dat2 <-  do.call(paste("dgpquad",sep=""),list(n=1000,typeg=1,typelift=2,typex2=2))
datTest2 <- do.call(paste("dgpquad",sep=""),list(n=100,typeg=1,typelift=2,typex2=2))


dat2$y[1:5]

v.obj1_CResponse <- irfsrc(indMAXQC(y, g) ~ x1+x2+x3+x4+x5, data = dat2,incremental.CenterResponse = TRUE,
                  ntree = 100, block.size = 1,nodesize = 50, do.trace = FALSE, membership = TRUE, statistics = TRUE)

v.obj1_CResponse$yvar


v.obj1_No_CResponse <- irfsrc(indMAXQC(y, g) ~ x1+x2+x3+x4+x5, data = dat2,incremental.CenterResponse = FALSE,
                  ntree = 100, block.size = 1,nodesize = 50, do.trace = FALSE, membership = TRUE, statistics = TRUE)

v.obj1_No_CResponse$yvar




v.obj1_CResponseTreatment <- irfsrc(indMAXQC(y, g) ~ x1+x2+x3+x4+x5, data = dat2,incremental.CenterResponse = TRUE,incremental.CenterTreatment = TRUE,
                  ntree = 100, block.size = 1,nodesize = 50, do.trace = FALSE, membership = TRUE, statistics = TRUE)
dat2$y
dat2$g
v.obj1_CResponse$yvar


v.obj1_No_CResponse <- irfsrc(indMAXQC(y, g) ~ x1+x2+x3+x4+x5, data = dat2,incremental.CenterResponse = FALSE,
                  ntree = 100, block.size = 1,nodesize = 50, do.trace = FALSE, membership = TRUE, statistics = TRUE)

v.obj1_No_CResponse$yvar













randomForestSRC::rfsrc(y ~ "x1 + x2", data = dat2)

vary <- "y"
xvar.names <- c("a","b","c")
varx <- paste0(paste(xvar.names,'+'))
xnames <- ""
for (i in 1 : length(xvar.names))
{
  xnames <- paste0(xnames, xvar.names[i], "+")
}
xnames <- substr(xnames,1,nchar(xnames)-1)  

str(xvar.names)
  
eval(parse(text = paste0("randomForestSRC::rfsrc(",vary,"~ " , varx , ", data = dat2)")))

v.obj122$yvar[1:5,]


incremental.CenterResponse