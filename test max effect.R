dgp4107variation=function(n,typeg,typelift,typex2)
{
  # n =sample size
  # typeg = 1 or 2 (see below)
  # typelift = 1,2 or 3 (see below)
  # typex2 = 1 or 2 (see below) 
  
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
  if(typex2==2){x2=(x1<=5)*(x3<=5)*1 + (x1<=5)*(x3>5)*2+(x1>5)*(x4<=5)*3 + (x1>5)*(x4>5)*4 + g*lift +  rnorm(n,0,.5) }
  
  y=(x1<=5)*(x3<=5)*1 + (x1<=5)*(x3>5)*2+(x1>5)*(x4<=5)*3 + (x1>5)*(x4>5)*4 + g*lift + rnorm(n,0,1)
  
  dat <- data.frame(y , g , x5,x1,x2,x3,x4, b1 = lift,b0=0,b2=0,truemax=0,maxeffect=0)
  
  dat$p=5
  dat$scenario=4107
  dat$idobs=1:n
  dat
}

dat <-  do.call(paste("dgp4107variation",sep=""),list(n=1000,typeg=1,typelift=1,typex2=1))
datTest <- do.call(paste("dgp4107variation",sep=""),list(n=100,typeg=1,typelift=1,typex2=1))


v.obj12 <- irfsrc(indMAXQC(y, g) ~ x1+x2+x3+x4+x5, data = dat,
                  ntree = 100, block.size = 1,nodesize = 50, do.trace = FALSE, membership = TRUE, statistics = TRUE)

pr12 <- predict(v.obj12, newdata = datTest, membership = TRUE)

pr12$predicted[90,]

getMaxEffect(pr12)

getMaxEffect2(pr12)

getMaxEffect <- function(predictResult)
{
  apply(
    predictResult$predicted, 1, function(x) 
  {
    Beta0 = x[4]
    Beta1 = x[5]
    Beta2 = x[6]
    if (Beta2 < 0 && (-Beta1 / (2 * Beta2)) >= 0 && (-Beta1 / (2 * Beta2)) <= 1)
    {
      tau = Beta0 - ((Beta1*Beta1) / (4 * Beta2));
    }
    else
    {
      if (Beta0 > (Beta0 + Beta1 + Beta2))
      {
        tau = Beta0;
      }
      else
      {
        tau = (Beta0 + Beta1 + Beta2);
      }
    }
    tau
  }
  )
}

getMaxEffect2 <- function(predictResult)
{
  apply(
    predictResult$predicted, 1, function(x) 
    {
      Beta0 = x[4]
      Beta1 = x[5]
      Beta2 = x[6]
      if (Beta2 < 0 && (-Beta1 / (2 * Beta2)) >= 0 && (-Beta1 / (2 * Beta2)) <= 1)
      {
        tau = 0
      }
      else
      {
        if (Beta0 > (Beta0 + Beta1 + Beta2))
        {
          tau = 1
        }
        else
        {
          tau = 2
        }
      }
      tau
    }
  )
}

function(x) 
{
  Beta0 = x[4]
  Beta1 = x[5]
  Beta2 = x[6]
  if (Beta2 < 0 && (-Beta1 / (2 * Beta2)) >= 0 && (-Beta1 / (2 * Beta2)) <= 1)
  {
    tau = Beta0 - ((Beta1*Beta1) / (4 * Beta2));
  }
  else
  {
    if (Beta0 > (Beta0 + Beta1 + Beta2))
    {
      tau = Beta0;
    }
    else
    {
      tau = (Beta0 + Beta1 + Beta2);
    }
  }
  tau
}



##################


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


v.obj122 <- irfsrc(indMAXQC(y, g) ~ x1+x2+x3+x4+x5, data = dat2,
                  ntree = 100, block.size = 1,nodesize = 50, do.trace = FALSE, membership = TRUE, statistics = TRUE)

pr122 <- predict(v.obj122, newdata = datTest2, membership = TRUE)


v.mobq <- irfsrc(indMOBQC(y, g) ~ x1+x2+x3+x4+x5, data = dat2,
                   ntree = 100, block.size = 1,nodesize = 50, do.trace = FALSE, membership = TRUE, statistics = TRUE)

prmobq <- predict(v.mobq, newdata = datTest2, membership = TRUE)


pr122$predicted[21,]

1.115019 + 8.034246 -3.858281 #5.290984
8.034246 / (2 * 3.858281) #1.041169

datTest2$maxeffect
getMaxEffect(pr122)
get.max.effect(pr122)

plot(get.max.effect(pr122),datTest2$maxeffect)
plot(datTest2$maxeffect,get.max.effect(pr122))

plot(datTest2$maxeffect,get.max.effect(prmobq))

mean((datTest2$maxeffect - get.max.effect(prmobq))^2)
mean((datTest2$maxeffect - get.max.effect(pr122))^2)

getMaxEffect2(pr122)


subset(iextractBBOPid(v.obj122$membership,pr122$membership,v.obj122$inbag,1,1),TestId == 21)

subset(iextractBBOP(v.obj122$yvar,v.obj122$membership,pr122$membership,v.obj122$inbag,1,1),TestId == 21)

iextractBeta(v.obj122$yvar,v.obj122$membership,pr122$membership,v.obj122$inbag,1,1)[21,]
