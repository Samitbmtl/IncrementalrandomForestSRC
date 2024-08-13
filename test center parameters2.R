dgp5107variation=function(n,typeg,typelift,typex2,sigmay=1,sigmacollider=.5,sigmag=.3)
  
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
    g=-x1/50 + .6 + rnorm(n,0,sigmag)
    g=.01*(g<=.01) + .99*(g>=.99) + g*(g>.01)*(g<.99) }
  
  if(typelift==1){lift=1*(x5<5) + 5*(x5>=5)}
  if(typelift==2){lift=x5/2}
  if(typelift==3){lift=((x5-5)^2)/5}
  
  if(typex2==1){x2=runif(n)*10}
  if(typex2==2){x2=(x1<=5)*(x3<=5)*1 + (x1<=5)*(x3>5)*2+(x1>5)*(x4<=5)*3 + (x1>5)*(x4>5)*4 + g*lift +  rnorm(n,0,sigmacollider) }
  
  y=(x1<=5)*(x3<=5)*1 + (x1<=5)*(x3>5)*2+(x1>5)*(x4<=5)*3 + (x1>5)*(x4>5)*4 + g*lift + rnorm(n,0,sigmay)
  
  
  #### Conversion pour que les variables correspondent à la description de l'article
  xt1=x1
  xt3=x3
  xt2=x4
  xt4=x5
  xt5=x2
  
  x1=xt1
  x2=xt2
  x3=xt3
  x4=xt4
  x5=xt5
  
  #### Fin conversion
  
  # print("SNR")
  # print(var(lift)/sigmay^2)
  # print("cor(cbind(x1,g))")
  # print(cor(cbind(x1,g))[1,2])
  # print("cor(cbind(y,x5))")
  # print(cor(cbind(y,x5))[1,2])
  
  dat <- data.frame(y , g , x5,x1,x2,x3,x4, b1 = lift,b0=0,b2=0,truemax=0,maxeffect=0)
  
  dat$p=5
  dat$scenario=5107
  dat$idobs=1:n
  dat
  
}


dat <-  do.call(paste("dgp5107variation",sep=""),list(n=1000,typeg=2,typelift=2,typex2=1,sigmacollider=0,sigmag=.07))
datTest <-  do.call(paste("dgp5107variation",sep=""),list(n=500,typeg=2,typelift=2,typex2=1,sigmacollider=0,sigmag=.07))


v.LinearC_NoC <- irfsrc(indLinearC(y, g) ~ x5+x1+x2+x3+x4, data = dat, seed=-1, mtry=5, bootstrap = "by.root",samptype="swr",
                          ntree = 100, block.size = 1,nodesize = 30, do.trace = FALSE, membership = TRUE, statistics = TRUE)


v.LinearC_Cy <- irfsrc(indLinearC(y, g) ~ x5+x1+x2+x3+x4, data = dat, seed=-1, mtry=5, bootstrap = "by.root",samptype="swr",incremental.CenterResponse = TRUE,
                          ntree = 100, block.size = 1,nodesize = 30, do.trace = FALSE, membership = TRUE, statistics = TRUE)


v.LinearC_Cg <- irfsrc(indLinearC(y, g) ~ x5+x1+x2+x3+x4, data = dat, seed=-1, mtry=5, bootstrap = "by.root",samptype="swr",incremental.CenterTreatment = TRUE,
                          ntree = 100, block.size = 1,nodesize = 30, do.trace = FALSE, membership = TRUE, statistics = TRUE)


v.LinearC_Cyg <- irfsrc(indLinearC(y, g) ~ x5+x1+x2+x3+x4, data = dat, seed=-1, mtry=5, bootstrap = "by.root",samptype="swr",incremental.CenterResponse = TRUE, incremental.CenterTreatment = TRUE,
                          ntree = 100, block.size = 1,nodesize = 30, do.trace = FALSE, membership = TRUE, statistics = TRUE)


pr.LinearC_NoC <- predict(v.LinearC_NoC,newdata=datTest, membership = TRUE)
pr.LinearC_Cy <- predict(v.LinearC_Cy,newdata=datTest, membership = TRUE)
pr.LinearC_Cg <- predict(v.LinearC_Cg,newdata=datTest, membership = TRUE)
pr.LinearC_Cyg <- predict(v.LinearC_Cyg,newdata=datTest, membership = TRUE)

mean((pr.LinearC_NoC$predicted$b1lin - datTest$b1)^2)
mean((pr.LinearC_Cy$predicted$b1lin - datTest$b1)^2)
mean((pr.LinearC_Cg$predicted$b1lin - datTest$b1)^2)
mean((pr.LinearC_Cyg$predicted$b1lin - datTest$b1)^2)

irfsrc(HET(y, g) ~ x5+x1+x2+x3+x4, data = dat, seed=-1, mtry=5, bootstrap = "by.root",samptype="swr",incremental.CenterTreatment = TRUE,
                          ntree = 100, block.size = 1,nodesize = 30, do.trace = FALSE, membership = TRUE, statistics = TRUE)


irfsrc(CMB(y, g) ~ x5+x1+x2+x3+x4, data = dat, seed=-1, mtry=5, bootstrap = "by.root",samptype="swr",incremental.CenterTreatment = TRUE,
                          ntree = 100, block.size = 1,nodesize = 30, do.trace = FALSE, membership = TRUE, statistics = TRUE)

HET(f(y, g)~ x5+x1, membership = TRUE)

h <- HET(f(y, g) ~ x5+x1+x2+x3+x4, data = dat, seed=-1, mtry=5, bootstrap = "by.root",samptype="swr",incremental.CenterTreatment = TRUE,
                          ntree = 100, block.size = 1,nodesize = 30, do.trace = FALSE, membership = TRUE, statistics = TRUE)

ph <- predict(h,data=dat)

ph$predicted

h <- CMB(f(y, g) ~ x5+x1+x2+x3+x4, data = dat, seed=-1, mtry=5, bootstrap = "by.root",samptype="swr",incremental.CenterTreatment = TRUE,
                          ntree = 100, block.size = 1,nodesize = 30, do.trace = FALSE, membership = TRUE, statistics = TRUE)
print(h)

ph <- predict(h,data=dat)

dat$g

dat$g[1] <- -dat$g[1]

max(dat$g)
min(dat$g)

hm <- HETMaxQuad(f(y, g) ~ x5+x1+x2+x3+x4, data = dat, seed=-1, mtry=5, bootstrap = "by.root",samptype="swr",
                          ntree = 100, block.size = 1,nodesize = 30, do.trace = FALSE, membership = TRUE, statistics = TRUE)

dat3 <- dat[1:5,]
pm <- predict(hm,data=dat)

pm$predicted







formula2 <- as.formula("y ~ x5 + x1")
format(formula2)

all.names(formula, max.names = 1e7)[2]
all.vars(formula, max.names = 1e7)

inherits(formula, "formula")

CHAR    <- as.character(formula2)

format(formula)

formula <- as.formula("f(y, g) ~ x5+x1+x2+x3+x4")

gsub("f\\(", "HET(", "f(y, g) ~ x5 + x1 + x2 + x3 + x4")
gsub("f\\(", "HET(", formula)

# CMB(f(y, g)~ x5+x1 #centrag
# HET(f(y, g)~ x5+x1 #centrage
# HETsurv(f((y,c,g)~ x5+x1)   # param esurv
# HETmaxquad (f(y, g)~ x5+x1 #limiter g entre 0 et 1
    