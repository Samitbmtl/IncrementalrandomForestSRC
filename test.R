library(survival)
library(IncrementalrandomForestSRC)
data(veteran, package = "IncrementalrandomForestSRC")

v.obj <- rfsrc(time ~ karno+diagtime+age+prior, data = veteran, 
               ntree = 100, block.size = 1,nsplit = 0,mtry = 4, do.trace = TRUE, membership = TRUE, statistics = TRUE)

v.obj <- rfsrc(Surv(time, status) ~ karno+diagtime+age+prior, data = veteran, 
                ntree = 100, block.size = 1,nsplit = 0,mtry = 4, do.trace = TRUE, membership = TRUE, statistics = TRUE)

v.obj <- irfsrc(time ~ karno+diagtime+age+prior, data = veteran, 
               ntree = 100, block.size = 1,nsplit = 0,mtry = 4, do.trace = TRUE, membership = TRUE, statistics = TRUE)

v.obj <- irfsrc(Surv(time, status) ~ karno+diagtime+age+prior, data = veteran, 
               ntree = 100, block.size = 1,nsplit = 0,mtry = 4, do.trace = TRUE, membership = TRUE, statistics = TRUE)

t <- c(72,97)
g <- c(69,67)
d <- data.frame(t,g)
lm(t~g, data=d)

pr <- predict(v.obj, membership = TRUE)
pr <- predict(v.obj, membership = TRUE, newdata=veteran)
pr$predicted

typeof(v.obj$xvar)

v.obj <- irfsrc(indLSC(time, age) ~ karno+diagtime+prior, data = veteran,
                ntree = 100, block.size = 1,nsplit = 0,mtry = 3, do.trace = TRUE, membership = TRUE, statistics = TRUE)


v.obj <- irfsrc(indLSC(time, status) ~ karno+diagtime+age+prior, data = veteran, bootstrap = "none",
                ntree = 1, block.size = 1,nsplit = 0,mtry = 4, do.trace = TRUE, membership = TRUE, statistics = TRUE)

v.obj2 <- rfsrc(time ~ karno+diagtime+age+prior, data = veteran, bootstrap = "none",
                ntree = 1, block.size = 1,nsplit = 0,mtry = 4, do.trace = TRUE, membership = TRUE, statistics = TRUE)


pr2 <- predict(v.obj2, membership = TRUE)
pr2$predicted

stat.split.irfsrc(v.obj)[[1]]
stat.split(v.obj)[[1]]
stat.split(v.obj2)[[1]]

library(IncrementalrandomForestSRC)
library(parallel)
d <- veteran
d$trt <- ifelse(d$trt == 2,1,0)

v.obj2 <- irfsrc(indSurvB(time, status, trt) ~ karno+diagtime+age+prior, data = d, 
                 ntree = 3, block.size = 1,nsplit = 0,nodesize = 50,mtry = 4 ,bootstrap = "none", do.trace = TRUE, membership = TRUE, statistics = TRUE)

pr <- predict(v.obj2, membership = TRUE)

pr$predicted


d$trt <- ifelse(d$trt == 2,.5,.9)


v.obj2 <- irfsrc(indLinearC(time, trt) ~ karno+diagtime+age+prior, data = d, 
                 ntree = 3, block.size = 1,nsplit = 0,nodesize = 50,mtry = 4 ,bootstrap = "none", do.trace = TRUE, membership = TRUE, statistics = TRUE)

pr <- predict(v.obj2, membership = TRUE)

pr$predicted
pr$membership

# indSurvB(time, status, trt) : trt 0/1
# indSurvLinearC(time, status, trt) : trt continue
# indSurvQuadC(time, status, trt) : trt continue
# indSurvQuad2C(time, status, trt) : trt continue

# indLinearC(time, trt) : trt continue
# indQuadC(time, trt) : trt continue
# indQuad2C(time, trt) : trt continue


p_Path <- 'C:\\Article2b\\SimData\\Simulations\\1\\'

p_TrainingFile <- 'dataTrainingFplusplus.csv'
p_TestFile <- 'dataTestR.csv'
#p_ResultFile <- 'Reg_Results.txt'

getwd()

d_training <-read.csv(paste(p_Path,p_TrainingFile,sep=""), header = TRUE, sep = ";")
#d_results <- read.csv(paste(p_Path,p_ResultFile,sep=""), header = TRUE, sep = ";")
d_test <- read.csv(paste(p_Path,p_TestFile,sep=""), header = TRUE, sep = ";")
#d_results$X <- NULL

d_training$logt = log(d_training$T)
d_test$logt = log(d_test$T)

names(d_training)
names(d_test)

v.obj2 <- irfsrc(indLinearC(T, G) ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data = d_training, 
                 ntree = 100,nodesize = 50, do.trace = FALSE, membership = TRUE, statistics = TRUE)


v.obj2 <- irfsrc(indLinearC(T, G) ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data = d_training, 
                 ntree = 1, block.size = 1,nsplit = 0,nodesize = 50,mtry = 10 ,bootstrap = "none", do.trace = TRUE, membership = TRUE, statistics = TRUE)


v.obj2 <- irfsrc(indLinearC(logt, G) ~ x1, data = d_training, 
                 ntree = 1, block.size = 1,nsplit = 0,nodesize = 50,mtry = 10 ,bootstrap = "none", do.trace = TRUE, membership = TRUE, statistics = TRUE)

v.obj2 <- irfsrc(indQuadC(logt, G) ~ x1, data = d_training, 
                 ntree = 1, block.size = 1,nsplit = 0,nodesize = 50,mtry = 10 ,bootstrap = "none", do.trace = TRUE, membership = TRUE, statistics = TRUE)

v.obj2 <- irfsrc(indQuad2C(logt, G) ~ x1, data = d_training, 
                 ntree = 1, block.size = 1,nsplit = 0,nodesize = 50,mtry = 10 ,bootstrap = "none", do.trace = TRUE, membership = TRUE, statistics = TRUE)



pr <- predict(v.obj2, membership = TRUE)
pr <- predict(v.obj2, newdata = d_test, membership = TRUE)
pr <- predict(v.obj2, newdata= d_test[1,], membership = TRUE)

pr$predicted
v.obj2 <- NULL
pr$predicted[,2]

v.obj2 <- irfsrc(indLinearC(logt, G) ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data = d_training, 
                 ntree = 1, block.size = 1,nsplit = 0,nodesize = 50,mtry = 10 ,bootstrap = "none", do.trace = TRUE, membership = TRUE, statistics = TRUE)


v.obj2 <- irfsrc(indLinearC(logt, G) ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data = d_training, 
                 ntree = 1, block.size = 1,nsplit = 0,nodesize = 50,mtry = 10 ,bootstrap = "none", do.trace = FALSE, membership = TRUE, statistics = TRUE)

v.obj2 <- irfsrc(indLinearC(logt, G) ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data = d_training, 
                 ntree = 1, block.size = 1,nsplit = 0,nodesize = 50,mtry = 10,bootstrap = "by.root",samptype="swr", do.trace = FALSE, membership = TRUE, statistics = TRUE)


v.obj2 <- irfsrc(indQuadC(logt, G) ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data = d_training, 
                 ntree = 1, block.size = 1,nsplit = 0,nodesize = 50,mtry = 10 ,bootstrap = "none", do.trace = TRUE, membership = TRUE, statistics = TRUE)

v.obj2 <- irfsrc(indQuadC(logt, G) ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data = d_training, 
                 ntree = 1, block.size = 1,nsplit = 0,nodesize = 50,mtry = 10 ,bootstrap = "none", do.trace = FALSE, membership = TRUE, statistics = TRUE)

v.obj2 <- irfsrc(indQuadC(logt, G) ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data = d_training, 
                 ntree = 1, block.size = 1,nsplit = 0,nodesize = 50,mtry = 10,bootstrap = "by.root",samptype="swr", do.trace = FALSE, membership = TRUE, statistics = TRUE)



v.obj2 <- irfsrc(indQuad2C(logt, G) ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data = d_training, 
                 ntree = 1, block.size = 1,nsplit = 0,nodesize = 50,mtry = 10 ,bootstrap = "none", do.trace = TRUE, membership = TRUE, statistics = TRUE)

v.obj2 <- irfsrc(indQuad2C(logt, G) ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data = d_training, 
                 ntree = 1, block.size = 1,nsplit = 0,nodesize = 50,mtry = 10 ,bootstrap = "none", do.trace = FALSE, membership = TRUE, statistics = TRUE)

v.obj2 <- irfsrc(indQuad2C(logt, G) ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data = d_training, 
                 ntree = 1, block.size = 1,nsplit = 0,nodesize = 50,mtry = 10,bootstrap = "by.root",samptype="swr", do.trace = FALSE, membership = TRUE, statistics = TRUE)

d_test[106,]

pr$predicted[1]
pr$predicted[106]
pr$predicted[324]
pr$predicted[338]
pr$predicted[1000]

d_test[1,]

v.obj2$time.interest

pr$inbag

f_Withbag <- function(x,rowID, RFobj, PrObj,Objinbag) {
  t <- which(RFobj[,x] == PrObj[rowID,x])
  t[Objinbag[t,x] != 0]
}

#############################

p_PathSim <- 'C:\\Article2b\\SimData\\Simulations\\'

p_TrainingFile <- 'dataTrainingFplusplus.csv'
p_TestFile <- 'dataTestR.csv'
p_ResultatFile <- 'ResultatFinal_RegB_R.csv'

for (i in 1 : 10)
{
  d_training <-read.csv(paste(p_PathSim,i,"\\",p_TrainingFile,sep=""), header = TRUE, sep = ";")
  d_test <- read.csv(paste(p_PathSim,i,"\\",p_TestFile,sep=""), header = TRUE, sep = ";")

  d_training$logt = log(d_training$T)
  d_test$logt = log(d_test$T)  
  
  # v.obj2 <- irfsrc(indLinearC(logt, G) ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data = d_training, 
  #                  ntree = 100, block.size = 1,nsplit = 0,nodesize = 50,mtry = 10 ,bootstrap = "by.root",samptype="swr", do.trace = FALSE, membership = TRUE, statistics = TRUE)
  
  v.obj2 <- irfsrc(indLinearC(logt, G) ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data = d_training, 
                   ntree = 100, block.size = 1,nodesize = 50, do.trace = FALSE, membership = TRUE, statistics = TRUE)
  
  pr <- predict(v.obj2, newdata = d_test, membership = TRUE)
  
  #ds <- data.frame(Reallift = d_test$lift ,Beta= pr$predicted ,BetaMSE_RegB_R= (pr$predicted - d_test$lift)^2 )
  ds <- data.frame(Reallift = d_test$lift ,Beta= pr$predicted ,BetaMSE_RegB_R= (pr$predicted[,2] - d_test$lift)^2 )
  
  write.table(ds, file = paste(p_PathSim,i,"\\",p_ResultatFile,sep=""), sep = ";",quote = FALSE,row.names = FALSE, qmethod = "escape")
  
  print(i)
}



ResultatCompile <- matrix(ncol=3, nrow=10)
colnames(ResultatCompile) <-c("Simulation","Nbr","MSE_RegB_R")

for (i in 1 : 10)
{
  
  ResultatFinalPath <- paste(p_PathSim,'\\',i,'\\ResultatFinal_RegB_R.csv' ,sep="")
  if (file.exists(ResultatFinalPath)) 
  {
    tab <- read.csv(ResultatFinalPath, header = TRUE, sep = ";")
    ResultatCompile[i,] <-  cbind(i,nrow(tab),mean(tab$b1lin))
  }
}
write.table(ResultatCompile, file = paste(p_PathSim,"\\ResultatCompile_RegB_R.csv",sep=""), sep = ";",quote = FALSE,row.names = FALSE, qmethod = "escape")



#############################

p_PathSim <- 'C:\\Article2b\\SimData\\Simulations\\'

p_TrainingFile <- 'dataTrainingFplusplus.csv'
p_TestFile <- 'dataTestR.csv'
p_ResultatFile <- 'ResultatFinal_RegQuad_R.csv'

for (i in 1 : 10)
{
  d_training <-read.csv(paste(p_PathSim,i,"\\",p_TrainingFile,sep=""), header = TRUE, sep = ";")
  d_test <- read.csv(paste(p_PathSim,i,"\\",p_TestFile,sep=""), header = TRUE, sep = ";")
  
  d_training$logt = log(d_training$T)
  d_test$logt = log(d_test$T)  
  
  # v.obj2 <- irfsrc(indQuadC(logt, G) ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data = d_training,
  #                  ntree = 100, block.size = 1,nsplit = 0,nodesize = 50,mtry = 10 ,bootstrap = "by.root",samptype="swr", do.trace = FALSE, membership = TRUE, statistics = TRUE)
  
  v.obj2 <- irfsrc(indQuadC(logt, G) ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data = d_training,
                   ntree = 100, block.size = 1,nodesize = 50, do.trace = FALSE, membership = TRUE, statistics = TRUE)
  
  pr <- predict(v.obj2, newdata = d_test, membership = TRUE)
  
  ds <- data.frame(Reallift = d_test$lift ,Beta= pr$predicted ,BetaMSE_RegQuad_R= (pr$predicted[,2] - d_test$lift)^2 )
  
  write.table(ds, file = paste(p_PathSim,i,"\\",p_ResultatFile,sep=""), sep = ";",quote = FALSE,row.names = FALSE, qmethod = "escape")
  
  print(i)
}



ResultatCompile <- matrix(ncol=3, nrow=10)
colnames(ResultatCompile) <-c("Simulation","Nbr","MSE_RegQuad_R")

for (i in 1 : 10)
{
  
  ResultatFinalPath <- paste(p_PathSim,'\\',i,'\\ResultatFinal_RegQuad_R.csv' ,sep="")
  if (file.exists(ResultatFinalPath)) 
  {
    tab <- read.csv(ResultatFinalPath, header = TRUE, sep = ";")
    ResultatCompile[i,] <-  cbind(i,nrow(tab),mean(tab$BetaMSE_RegQuad_R))
  }
}
write.table(ResultatCompile, file = paste(p_PathSim,"\\ResultatCompile_RegQuad_R.csv",sep=""), sep = ";",quote = FALSE,row.names = FALSE, qmethod = "escape")



#############################

p_PathSim <- 'C:\\Article2b\\SimData\\Simulations\\'

p_TrainingFile <- 'dataTrainingFplusplus.csv'
p_TestFile <- 'dataTestR.csv'
p_ResultatFile <- 'ResultatFinal_RegQuad2_R.csv'

for (i in 1 : 10)
{
  d_training <-read.csv(paste(p_PathSim,i,"\\",p_TrainingFile,sep=""), header = TRUE, sep = ";")
  d_test <- read.csv(paste(p_PathSim,i,"\\",p_TestFile,sep=""), header = TRUE, sep = ";")
  
  d_training$logt = log(d_training$T)
  d_test$logt = log(d_test$T)  
  
  # v.obj2 <- irfsrc(indQuad2C(logt, G) ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data = d_training,
  #                  ntree = 100, block.size = 1,nsplit = 0,nodesize = 50,mtry = 10 ,bootstrap = "by.root",samptype="swr", do.trace = FALSE, membership = TRUE, statistics = TRUE)
  
  v.obj2 <- irfsrc(indQuad2C(logt, G) ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data = d_training,
                   ntree = 100, block.size = 1,nodesize = 50, do.trace = FALSE, membership = TRUE, statistics = TRUE)
  
  pr <- predict(v.obj2, newdata = d_test, membership = TRUE)
  
  #ds <- data.frame(Reallift = d_test$lift ,Beta= pr$predicted ,BetaMSE_RegQuad2_R= (pr$predicted - d_test$lift)^2 )
  ds <- data.frame(Reallift = d_test$lift ,Beta= c(pr$predicted) ,BetaMSE_RegQuad2_R= c((pr$predicted[,2] - d_test$lift)^2) )
  names(ds) = c("Reallift","Beta","BetaMSE_RegQuad2_R")
  
  write.table(ds, file = paste(p_PathSim,i,"\\",p_ResultatFile,sep=""), sep = ";",quote = FALSE,row.names = FALSE, qmethod = "escape")
  
  print(i)
}


ResultatCompile <- matrix(ncol=3, nrow=10)
colnames(ResultatCompile) <-c("Simulation","Nbr","MSE_RegQuad2_R")

for (i in 1 : 10)
{
  
  ResultatFinalPath <- paste(p_PathSim,'\\',i,'\\ResultatFinal_RegQuad2_R.csv' ,sep="")
  if (file.exists(ResultatFinalPath)) 
  {
    tab <- read.csv(ResultatFinalPath, header = TRUE, sep = ";")
    ResultatCompile[i,] <-  cbind(i,nrow(tab),mean(tab$BetaMSE_RegQuad2_R))
  }
}
write.table(ResultatCompile, file = paste(p_PathSim,"\\ResultatCompile_RegQuad2_R.csv",sep=""), sep = ";",quote = FALSE,row.names = FALSE, qmethod = "escape")


#############################

p_PathSim <- 'C:\\Article2b\\SimData\\Simulations\\'

p_TrainingFile <- 'dataTrainingFplusplus.csv'
p_TestFile <- 'dataTestR.csv'
p_ResultatFile <- 'ResultatFinal_RegTST_R.csv'

for (i in 1 : 10)
{
  d_training <-read.csv(paste(p_PathSim,i,"\\",p_TrainingFile,sep=""), header = TRUE, sep = ";")
  d_test <- read.csv(paste(p_PathSim,i,"\\",p_TestFile,sep=""), header = TRUE, sep = ";")
  
  d_training$logt = log(d_training$T)
  d_test$logt = log(d_test$T)  
  
  # v.obj2 <- irfsrc(indTSTC(logt, G) ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data = d_training,
  #                  ntree = 100, block.size = 1,nsplit = 0,nodesize = 50,mtry = 10 ,bootstrap = "by.root",samptype="swr", do.trace = FALSE, membership = TRUE, statistics = TRUE)
  
  v.obj2 <- irfsrc(indTSTC(logt, G) ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data = d_training,
                   ntree = 100, block.size = 1,nodesize = 50, do.trace = FALSE, membership = TRUE, statistics = TRUE)
  
  pr <- predict(v.obj2, newdata = d_test, membership = TRUE)
  
  ds <- data.frame(Reallift = d_test$lift ,Beta= c(pr$predicted[,2]) ,BetaMSE_RegTST_R= c((pr$predicted[,2] - d_test$lift)^2) )
  names(ds) = c("Reallift","Beta","BetaMSE_RegTST_R")
  
  write.table(ds, file = paste(p_PathSim,i,"\\",p_ResultatFile,sep=""), sep = ";",quote = FALSE,row.names = FALSE, qmethod = "escape")
  
  print(i)
}


ResultatCompile <- matrix(ncol=3, nrow=10)
colnames(ResultatCompile) <-c("Simulation","Nbr","MSE_RegTST_R")

for (i in 1 : 10)
{
  
  ResultatFinalPath <- paste(p_PathSim,'\\',i,'\\ResultatFinal_RegTST_R.csv' ,sep="")
  if (file.exists(ResultatFinalPath)) 
  {
    tab <- read.csv(ResultatFinalPath, header = TRUE, sep = ";")
    ResultatCompile[i,] <-  cbind(i,nrow(tab),mean(tab$BetaMSE_RegTST_R))
  }
}
write.table(ResultatCompile, file = paste(p_PathSim,"\\ResultatCompile_RegTST_R.csv",sep=""), sep = ";",quote = FALSE,row.names = FALSE, qmethod = "escape")


#############################


p_PathSim <- 'C:\\Article2b\\SimData\\Simulations\\'

p_TrainingFile <- 'dataTrainingFplusplus.csv'
p_TestFile <- 'dataTestR.csv'
p_ResultatFile <- 'ResultatFinal_RegNSR_R.csv'
i = 1
d_training <-read.csv(paste(p_PathSim,i,"\\",p_TrainingFile,sep=""), header = TRUE, sep = ";")
d_test <- read.csv(paste(p_PathSim,i,"\\",p_TestFile,sep=""), header = TRUE, sep = ";")

d_training$logt = log(d_training$T)
d_test$logt = log(d_test$T)  

v.obj2 <- irfsrc(indNSRC(logt, G) ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data = d_training,nsplit = 0,mtry = 10,bootstrap = "none",
                 ntree = 1, block.size = 1,nodesize = 400, do.trace = TRUE, membership = TRUE, statistics = TRUE)

v.obj2 <- irfsrc(indTSTC(logt, G) ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data = d_training,nsplit = 0,mtry = 10,bootstrap = "none",
                 ntree = 1, block.size = 1,nodesize = 400, do.trace = TRUE, membership = TRUE, statistics = TRUE)

v.obj2 <- irfsrc(indTST2C(logt, G) ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data = d_training,nsplit = 0,mtry = 10,bootstrap = "none",
                 ntree = 1, block.size = 1,nodesize = 400, do.trace = TRUE, membership = TRUE, statistics = TRUE)


#############################


p_PathSim <- 'C:\\Article2b\\SimData\\Simulations\\'

p_TrainingFile <- 'dataTrainingFplusplus.csv'
p_TestFile <- 'dataTestR.csv'
p_ResultatFile <- 'ResultatFinal_RegNSR_R.csv'

for (i in 1 : 10)
{
  d_training <-read.csv(paste(p_PathSim,i,"\\",p_TrainingFile,sep=""), header = TRUE, sep = ";")
  d_test <- read.csv(paste(p_PathSim,i,"\\",p_TestFile,sep=""), header = TRUE, sep = ";")
  
  d_training$logt = log(d_training$T)
  d_test$logt = log(d_test$T)  
  
  # v.obj2 <- irfsrc(indNSRC(logt, G) ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data = d_training,
  #                  ntree = 100, block.size = 1,nsplit = 0,nodesize = 50,mtry = 10 ,bootstrap = "by.root",samptype="swr", do.trace = FALSE, membership = TRUE, statistics = TRUE)
  
  v.obj2 <- irfsrc(indNSRC(logt, G) ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data = d_training,
                   ntree = 100, block.size = 1,nodesize = 50, do.trace = FALSE, membership = TRUE, statistics = TRUE)
  
  pr <- predict(v.obj2, newdata = d_test, membership = TRUE)
  
  ds <- data.frame(Reallift = d_test$lift ,Beta= c(pr$predicted[,2]) ,BetaMSE_RegNSR_R= c((pr$predicted[,2] - d_test$lift)^2) )
  names(ds) = c("Reallift","Beta","BetaMSE_RegNSR_R")
  
  write.table(ds, file = paste(p_PathSim,i,"\\",p_ResultatFile,sep=""), sep = ";",quote = FALSE,row.names = FALSE, qmethod = "escape")
  
  print(i)
}


ResultatCompile <- matrix(ncol=3, nrow=10)
colnames(ResultatCompile) <-c("Simulation","Nbr","MSE_RegNSR_R")

for (i in 1 : 10)
{
  
  ResultatFinalPath <- paste(p_PathSim,'\\',i,'\\ResultatFinal_RegNSR_R.csv' ,sep="")
  if (file.exists(ResultatFinalPath)) 
  {
    tab <- read.csv(ResultatFinalPath, header = TRUE, sep = ";")
    ResultatCompile[i,] <-  cbind(i,nrow(tab),mean(tab$BetaMSE_RegNSR_R))
  }
}
write.table(ResultatCompile, file = paste(p_PathSim,"\\ResultatCompile_RegNSR_R.csv",sep=""), sep = ";",quote = FALSE,row.names = FALSE, qmethod = "escape")

################ Find BOP

v.obj2

pr$predicted

###### Version Simple estim: 
f_Withbag <- function(x,rowID, RFobj, PrObj,Objinbag) {
  t <- which(RFobj[,x] == PrObj[rowID,x])
  t[Objinbag[t,x] != 0]
}

fun_Test <- function(rowID,dat, RFobj, PrObj,Objinbag) {
  # print(rowID)  
  BOPdata <- dat[c(Reduce(c,parallel::mclapply(1:length(RFobj[1,]), f_Withbag,rowID,RFobj,PrObj,Objinbag))),]
  colnames(BOPdata) <- c("T","G")
  mod <- lm(T ~ G, data=BOPdata)
  return (coef(mod)["G"])
}

sapply(1:nrow(pr$membership), fun_Test,v.obj2$yvar,v.obj2$membership,pr$membership,v.obj2$inbag)
######

###### Version complique estim: 
sapply(1:nrow(pr$membership), function(rowID,dat, RFobj, PrObj,Objinbag) {
  BOPdata <- dat[c(Reduce(c,parallel::mclapply(1:length(RFobj[1,]), function(x,rowID, RFobj, PrObj,Objinbag) {
    t <- which(RFobj[,x] == PrObj[rowID,x])
    t[Objinbag[t,x] != 0]
  },rowID,RFobj,PrObj,Objinbag))),]
  colnames(BOPdata) <- c("T","G")
  mod <- lm(T ~ G, data=BOPdata)
  return (coef(mod)["G"])
},v.obj2$yvar,v.obj2$membership,pr$membership,v.obj2$inbag)
######

###### Version complique BOP  : 
sapply(1:nrow(pr$membership), function(rowID,dat, RFobj, PrObj,Objinbag) {
  return (paste(c(Reduce(c,parallel::mclapply(1:length(RFobj[1,]), function(x,rowID, RFobj, PrObj,Objinbag) {
    t <- which(RFobj[,x] == PrObj[rowID,x])
    t[Objinbag[t,x] != 0]
  },rowID,RFobj,PrObj,Objinbag)))
  ,collapse= ";"))
},v.obj2$yvar,v.obj2$membership,pr$membership,v.obj2$inbag)
######

ObjMembership <- v.obj2$membership
Objinbag <- v.obj2$inbag
Objyvar <- v.obj2$yvar

no_cores <-parallel::detectCores()
cl_compile<- parallel::makeCluster(no_cores)
predicted <- parallel::parSapply(cl_compile,1:nrow(irfsrcOutput$membership), fun_pred_indQuadC,Objyvar,ObjMembership,irfsrcOutput$membership,Objinbag)
parallel::stopCluster(cl_compile)

fun_pred_indQuadC <- function(rowID,dat, RFobj, PrObj,Objinbag) {
  BOPdata <- dat[c(Reduce(c,parallel::mclapply(1:length(RFobj[1,]), f_Withbag,rowID,RFobj,PrObj,Objinbag))),]
  colnames(BOPdata) <- c("T","G")
  mod <- lm(T ~ G, data=BOPdata)
  return (coef(mod)["G"])
  # return (paste(c(Reduce(c,parallel::mclapply(1:length(RFobj[1,]), f_Withbag,rowID,RFobj,PrObj,Objinbag))),collapse= ";"))
}

###### Version complique betas: 
r <- data.frame(t(sapply(1:nrow(pr$membership), function(rowID,dat, RFobj, PrObj,Objinbag) {
  BOPdata <- dat[c(Reduce(c,parallel::mclapply(1:length(RFobj[1,]), function(x,rowID, RFobj, PrObj,Objinbag) {
    t <- which(RFobj[,x] == PrObj[rowID,x])
    t[Objinbag[t,x] != 0]
  },rowID,RFobj,PrObj,Objinbag))),]
  colnames(BOPdata) <- c("T","G")
  out=c(as.numeric(c(coef(lm(T ~ G, data=BOPdata)),coef(lm(T ~ G+I(G^2) , data=BOPdata)))))
  return (out)
},v.obj2$yvar,v.obj2$membership,pr$membership,v.obj2$inbag)))
names(r) <- c("b0lin","b1lin","b0quad","b1quad","b2quad")
r
######

r[2,]

mod <- lm(T ~ G, data=d_test)

c(as.numeric(c(coef(lm(T ~ G, data=d_test)))))

#############################


library(IncrementalrandomForestSRC)
library(randomForestSRC)

###############################

d <- veteran
d$trt <- ifelse(d$trt == 2,.5,.9)

v.obj <- irfsrc(indLSC(time, age) ~ karno+diagtime+prior, data = veteran, bootstrap = "none",
                ntree = 1, block.size = 1,nsplit = 0,mtry = 3, do.trace = TRUE, membership = TRUE, statistics = TRUE)

pr <- predict(v.obj, membership = TRUE)
pr <- predict(v.obj, membership = TRUE, newdata=veteran)
pr$predicted

v.obj2 <- irfsrc(indTSTC(time, status) ~ karno+diagtime+age+prior, data = d, bootstrap = "none",
                ntree = 1, block.size = 1,nsplit = 0,mtry = 4, do.trace = TRUE, membership = TRUE, statistics = TRUE)

v.obj1 <- irfsrc(indNSRC(time, status) ~ karno+diagtime+age+prior, data = d, bootstrap = "none",
                 ntree = 1, block.size = 1,nsplit = 0,mtry = 4, do.trace = TRUE, membership = TRUE, statistics = TRUE)

v.obj3 <- irfsrc(indLinearC(time, status) ~ karno+diagtime+age+prior, data = d, bootstrap = "none",
                 ntree = 1, block.size = 1,nsplit = 0,mtry = 4, do.trace = TRUE, membership = TRUE, statistics = TRUE)

v.obj4 <- irfsrc(indLSC(time, status) ~ karno+diagtime+age+prior, data = d, bootstrap = "none",
                 ntree = 1, block.size = 1,nsplit = 0,mtry = 4, do.trace = TRUE, membership = TRUE, statistics = TRUE)



stat.split.irfsrc(v.obj1)[[1]]
stat.split.irfsrc(v.obj2)[[1]]
stat.split.irfsrc(v.obj4)[[1]]
stat.split.irfsrc(v.obj3)[[1]]

pr2 <- predict(v.obj2, membership = TRUE)
pr4 <- predict(v.obj4, membership = TRUE)

pr2$predicted
pr4$predicted

################################## bug plantage ###############


library(IncrementalrandomForestSRC)


v.obj2 <- irfsrc(indTSTC(time, status) ~ karno+diagtime+age+prior, data = d, bootstrap = "none",
                 ntree = 1, block.size = 1,nsplit = 0,mtry = 4, do.trace = TRUE, membership = TRUE, statistics = TRUE)

v.obj1 <- irfsrc(indNSRC(time, status) ~ karno+diagtime+age+prior, data = d, bootstrap = "none",
                 ntree = 1, block.size = 1,nsplit = 0,mtry = 4, do.trace = TRUE, membership = TRUE, statistics = TRUE)

v.obj3 <- irfsrc(indLinearC(time, status) ~ karno+diagtime+age+prior, data = d, bootstrap = "none",
                 ntree = 1, block.size = 1,nsplit = 0,mtry = 4, do.trace = TRUE, membership = TRUE, statistics = TRUE)

v.obj4 <- irfsrc(indLSC(time, status) ~ karno+diagtime+age+prior, data = d, bootstrap = "none",
                 ntree = 1, block.size = 1,nsplit = 0,mtry = 4, do.trace = TRUE, membership = TRUE, statistics = TRUE)



pr2 <- predict(v.obj2, membership = TRUE)
pr4 <- predict(v.obj4, membership = TRUE)

a <- showConnections()

pr2 <- predict(v.obj2, membership = TRUE)
pr4 <- predict(v.obj4, membership = TRUE)

###############################


library(IncrementalrandomForestSRC)


v.obj4 <- irfsrc(indLSC(time, status) ~ karno+diagtime+age+prior, data = d, bootstrap = "none",
                 ntree = 1, block.size = 1,nsplit = 0,mtry = 4, do.trace = TRUE, membership = TRUE, statistics = TRUE)



v.obj2 <- irfsrc(indTSTC(time, status) ~ karno+diagtime+age+prior, data = d, bootstrap = "none",
                 ntree = 1, block.size = 1,nsplit = 0,mtry = 4, do.trace = TRUE, membership = TRUE, statistics = TRUE)


pr2 <- predict(v.obj2, membership = TRUE)

pr2$predicted[1:10,]

pr2 <- predict(v.obj2, membership = TRUE, newdata=d )


pr4 <- predict(v.obj4, membership = TRUE)
pr2 <- predict(v.obj2, membership = TRUE)

######################################

v.obj2 <- irfsrc(indLSC(time, status) ~ karno+diagtime+age+prior, data = d, bootstrap = "none",nodesize = 50,
                 ntree = 1, block.size = 1,nsplit = 0,mtry = 4, do.trace = TRUE, membership = TRUE, statistics = TRUE)


v.obj2 <- irfsrc(indTSTC(time, status) ~ karno+diagtime+age+prior, data = d, bootstrap = "none",nodesize = 50,
                 ntree = 1, block.size = 1,nsplit = 0,mtry = 4, do.trace = TRUE, membership = TRUE, statistics = TRUE)

stat.split.irfsrc(v.obj2)[[1]]

v.obj2 <- irfsrc(indTST2C(time, status) ~ karno+diagtime+age+prior, data = d, bootstrap = "none",nodesize = 50,
                 ntree = 1, block.size = 1,nsplit = 0,mtry = 4, do.trace = TRUE, membership = TRUE, statistics = TRUE)

v.obj2 <- irfsrc(indNSRC(time, status) ~ karno+diagtime+age+prior, data = d, bootstrap = "none",nodesize = 50,
                 ntree = 1, block.size = 1,nsplit = 0,mtry = 4, do.trace = TRUE, membership = TRUE, statistics = TRUE)

v.obj2 <- irfsrc(indMOB0C(time, status) ~ karno+diagtime+age+prior, data = d, bootstrap = "none",
                 ntree = 1, block.size = 1,nsplit = 0,mtry = 4, do.trace = TRUE, membership = TRUE, statistics = TRUE)

pr2 <- predict(v.obj2, membership = TRUE)
pr2$predicted[1:10,]

v.obj2 <- irfsrc(indMOBLC(time, status) ~ karno+diagtime+age+prior, data = d, bootstrap = "none",
                 ntree = 1, block.size = 1,nsplit = 0,mtry = 4, do.trace = TRUE, membership = TRUE, statistics = TRUE)

pr2 <- predict(v.obj2, membership = TRUE)
pr2$predicted[1:10,]

v.obj2 <- irfsrc(indMOBQC(time, status) ~ karno+diagtime+age+prior, data = d, bootstrap = "none",
                 ntree = 1, block.size = 1,nsplit = 0,mtry = 4, do.trace = TRUE, membership = TRUE, statistics = TRUE)

v.obj2 <- irfsrc(indMOB0C(logt, G) ~ x1+x2, data = d_training, bootstrap = "none",
                 ntree = 1, block.size = 1,nsplit = 0,mtry = 4, do.trace = TRUE, membership = TRUE, statistics = TRUE)

pr2 <- predict(v.obj2, membership = TRUE)
pr2$predicted[1:10,]

v.obj2 <- irfsrc(indMAXLC(time, status) ~ karno+diagtime+age+prior, data = d, bootstrap = "none",
                 ntree = 1, block.size = 1,nsplit = 0,mtry = 4, do.trace = TRUE, membership = TRUE, statistics = TRUE)

pr2 <- predict(v.obj2, membership = TRUE)
pr2$predicted[1:10,]

v.obj2 <- irfsrc(indMAXQC(time, status) ~ karno+diagtime+age+prior, data = d, bootstrap = "none",
                 ntree = 1, block.size = 1,nsplit = 0,mtry = 4, do.trace = TRUE, membership = TRUE, statistics = TRUE)

v.obj2 <- irfsrc(indMAXQC(logt, G) ~ x1+x2, data = d_training, bootstrap = "none",
                 ntree = 1, block.size = 1,nsplit = 0,mtry = 4, do.trace = TRUE, membership = TRUE, statistics = TRUE)


pr2 <- predict(v.obj2, membership = TRUE)
pr2$predicted[1:10,]

v.obj2 <- irfsrc(indMAXMAXQC(time, status) ~ karno+diagtime+age+prior, data = d, bootstrap = "none",
                 ntree = 1, block.size = 1,nsplit = 0,mtry = 4, do.trace = TRUE, membership = TRUE, statistics = TRUE)

v.obj2 <- irfsrc(indMAXMAXQC(logt, G) ~ x1+x2, data = d_training, bootstrap = "none",
                 ntree = 1, block.size = 1,nsplit = 0,mtry = 4, do.trace = TRUE, membership = TRUE, statistics = TRUE)

pr2 <- predict(v.obj2, membership = TRUE)
pr2$predicted[1:10,]

##############################

library(partykit)

reglin <- function(y, x, start = NULL, weights = NULL, offset = NULL, ...) 
{ 
  glm(y ~ 0 + x, family = gaussian, start = start, ...)
}


m <- mob(logt ~ G | x1, data=d_training,fit=reglin)


v.obj2 <- irfsrc(indMOBLC(logt, G) ~ x1, data = d_training, bootstrap = "none",
                 ntree = 1, block.size = 1,nsplit = 0,mtry = 1, do.trace = FALSE, membership = TRUE, statistics = TRUE)


stat.split.irfsrc(v.obj2)[[1]]


####################################

library(IncrementalrandomForestSRC)

p_PathSim <- 'C:\\Article2b\\SimData\\Simulations\\'

p_TrainingFile <- 'dataTrainingFplusplus.csv'
p_TestFile <- 'dataTestR.csv'
p_ResultatFile <- 'ResultatFinal_RegTST_R.csv'

i =1
d_training <-read.csv(paste(p_PathSim,i,"\\",p_TrainingFile,sep=""), header = TRUE, sep = ";")
d_test <- read.csv(paste(p_PathSim,i,"\\",p_TestFile,sep=""), header = TRUE, sep = ";")


v.obj2 <- irfsrc(indLinearC(T, G) ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data = d_training, seed= 1,
                 ntree = 100,nodesize = 400, do.trace = FALSE, membership = TRUE, statistics = TRUE)

system.time(
pr <- predict(v.obj2, newdata = d_test, membership = TRUE)
)
#13.84

mean(pr$predicted$b1lin)#43.69277

#copy all debug variable to main env.
#lapply(ls(), function(o) assign(x = o, value = get(o), envir = .GlobalEnv))

#to save env.
#save.image("samiTest3.rdata")



###



p_PathSim <- 'C:\\Article2b\\SimData\\Simulations\\'

p_TrainingFile <- 'dataTrainingFplusplus.csv'
p_TestFile <- 'dataTestR.csv'
p_ResultatFile <- 'ResultatFinal_RegTST_R.csv'

i =1
d_training <-read.csv(paste(p_PathSim,i,"\\",p_TrainingFile,sep=""), header = TRUE, sep = ";")
d_test <- read.csv(paste(p_PathSim,i,"\\",p_TestFile,sep=""), header = TRUE, sep = ";")


v.obj2 <- irfsrc(indLinearC(T, G) ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data = d_training, seed= 1,
                 ntree = 100,nodesize = 5, do.trace = FALSE, membership = TRUE, statistics = TRUE)

pr <- predict(v.obj2, newdata = d_test, membership = TRUE)

pr_old <- predict.irfsrc_old(v.obj2, newdata = d_test, membership = TRUE)

sum(pr_old$predicted) - sum(pr$predicted)
max(pr_old$predicted - pr$predicted)

mean(pr$predicted$b1lin)#43.69277
mean(out$b1lin)#43.69277

v.obj2$membership[1:10,]
v.obj2$inbag[1:10,]
v.obj2$yvar[1:10,]
pr$membership[1:1,]

out <- iextractBeta(v.obj2$yvar,v.obj2$membership,pr$membership,v.obj2$inbag,0,1)

pr$predicted
out

out2 <- NULL
gc() 

out2 <- iextractBBOP(v.obj2$yvar,v.obj2$membership,pr$membership,v.obj2$inbag,0,1)
su <- subset(out2,TreeId==2)
nrow(su)
lm(T~G, data=su)

f_Withbag <- function(x,rowID, RFobj, PrObj,Objinbag) {
  t <- which(RFobj[,x] == PrObj[rowID,x])
  t[Objinbag[t,x] == 0]
}

v.obj2$yvar[c(Reduce(c,parallel::mclapply(1:2, f_Withbag,1,v.obj2$membership,pr$membership,v.obj2$inbag))),]

c(Reduce(c,parallel::mclapply(1:1, f_Withbag,1,v.obj2$membership,pr$membership,v.obj2$inbag)))
c(Reduce(c,parallel::mclapply(2:2, f_Withbag,1,v.obj2$membership,pr$membership,v.obj2$inbag)))

out2[out2$TreeId == 1,]$TrainId
out2[out2$TreeId == 2,]$TrainId

pr$membership

v.obj2$membership[1 == 7,]
v.obj2$inbag[479:480,]