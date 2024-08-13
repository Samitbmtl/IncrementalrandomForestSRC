
p_PathSim <- 'C:\\Article2b\\SimData\\Simulations\\'

p_TrainingFile <- 'dataTrainingFplusplus.csv'
p_TestFile <- 'dataTestR.csv'
p_ResultatFile <- 'ResultatFinal_All_R.csv'

for (i in 1 : 10)
{
  d_training <-read.csv(paste(p_PathSim,i,"\\",p_TrainingFile,sep=""), header = TRUE, sep = ";")
  d_test <- read.csv(paste(p_PathSim,i,"\\",p_TestFile,sep=""), header = TRUE, sep = ";")
  
  d_training$logt = log(d_training$T)
  d_test$logt = log(d_test$T)  
  
  v.obj1 <- irfsrc(indLSC(logt, G) ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data = d_training, 
                   ntree = 100, block.size = 1,nodesize = 50, do.trace = FALSE, membership = TRUE, statistics = TRUE)
  
  v.obj2 <- irfsrc(indLinearC(logt, G) ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data = d_training, 
                   ntree = 100, block.size = 1,nodesize = 50, do.trace = FALSE, membership = TRUE, statistics = TRUE)
  
  v.obj3 <- irfsrc(indQuadC(logt, G) ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data = d_training,
                   ntree = 100, block.size = 1,nodesize = 50, do.trace = FALSE, membership = TRUE, statistics = TRUE)
  
  v.obj4 <- irfsrc(indQuad2C(logt, G) ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data = d_training,
                   ntree = 100, block.size = 1,nodesize = 50, do.trace = FALSE, membership = TRUE, statistics = TRUE)
  
  v.obj5 <- irfsrc(indTSTC(logt, G) ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data = d_training,
                   ntree = 100, block.size = 1,nodesize = 50, do.trace = FALSE, membership = TRUE, statistics = TRUE)
  
  v.obj6 <- irfsrc(indNSRC(logt, G) ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data = d_training,
                   ntree = 100, block.size = 1,nodesize = 50, do.trace = FALSE, membership = TRUE, statistics = TRUE)
  
  v.obj7 <- irfsrc(indTST2C(logt, G) ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data = d_training,
                   ntree = 100, block.size = 1,nodesize = 50, do.trace = FALSE, membership = TRUE, statistics = TRUE)
  
  v.obj8 <- irfsrc(indMOB0C(logt, G) ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data = d_training,
                   ntree = 100, block.size = 1,nodesize = 50, do.trace = FALSE, membership = TRUE, statistics = TRUE)
  
  v.obj9 <- irfsrc(indMOBLC(logt, G) ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data = d_training,
                   ntree = 100, block.size = 1,nodesize = 50, do.trace = FALSE, membership = TRUE, statistics = TRUE)
  
  v.obj10 <- irfsrc(indMOBQC(logt, G) ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data = d_training,
                   ntree = 100, block.size = 1,nodesize = 50, do.trace = FALSE, membership = TRUE, statistics = TRUE)
  
  v.obj11 <- irfsrc(indMAXLC(logt, G) ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data = d_training,
                   ntree = 100, block.size = 1,nodesize = 50, do.trace = FALSE, membership = TRUE, statistics = TRUE)
  
  v.obj12 <- irfsrc(indMAXQC(logt, G) ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data = d_training,
                   ntree = 100, block.size = 1,nodesize = 50, do.trace = FALSE, membership = TRUE, statistics = TRUE)
  
  v.obj13 <- irfsrc(indMAXMAXQC(logt, G) ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data = d_training,
                   ntree = 100, block.size = 1,nodesize = 50, do.trace = FALSE, membership = TRUE, statistics = TRUE)
  
  pr1 <- predict(v.obj1, newdata = d_test, membership = TRUE)
  pr2 <- predict(v.obj2, newdata = d_test, membership = TRUE)
  pr3 <- predict(v.obj3, newdata = d_test, membership = TRUE)
  pr4 <- predict(v.obj4, newdata = d_test, membership = TRUE)
  pr5 <- predict(v.obj5, newdata = d_test, membership = TRUE)
  pr6 <- predict(v.obj6, newdata = d_test, membership = TRUE)
  pr7 <- predict(v.obj7, newdata = d_test, membership = TRUE)
  
  pr8 <- predict(v.obj3, newdata = d_test, membership = TRUE)
  pr9 <- predict(v.obj4, newdata = d_test, membership = TRUE)
  pr10 <- predict(v.obj5, newdata = d_test, membership = TRUE)
  pr11 <- predict(v.obj6, newdata = d_test, membership = TRUE)
  pr12 <- predict(v.obj7, newdata = d_test, membership = TRUE)
  pr13 <- predict(v.obj7, newdata = d_test, membership = TRUE)

  ds <- data.frame(Reallift = d_test$lift ,BetaLSC= pr1$predicted[,2] ,BetaMSE_RegLSC_R= (pr1$predicted[,2] - d_test$lift)^2 
                   ,BetaLinearC= pr2$predicted[,2] ,BetaMSE_RegLinearC_R= (pr2$predicted[,2] - d_test$lift)^2
                   ,BetaQuadC= pr3$predicted[,2] ,BetaMSE_RegQuadC_R= (pr3$predicted[,2] - d_test$lift)^2 
                   ,BetaQuad2C= pr4$predicted[,2] ,BetaMSE_RegQuad2C_R= (pr4$predicted[,2] - d_test$lift)^2 
                   ,BetaTSTC= pr5$predicted[,2] ,BetaMSE_RegTSTC_R= (pr5$predicted[,2] - d_test$lift)^2 
                   ,BetaNSRC= pr6$predicted[,2] ,BetaMSE_RegNSRC_R= (pr6$predicted[,2] - d_test$lift)^2 
                   ,BetaTST2C= pr7$predicted[,2] ,BetaMSE_RegTST2C_R= (pr7$predicted[,2] - d_test$lift)^2 
                   
                   ,BetaMOB0C= pr2$predicted[,2] ,BetaMSE_RegMOB0C_R= (pr2$predicted[,2] - d_test$lift)^2
                   ,BetaMOBLC= pr3$predicted[,2] ,BetaMSE_RegMOBLC_R= (pr3$predicted[,2] - d_test$lift)^2 
                   ,BetaMOBQC= pr4$predicted[,2] ,BetaMSE_RegMOBQC_R= (pr4$predicted[,2] - d_test$lift)^2 
                   ,BetaMAXLC= pr5$predicted[,2] ,BetaMSE_RegMAXLC_R= (pr5$predicted[,2] - d_test$lift)^2 
                   ,BetaMAXQC= pr6$predicted[,2] ,BetaMSE_RegMAXQC_R= (pr6$predicted[,2] - d_test$lift)^2 
                   ,BetaMAXMAXQC= pr7$predicted[,2] ,BetaMSE_RegMAXMAXQC_R= (pr7$predicted[,2] - d_test$lift)^2
                   )
  
  write.table(ds, file = paste(p_PathSim,i,"\\",p_ResultatFile,sep=""), sep = ";",quote = FALSE,row.names = FALSE, qmethod = "escape")
  
  print(i)
}




ResultatCompile <- matrix(ncol=15, nrow=10)
colnames(ResultatCompile) <-c("Simulation","Nbr","BetaMSE_RegLSC_R","BetaMSE_RegLinearC_R","BetaMSE_RegQuadC_R","BetaMSE_RegQuad2C_R","BetaMSE_RegTSTC_R","BetaMSE_RegNSRC_R","BetaMSE_RegTST2C_R","BetaMSE_RegMOB0C_R","BetaMSE_RegMOBLC_R","BetaMSE_RegMOBQC_R","BetaMSE_RegMAXLC_R","BetaMSE_RegMAXQC_R","BetaMSE_RegMAXMAXQC_R")

for (i in 1 : 10)
{
  
  ResultatFinalPath <- paste(p_PathSim,'\\',i,'\\ResultatFinal_All_R.csv' ,sep="")
  if (file.exists(ResultatFinalPath)) 
  {
    tab <- read.csv(ResultatFinalPath, header = TRUE, sep = ";")
    ResultatCompile[i,] <-  cbind(i,nrow(tab),mean(tab$BetaMSE_RegLSC_R),mean(tab$BetaMSE_RegLinearC_R),mean(tab$BetaMSE_RegQuadC_R),mean(tab$BetaMSE_RegQuad2C_R),mean(tab$BetaMSE_RegTSTC_R),mean(tab$BetaMSE_RegNSRC_R),mean(tab$BetaMSE_RegTST2C_R)
                                  ,mean(tab$BetaMSE_RegMOB0C_R),mean(tab$BetaMSE_RegMOBLC_R),mean(tab$BetaMSE_RegMOBQC_R),mean(tab$BetaMSE_RegMAXLC_R),mean(tab$BetaMSE_RegMAXQC_R),mean(tab$BetaMSE_RegMAXMAXQC_R))
  }
}
write.table(ResultatCompile, file = paste(p_PathSim,"\\ResultatCompile_All_R.csv",sep=""), sep = ";",quote = FALSE,row.names = FALSE, qmethod = "escape")

