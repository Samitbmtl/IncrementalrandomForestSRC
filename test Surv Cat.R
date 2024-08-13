library(survival)

d <- data(colon)
d <- colon

DsEvnt1 <- subset(d, etype==2)
DsEvnt1 <- subset(DsEvnt1, rx != 'Lev+5FU')

#0 - j'ai 625 enregistrement (obs et taitement)
DsEvnt1_1 <- na.omit(DsEvnt1)

#1 - On cre une variable Control = 0 pour obsevation sinon 1
DsEvnt1_1$G <- ifelse(DsEvnt1_1$rx == "Obs",0,1)
DsEvnt1_1$T <- DsEvnt1_1$time 
DsEvnt1_1$C <- DsEvnt1_1$status

names(DsEvnt1_1)



v.obj <- irfsrc(indSurvB(T, C, G) ~ sex+age+obstruct+perfor+adhere+nodes+differ+extent+surg+node4, data = DsEvnt1_1, , incremental.minObsControl = 10 , incremental.minObsTreatment = 10 , incremental.TreatmentLowBound = 0.25,incremental.TreatmentUpperBound = 0.75,
                 ntree = 1, block.size = 1,nsplit = 0,nodesize = 50,mtry =10,bootstrap = "none", do.trace = FALSE, membership = TRUE, statistics = TRUE)

pr <- predict(v.obj, membership = TRUE, OOB = 1)
pr <- predict(v.obj, membership = TRUE, OOB = 0)

pr$predicted

### 
v.objC <- irfsrc(indSurvCat(T, C, G) ~ sex+age+obstruct+perfor+adhere+nodes+differ+extent+surg+node4, data = DsEvnt1_1, , incremental.minObsControl = 10 , incremental.minObsTreatment = 10 , incremental.TreatmentLowBound = 0.25,incremental.TreatmentUpperBound = 0.75,
                ntree = 1, block.size = 1,nsplit = 0,nodesize = 50,mtry =10,bootstrap = "none", do.trace = FALSE, membership = TRUE, statistics = TRUE)


prC <- predict(v.objC, membership = TRUE, OOB = 1)
prC <- predict(v.objC, membership = TRUE, OOB = 0)

prC$predicted

#prC$predicted[1:10,] - pr$predicted[1:10] 

iextractBBOP(v.obj$yvar, v.obj$membership,  )

d <- data(veteran)
d <- veteran
d$trt <- ifelse(d$trt == 2,1,0)

v.obj2 <- irfsrc(indSurvB(time, status, trt) ~ karno+diagtime+age+prior, data = d, incremental.minObsControl = 10 , 
                 ntree = 1, block.size = 1,nsplit = 0,nodesize = 20,mtry = 4 ,bootstrap = "none", do.trace = FALSE, membership = TRUE, statistics = TRUE)

pr <- predict(v.obj2, membership = TRUE)

pr$predicted


v.obj2C <- irfsrc(indSurvCat(time, status, trt) ~ karno+diagtime+age+prior, data = d, incremental.minObsControl = 10 , 
                 ntree = 1, block.size = 1,nsplit = 0,nodesize = 20,mtry = 4 ,bootstrap = "none", do.trace = FALSE, membership = TRUE, statistics = TRUE)

prC2 <- predict(v.obj2C, membership = TRUE)

prC2$predicted

pr$predicted

#prC2$predicted[1:10] - pr$predicted[1:10] 

warnings()


#################################

ds <- dgpCat(10)

fun_test<-function(x){
  print(x)
  #val = lm(y ~ idobs, data=as.data.frame(x))
  return(1)
}
ds[1,"x3"]

subset(ds, g == 0 & idobs == 3)

as.data.frame(ds)
aggregate(ds, list(ds$g), FUN=mean) 
aggregate(ds[,c("y","idobs")], list(ds$g), FUN=fun_test)
aggregate(y ~ idobs, data=ds, list(ds$g), FUN=mean) 