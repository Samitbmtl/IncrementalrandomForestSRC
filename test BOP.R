getwd()
setwd("C:\\Article3a")

#"C:/Article3a/incrementalRandomForestSRC"

library(tools)

tools::package_native_routine_registration_skeleton("IncrementalrandomForestSRC")

library(IncrementalrandomForestSRC)

IncrementalrandomForestSRC::iex
out <- iextractBBOPid(v.obj2$membership,pr$membership,v.obj2$inbag,0,1)
out <- iextractBBOP(v.obj2$yvar,v.obj2$membership,pr$membership,v.obj2$inbag,0,1)
out <- iextractBeta(v.obj2$yvar,v.obj2$membership,pr$membership,v.obj2$inbag,0,1)

pr$predicted

out <- iextractBBOPid(myobjmem,myprmem,myobjinbag,1,1)
out <- iextractBBOPid(myobjmem,myprmem,myobjinbag,1,0)
out <- iextractBBOPid(myobjmem,myprmem,myobjinbag,0,1)

object.size(out) / 1024^2
out <- iextractBBOP(myobjyvar,myobjmem,myprmem,myobjinbag,0,1)

out <- iextractBeta(myobjyvar,myobjmem,myprmem,myobjinbag,0,1)

out[1:10,]
)
subset(subset(out, TestId==1) , TreeId ==1)$TrainId

c("a",names(d))

gc()

1000 * 1000 * 3
a <- out$sami

out <- NULL
gc() 

matrix(c(1,2,4,3,3,3,1,2,2,1,2,1,1,4,2,2,1,1,2,4), c(10, 2))

t(matrix(c(1,2,4,3,3,3,1,2,2,1,2,1,1,4,2,2,1,1,2,4), c(2, 10)))

as.integer(as.vector(myprmem))

myobjmem <-  matrix(c(1,2,4,3,3,3,1,2,2,1,2,1,1,4,2,2,1,1,2,4), nrow = 10, ncol = 2)

myprmem <- matrix(c(4,3,3,2,3,1), nrow = 3, ncol = 2)

myobjinbag <- matrix(c(0,1,0,1,0,1,0,1,0,3,1,0,1,0,1,0,1,0,1,0), nrow = 10, ncol = 2)

myobjyvar <- matrix(c(1.1,1,1,2.2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9,10,10,10), nrow = 10, ncol = 3)
myobjyvar <- matrix(c(1.1,1,2,2.2,3,3.7,4,4.3,5,5,6.1,6,7,.87,8,.8,9,1.9,10.1,10), nrow = 10, ncol = 2)

as.integer(as.vector(v.obj2$inbag))
as.double(as.vector(as.matrix(v.obj2$yvar)))

as.double(as.vector(t(as.matrix(v.obj2$yvar))))