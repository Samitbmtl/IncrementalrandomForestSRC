f_Withbag <- function(x,rowID, RFobj, PrObj,Objinbag, getInBag) {
  t <- which(RFobj[,x] == PrObj[rowID,x])
  t[Objinbag[t,x] == getInBag]
}
fun_IntSurvivalNoTail<-function(KM){
  if (length(KM$time)==0) {return(0)}
  if (length(KM$time)==1) {return(KM$time[1])}
  return(KM$time[1] + sum(sapply(1:(length(KM$time)-1), function(i,s) (s$surv[i] * (s$time[i+1] - s$time[i])),s=KM)))
}
fun_IntSurvival<-function(KM){
  val <- fun_IntSurvivalNoTail(KM)
  if (KM$surv[length(KM$time)] !=0)
  {
    val = val - (KM$time[length(KM$time)] * KM$surv[length(KM$time)])/log(KM$surv[length(KM$time)])
  }
  return(val)
}
fun_pred <- function(rowID,dat, RFobj, PrObj,Objinbag, getInBag) {
  BOPdata <- dat[c(Reduce(c,parallel::mclapply(1:length(RFobj[1,]), f_Withbag,rowID,RFobj,PrObj,Objinbag,getInBag))),]
  colnames(BOPdata) <- c("T","C","G")
  
  if (nrow(BOPdata) == 0)
  {
    stop("No rows in BBOP. Check OOB parameter/bootstrap !")
  }
  KM_T <- survival::survfit(survival::Surv(T, C) ~ 1, type="kaplan-meier", data=BOPdata[which(BOPdata$G==1),])
  KM_C <- survival::survfit(survival::Surv(T, C) ~ 1, type="kaplan-meier", data=BOPdata[which(BOPdata$G==0),])
  valTreeT <- fun_IntSurvival(KM_T)
  valTreeC <- fun_IntSurvival(KM_C)
  return (valTreeT - valTreeC)
}

fun_pred_indLinearC <- function(rowID,dat, RFobj, PrObj,Objinbag) {
  BOPdata <- dat[c(Reduce(c,parallel::mclapply(1:length(RFobj[1,]), f_Withbag,rowID,RFobj,PrObj,Objinbag))),]
  colnames(BOPdata) <- c("T","G")
  mod <- lm(T ~ G, data=BOPdata)
  return (coef(mod)["G"])
  # return (paste(c(Reduce(c,parallel::mclapply(1:length(RFobj[1,]), f_Withbag,rowID,RFobj,PrObj,Objinbag))),collapse= ";"))
}

fun_pred_indQuadC <- function(rowID,dat, RFobj, PrObj,Objinbag) {
  BOPdata <- dat[c(Reduce(c,parallel::mclapply(1:length(RFobj[1,]), f_Withbag,rowID,RFobj,PrObj,Objinbag))),]
  colnames(BOPdata) <- c("T","G")
  mod <- lm(T ~ G, data=BOPdata)
  return (coef(mod)["G"])
  # return (paste(c(Reduce(c,parallel::mclapply(1:length(RFobj[1,]), f_Withbag,rowID,RFobj,PrObj,Objinbag))),collapse= ";"))
}

fun_pred_indQuad2C <- function(rowID,dat, RFobj, PrObj,Objinbag) {
  BOPdata <- dat[c(Reduce(c,parallel::mclapply(1:length(RFobj[1,]), f_Withbag,rowID,RFobj,PrObj,Objinbag))),]
  colnames(BOPdata) <- c("T","G")
  mod <- lm(T ~ G, data=BOPdata)
  return (coef(mod)["G"])
  # return (paste(c(Reduce(c,parallel::mclapply(1:length(RFobj[1,]), f_Withbag,rowID,RFobj,PrObj,Objinbag))),collapse= ";"))
}

fun_pred_5beta <- function(rowID,dat, RFobj, PrObj,Objinbag) {
  BOPdata <- as.data.frame(dat[c(Reduce(c,parallel::mclapply(1:length(RFobj[1,]), f_Withbag,rowID,RFobj,PrObj,Objinbag))),])
  #BOPdata <- dat[c(Reduce(c,parallel::mclapply(1:length(RFobj[1,]), f_Withbag,rowID,RFobj,PrObj,Objinbag))),]
  colnames(BOPdata) <- c("T","G")
  return (c(as.numeric(c(mean(BOPdata$T),coef(lm(T ~ G, data=BOPdata)),coef(lm(T ~ G+I(G^2) , data=BOPdata))))))
}

###Ajout du R2
fun_pred_5beta_withR2 <- function(rowID,dat, RFobj, PrObj,Objinbag) {
  BOPdata <- as.data.frame(dat[c(Reduce(c,parallel::mclapply(1:length(RFobj[1,]), f_Withbag,rowID,RFobj,PrObj,Objinbag,TRUE))),])
  #BOPdata <- dat[c(Reduce(c,parallel::mclapply(1:length(RFobj[1,]), f_Withbag,rowID,RFobj,PrObj,Objinbag))),]
  colnames(BOPdata) <- c("T","G")
  lm1 <- lm(T ~ G, data=BOPdata)
  lm2 <- lm(T ~ G+I(G^2) , data=BOPdata)
  return (c(as.numeric(c(mean(BOPdata$T),coef(lm1),coef(lm2),summary(lm1)$r.squared,summary(lm2)$r.squared))))
}

fun_pred_Categorical <- function(predicted) {
  #predicted rece a partir de iextractBBOP
  agg <- aggregate(predicted$y, list(predicted$TestId, predicted$g), FUN=mean) 
  names(agg) <- c("TestId","g","y")
  
  agg0 <- subset(agg,g==0)
  agg0$y0 <- agg0$y
  agg0$y <- NULL
  agg0$g <- NULL
  
  resAgg <- merge(agg,agg0,by=c("TestId")) 
  resAgg$tau <- resAgg$y - resAgg$y0
  resAgg <- subset(resAgg,g!=0)
  
  resAgg$y <- NULL
  resAgg$y0 <- NULL
  
  resAgg <- resAgg[order(resAgg$g),]
  
  res <- reshape(resAgg, direction = "wide", idvar = c("TestId"), timevar = "g")
  return (res)
}

fun_CalculKM <- function(predicted,pTestId,pg) {
  
  dat <- subset(predicted , TestId == pTestId & G == pg)

  if (nrow(dat) == 0)
  {
    stop("No rows to compute KM")
  }
  KM <- survival::survfit(survival::Surv(T, C) ~ 1, type="kaplan-meier", data=dat)
  return (fun_IntSurvival(KM))
}

fun_pred_SurvCat <- function(predicted) {
  #predicted rece a partir de iextractBBOP
  colnames(predicted) <- c("TestId","TreeId","TrainId","T","C","G")
  
  agg <- aggregate(predicted$T, list(predicted$TestId, predicted$G), FUN=mean) 
  names(agg) <- c("TestId","G","y")
  
  for (i in 1:nrow(agg))
  {
    agg[i,"y"] =  fun_CalculKM(predicted, agg[i,"TestId"] , agg[i,"G"])
  }
  
  agg0 <- subset(agg,G==0)
  agg0$y0 <- agg0$y
  agg0$y <- NULL
  agg0$G <- NULL
  
  resAgg <- merge(agg,agg0,by=c("TestId")) 
  resAgg$tau <- resAgg$y - resAgg$y0
  resAgg <- subset(resAgg,G!=0)
  
  resAgg$y <- NULL
  resAgg$y0 <- NULL
  
  resAgg <- resAgg[order(resAgg$G),]
  
  res <- reshape(resAgg, direction = "wide", idvar = c("TestId"), timevar = "G")
  return (res)
}

# f_Withbag_regression <- function(x,rowID, RFobj, PrObj,Objinbag,dat) {
#   t <- which(RFobj[,x] == PrObj[rowID,x])
#   d <- dat[t[Objinbag[t,x] != 0],]
#   #coef(lm(T ~ G, data=d))
#   return (c(as.numeric(c(mean(d$T),coef(lm(T ~ G, data=d)),coef(lm(T ~ G+I(G^2) , data=d))))))
# }
# 
# fun_pred_5beta <- function(rowID,dat, RFobj, PrObj,Objinbag) {
#   colnames(dat) <- c("T","G")
#   bops <- parallel::mclapply(1:length(RFobj[1,]), f_Withbag_regression,rowID,RFobj,PrObj,Objinbag,dat)
#   # return (bops)
#   cof  <-  as.data.frame(do.call(rbind, bops))
#   #return (cof)
#   names(cof) <- c("N_b0int","N_b0lin","N_b1lin","N_b0quad","N_b1quad","N_b2quad")
#   # return (bops)
#   BOPdata <- dat[c(Reduce(c,parallel::mclapply(1:length(RFobj[1,]), f_Withbag,rowID,RFobj,PrObj,Objinbag))),]
#   #colnames(BOPdata) <- c("T","G")
#   return (c(as.numeric(c(mean(BOPdata$T),coef(lm(T ~ G, data=BOPdata)),coef(lm(T ~ G+I(G^2) , data=BOPdata)),mean(cof$N_b0int),mean(cof$N_b0lin),mean(cof$N_b1lin),mean(cof$N_b0quad),mean(cof$N_b1quad),mean(cof$N_b2quad)))))
# }



