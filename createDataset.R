library(causaldrf)
library(survival)
library(TH.data)

##############################


getwd()
#setwd("D:/r-temp")

save(nmes_data, file = "nmes_data.rda")



dat <- nmes_data
dat <- dat[dat$TOTALEXP != 0,]
dat$uid <- seq.int(nrow(dat))
dat$beltuse <- as.numeric(unclass(dat$beltuse))
dat$POVSTALB <- as.numeric(unclass(dat$POVSTALB))
dat$Y = log(dat$TOTALEXP+1)
dat$G = dat$packyears
dat$G = dat$G - min(dat$packyears)
dat$G = dat$G / max(dat$G)

save(dat, file = "nmes_data_modified.rda")


save(colon, file = "colon_data.rda")


d <- data(colon)

d <- colon

DsEvnt1 <- subset(d,etype==2 & rx != 'Lev+5FU')

#0 - j'ai 625 enregistrement (obs et taitement)
DsEvnt1_1 <- na.omit(DsEvnt1)

#1 - On cre une variable Control = 0 pour obsevation sinon 1
DsEvnt1_1$G <- ifelse(DsEvnt1_1$rx == "Obs",0,1)
DsEvnt1_1$T <- DsEvnt1_1$time 
DsEvnt1_1$C <- DsEvnt1_1$status
#DsEvnt1_1$C <- 1

save(DsEvnt1_1, file = "colon_data_modified.rda")


save(colon, file = "colon_data.rda")

save(GBSG2, file = "GBSG2_data.rda")

d <- data(GBSG2)

d <- GBSG2

DsEvnt1 <- d

#0 - j'ai 625 enregistrement (obs et taitement)
DsEvnt1_1 <- na.omit(DsEvnt1)

#1 - On cre une variable Control = 0 pour obsevation sinon 1
DsEvnt1_1$G <- ifelse(DsEvnt1_1$horTh == "no",0,1)
DsEvnt1_1$T <- DsEvnt1_1$time 
DsEvnt1_1$C <- DsEvnt1_1$cens

DsEvnt1_1$menostat <- ifelse(DsEvnt1_1$menostat == "Pre",0,1)
DsEvnt1_1$tgrade <-  ifelse(DsEvnt1_1$tgrade == "I",1, ifelse(DsEvnt1_1$tgrade == "II",2,3))

save(DsEvnt1_1, file = "GBSG2_data_modified.rda")
