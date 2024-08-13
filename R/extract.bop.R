iextractBBOPid <- function(ObjMembership, PrMembership, ObjInbag,getInBag,obsUniqueInTree)
{
  ntrain = nrow(ObjMembership)
  ntrees = ncol(ObjMembership)
  ntest = nrow(PrMembership)
  
  if (ncol(ObjInbag) != ntrees )
  {
    stop("trees number in ObjMembership not the same as ObjInbag !")
  }
  if (ncol(PrMembership) != ntrees )
  {
    stop("trees number in ObjMembership not the same as PrMembership !")
  }
  if (ntrain != nrow(ObjInbag))
  {
    stop("observation number not the same in ObjMembership and ObjInbag !")
  }
  if (getInBag != 0 & getInBag != 1)
  {
    stop("getInBag should be 0 or 1 !")
  }
  if (obsUniqueInTree != 0 & obsUniqueInTree != 1)
  {
    stop("obsUniqueInTree should be 0 or 1 !")
  }
  
  nativeOutput <- tryCatch({.Call("extractBBOPid",
                                  ntrain,
                                  ntest,
                                  ntrees,
                                  as.integer(as.vector(t(ObjMembership))),
                                  as.integer(as.vector(t(ObjInbag))),
                                  as.integer(as.vector(t(PrMembership))),
                                  as.integer(getInBag),           
                                  as.integer(obsUniqueInTree)
                                  )
  
  })
  
  if (is.null(nativeOutput)) {
    stop("An error has occurred in extractBBOPid.")
  }
  
  dat <- data.frame(matrix(nativeOutput$res, byrow = TRUE,  ncol = 3))
  names(dat) = c("TestId","TreeId","TrainId")
  
  return (dat)
}

iextractBBOP <- function(Objyvar,ObjMembership, PrMembership, ObjInbag,getInBag,obsUniqueInTree)
{
  ntrain = nrow(ObjMembership)
  ntrees = ncol(ObjMembership)
  ntest = nrow(PrMembership)
  nYcols = ncol(Objyvar)
  
  if (ncol(ObjInbag) != ntrees )
  {
    stop("trees number in ObjMembership not the same as ObjInbag !")
  }
  if (ncol(PrMembership) != ntrees )
  {
    stop("trees number in ObjMembership not the same as PrMembership !")
  }
  if (ntrain != nrow(ObjInbag))
  {
    stop("observation number not the same in ObjMembership and ObjInbag !")
  }
  if (ntrain != nrow(Objyvar))
  {
    stop("observation number not the same in ObjMembership and Objyvar !")
  }
  if (getInBag != 0 & getInBag != 1)
  {
    stop("getInBag should be 0 or 1 !")
  }
  if (obsUniqueInTree != 0 & obsUniqueInTree != 1)
  {
    stop("obsUniqueInTree should be 0 or 1 !")
  }
  
  nativeOutput <- tryCatch({.Call("extractBBOP",
                                  ntrain,
                                  ntest,
                                  ntrees,
                                  nYcols,
                                  as.integer(as.vector(t(ObjMembership))),
                                  as.integer(as.vector(t(ObjInbag))),
                                  as.integer(as.vector(t(PrMembership))),
                                  as.double(as.vector(t(as.matrix(Objyvar)))),
                                  as.integer(getInBag),           
                                  as.integer(obsUniqueInTree)
  )
    
  })
  
  if (is.null(nativeOutput)) {
    stop("An error has occurred in extractBBOPid.")
  }
  
  dat <- data.frame(matrix(nativeOutput$res, byrow = TRUE,  ncol = (3+nYcols)))
  names(dat) = c("TestId","TreeId","TrainId",names(Objyvar))
  
  return (dat)
}

iextractBeta <- function(Objyvar,ObjMembership, PrMembership, ObjInbag,getInBag,obsUniqueInTree)
{
  ntrain = nrow(ObjMembership)
  ntrees = ncol(ObjMembership)
  ntest = nrow(PrMembership)
  nYcols = ncol(Objyvar)
  
  if (ncol(ObjInbag) != ntrees )
  {
    stop("trees number in ObjMembership not the same as ObjInbag !")
  }
  if (ncol(PrMembership) != ntrees )
  {
    stop("trees number in ObjMembership not the same as PrMembership !")
  }
  if (ntrain != nrow(ObjInbag))
  {
    stop("observation number not the same in ObjMembership and ObjInbag !")
  }
  if (ntrain != nrow(Objyvar))
  {
    stop("observation number not the same in ObjMembership and Objyvar !")
  }
  if (getInBag != 0 & getInBag != 1)
  {
    stop("getInBag should be 0 or 1 !")
  }
  if (obsUniqueInTree != 0 & obsUniqueInTree != 1)
  {
    stop("obsUniqueInTree should be 0 or 1 !")
  }
  if (nYcols != 2 )
  {
    stop("Objyvar should have T and G columns !")
  }
  
  nativeOutput <- tryCatch({.Call("extractBeta",
                                  ntrain,
                                  ntest,
                                  ntrees,
                                  nYcols,
                                  as.integer(as.vector(t(ObjMembership))),
                                  as.integer(as.vector(t(ObjInbag))),
                                  as.integer(as.vector(t(PrMembership))),
                                  as.double(as.vector(t(as.matrix(Objyvar)))),
                                  as.integer(getInBag),           
                                  as.integer(obsUniqueInTree)
  )
    
  })
  
  if (is.null(nativeOutput)) {
    stop("An error has occurred in extractBBOPid.")
  }
  
  dat <- data.frame(matrix(nativeOutput$res, byrow = TRUE,  ncol = 6))
  names(dat) = c("b0int","b0lin","b1lin","b0quad","b1quad","b2quad")
  
  return (dat)
}

