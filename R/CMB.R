CMB <- function(formula, data, ntree = 1000,
                  mtry = NULL, ytry = NULL,
                  nodesize = NULL, nodedepth = NULL,
                  splitrule = NULL, nsplit = 10,
                  importance = c(FALSE, TRUE, "none", "permute", "random", "anti"),
                  block.size = if (any(is.element(as.character(importance), c("none", "FALSE")))) NULL else 10,
                  ensemble = c("all", "oob", "inbag"),
                  bootstrap = c("by.root", "by.node", "none", "by.user"),
                  samptype = c("swor", "swr"),  samp = NULL, membership = FALSE,
                  sampsize = if (samptype == "swor") function(x){x * .632} else function(x){x},
                  na.action = c("na.omit", "na.impute"), nimpute = 1,
                  ntime, cause,
                  proximity = FALSE, distance = FALSE, forest.wt = FALSE,
                  xvar.wt = NULL, yvar.wt = NULL, split.wt = NULL, case.wt = NULL, 
                  forest = TRUE,
                  var.used = c(FALSE, "all.trees", "by.tree"),
                  split.depth = c(FALSE, "all.trees", "by.tree"),
                  seed = NULL,
                  do.trace = FALSE,
                  statistics = FALSE,
				  #Sami
				  incremental.CenterResponse = FALSE,
				  incremental.CenterTreatment = FALSE)
{
  #print (formula)
  if (!inherits(formula, "formula")) {
    stop("'formula' is not a formula object.")
  }
  
  var.to.replace <- "f"
  new.terms1 <-  "CMB"
  formula <- as.formula(do.call("substitute", list(formula, setNames(list(str2lang(new.terms1)), var.to.replace))) )
  
  #print (formula)
  
  irfsrc(formula, data, ntree,mtry , ytry,
    nodesize, nodedepth,
    splitrule, nsplit ,
    importance,
    block.size ,
    ensemble,
    bootstrap ,
    samptype,  samp, membership ,
    sampsize,
    na.action, nimpute,
    ntime, cause,
    proximity , distance , forest.wt ,
    xvar.wt, yvar.wt, split.wt, case.wt,
    forest ,
    var.used ,
    split.depth ,
    seed,
    do.trace ,
    statistics ,
	  incremental.CenterResponse ,
	  incremental.CenterTreatment
    )
}