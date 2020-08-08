### Lifespan-analysis for Biologists ###

## Chapter 5.2: BaSTA ##

#Load packages and data
suppressMessages(library(dplyr))
suppressMessages(library(BaSTA))
#library(parallel)  #for better performance on your local machine
#library(snowfall)  #for better performance on your local machine
#script for function 'plotFancyBaSTA' is at bottom of this script file

data1 <- read.csv("https://github.com/ZajitschekTeam/lifespananalysis/raw/master/binder/data/expevol_male_flies.csv")
data1 <- data1 %>% mutate(across(where(is.integer), as.factor))
data1$status <- 1

# Create a subset of three treatment groups/cohorts 
# (3 assay diets (= all that are available) from 1 cagediet)
group1 <- subset(data1, cagediet== 1)

## Setup data for BaSTA ##

basta_v2 <- matrix(0, nrow = nrow(group1), ncol = max(group1$lifespan))
basta_v2 <- as.data.frame(basta_v2)
vbasta3 <- cbind(group1[,c(1,2,8)], basta_v2)
names(vbasta3) <- c("AD", "ED", "age", seq(1, max(group1$lifespan)) )
vbasta3$age <- vbasta3$age -1
#str(vbasta3)
idx <- do.call(rbind, Map(function(a,b) 
    cbind(a, match(1:b, colnames(vbasta3))), 
    seq_along(vbasta3$age), vbasta3$age))
vbasta3[idx]<-1
head(vbasta3)
##below: split ED into 3 separate columns, containing 0 and 1
basta_v <- cbind(seq(1,nrow(vbasta3)), rep(1, nrow(vbasta3)), vbasta3$age+1,  vbasta3[,4:(max(V_surv$age)+3)], vbasta3$ED,vbasta3$ED,vbasta3$ED)
basta_v$"1" <- rep(0, nrow(vbasta3))
str(basta_v)
names(basta_v) <- c("ID", "BIRTH", "DEATH", seq(1, max(V_surv$age)), "ED_LD", "ED_SD", "ED_HD")
str(basta_v)
head(basta_v)
##recode LD,SD,HD into 0 and 1
basta_v2 <- basta_v
basta_v2$ED_LD <- as.character(basta_v2$ED_LD); basta_v2$ED_SD <- as.character(basta_v2$ED_SD); basta_v2$ED_HD <- as.character(basta_v2$ED_HD)
basta_v2$ED_LD[basta_v2$ED_LD == "LD"] = 1
basta_v2$ED_LD[basta_v2$ED_LD != "1"] = 0
basta_v2$ED_SD[basta_v2$ED_SD == "SD"] = 1
basta_v2$ED_SD[basta_v2$ED_SD != "1"] = 0
basta_v2$ED_HD[basta_v2$ED_HD == "HD"] = 1
basta_v2$ED_HD[basta_v2$ED_HD != "1"] = 0
basta_v2$ED_LD <- as.numeric(basta_v2$ED_LD); basta_v2$ED_SD <- as.numeric(basta_v2$ED_SD); basta_v2$ED_HD <- as.numeric(basta_v2$ED_HD)
str(basta_v2)
head(basta_v2)
basta_ADLD <- basta_v2

## V_ADSD ##
#fuer basta: use AD subsets, i.e. only 3 levels of ED
basta_v2 <- matrix(0, nrow = nrow(V_ADSD), ncol = max(V_ADSD$age))
basta_v2 <- as.data.frame(basta_v2)
vbasta3 <- cbind(V_ADSD[,c(1,2,8)], basta_v2)
names(vbasta3) <- c("AD", "ED", "age", seq(1, max(V_ADSD$age)) )
vbasta3$age <- vbasta3$age -1
str(vbasta3)
idx <- do.call(rbind, Map(function(a,b) 
    cbind(a, match(1:b, colnames(vbasta3))), 
    seq_along(vbasta3$age), vbasta3$age))
vbasta3[idx]<-1
head(vbasta3)
##below: split ED into 3 separate columns, containing 0 and 1
basta_v <- cbind(seq(1,nrow(vbasta3)), rep(1, nrow(vbasta3)), vbasta3$age+1,  vbasta3[,4:(max(V_surv$age)+3)], vbasta3$ED,vbasta3$ED,vbasta3$ED)
basta_v$"1" <- rep(0, nrow(vbasta3))
str(basta_v)
names(basta_v) <- c("ID", "BIRTH", "DEATH", seq(1, max(V_surv$age)), "ED_LD", "ED_SD", "ED_HD")
str(basta_v)
head(basta_v)
##recode LD,SD,HD into 0 and 1
basta_v2 <- basta_v
basta_v2$ED_LD <- as.character(basta_v2$ED_LD); basta_v2$ED_SD <- as.character(basta_v2$ED_SD); basta_v2$ED_HD <- as.character(basta_v2$ED_HD)
basta_v2$ED_LD[basta_v2$ED_LD == "LD"] = 1
basta_v2$ED_LD[basta_v2$ED_LD != "1"] = 0
basta_v2$ED_SD[basta_v2$ED_SD == "SD"] = 1
basta_v2$ED_SD[basta_v2$ED_SD != "1"] = 0
basta_v2$ED_HD[basta_v2$ED_HD == "HD"] = 1
basta_v2$ED_HD[basta_v2$ED_HD != "1"] = 0
basta_v2$ED_LD <- as.numeric(basta_v2$ED_LD); basta_v2$ED_SD <- as.numeric(basta_v2$ED_SD); basta_v2$ED_HD <- as.numeric(basta_v2$ED_HD)
str(basta_v2)
head(basta_v2)
basta_ADSD <- basta_v2

## V_ADHD ##
#fuer basta: use AD subsets, i.e. only 3 levels of ED
basta_v2 <- matrix(0, nrow = nrow(V_ADHD), ncol = max(V_ADHD$age))
basta_v2 <- as.data.frame(basta_v2)
vbasta3 <- cbind(V_ADHD[,c(1,2,8)], basta_v2)
names(vbasta3) <- c("AD", "ED", "age", seq(1, max(V_ADHD$age)) )
vbasta3$age <- vbasta3$age -1
str(vbasta3)
idx <- do.call(rbind, Map(function(a,b) 
    cbind(a, match(1:b, colnames(vbasta3))), 
    seq_along(vbasta3$age), vbasta3$age))
vbasta3[idx]<-1
head(vbasta3)
##below: split ED into 3 separate columns, containing 0 and 1
basta_v <- cbind(seq(1,nrow(vbasta3)), rep(1, nrow(vbasta3)), vbasta3$age+1,  vbasta3[,4:(max(V_surv$age)+3)], vbasta3$ED,vbasta3$ED,vbasta3$ED)
basta_v$"1" <- rep(0, nrow(vbasta3))
str(basta_v)
names(basta_v) <- c("ID", "BIRTH", "DEATH", seq(1, max(V_surv$age)), "ED_LD", "ED_SD", "ED_HD")
str(basta_v)
head(basta_v)
##recode LD,SD,HD into 0 and 1
basta_v2 <- basta_v
basta_v2$ED_LD <- as.character(basta_v2$ED_LD); basta_v2$ED_SD <- as.character(basta_v2$ED_SD); basta_v2$ED_HD <- as.character(basta_v2$ED_HD)
basta_v2$ED_LD[basta_v2$ED_LD == "LD"] = 1
basta_v2$ED_LD[basta_v2$ED_LD != "1"] = 0
basta_v2$ED_SD[basta_v2$ED_SD == "SD"] = 1
basta_v2$ED_SD[basta_v2$ED_SD != "1"] = 0
basta_v2$ED_HD[basta_v2$ED_HD == "HD"] = 1
basta_v2$ED_HD[basta_v2$ED_HD != "1"] = 0
basta_v2$ED_LD <- as.numeric(basta_v2$ED_LD); basta_v2$ED_SD <- as.numeric(basta_v2$ED_SD); basta_v2$ED_HD <- as.numeric(basta_v2$ED_HD)
str(basta_v2)
head(basta_v2)
basta_ADHD <- basta_v2

#### BASTA RUNS ###

DataCheck(basta_ADHD, studyStart = 1, studyEnd = 86, silent=FALSE)

#out_basta_ADHDmulti <- multibasta(object= basta_ADHD, studyStart= 1, studyEnd= 86, niter=50000, burnin=5001, thinning=50,  
        models = c("EX", "GO", "LO"), shape = "simple", covarsStruct = "fused", nsim = 4, ncpus = 4, parallel = TRUE, updateJumps=TRUE)

out_basta_ADHD_exp <- basta(object= basta_ADHD, studyStart= 1, studyEnd= 86, niter=50000, burnin=5001, thinning=50,  
        model = "EX", shape = "simple", covarsStruct = "fused", nsim = 4, ncpus = 4, parallel = TRUE, updateJumps=TRUE)
out_basta_ADHD_exp$DIC[5] #15517.83
out_basta_ADHD$DIC[5] #13019.75
out_basta_ADHD_exp$DIC[5] - out_basta_ADHD$DIC[5] #2498.078 : Gompertz much better!

out_basta_ADHD <- basta(object= basta_ADHD, studyStart= 1, studyEnd= 86, niter=50000, burnin=5001, thinning=50,  
        model = "GO", shape = "simple", covarsStruct = "fused", nsim = 4, ncpus = 4, parallel = TRUE, updateJumps=TRUE)

#plot(out_basta_ADHD)
#plot(out_basta_ADHD, plot.trace = FALSE)
plotFancyBaSTA(out_basta_ADHD)
summary(out_basta_ADHD)


#######################

plotFancyBaSTA <-
function(x, open.new = FALSE, ...){
  # This function creates denstity plots for mortality parameters
  # and plots the predicted mortality rates and survival probabilities
  # from a 'basta' object.
  Palette <- c('#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#FF7F00', 
               '#FFFF33', '#A65628', '#F781BF', '#999999')
  if (.Platform$OS.type == "unix") {
    devtype <- quartz 
  } else {
    devtype <- windows
  }
  allThetaNames <- c("a0", "a1", "c", "b0", "b1", "b2")
  allThetaExpr  <- expression(a[0], a[1], c, b[0], b[1], b[2])
  catNames <- names(x$survQuant)
  lenCat <- length(catNames)
  varNames <- substr(colnames(x$params), 1, 2)
  idTh <- is.element(substr(varNames, 1, 1), c("a", "b", "c"))
  thetaMat <- x$params[,idTh]
  model <- as.character(x$modelSpecs['model'])
  shape <- as.character(x$modelSpecs['shape'])
  if ("col" %in% names(list(...))) {
  	if (length(list(...)$col) < lenCat) {
  		warning("Insufficient number of colors. Not all covariates will be displayed.",
  		        call. = FALSE)
  	} else {
	  	Bord <- list(...)$col[1:lenCat]
  	}
  } else {
	Bord <- Palette[1:lenCat]
  }
  Cols <- adjustcolor(Bord, alpha.f = 0.5)
  if (model %in% c("EX", "1")) {
    nTheta          <- 1
    idAllTheta      <- 4
  } else if (model %in% c("GO", "WE", "2", "3")) {
    nTheta          <- 2
    idAllTheta      <- c(4, 5)
  }  else {
    nTheta          <- 3
    idAllTheta      <- c(4, 5, 6)
  }
  if (shape %in% c("Makeham", "2")) {
    nTheta          <- nTheta + 1
    idAllTheta      <- c(3, idAllTheta)
  } else if(shape %in% c("bathtub", "3")) {
    nTheta          <- nTheta + 3
    idAllTheta      <- c(1:3, idAllTheta)
  }

  # Build plots:	
  if (open.new) devtype(width = 6, height = 6)
  op <- par(no.readonly = TRUE)
  layout(mat = matrix(data = c(rep(1:nTheta, each = 2), 
                              rep(rep(c(nTheta + 1, nTheta + 2), 
                              each = nTheta), 2)), 
                          nrow  = nTheta * 2, 
                          ncol  = 3), 
         widths  = rep(2, 3), 
         heights = rep(1, nTheta))
  par(mar=c(3, 3, 0.5, 0.5))
  for(i in 1:nTheta){
    dez <- list()
    ylz <- rep(NA, lenCat)
    xlz <- matrix(0, lenCat, 2)
    for(j in 1:lenCat){
      dez[[catNames[j]]] <- density(thetaMat[, 
                                 paste(allThetaNames[idAllTheta[i]],
                                 catNames[j], sep = "")])
      ylz[j] <- max(dez[[catNames[j]]]$y)
      xlz[j,] <- range(dez[[j]]$x)
    }
    xr <- range(xlz)
    xl <- c(floor(xr[1] * 10) / 10, ceiling(xr[2] * 10) / 10)
    xd <- ceiling(diff(xl) * 10) / 10
    plot(x = dez[[1]], xlab = "", ylab = "", xlim = xl, ylim = c(0, max(ylz)), 
         lwd = 3, axes = FALSE, main = "", col = NA)
    for(j in 1:lenCat) {
      polygon(x = c(dez[[j]]$x, dez[[j]]$x[1]), y = c(dez[[j]]$y, dez[[j]]$y[1]), 
              col = Cols[j], border = Bord[j], lwd = 1.5)
    }
    axis(side = 1, at = seq(xl[1], xl[2], length = 5), 
         line = 0.5, labels = NA, tcl = 0.4)
    axis(side = 1, at = seq(xl[1], xl[2], length = 3), lwd = NA)
    mtext(text = allThetaExpr[idAllTheta[i]], side = 2, line = 0, 
          at = max(ylz) * 0.8, las = 2, cex = 1.25)
  }

  # Plot survival probability:
  par(mar = c(4, 7, 0.5, 0.5))
  xv <- lapply(1:lenCat, function(idcovs) as.numeric(colnames(x$mortQuant[[idcovs]])))
  mxv <- ceiling(max(unlist(xv)) / 5) * 5
  plot(x = c(0, mxv), y = range(0, 1), col = NA, axes = FALSE, xlab = "", ylab = "")
  for(i in 1:lenCat){
    polygon(x = c(xv[[i]], rev(xv[[i]])), 
            y = c(x$survQuant[[i]][2, ], rev(x$survQuant[[i]][3, ])), 
            col = Cols[i], border = Bord[i])
    lines(x = xv[[i]], y = x$survQuant[[i]][1, ], col = Bord[i], lty = 3)
  }
  if (lenCat > 1) {
    legend(x = 'topright', legend = substr(catNames, 2, nchar(catNames)), 
           pch = 15, pt.cex = 3, cex = 1.5, col = Cols, bty = 'n')
  }

  axis(side = 1, at = seq(0, mxv, 5), labels = NA, tcl = 0.4, line = 0.5)
  axis(side = 2, at = seq(0, 1, 0.2), tcl = 0.4, las = 2, cex.axis = 1.2)
  mtext(text = expression(paste("Survival ", italic(S(x)))), 
        side = 2, line = 3.5, cex = 1.25)

  # Plot mortality rates:
  ylmx            <- c(0, round(max(unlist(x$mortQuant)), 1))
  plot(x = c(0, mxv), y = ylmx, col = NA, ylim = ylmx, xlab = "", ylab = "",
       axes = FALSE)
  for(i in 1:lenCat){
    polygon(x = c(xv[[i]], rev(xv[[i]])), 
            y = c(x$mortQuant[[i]][2, ], rev(x$mortQuant[[i]][3, ])), 
            col = Cols[i], border = Bord[i])
    lines(x = xv[[i]], y = x$mortQuant[[i]][1, ], col = Bord[i], lty = 3)
  }
  axis(side = 1, at = seq(0, mxv, 5), labels = NA, tcl = 0.4, line = 0.5)
  axis(side = 1, at = seq(0, mxv, 5), lwd = NA, cex.axis = 1.2)
  axis(side = 2, at = seq(0, ylmx[2], length = 5), 
       tcl = 0.4, las = 2, cex.axis = 1.2)
  mtext(text = expression(paste("Mortality ", italic(mu(x)))), 
        side = 2, line = 3.5, cex = 1.25)
  mtext(text = expression(paste("Age ", italic(x), " (years)")), 
        side = 1, cex = 1.25, line = 3)
  par(op)
}