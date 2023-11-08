## ----opts, echo=FALSE-------------------------------------------------------------------------
knitr::opts_chunk$set(fig.width=8, fig.height=8)
options(width=96)
.format <- knitr::opts_knit$get("rmarkdown.pandoc.to")
.tag <- function(N, cap ) ifelse(.format == "html",
                                 paste("Figure", N, ":",  cap),
                                 cap)

## ----mycss, results="asis", echo=FALSE--------------------------------------------------------
cat('
<style type="text/css">
.figure { text-align: center; }
.caption { font-weight: bold; }
</style>
')

## ----Setup------------------------------------------------------------------------------------
library(SillyPutty)
library(Umpire)
suppressMessages( library(Mercator) )
suppressMessages( library(mclust) ) # for adjusted rand index

## ----genData----------------------------------------------------------------------------------
set.seed(21315)
trueK <- 4
## Set up survival outcome; baseline is exponential
sm <- SurvivalModel(baseHazard=1/5, accrual=5, followUp=1)
## Build a CancerModel with four subtypes
nBlocks <- 20    
cm <- CancerModel(name="cansim",
                  nPossible=nBlocks,
                  nPattern=trueK,
                  OUT = function(n) rnorm(n, 0, 1), 
                  SURV= function(n) rnorm(n, 0, 1),
                  survivalModel=sm)
## Include 100 blocks/pathways that are not hit by cancer
nTotalBlocks <- nBlocks + 100
## Assign values to hyperparameters
## block size
blockSize <- round(rnorm(nTotalBlocks, 100, 30))
## log normal mean hyperparameters
mu0    <- 6
sigma0 <- 1.5
## log normal sigma hyperparameters
rate   <- 28.11
shape  <- 44.25
## block correlation
p <- 0.6
w <- 5
## Set up the baseline Engine
rho <- rbeta(nTotalBlocks, p*w, (1-p)*w)
base <- lapply(1:nTotalBlocks,
               function(i) {
                 bs <- blockSize[i]
                 co <- matrix(rho[i], nrow=bs, ncol=bs)
                 diag(co) <- 1
                 mu <- rnorm(bs, mu0, sigma0)
                 sigma <- matrix(1/rgamma(bs, rate=rate, shape=shape), nrow=1)
                 covo <- co *(t(sigma) %*% sigma)
                 MVN(mu, covo)
               })
eng <- Engine(base)
## Alter the means if there is a hit
altered <- alterMean(eng, normalOffset, delta=0, sigma=1)
## Build the CancerEngine using character strings
object <- CancerEngine(cm, "eng", "altered")
rm(sm, nBlocks, cm, nTotalBlocks, blockSize, mu0, sigma0, rate, shape, p, w, rho, base, eng, altered)

## ----simData----------------------------------------------------------------------------------
trueN <- 144
dset <- rand(object, trueN, keepall = TRUE) # contains two objects
labels <- dset$clinical$CancerSubType # the true clusters/types
d1 <- dset$data # the noise-free simulated data

## ----noiseModel-------------------------------------------------------------------------------
SpecialNoise <- function(nFeat, nu = 0.1, shape = 1.02, scale = 0.05/shape) {
  NoiseModel(nu = nu,
             tau = rgamma(nFeat, shape = shape, scale = scale),
             phi = 0)
}
nm <- SpecialNoise(nrow(d1), nu = 0)
d1 <- blur(nm, d1)
dim(d1)

## ----distancematrix---------------------------------------------------------------------------
tdis <- t(d1)
dimnames(tdis) <- list(paste("Sample", 1:nrow(tdis), sep=''),
                     paste("Feature", 1:ncol(tdis), sep=''))
dis <- dist(tdis)   ## This step is the rate-liomiting factor. Only way to speed up is to use fewerw samples
names(labels) <- rownames(tdis)

## ----eval=FALSE, echo=FALSE, results='hide'---------------------------------------------------
#  dataset <- tdis
#  eucdist <- dis
#  trueGroups <- labels
#  save(eucdist, trueGroups, file="../data/eucdist.rda")

## ----mercViews--------------------------------------------------------------------------------
mercViews <- function(object, main, tag = NULL) {
  opar <- par(mfrow = c(2, 2))
  on.exit(par(opar))
  pts <- barplot(object, main = main)
  if (!is.null(tag)) {
    gt <- as.vector(as.matrix(table(getClusters(object))))
    loc <- pts[round((c(0, cumsum(gt))[-(1 + length(gt))] + cumsum(gt))/2)]
    mtext(tag, side =1, line = 0, at = loc, col = object@palette, font = 2)
  }
  plot(object, view = "tsne", main = "t-SNE")
  plot(object, view = "hclust")
  plot(object, view = "mds", main = "MDS")
}

## ----fig01, fig.cap = .tag(1, "Hierachical Clustering, with four clusters.")------------------
set.seed(1987)
vis <- Mercator(dis, "euclid", "hclust", K = trueK)
palette <- vis@palette[c(1:3, 7, 8, 6, 10, 4, 11, 5, 15, 14, 17:18, 9, 12, 16, 19:24)]
vis@palette <- palette
vis <- addVisualization(vis, "mds")
vis <- addVisualization(vis, "tsne")
mercViews(vis, "Hierarchical Clustering, Five Clusters")

## ----ari--------------------------------------------------------------------------------------
ari.hier <- adjustedRandIndex(labels, vis@clusters)
ari.hier

## ----fig02, fig.cap = .tag(2, "Visualization of true cancer clusters.")-----------------------
truebin <- remapColors(vis, setClusters(vis, labels))
mercViews(truebin, main = "True Cluster Types", 
          tag = unique(sort(labels)))

## ----fig03, fig.cap = .tag(3, "PAM Clustering, K = 4.")---------------------------------------
pc <- pam(dis, k = trueK, diss=TRUE)
pamc <- remapColors(vis, setClusters(vis, pc$clustering))
mercViews(pamc, main = "PAM, K = 4", 
          tag = paste("P", 1:trueK, sep = ""))
ari.pam <- adjustedRandIndex(labels, pamc@clusters)
ari.pam

## ----RandomSilly------------------------------------------------------------------------------
set.seed(12)
y2 <- suppressWarnings(RandomSillyPutty(dis, trueK, N = 100)) ## this is also slow
ari.max <- adjustedRandIndex(truebin@clusters, y2@MX)
ari.min <- adjustedRandIndex(truebin@clusters, y2@MN)
ari.max
ari.min

## ----fig04, fig.cap = .tag(4, "Random SillyPutty clustering,  K = 4.")------------------------
randSillyBin <- remapColors(vis, setClusters(vis, y2@MX))
mercViews(randSillyBin, main = "SillyPutty Cluster Types, K = 4", 
          tag = paste("C", 1:trueK, sep = ""))

## ----fig05, fig.cap = .tag(5, "Cluster assignements with best and worst silhouette widths after random starts.")----
plot(y2, randSillyBin@view[["mds"]], distobj = dis, col = randSillyBin@palette)
summary(y2)

## ----fig06, fig.cap = .tag(6, "Hierarchical Clustering + SillyPutty, K = 4.")-----------------
hierSilly <- SillyPutty(vis@clusters, dis)
hierSillyBin <- remapColors(vis, setClusters(vis, hierSilly@cluster))
mercViews(hierSillyBin, main = "HClust + Silly, k = 4",tag = paste("C", 1:trueK, sep = ""))
ari.Sillyhier <- adjustedRandIndex(labels, hierSillyBin@clusters)
ari.Sillyhier

## ----fig07, fig.width=6, fig.height=5, fig.cap=.tag(7, "Best mean siilhouette width, by number of clusters, found by combining huierarchical clustering with Silly Putty.")----
y <- findClusterNumber(dis, start = 2, end = 12, method = "HCSP")
plot(names(y), y, xlab = "K", ylab = "Silhouette Width", type = "b", lwd = 2, pch = 16)

## ---------------------------------------------------------------------------------------------
sessionInfo()

