#' calculate distance between clustering results
#' @param clustering1 result of some clustering, for example output from hclust(). A clustering
#' can also be an n by m  matrix, where n is the number of data points and m is the number
#' of levels in the clustering hierarchy. 
#' @param clustering2 result of a second clustering, to be combined with the first.
#' @param ... results of other clustering methods, to be combined with the first two.
#' @param diag whether to plot the diagonal of distance matrix
#' @param upper whether to plot the upper half of diagonal matrix
#' @examples 
#' data <- rbind(matrix(rnorm(100, mean = 10, sd = 2), nrow = 50),
#'               matrix(rnorm(100, mean = 0, sd = 1), nrow = 50),
#'               matrix(rnorm(100, mean = -10, sd = 3), nrow = 50)
#'               )
#' clustering1 <- stats::hclust(dist(data),method='single')
#' clustering2 <- kmeans(data,centers=3)
#' clustering3 <- dbscan::dbscan(data,eps=.8)
#' res <- clusDist(clustering1,clustering2,clustering3)
#' @export
clusDist <- function(clustering1, clustering2, ..., diag = FALSE, upper = FALSE){
  if (missing(clustering1)) stop("Must provide output of clustering for first argument")
  if (missing(clustering2)) stop("Must provide output of clustering for second argument")
  clusterTrees <- list(clustering1, clustering2, ...)
  n <- length(clusterTrees)
  #for (ii in c(1:n)) {
  #  if(!(class(clusterTrees[[ii]])=='clusterTree')[1])
  #    stop(ii,'th argument is not clusterTree object')
  #}
  distance <- rep(1, length.out = n*(n-1)/2L)
  for (ii in c(1:(n-1))) {
    count <- (2*n - ii)*(ii-1)/2L
    for (jj in c((ii+1):n)) {
      distance[count + jj - ii] <- getClusDist(clusterTrees[[ii]], clusterTrees[[jj]])
    }
  }
  # rownames(distance) <- c(1:n)
  # colnames(distance) <- c(1:n)
  attr(distance,'Size') <- n
  attr(distance,'Labels') <- dimnames(distance)[[1L]]
  attr(distance,'Diag') <- diag
  attr(distance,'Upper') <- upper
  attr(distance, 'method') <- 'clustering distance'
  attr(distance, 'call') <- match.call()
  attr(distance, 'class') <- 'dist'
  class(distance) <- 'dist'
  distance
}


#' transform dist object to Gram matrix
#' @param x a dist object
#' @examples 
#' data <- matrix(data = rnorm(300),nrow = 100)
#' distance <- dist(data)
#' g <- distToGram(distance)
#' m <- as.matrix(scale(data, center = TRUE, scale = FALSE) )
#' g2 <- m %*% t(m)
#' round(sum(g - g2), 7)
#' 
#' # Using the between cluster distances
#' set.seed(3141593)
#' data <- rbind(matrix(rnorm(100, mean = 10, sd = 2), nrow = 50),
#'               matrix(rnorm(100, mean = 0, sd = 1), nrow = 50),
#'               matrix(rnorm(100, mean = -10, sd = 3), nrow = 50))
#'               
#' clustering1 <- stats::hclust(dist(data),method='single')
#' clustering2 <- kmeans(data,centers=3)
#' clustering3 <- dbscan::dbscan(data,eps=.8)
#' library(mclust)
#' clustering4 <- Mclust(data)
#' distance <- clusDist(clustering1,clustering2,
#'                      clustering3, clustering4)
#' distance
#' gram <- distToGram(distance)
#' decomp <- eigen(gram)
#' evals <- eigen(gram)$values
#' coords <- eigen(gram)$vectors
#' savePar <- par(mfrow = c(1,2))
#' plot(evals/max(evals), type ="b", ylab = "contribution",
#'      main = "Contributions to dimensionality",
#'      sub = "Only two dimensions needed")
#' plot(coords[, 1:2], pch = 0:3, cex = 3,
#'      xlab = "Var 1", ylab = "Var 2",
#'      main = "Comparing clusters in cluster space",
#'      sub = "kmeans and model based agree")
#' par(savePar)
#' @export
distToGram <- function(x){
  if(!(class(x)=='dist'))
      stop('input argument is not a clusDist object')
  n <- attr(x,'Size')
  # transform_matrix <- diag(nrow = n, ncol = n) - 1/n * matrix(data = 1, nrow = n, ncol = n)
  x <- as.matrix(x)
  res <- x * x 
  rowAverage <- colSums(res)/n
  res <- sweep(res,2,rowAverage)
  colAverage <- rowSums(res)/n
  res <- sweep(res,1,colAverage)
  res <- - res / 2
  res
}

#' Calculate distance of two hierarchical clustering matrix
#' @param clustering1 The first clustering as a branch component
#' @param clustering2 The second clustering as a branch component
#' @return distance between 0 and square root of 2
#' @examples
#' data <- rbind(matrix(rnorm(100, mean = 10, sd = 2), nrow = 50),
#'               matrix(rnorm(100, mean = 0, sd = 1), nrow = 50),
#'               matrix(rnorm(100, mean = -10, sd = 3), nrow = 50)
#'               )
#' clust1 <- stats::hclust(dist(data),method='complete')
#' clust2 <- stats::hclust(dist(data),method='single')
#' clusDist(clust1,clust2)
#' @export
#' 
getClusDist <- function(clustering1, clustering2)
{
  branchComponent1 <- getClusterTree(clustering1)$treeMatrix
  branchComponent2 <- getClusterTree(clustering2)$treeMatrix
  branchComponent1[is.na(branchComponent1)] <- 0
  branchComponent2[is.na(branchComponent2)] <- 0
  n <- nrow(branchComponent1)
  w1 <- vector(mode = 'integer', length = n*(n-1)/2)
  for(j in 2:n)
  {
    for(k in 1:(j-1))
    {
      w1[(j-1)*(j-2)/2+k] <- sum((branchComponent1[j,]!=0)&
                                   (branchComponent1[k,]!=0)&
                                   (branchComponent1[j,]==branchComponent1[k,]))
    }
  }
  w2 <- vector(mode = 'integer', length = n*(n-1)/2)
  for(j in 2:n)
  {
    for(k in 1:(j-1))
    {
      w2[(j-1)*(j-2)/2+k] <- sum((branchComponent2[j,]!=0)&
                                   (branchComponent2[k,]!=0)&
                                   (branchComponent2[j,]==branchComponent2[k,]))
    }
  }
  #w1<- w1 - 1
  if (sum(w1^2) > 0) {
    w1 <- w1/sqrt(sum(w1^2))
  }
  
  #w2 <- w2 - 1
  if (sum(w2^2) > 0) {
    w2 <- w2/sqrt(sum(w2^2))
  }
  sqrt(sum((w1-w2)^2))
  # class(result) <- c("clusDist", "dist")
}

