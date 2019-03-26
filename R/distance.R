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


# #' print distance(a clusDist object) between clustering results
# #' @param x a clusDist object
# #' @param diag whether to plot diagonal of distance matrix
# #' @param upper whether to plot upper half of distance matrix
# #' @param digits precision digits to retain
# #' @param justify should a character vector be left-justified (the default), right-justified, centred or left alone. Can be abbreviated.
# #' @param right logical, indicating whether or not strings should be right aligned.
# #' @param ... remain to be processed
# #' @export
# print.clusDist <- function(x, diag = NULL, upper = NULL,
# 	     digits = getOption("digits"), justify = "none", right = TRUE, ...){
#   if(length(x)) {
# 	  if(is.null(diag))
# 	    diag <- if(is.null(a <- x$Diag)) FALSE else a
# 	if(is.null(upper))
# 	    upper <- if(is.null(a <- x$Upper)) FALSE else a
# 
# 	m <- as.matrix(x$distance)
# 	cf <- format(m, digits = digits, justify = justify)
# 	if(!upper)
# 	    cf[row(cf) < col(cf)] <- ""
# 	if(!diag)
# 	    cf[row(cf) == col(cf)] <- ""
# 
# 	## Better: use an improved prettyNum() function -> ../../base/R/format.R
# 	##-	if(any((i <- m == floor(m))))
# 	##-	    cf[i] <- sub("0+$", "", cf[i])
# 	print(if(diag || upper) cf else cf[-1, -x$Size, drop = FALSE],
# 	      quote = FALSE, right = right, ...)
#   } else {
# 	  cat(data.class(x),"(0)\n", sep = "")
#   }
#   invisible(x)
# }

#' transform clusDist object to Gram matrix
#' @param x a clusDist object
#' @examples 
#' data <- rbind(matrix(rnorm(100, mean = 10, sd = 2), nrow = 50),
#'               matrix(rnorm(100, mean = 0, sd = 1), nrow = 50),
#'               matrix(rnorm(100, mean = -10, sd = 3), nrow = 50)
#'               )
#' clustering1 <- stats::hclust(dist(data),method='single')
#' clustering2 <- kmeans(data,centers=3)
#' clustering3 <- dbscan::dbscan(data,eps=.8)
#' distance <- clusDist(clustering1,clustering2,clustering3)
#' res <- distToGram(distance)
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
  branchComponent1 <- getClusterTree(clustering1)$tree
  branchComponent2 <- getClusterTree(clustering2)$tree
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

