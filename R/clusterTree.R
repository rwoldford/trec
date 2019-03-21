# transfrom a matrix into a clusterTree class object
matrixToClusterTree <- function(x, labels = NULL){
  if (!is.matrix(x)) {
    # transform a vector if it is not a matrix
    if (is.vector(x)) {
      x <- matrix(x, ncol = 1)
    } else stop("Argument must be a matrix or vector")
  }
  # assign column names and row names to this matrix
  colnames(x) <- paste("Level", 1:ncol(x))
  if (is.null(labels)) labels <- paste("object", 1:nrow(x))
  if (length(labels)!=nrow(x))
    stop("length of labels and rows of matrix do not match!")
  rownames(x) <- labels
  ####
  # map all zeros to NA
  mapZerosToNA <- function(x){
    for (i in 1:length(x)) {
      if(x[i]==0){
        x[i:length(x)] <- NA
        break
      }
    }
    x
  }
  if(ncol(x)==1){
    x <- as.matrix(apply(x, 1, mapZerosToNA))
  }
  else{
    x <- t(apply(x, 1, mapZerosToNA))
  }
  ####
  ####
  tree <- list(tree = x, labels = labels)
  class(tree) <- c("clusterTree", class(tree))
  tree
}

#' Get the clusterTree from the result of (single/average/complete) linkage clustering.
#' @param merge The n-1 by 2 matrix that is the merge component output from
#' a hierarchical clustering which describes how the n points are merged in the cluster
#' hierarchy.
#' @return A matrix which represent the clustering result. This matrix is a n by nlevels matrix
#' where each row corresponds to a data point and each column identifies a cluster number for that data point
#' in the hierarchy.
#' @examples
#' data <- rbind(matrix(rnorm(100, mean = 10, sd = 2), nrow = 50),
#'               matrix(rnorm(100, mean = 0, sd = 1), nrow = 50),
#'               matrix(rnorm(100, mean = -10, sd = 3), nrow = 50)
#'               )
#'               
#' clust <- stats::hclust(dist(data),method='single')
#' mergeToMatrix(clust$merge)
#' 
#' @export
mergeToMatrix <- function( merge ) {
  
  # m is number of points minus one
  # n is number of points
  m <- dim(merge)[1]
  n <- m + 1
  
  # verticeSets record data points in each cluster by
  # decoding input merge
  verticesSets <- array(0,dim=c(m,n))
  # assign clusterId based on verticeSets
  componentSets <- array(0,dim=c(m,n))
  # record which layer each data point belongs to 
  layerSet <- array(0,dim=c(m,1))
  
  # get the coorresponding vertices set for each row in merge
  for (i in 1:m) {
    location <- 1
    for (j in 1:2) {
      getValue <- merge[i,j]
      if (getValue < 0) {
        verticesSets[i,location] <- getValue * -1
        location <- location + 1
      }
      if (getValue > 0) {
        k <- 1
        while (verticesSets[getValue,k] != 0) {
          verticesSets[i,location] <- verticesSets[getValue,k]
          location <- location + 1
          k <- k + 1
        }
      }
    }
  }
  
  # get the coorresponding components for each row in merge
  for (i in m:1) {
    index <- i
    if (i==m) {
      layerSet[i] <- 1
    }
    componentId <- 1
    
    for (j in 1:2) {
      getValue <- merge[i,j]
      if (getValue < 0) {
        index2 <- -1 * getValue
        componentSets[index,index2] <- componentId
        componentId <- componentId + 1
      }
      if (getValue > 0) {
        layerSet[getValue] <- layerSet[i] + 1
        k <- 1
        while (verticesSets[getValue,k] != 0) {
          index2 <- verticesSets[getValue,k]
          componentSets[index,index2] <- componentId
          k <- k + 1
          
        }
        componentId <- componentId + 1
      }
      
    }
  }
  
  branchComponentFamily <- array(0, dim=c(n,max(layerSet)+1))
  branchComponentFamily[,1] <- rep(1,n)
  
  # get all the components for each layer of the dendrogram, and then form the branch component family
  for (i in 1:max(layerSet)) {
    clusterId <- 0
    for (j in 1:m) {
      if (layerSet[j] == i) {
        for (k in 1:n) {
          if (componentSets[j,k] > 0) {
            branchComponentFamily[k,i+1] <- clusterId + componentSets[j,k]
          }
        }
        clusterId <- clusterId + 2
      }
    }
  }
  branchComponentFamily
}

#' An ensemble method for combining multiple clustering outcomes based
#' on monotonic graph families of Zhou and Oldford
#' @param clustering1 result of some clustering, for example output from hclust(). A clustering
#' can also be an n by m  matrix, where n is the number of data points and m is the number
#' of levels in the clustering hierarchy. 
#' @param clustering2 result of a second clustering, to be combined with the first.
#' @param ... results of other clustering methods, to be combined with the first two.
#' @param labels labels of data points in clustering results
#' @return a clusterTree object, which is the final clustering result from combining
#' all input clustering results
#' @examples
#' data <- rbind(matrix(rnorm(100, mean = 10, sd = 2), nrow = 50),
#'               matrix(rnorm(100, mean = 0, sd = 1), nrow = 50),
#'               matrix(rnorm(100, mean = -10, sd = 3), nrow = 50)
#'               )
#' clustering1 <- stats::hclust(dist(data),method='single')
#' clustering2 <- kmeans(data,centers=3)
#' clustering3 <- dbscan::dbscan(data,eps=.1)
#' res <- combineClusterings(clustering1,clustering2,clustering3)
#' 
#' @export
combineClusterings <- function(clustering1, clustering2, 
                               ..., labels = NULL) {
  if (missing(clustering1)) stop("Must provide output of clustering for first argument")
  if (missing(clustering2)) stop("Must provide output of clustering for second argument")
  clusterTrees <- Map(getClusterTree, list(clustering1, clustering2, ...))
  
  # n is number of data points
  n <- nrow(clusterTrees[[1]]$tree)
  # combine clusterTrees
  clustering<-clusterTrees[[1]]$tree
  for(i in 2:length(clusterTrees))
  {
    clustering<-cbind(clustering,clusterTrees[[2]]$tree)  
  }
  clustering[is.na(clustering)] <- 0
  
  clusteringsum <- array(0,dim = c(n,n))
  for(j in 2:n)
  {
    for(k in 1:(j-1))
    {
      clusteringsum[j,k] <- sum((clustering[j,]!=0)&
                                  (clustering[k,]!=0)&
                                  (clustering[j,]==clustering[k,]))
    }
  }
  
  distance <- stats::as.dist(clusteringsum)
  singlelinkage <- stats::hclust(-distance, method = "single")
  merge <- singlelinkage$merge
  height <- singlelinkage$height
  
  m <- dim(merge)[1]
  n <- m + 1
  
  verticesSets <- array(0,dim = c(m,n))
  layerSet <- array(0, dim = c(m,1))
  
  shrink <- array(0, dim = c(m,1))
  for(i in m:1)
  {
    if(merge[i,1]>0 && height[i]==height[merge[i,1]])
    {
      shrink[merge[i,1]]=TRUE
    }
    if(merge[i,2]>0 && height[i]==height[merge[i,2]])
    {
      shrink[merge[i,2]]=TRUE
    }
    if(merge[i,1]<0 && merge[i,2]>0)
    {
      shrink[merge[i,2]]=TRUE
    }
    #if(merge[i,1]>0 && merge[i,2]<0)
    #{
    #  shrink[merge[i,1]]=TRUE
    #}
  }
  
  for (i in 1:m) {
    location <- 1
    for (j in 1:2) {
      getValue <- merge[i,j]
      if (getValue < 0) {
        verticesSets[i,location] <- getValue * -1
        location <- location + 1
      }
      if (getValue > 0) {
        k <- 1
        while (verticesSets[getValue,k] != 0) {
          verticesSets[i,location] <- verticesSets[getValue,k]
          location <- location + 1
          k <- k + 1
        }
      }
    }
  }
  
  for (i in m:1) {
    index <- i
    if (i==m) {
      layerSet[i] <- 1
    }
    
    for (j in 1:2) {
      getValue <- merge[i,j]
      if (getValue < 0) {
      }
      if (getValue > 0) {
        if(shrink[getValue])
        {
          layerSet[getValue] <- layerSet[i]
        }
        else
        {
          layerSet[getValue] <- layerSet[i] + 1
        }
      }
    }
  }
  
  
  branchComponentFamily <- array(0, dim=c(n,max(layerSet)))
  for (i in 1:max(layerSet)) {
    clusterId <- 1
    for (j in 1:m) {
      if (layerSet[j] == i && shrink[j]==FALSE) {
        for(k in 1:n)
        {
          if(verticesSets[j,k]>0)
          {
            branchComponentFamily[verticesSets[j,k],i] <- clusterId
          }
        }
        clusterId <- clusterId + 1
      }
    }
  }
  tree <- matrixToClusterTree(branchComponentFamily[,2:max(layerSet)],labels = labels)
  tree
}


#' reorder rows of tree attribute of clusterTree object
#' @param x is the tree attribute of clusterTree object
#' @param labels labels is the order of rows of x
#' @return an order which simplies the process of plotting dendogram/density plot
reOrderClusterTreeMatrix <- function(x,labels=NULL)
{
  if(!is.matrix(x))
    stop('input must be a matrix!')
  n <- dim(x)[1]
  if(is.null(labels)){ labels <- c(1:n) }
  if(!is.vector(labels)){ stop("labels are not a vector!") }
  if(!is.matrix(x)){ stop("input is not a matrix!") }
  if(dim(x)[2]==1){
    IDs <- unique(x)
    IDs <- IDs[!is.na(IDs)]
    res <- c()
    for (id in IDs) {
      tmp <- (x[,1]==id) & (!is.na(x[,1]))
      res<-c(res,labels[tmp])
    }
    res <- c(res,labels[is.na(x[,1])])
    res
  }
  else{
    IDs <- unique(x[,1])
    IDs <- IDs[!is.na(IDs)]
    res <- c()
    for (id in IDs) {
      tmp <- (x[,1]==id) & (!is.na(x[,1]))
      tmpx <- as.matrix(x[tmp,][,2:dim(x)[2]])
      tmplabels <- labels[tmp]
      res <- c(res,reOrderClusterTreeMatrix(tmpx,tmplabels))
    }
    res <- c(res,labels[is.na(x[,1])])
    res
  }
}



