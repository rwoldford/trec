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
  tree <- list(treeMatrix = x, labels = labels)
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
  layerHeight <- array(0,dim=c(m,1))
  
  # get the corresponding vertices set for each row in merge
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
      layerHeight[i] <- 1
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
        layerHeight[getValue] <- layerHeight[i] + 1
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
  
  branchComponentFamily <- array(0, dim=c(n,max(layerHeight)+1))
  branchComponentFamily[,1] <- rep(1,n)
  
  # get all the components for each layer of the dendrogram, and then form the branch component family
  for (i in 1:max(layerHeight)) {
    clusterId <- 0
    for (j in 1:m) {
      if (layerHeight[j] == i) {
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
#' @param pruneNumber set number for pruning trivial components
#' @return a clusterTree object, which is the final clustering result from combining
#' all input clustering results
#' @examples
#' data <- rbind(matrix(rnorm(100, mean = 10, sd = 2), nrow = 50),
#'               matrix(rnorm(100, mean = 0, sd = 1), nrow = 50),
#'               matrix(rnorm(100, mean = -10, sd = 3), nrow = 50)
#'               )
#' clustering1 <- stats::hclust(dist(data),method='single')
#' clustering2 <- kmeans(data,centers=3)
#' clustering3 <- dbscan::dbscan(data,eps=.8)
#' res <- combineClusterings(clustering1,clustering2,clustering3)
#' 
#' @export
combineClusterings <- function(clustering1, clustering2, 
                               ..., labels = NULL, weights = NULL, pruneNumber = 1) {
  if (missing(clustering1)) stop("Must provide output of clustering for first argument")
  if (missing(clustering2)){
    if (!methods::is(clustering1, "clusterings")) {
      stop(paste("Either clustering1 is a list of clusterings of class `clusterings`",
                 "or at least one more clustering must be provided."))
    } else {
      clusterTrees <- Map(getClusterTree, clustering1)
    }
  } else {
    # map clusterings to cluster trees
    clusterTrees <- Map(getClusterTree, list(clustering1, clustering2, ...))
  }
  if(is.null(weights)){
    # initialize weights to all one vector if it is null
    weights <- rep(1,length(clusterTrees))
  }
  if(length(weights)!=length(clusterTrees)){
    stop('length of weights and length of clusterings do not match!')
  }
  # n is number of data points
  n <- nrow(clusterTrees[[1]]$treeMatrix)
  extendedWeights <- c()
  # combine clusterTrees
  clustering <- clusterTrees[[1]]$treeMatrix
  extendedWeights <- c(extendedWeights,rep(weights[1],ncol(clusterTrees[[1]]$treeMatrix)))
  for(i in 2:length(clusterTrees))
  {
    # add each treeMatrix into final matrix
    clustering<-cbind(clustering,clusterTrees[[i]]$treeMatrix)
    # extend weights because hierarchical clustering's treeMatrix may have more than one column/layer
    # we need to assign weights to each layer of hierarchical clustering
    extendedWeights <- c(extendedWeights,ncol(clusterTrees[[i]]$treeMatrix))
  }
  # change NA to zeros
  clustering[is.na(clustering)] <- 0
  if(sum(clustering!=0) == 0)
    stop('at least one entry of clustering result must be nonzero')
  clusteringsum <- array(0,dim = c(n,n))
  for(j in 2:n)
  {
    for(k in 1:(j-1))
    {
      # use extendedWeights to dot product 
      clusteringsum[j,k] <- sum((extendedWeights)*((clustering[j,]!=0)&
                                  (clustering[k,]!=0)&
                                  (clustering[j,]==clustering[k,])))
    }
  }
  distance <- -(stats::as.dist(clusteringsum))  # clustering sum is a set of similarities
  # apply single linkage on negative similarities
  singlelinkage <- stats::hclust(distance, method = "single")
  merge <- singlelinkage$merge
  height <- singlelinkage$height
  
  m <- dim(merge)[1]
  n <- m + 1
  
  #count number of data points for each row in merge
  numberOfPointsInMerge <- rep(0,m)
  for (i in 1:m) {
    currentNumber <- 0
    if(merge[i,1]<0){
      # less than zero means size of cluster is one
      currentNumber <- currentNumber + 1
    }
    else{
      # greater than zero means size of cluster can be retrieved recursively
      currentNumber <- currentNumber + numberOfPointsInMerge[merge[i,1]]
    }
    
    if(merge[i,2]<0){
      # less than zero means size of cluster is one
      currentNumber <- currentNumber + 1
    }
    else{
      # greater than zero means size of cluster can be retrieved recursively
      currentNumber <- currentNumber + numberOfPointsInMerge[merge[i,2]]
    }
    # assign number of points for each row in merge
    numberOfPointsInMerge[i] <- currentNumber
  }
  
  # Final cluster tree will be a collection of layers
  # verticeSets record all vertices ID for each row in merge
  verticesSets <- array(0,dim = c(m,n))
  # layerHeight record layer number for each row in merge
  layerHeight <- array(0, dim = c(m,1))
  
  # decide whether to shrink for each row in merge
  # if shrink is true, then we do not show 
  shrink <- array(0, dim = c(m,1))
  for(i in m:1)
  {
    if(merge[i,1]>0 && height[i]==height[merge[i,1]])
      # Left branch: Split distance is identical ... more than binary split 
    {
      shrink[merge[i,1]]=TRUE
    }
    if(merge[i,2]>0 && height[i]==height[merge[i,2]])
      # Right branch: Split distance is identical ... more than binary split 
    {
      shrink[merge[i,2]]=TRUE
    }
    if(merge[i,1]<0 && merge[i,2]>0)
      # Trivial pruning
    {
      # if merge[i,1] is trivial components, merge[i,2] is unnecessary
      shrink[merge[i,2]]=TRUE
    }
    if(merge[i,1]>0)
    {
      if(numberOfPointsInMerge[merge[i,1]]<=pruneNumber)
      # if number of points is less than pruneNumber 
      {
        # shrink this row in merge
        shrink[merge[i,1]]=TRUE
        if(merge[i,2]>0){
          # shrink this row in merge ... telescopic contraction
          shrink[merge[i,2]]=TRUE
        }
      }
    }
    if(merge[i,2]>0)
    {
      if(numberOfPointsInMerge[merge[i,2]]<=pruneNumber)
      # if number of points is less than pruneNumber
      {
        shrink[merge[i,2]]=TRUE
        if(merge[i,1]>0){
          # shrink this row in merge ... telescopic contraction
          shrink[merge[i,1]]=TRUE
        }
      }
    }
  }
  
  # construct verticeSets for each row
  # note verticeSets represent vertices correspond to each row in merge
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
  
  # record layer height for each row in merge
  for (i in m:1) {
    if (i==m) {
      layerHeight[i] <- 1
    }
    
    for (j in 1:2) {
      getValue <- merge[i,j]
      if (getValue > 0) { # have a branch
        if(shrink[getValue])
        {
          layerHeight[getValue] <- layerHeight[i]
        }
        else
        {
          layerHeight[getValue] <- layerHeight[i] + 1
        }
      }
    }
  }
  
  # branchComponentFamily is actually a treeMatrix
  # put verticeSets in branchComponentFamily if shrink if false for each row
  branchComponentFamily <- array(0, dim=c(n,max(layerHeight)))
  for (i in 1:max(layerHeight)) {
    clusterId <- 1
    for (j in 1:m) {
      if (layerHeight[j] == i && shrink[j]==FALSE) {  # effects the telescoping
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
  if(ncol(branchComponentFamily)>1){
    branchComponentFamily <- branchComponentFamily[,-1]
  }
  tree <- matrixToClusterTree(branchComponentFamily,labels = labels)
  tree
}


#' reorder rows of treeMatrix attribute of clusterTree so that data points 
#' that belong to the same cluster are always next to each other. This will
#' help us process the plot function
#' @param treeMatrix is the treeMatrix attribute of clusterTree object
#' @param order order is the inital order of rows of x, default is c(1:n)
#' @return an order which data points belong to same cluster are
#' always next to each other
#' @examples
#' data <- iris[,1:4]
#' clustering1 <- stats::hclust(dist(data),method='single')
#' clustering2 <- kmeans(data,centers=3)
#' clustering3 <- dbscan::dbscan(data,eps=.78)
#' res <- combineClusterings(clustering1,clustering2,clustering3) 
#' # obtain the order for reordering cluster tree
#' order <- reorderClusterTreeMatrix(res$treeMatrix)
#' # After reordering, data points with same cluster id are always next
#' # to each other
#' reordered_treeMatrix <- res$treeMatrix[order,]
#' # check reordered_treeMatrix
#' reordered_treeMatrix
#' @export
reOrderClusterTreeMatrix <- function(treeMatrix,order=NULL)
{
  # reordering is implemented recursively
  # for each function call, the first column of treeMatrix is reordered correctly
  # then a recursive call is applied to treeMatrix except for the first column to 
  # reorder the rest of the matrix 
  if(!is.matrix(treeMatrix))
    stop('input must be a matrix!')
  n <- nrow(treeMatrix)
  if(is.null(order)){ order <- c(1:n) }
  if(!is.vector(order)){ stop("order are not a vector!") }
  if(!is.matrix(treeMatrix)){ stop("input is not a matrix!") }
  if(ncol(treeMatrix)==1){
    # check all unique ids in treeMatrix
    IDs <- unique(treeMatrix)
    # find all unique ids that is not NA
    IDs <- IDs[!is.na(IDs)]
    groupIndexWithSameId(treeMatrix = treeMatrix, IDs = IDs,initialOrder = order)
  }
  else{
    # find all IDs in first column
    IDs <- unique(treeMatrix[,1])
    # find all non-NA IDs
    IDs <- IDs[!is.na(IDs)]
    res <- c()
    for (id in IDs) {
      # find Indexes that belong to same cluster
      IndexesWithSameIds <- (treeMatrix[,1]==id) & (!is.na(treeMatrix[,1]))
      # form treeMatrix and order for recursive calls
      recursiveTreeMatrix <- as.matrix(treeMatrix[IndexesWithSameIds,-1])
      recursiveOrder <- order[IndexesWithSameIds]
      # call this function recursively to get order for next level
      # add recursive answers into the result
      res <- c(res,reOrderClusterTreeMatrix(recursiveTreeMatrix,recursiveOrder))
    }
    # add indexes that corresponds to NAs in the end
    res <- c(res,order[is.na(treeMatrix[,1])])
    res
  }
}

#' group index with same id in first column of treeMatrix
#' @param treeMatrix treeMatrix attribute of cluster tree
#' @param IDs cluster id in treeMatrix
#' @param initialOrder  initialOrder is the initial order of treeMatrix
groupIndexWithSameId <- function(treeMatrix,IDs, initialOrder){
  # group index with same id in first column of treeMatrix
  # together and form a vector
  indexes <- c()
  for (id in IDs) {
    # group all indexes which corresponds to the current id together
    # note NA values will be processed at the end to avoid error
    indexesGroup <- (treeMatrix[,1]==id) & (!is.na(treeMatrix[,1]))
    indexes <- c(indexes, initialOrder[indexesGroup])
  }
  indexes <- c(indexes, initialOrder[is.na(treeMatrix[,1])])
  indexes
}

