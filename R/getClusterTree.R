#' Transform all data structures into BranchComponent
#'
#' @param x data structure output by some clustering method 
#' (e.g. hclust, kmeans, dbscan, etc.)
#' @return a matrix providing the mapping 
#' between data points and cluster id.
#' @examples
#' x <- kmeans(matrix(rnorm(100),nrow=50),centers=3)
#' getClusterTree(x)
#' @export
getClusterTree <- function(x)
{
    UseMethod("getClusterTree")
}

#' @export
getClusterTree.clusterTree <- function(x) {
    class(x) <- unique(c("clusterTree", class(x)))
    x
}

#' @export
getClusterTree.default <- function(x)
{
  if(is.vector(x)) {
    getClusterTree(matrix(x, ncol = 1))
  } else stop("Not a clustering")
}

#' @export
getClusterTree.factor <- function(x)
{
  if(is.numeric(x)) {
    getClusterTree(matrix(x, ncol = 1))
  } else stop("Not a clustering")
}

#' @export
getClusterTree.matrix <- function(x)
{
    matrixToClusterTree(x) 
}

## adjclust package
#' @export
getClusterTree.chac <- function(x)
{
    matrixToClusterTree(mergeToMatrix(x$merge),labels = x$labels)
}

## adpclust package
#' @export
getClusterTree.adpclust <- function(x)
{
    x <- as.matrix(x$clusters)
    tree <- matrixToClusterTree(x)
    tree
}

## afCEC package
#' @export
getClusterTree.afCEC <- function(x)
{
    x<-as.matrix(x$labels)
    tree <- matrixToClusterTree(x)
    tree
}

## apcluster package
#' @export
getClusterTree.apcluster <- function(x)
{
    x<-as.matrix(x)
    tree <- matrixToClusterTree(x)
    tree
}
#' @export
getClusterTree.AggExResul <- function(x)
{
    matrixToClusterTree(mergeToMatrix(x$merge),x$labels)
}

## bclust package
#' @export
getClusterTree.bclust <- function(x)
{
    matrixToClusterTree(mergeToMatrix(x$merge),x$labels)
}

#' @export
getClusterTree.bclustvs <- function(x)
{
    matrixToClusterTree(mergeToMatrix(x$merge),x$labels)
}

## biclust package  #Nothing needed here

## package cba
#' @export
getClusterTree.ccfkms <- function(x)
{
    x<-as.matrix(x$cl)
    tree <- matrixToClusterTree(x)
    tree
}

## package cclust
#' @export
getClusterTree.cclust <- function(x)
{
    x <- as.matrix(x$cluster)
    tree <- matrixToClusterTree(x)
    tree
}

##  package CEC
#' @export
getClusterTree.cec <- function(x)
{
    x<-as.matrix(x$cluster)
    tree <- matrixToClusterTree(x)
    tree
}

####  package Ckmeans.1d.dp  #### TO DO
#' @export
getClusterTree.Ckmeans.1d.dp <- function(x)
{
    x<-as.matrix(x$cluster)
    tree <- matrixToClusterTree(x)
    tree
}

##  package clues
#' @export
getClusterTree.clues <- function(x)
{
    x<-as.matrix(x$mem)
    tree <- matrixToClusterTree(x)
    tree
}

## package cluster
#' @export
getClusterTree.agnes <- function(x)
{
   matrixToClusterTree(mergeToMatrix(x$merge),x$labels)
}

#' @export
getClusterTree.clara <- function(x)
{
    x<-as.matrix(x$clustering)
    tree <- matrixToClusterTree(x)
    tree
}

#' @export
getClusterTree.diana <- function(x)
{
    matrixToClusterTree(mergeToMatrix(x$merge),x$labels)
}

#' @export
getClusterTree.fanny <- function(x)
{
    x<-as.matrix(x$clustering)
    tree <- matrixToClusterTree(x)
    tree
}

#' @export
getClusterTree.pam <- function(x)
{
    x<-as.matrix(x$clustering)
    tree <- matrixToClusterTree(x)
    tree
}

## ClusterR package  TODO
 
## clustMD package   TODO

## CoClust package  TODO
#' @export
getClusterTree.CoClust <- function(x)
{
    stop("don't know how to handle result from CoClust package")
}

##  compHclust package
#' @export
#getClusterTree.compHclust <- function(x)
#{
#}

## conclust package
## only need default processing

## contaminatedmixt package
#' @export
#getClusterTree.CNmixt <- function(x)
#{
#don't know how to do this!!!
#}

## CORM package  ###
## don't know how to do this!
## don't need to deal with CREAM package ##

## package CrossClustering returns list??????????????????????

## Package 'CRPClustering' only requires default processing ##
#' @export
getClusterTree.dbscan_fast <- function(x)
{
    tmp <- x$cluster
    cl<-max(tmp)+1
    for (ii in 1:length(tmp)) {
        if(tmp[ii]==0){
            tmp[ii]<-cl
            cl <- cl+1
        }
    }
    x<-as.matrix(tmp)
    tree <- matrixToClusterTree(x)
    tree
}


## densityClust package ##
#' @export
getClusterTree.densityCluster <- function(x){
  tree <- matrixToClusterTree(x$clusters)
}

## Fclust package  ##
#' @export
getClusterTree.fclust <- function(x)
{
    x <- as.matrix(x$clus)
    tree <- matrixToClusterTree(x)
    tree
}

##  hierarchical based clustering
#' @export
getClusterTree.hclust <- function(x)
{
    matrixToClusterTree(mergeToMatrix(x$merge),x$labels)
}

## package kernlab
#' @export
getClusterTree.specc <- function(x)
{
    x<-as.matrix(x@.Data)
    tree <- matrixToClusterTree(x)
    tree
}

## kmeans package
#' @export
getClusterTree.kmeans <- function(x)
{
    x<-as.matrix(x$cluster)
    tree <- matrixToClusterTree(x)
    tree
}

##  mixture based clustering
#' @export
getClusterTree.Mclust <- function(x)
{
    x<-as.matrix(x$classification)
    tree <- matrixToClusterTree(x)
    tree
}

