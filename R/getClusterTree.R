#' Transform all data structures into clusterTree
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
    matrixToClusterTree(x)
}

## afCEC package
#' @export
getClusterTree.afCEC <- function(x)
{
    x<-as.matrix(x$labels)
    matrixToClusterTree(x)
}

## apcluster package
#' @export
getClusterTree.apcluster <- function(x)
{
    x<-as.matrix(x)
    matrixToClusterTree(x)
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
    matrixToClusterTree(x)
}

## package cclust
#' @export
getClusterTree.cclust <- function(x)
{
    x <- as.matrix(x$cluster)
    matrixToClusterTree(x)
}

##  package CEC
#' @export
getClusterTree.cec <- function(x)
{
    x<-as.matrix(x$cluster)
    matrixToClusterTree(x)
}

####  package Ckmeans.1d.dp  #### TO DO
#' @export
getClusterTree.Ckmeans.1d.dp <- function(x)
{
    x<-as.matrix(x$cluster)
    matrixToClusterTree(x)
}

##  package clues
#' @export
getClusterTree.clues <- function(x)
{
    x<-as.matrix(x$mem)
    matrixToClusterTree(x)
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
    matrixToClusterTree(x)
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
    matrixToClusterTree(x)
}

#' @export
getClusterTree.pam <- function(x)
{
    x<-as.matrix(x$clustering)
    matrixToClusterTree(x)
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
  matrixToClusterTree(as.matrix(x$cluster))
}

getClusterTree.dbscan <- function(x)
{
  matrixToClusterTree(as.matrix(x$cluster))
}

## densityClust package ##
#' @export
getClusterTree.densityCluster <- function(x){
  matrixToClusterTree(x$clusters)
}

## Fclust package  ##
#' @export
getClusterTree.fclust <- function(x)
{
    matrixToClusterTree(as.matrix(x$clus))
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
    matrixToClusterTree(as.matrix(x@.Data))
}

## kmeans package
#' @export
getClusterTree.kmeans <- function(x)
{
    matrixToClusterTree(as.matrix(x$cluster))
}

##  mixture based clustering
#' @export
getClusterTree.Mclust <- function(x)
{
    matrixToClusterTree(as.matrix(x$classification))
}

