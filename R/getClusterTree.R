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
getClusterTree.default<-function(x)
{
    stop("Not a clustering")
}

## adjclust package
#' @export
getClusterTree.chac <- function(x)
{
    mergeToTree(x$merge)
}

## adpclust package
#' @export
getClusterTree.adpclust <- function(x)
{
    x <- as.matrix(x$clusters)
    tree <- list(tree=x)
    class(tree)<-c("clusterTree",class(tree))
    tree
}

## afCEC package
#' @export
getClusterTree.afCEC <- function(x)
{
    x<-as.matrix(x$labels)
    tree<-list(tree=x)
    class(tree)<-c("clusterTree",class(tree))
    tree
}

## apcluster package
#' @export
getClusterTree.apcluster <- function(x)
{
    x<-as.matrix(x)
    tree<-list(tree=x)
    class(tree)<-c("clusterTree",class(tree))
    tree
}
#' @export
getClusterTree.AggExResul <- function(x)
{
    mergeToTree(x$merge)
}

## bclust package
#' @export
getClusterTree.bclust <- function(x)
{
    mergeToTree(x$merge)
}

#' @export
getClusterTree.bclustvs <- function(x)
{
    mergeToTree(x$merge)
}

## biclust package  #Nothing needed here

## package cba
#' @export
getClusterTree.ccfkms <- function(x)
{
    x<-as.matrix(x$cl)
    tree<-list(tree=x)
    class(tree)<-c("clusterTree",class(tree))
    tree
}

## package cclust
#' @export
getClusterTree.cclust <- function(x)
{
    x$cluster
}

##  package CEC
#' @export
getClusterTree.cec <- function(x)
{
    x<-as.matrix(x$cluster)
    tree<-list(tree=x)
    class(tree)<-c("clusterTree",class(tree))
    tree
}

####  package Ckmeans.1d.dp  #### TO DO
#' @export
getClusterTree.Ckmeans.1d.dp <- function(x)
{
    x<-as.matrix(x$cluster)
    tree<-list(tree=x)
    class(tree)<-c("clusterTree",class(tree))
    tree
}

##  package clues
#' @export
getClusterTree.clues <- function(x)
{
    x<-as.matrix(x$mem)
    tree<-list(tree=x)
    class(tree)<-c("clusterTree",class(tree))
    tree
}

## package cluster
#' @export
getClusterTree.agnes <- function(x)
{
    mergeToTree(x$merge)
}

#' @export
getClusterTree.clara <- function(x)
{
    x<-as.matrix(x$clustering)
    tree<-list(tree=x)
    class(tree)<-c("clusterTree",class(tree))
    tree
}

#' @export
getClusterTree.diana <- function(x)
{
    mergeToTree(x$merge)
}

#' @export
getClusterTree.fanny <- function(x)
{
    x<-as.matrix(x$clustering)
    tree<-list(tree=x)
    class(tree)<-c("clusterTree",class(tree))
    tree
}

#' @export
getClusterTree.pam <- function(x)
{
    x<-as.matrix(x$clustering)
    tree<-list(tree=x)
    class(tree)<-c("clusterTree",class(tree))
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
    tree<-list(tree=x)
    class(tree)<-c("clusterTree",class(tree))
    tree
}

## Fclust package  ##
#' @export
getClusterTree.fclust <- function(x)
{
    x<-as.matrix(x$clus)
    tree<-list(tree=x)
    class(tree)<-c("clusterTree",class(tree))
    tree
}

##  hierarchical based clustering
#' @export
getClusterTree.hclust <- function(x)
{
    mergeToTree(x$merge)
}

## package kernlab
#' @export
getClusterTree.specc <- function(x)
{
    x<-as.matrix(x@.Data)
    tree<-list(tree=x)
    class(tree)<-c("clusterTree",class(tree))
    tree
}

## kmeans package
#' @export
getClusterTree.kmeans <- function(x)
{
    x<-as.matrix(x$cluster)
    tree<-list(tree=x)
    class(tree)<-c("clusterTree",class(tree))
    tree
}

##  mixture based clustering
#' @export
getClusterTree.Mclust <- function(x)
{
    x<-as.matrix(x$classification)
    tree<-list(tree=x)
    class(tree)<-c("clusterTree",class(tree))
    tree
}

