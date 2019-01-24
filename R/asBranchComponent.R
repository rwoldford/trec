#' Transform all data structures into BranchComponent
#'
#' @param x data structure output by some clustering method 
#' (e.g. hclust, kmeans, dbscan, etc.)
#' @return a matrix providing the mapping 
#' between data points and cluster id.
#' @examples
#' x <- kmeans(matrix(rnorm(100),nrow=50),centers=3)
#' asBranchComponent(x)
#' @export
asBranchComponent <- function(x)
{
    UseMethod("asBranchComponent")
}

#' @export
asBranchComponent.default<-function(x)
{
    as.matrix(x)
}

## adjclust package
#' @export
asBranchComponent.chac <- function(x)
{
    getComponentsfromMerge(x$merge)
}

## adpclust package
#' @export
asBranchComponent.adpclust <- function(x)
{
    x$clusters
}

## afCEC package
#' @export
asBranchComponent.afCEC <- function(x)
{
    x$labels
}

## apcluster package
#' @export
asBranchComponent.apcluster <- function(x)
{
    x
}
#' @export
asBranchComponent.AggExResul <- function(x)
{
    getComponentsfromMerge(x$merge)
}

## bclust package
#' @export
asBranchComponent.bclust <- function(x)
{
    getComponentsfromMerge(x$merge)
}

#' @export
asBranchComponent.bclustvs <- function(x)
{
    getComponentsfromMerge(x$merge)
}

## biclust package  #Nothing needed here

## package cba
#' @export
asBranchComponent.ccfkms <- function(x)
{
    x$cl
}

## package cclust
#' @export
asBranchComponent.cclust <- function(x)
{
    x$cluster
}

##  package CEC
#' @export
asBranchComponent.cec <- function(x)
{
    x$cluster
}

####  package Ckmeans.1d.dp  #### TO DO
#' @export
asBranchComponent.Ckmeans.1d.dp <- function(x)
{
    x$cluster
}

##  package clues
#' @export
asBranchComponent.clues <- function(x)
{
    x$mem
}

## package cluster
#' @export
asBranchComponent.agnes <- function(x)
{
    getComponentsfromMerge(x$merge)
}

#' @export
asBranchComponent.clara <- function(x)
{
    x$clustering
}

#' @export
asBranchComponent.diana <- function(x)
{
    getComponentsfromMerge(x$merge)
}

#' @export
asBranchComponent.fanny <- function(x)
{
    x$clustering
}

#' @export
asBranchComponent.pam <- function(x)
{
    x$clustering
}

## ClusterR package  TODO
 
## clustMD package   TODO

## CoClust package  TODO
#' @export
asBranchComponent.CoClust <- function(x)
{
    stop("don't know how to handle result from CoClust package")
}

##  compHclust package
#' @export
#asBranchComponent.compHclust <- function(x)
#{
#}

## conclust package
## only need default processing

## contaminatedmixt package
#' @export
#asBranchComponent.CNmixt <- function(x)
#{
#don't know how to do this!!!
#}

## CORM package  ###
## don't know how to do this!
## don't need to deal with CREAM package ##

## package CrossClustering returns list??????????????????????

## Package 'CRPClustering' only requires default processing ##
#' @export
asBranchComponent.dbscan_fast <- function(x)
{
    tmp <- x$cluster
    cl<-max(tmp)+1
    for (ii in 1:length(tmp)) {
        if(tmp[ii]==0){
            tmp[ii]<-cl
            cl <- cl+1
        }
    }
    tmp
}

## Fclust package  ##
#' @export
asBranchComponent.fclust <- function(x)
{
    x$clus
}

##  hierarchical based clustering
#' @export
asBranchComponent.hclust <- function(x)
{
    getComponentsfromMerge(x$merge)
}

## package kernlab
#' @export
asBranchComponent.specc <- function(x)
{
    x@.Data
}

## kmeans package
#' @export
asBranchComponent.kmeans <- function(x)
{
    x$cluster
}

##  mixture based clustering
#' @export
asBranchComponent.Mclust <- function(x)
{
    x$classification
}

