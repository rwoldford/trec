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

asBranchComponent.default<-function(x)
{
    as.matrix(x)
}

## adjclust package
asBranchComponent.chac <- function(x)
{
    getComponentsfromMerge(x$merge)
}

## adpclust package
asBranchComponent.adpclust <- function(x)
{
    x$clusters
}

## afCEC package
asBranchComponent.afCEC <- function(x)
{
    x$labels
}

## apcluster package
asBranchComponent.apcluster <- function(x)
{
    x
}
asBranchComponent.AggExResul <- function(x)
{
    getComponentsfromMerge(x$merge)
}

## bclust package
asBranchComponent.bclust <- function(x)
{
    getComponentsfromMerge(x$merge)
}

asBranchComponent.bclustvs <- function(x)
{
    getComponentsfromMerge(x$merge)
}
## bclust package

## biclust package

## package cba
asBranchComponent.ccfkms <- function(x)
{
    x$cl
}

## package cclust
asBranchComponent.cclust <- function(x)
{
    x$cluster
}

##  package CEC
asBranchComponent.cec <- function(x)
{
    x$cluster
}

####  package Ckmeans.1d.dp  ####???????????????????????????????????????????????????????????
asBranchComponent.Ckmeans.1d.dp <- function(x)
{
    x$cluster
}

##  package clues
asBranchComponent.clues <- function(x)
{
    x$mem
}

## package cluster
asBranchComponent.agnes <- function(x)
{
    getComponentsfromMerge(x$merge)
}

asBranchComponent.clara <- function(x)
{
    x$clustering
}

asBranchComponent.diana <- function(x)
{
    getComponentsfromMerge(x$merge)
}

asBranchComponent.fanny <- function(x)
{
    x$clustering
}

asBranchComponent.pam <- function(x)
{
    x$clustering
}

## ClusterR package

## clustMD package

## CoClust package
asBranchComponent.CoClust <- function(x)
{
    # don't know how to do this!!!!
}

##  compHclust package
#asBranchComponent.compHclust <- function(x)
#{
#}

## conclust package
## only need default processing

## contaminatedmixt package
#asBranchComponent.CNmixt <- function(x)
#{
#don't know how to do this!!!
#}

## CORM package  ###
## don't know how to do this!
## don't need to deal with CREAM package ##

## package CrossClustering returns list??????????????????????

## Package 'CRPClustering' only requires default processing ##
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
asBranchComponent.fclust <- function(x)
{
    x$clus
}

##  hierarchical based clustering
asBranchComponent.hclust <- function(x)
{
    getComponentsfromMerge(x$merge)
}

## package kernlab
asBranchComponent.specc <- function(x)
{
    x@.Data
}

## kmeans package
asBranchComponent.kmeans <- function(x)
{
    x$cluster
}

##  mixture based clustering
asBranchComponent.Mclust <- function(x)
{
    x$classification
}

