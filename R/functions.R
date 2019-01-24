#' Get the branch component family from the result of (single/average/complete) linkage clustering.
#' @param merge The n-1 by 2 matrix that is the merge component output from
#' a hierarchical clustering which describes how the n points are merged in the cluster
#' hierarchy.
#' @return An n by nlevels matrix where each row corresponds to a data point and each
#' column identifies a cluster number for that data point in the hierarchy.
#' @examples
#' data <- rbind(matrix(rnorm(100, mean = 10, sd = 2), nrow = 50),
#'               matrix(rnorm(100, mean = 0, sd = 1), nrow = 50),
#'               matrix(rnorm(100, mean = -10, sd = 3), nrow = 50)
#'               )
#'               
#' cl <- stats::hclust(dist(data),method='single')
#' getComponentsfromMerge(cl$merge)
#' 
#' @export
#' 
getComponentsfromMerge <- function( merge ) {
    
    # m is number of points
    m <- dim(merge)[1]
    n <- m + 1
    
    verticesSets <- array(0,dim=c(m,n))
    componentSets <- array(0,dim=c(m,n))
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

#' This is a new ensemble method for combining multiple clustering outcomes
#' @param clustering a matrix which stands for multiple clustering outcomes each row represent a data point, each column is an assignment of cluster id
#' @return final one hierarchical clustering
#' @examples
#' data <- rbind(matrix(rnorm(100, mean = 10, sd = 2), nrow = 50),
#'               matrix(rnorm(100, mean = 0, sd = 1), nrow = 50),
#'               matrix(rnorm(100, mean = -10, sd = 3), nrow = 50)
#'               )
#' cl1<-asBranchComponent(stats::hclust(dist(data),method='single'))
#' cl2<-asBranchComponent(kmeans(data,centers=3))
#' TRECgetComponentsfromClustering(cbind(cl1,cl2))
#' @export
#' 
TRECgetComponentsfromClustering <- function( clustering ) {
    # n is number of data points
    n <- dim(clustering)[1]
    # clustering sum of multiple clustering outcomes
    clusteringsum <- array(0,dim = c(n,n))
    for(j in 1:n)
    {
        for(k in 1:j)
        {
            clusteringsum[j,k] <- sum((clustering[j,]!=0)&
                                          (clustering[k,]!=0)&
                                          (clustering[j,]==clustering[k,]))
        }
    }
    
    # use single linkage method to produce result for trec
    distance <- stats::as.dist(clusteringsum)
    singlelinkage <- stats::hclust(-distance, method = "single")
    merge <- singlelinkage$merge
    height <- singlelinkage$height
    
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
    
    branchComponentFamily[,2:max(layerSet)]
}


#' Calculate distance of two hierarchical clustering matrix
#' @param branchComponent1 The first clustering as a branch component
#' @param branchComponent2 The second clustering as a branch component
#' @return distance between 0 and square root of 2
#' @examples
#' data <- rbind(matrix(rnorm(100, mean = 10, sd = 2), nrow = 50),
#'               matrix(rnorm(100, mean = 0, sd = 1), nrow = 50),
#'               matrix(rnorm(100, mean = -10, sd = 3), nrow = 50)
#'               )
#' clu1<-asBranchComponent(stats::hclust(dist(data),method='complete'))
#' clu2<-asBranchComponent(stats::hclust(dist(data),method='single'))
#' getClusteringDistance(clu1,clu2)
#' @export
#' 
getClusteringDistance <- function(branchComponent1, branchComponent2)
{
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
}

########################################################################################################################
# require(igraph)
#
# setClass("componentTree",representation(components = "matrix",content = "matrix", mapping = "matrix",total="matrix"))
#
# setMethod("show","componentTree",
# 	# show the summary of a component tree
# 	function(object) {
# 		out <- paste("This component tree has", dim(object@components)[2],"layers with all together", max(object@content[1,])+1,"nodes.")
# 		cat(out)
# 		cat("\n\n")
#
# 		cat("The available slots:\n")
# 		theSlots <- c("components","content","mapping","total")
# 		show(theSlots)
# 		cat("\n")
#
# 		cat("The available methods:\n")
# 		theMethods <- c("nodes","plot","mapping","average")
# 		show(theMethods)
# 	}
# )
#
# setMethod("plot",signature(x="componentTree",y="missing"),
# 	# plot the component tree
# 	function(x,y,...) {
# 		nv <- max(x@content[1,]+1)
# 		if (dim(x@mapping)[2] == 2) {
# 			adjacency_matrix <- array(0, dim=c(nv,nv))
# 			for (i in 1:dim(x@mapping)[1]){
# 				adjacency_matrix[x@mapping[i,1]+1,x@mapping[i,2]+1] <- 1
# 			}
# 		}
# 		else {
# 			adjacency_matrix <- diag(0,nv)
# 		}
# 		g <- graph.adjacency(adjacency_matrix)
# 		if (dim(x@mapping)[2] == 2) {
# 			igraph.par("plot.layout", layout.reingold.tilford)
# 		}
# 		else {
# 			igraph.par("plot.layout", layout.circle)
# 		}
# 		nodes(x)
# 		plot(g, ...)
# 	}
# )
#
# nodes <- function(cTree) {
# 	# show the content of each node in a component tree
# 	if (data.class(cTree)=="componentTree") {
# 		cat("The list of the node content:\n")
# 		class(cTree) <- "componentTree"
# 		for (i in 0:max(cTree@content[1,])) {
# 			out <- paste("node",i," \t")
# 			cat(out)
# 			out <- paste(cTree@content[2,cTree@content[1,]==i])
# 			cat(out)
# 			cat("\n")
# 		}
# 	}
# 	else {
# 		cat("This method only takes an object from \"componentTree\" as the input.\n")
# 	}
# }
#
# getChildren <- function(node,cTree) {
# 	# get all the list of IDs from the child nodes of the given node
# 	if (data.class(cTree) != "componentTree") {
# 		stop("Check if the argument is valid.")
# 	}
# 	class(cTree) <- "componentTree"
# 	children <- cTree@mapping[cTree@mapping[,1] == node,2]
# 	children
# }
#
# getNodesContent <- function(nodesList,cTree) {
# 	# get the union of content/points of the given node list
# 	if (data.class(cTree) != "componentTree") {
# 		stop("Check if the argument is valid.")
# 	}
# 	class(cTree) <- "componentTree"
#
# 	pointsList <- NULL
#
# 	for (i in nodesList) {
# 		pointsList <- append(pointsList,cTree@content[2,cTree@content[1,] == i])
# 	}
# 	as.integer(levels(factor(pointsList)))
# }
#
# cutTree <- function(cTree,layer) {
# 	# cut the component tree to keep the top layers by the value of 'layer'
# 	if ((data.class(cTree) != "componentTree")||(layer < 2)) {
# 		stop("Check if the argument is valid.")
# 	}
# 	class(cTree) <- "componentTree"
#
# 	if (dim(cTree@components)[2] <= layer) {
# 		cTree
# 	}
# 	else {
# 		branchComponentFamily <- cTree@components[,1:layer]
# 		generateTreebyBranchComponents(branchComponentFamily)
# 	}
#
# }
#
# searchTree <- function(root,cTree,treeInfo,rootLevel) {
# 	# recursively searching the tree to get the similar information of "order","merge" and "height"
# 	# as those in a hclust class, the first component of the output is "order", the second component
# 	# of the output combinds "merge"(the second and third columns) and levels of "height" (the fourth column)
# 	if (data.class(cTree) != "componentTree") {
# 		stop("Check if the argument is valid.")
# 	}
# 	#print(paste("node:",root))
# 	#print(paste("level:",rootLevel))
# 	class(cTree) <- "componentTree"
# 	children <- getChildren(root,cTree)
# 	#print(paste("children:",children))
#
# 	# the current root node has children
# 	if (length(children) > 0) {
# 		#currentLevel <- max(treeInfo$mergeList[,4])
# 		dataInParent <- getNodesContent(root,cTree)
# 		dataInChildren <- getNodesContent(children,cTree)
# 		#print(dataInParent)
# 		#print(dataInChildren)
#
# 		# get the orders for the points which are in the root node but not in its child nodes
# 		for (i in dataInParent) {
# 			if (sum(dataInChildren == i) == 0) {
# 				treeInfo$dataOrder[i] <- max(treeInfo$dataOrder) + 1
# 			}
# 		}
#
# 		# mergeStages is used to get the highest row number in "merge" for each child node
# 		# if that node is not singular, if the node is singular, just get the Id of that single point
# 		mergeStages <- children
#
# 		# reorder the nodes in children, to make singular nodes appear last
# 		temp1 <- array(0,dim=c(length(children),1))
# 		temp2 <- array(0,dim=c(length(children),1))
# 		index1 <- 1
# 		index2 <- 1
# 		for (i in children) {
# 			if (length(getNodesContent(i,cTree)) == 1) {
# 				temp1[index1] <- i
# 				index1 <- index1 + 1
# 			}
# 			else {
# 				temp2[index2] <- i
# 				index2 <- index2 + 1
# 			}
# 		}
# 		if (index2 > 1) {
# 			for (i in 1:(index2-1)) {
# 				children[i] <- temp2[i]
# 			}
# 		}
# 		if (index1 > 1) {
# 			for (i in 1:(index1-1)) {
# 				children[index2-1+i] <- temp1[i]
# 			}
# 		}
#
# 		# get both the orders for each point in every node and the merge and height information
# 		for (i in 1:length(children)) {
# 			children2 <- getChildren(children[i],cTree)
# 			if (length(children2) == 0) { # the node is a leaf node
# 				dataInChild <- getNodesContent(children[i],cTree)
#
# 				# get the orders
# 				for (j in dataInChild) {
# 					treeInfo$dataOrder[j] <- max(treeInfo$dataOrder) + 1
# 				}
# 				# get the merge and height information
# 				if (length(dataInChild) > 1) {
# 					treeInfo$mergeList[max(treeInfo$mergeList[,1])+1,1] <- max(treeInfo$mergeList[,1])+1
# 					treeInfo$mergeList[max(treeInfo$mergeList[,1]),2] <- -1*dataInChild[1]
# 					treeInfo$mergeList[max(treeInfo$mergeList[,1]),3] <- -1*dataInChild[2]
# 					#treeInfo$mergeList[max(treeInfo$mergeList[,1]),4] <- rootLevel + 1
# 					treeInfo$mergeList[max(treeInfo$mergeList[,1]),4] <- rootLevel + 2
# 					for (j in dataInChild) {
# 						if ((j != dataInChild[1]) && (j != dataInChild[2])) {
# 							treeInfo$mergeList[max(treeInfo$mergeList[,1])+1,1] <- max(treeInfo$mergeList[,1])+1
# 							treeInfo$mergeList[max(treeInfo$mergeList[,1]),2] <- -1*j
# 							treeInfo$mergeList[max(treeInfo$mergeList[,1]),3] <- max(treeInfo$mergeList[,1]) - 1
# 							#treeInfo$mergeList[max(treeInfo$mergeList[,1]),4] <- rootLevel + 1
# 							treeInfo$mergeList[max(treeInfo$mergeList[,1]),4] <- rootLevel + 2
# 						}
# 					}
# 					mergeStages[i] <- max(treeInfo$mergeList[,1])
# 				}
# 				else {
# 					mergeStages[i] <- dataInChild
# 				}
# 			}
# 			else { # the node has child nodes
# 				treeInfo <- searchTree(children[i],cTree,treeInfo,rootLevel+1)
# 				mergeStages[i] <- max(treeInfo$mergeList[,1])
# 			}
# 		}
#
# 		# get the information to merge each node together
# 		if ((index2-1) == 1) { # only one node is not singular
# 			mergePoint <- mergeStages[1]
# 		}
# 		if ((index2-1) >= 2) { # at least two nodes are not singular
# 			treeInfo$mergeList[max(treeInfo$mergeList[,1])+1,1] <- max(treeInfo$mergeList[,1])+1
# 			treeInfo$mergeList[max(treeInfo$mergeList[,1]),2] <- mergeStages[1]
# 			treeInfo$mergeList[max(treeInfo$mergeList[,1]),3] <- mergeStages[2]
# 			treeInfo$mergeList[max(treeInfo$mergeList[,1]),4] <- rootLevel + 1
#
# 			mergePoint <- max(treeInfo$mergeList[,1])
# 		}
# 		if ((index2-1) > 2) { # more than two nodes are not singular
# 			for (i in 3:(index2-1)) {
# 				treeInfo$mergeList[max(treeInfo$mergeList[,1])+1,1] <- max(treeInfo$mergeList[,1])+1
# 				treeInfo$mergeList[max(treeInfo$mergeList[,1]),2] <- mergeStages[i]
# 				treeInfo$mergeList[max(treeInfo$mergeList[,1]),3] <- mergePoint
# 				treeInfo$mergeList[max(treeInfo$mergeList[,1]),4] <- rootLevel + 1
#
# 				mergePoint <- max(treeInfo$mergeList[,1])
# 			}
# 		}
#
# 		if (((index1-1) >= 1)&&((index2-1)>=1)) { # at least a node is not singular and there are nodes are singular
# 			for (i in 1:(index1-1)) {
# 				treeInfo$mergeList[max(treeInfo$mergeList[,1])+1,1] <- max(treeInfo$mergeList[,1])+1
# 				treeInfo$mergeList[max(treeInfo$mergeList[,1]),2] <- -1*mergeStages[index2-1+i]
# 				treeInfo$mergeList[max(treeInfo$mergeList[,1]),3] <- mergePoint
# 				treeInfo$mergeList[max(treeInfo$mergeList[,1]),4] <- rootLevel + 1
#
# 				mergePoint <- max(treeInfo$mergeList[,1])
# 			}
# 		}
# 		if (((index1-1) >= 2)&&((index2-1)==0)) { # every node is singular
# 			treeInfo$mergeList[max(treeInfo$mergeList[,1])+1,1] <- max(treeInfo$mergeList[,1])+1
# 			treeInfo$mergeList[max(treeInfo$mergeList[,1]),2] <- -1*mergeStages[1]
# 			treeInfo$mergeList[max(treeInfo$mergeList[,1]),3] <- -1*mergeStages[2]
# 			treeInfo$mergeList[max(treeInfo$mergeList[,1]),4] <- rootLevel + 1
#
# 			for (j in children) {
# 				if ((j != children[1]) && (j != children[2])) {
# 					treeInfo$mergeList[max(treeInfo$mergeList[,1])+1,1] <- max(treeInfo$mergeList[,1])+1
# 					treeInfo$mergeList[max(treeInfo$mergeList[,1]),2] <- -1*j
# 					treeInfo$mergeList[max(treeInfo$mergeList[,1]),3] <- max(treeInfo$mergeList[,1]) - 1
# 					treeInfo$mergeList[max(treeInfo$mergeList[,1]),4] <- rootLevel + 1
# 				}
# 			}
# 		}
#
# 		# get the information of merge and height for every points in the current root node but not in any of its child node
# 		for (i in dataInParent) {
# 			if (sum(dataInChildren == i) == 0) {
# 				treeInfo$mergeList[max(treeInfo$mergeList[,1])+1,1] <- max(treeInfo$mergeList[,1])+1
# 				treeInfo$mergeList[max(treeInfo$mergeList[,1]),2] <- -1*i
# 				treeInfo$mergeList[max(treeInfo$mergeList[,1]),3] <- max(treeInfo$mergeList[,1]) - 1
# 				#treeInfo$mergeList[max(treeInfo$mergeList[,1]),4] <- rootLevel
# 				treeInfo$mergeList[max(treeInfo$mergeList[,1]),4] <- rootLevel + 1
# 			}
# 		}
# 	}
# 	else { # the current root is a leaf node
# 		dataInParent <- getNodesContent(root,cTree)
# 		for (i in dataInParent) {
# 			treeInfo$dataOrder[i] <- max(treeInfo$dataOrder) + 1
# 		}
# 		if (max(treeInfo$mergeList[,1]) == 0) {
# 			if (length(dataInParent) > 1) {
# 				treeInfo$mergeList[1,1] <- 1
# 				treeInfo$mergeList[1,2] <- -1*dataInParent[1]
# 				treeInfo$mergeList[1,3] <- -1*dataInParent[2]
# 				treeInfo$mergeList[1,4] <- max(treeInfo$mergeList[,4]) + 1
# 				#currentLevel <- max(treeInfo$mergeList[,4])
# 				for (i in dataInParent) {
# 					if ((i != dataInParent[1]) && (i != dataInParent[2])) {
# 						treeInfo$mergeList[max(treeInfo$mergeList[,1])+1,1] <- max(treeInfo$mergeList[,1])+1
# 						treeInfo$mergeList[max(treeInfo$mergeList[,1]),2] <- -1*i
# 						treeInfo$mergeList[max(treeInfo$mergeList[,1]),3] <- max(treeInfo$mergeList[,1]) - 1
# 						treeInfo$mergeList[max(treeInfo$mergeList[,1]),4] <- rootLevel + 1
# 					}
# 				}
# 			}
# 			else {
# 				stop("This data set is too small.")
# 			}
# 		}
# 	}
# 	treeInfo
# }
#
# dendr <- function(cTree) {
# 	# plot a dendrogram-like structure for a component tree
# 	if (data.class(cTree) != "componentTree") {
# 		stop("Check if the argument is valid.")
# 	}
# 	class(cTree) <- "componentTree"
#
# 	n <- dim(cTree@components)[1]
# 	treeInfo <- NULL
# 	treeInfo$dataOrder <- array(0,dim=c(n,1))
# 	treeInfo$mergeList <- array(0,dim=c(n-1,4))
#
# 	root <- min(cTree@mapping)
# 	startLevel <- 1
# 	getDendrogramInfo <- searchTree(root,cTree,treeInfo,startLevel)
#
# 	heightFactor <- 1
# 	heightList <- getDendrogramInfo$mergeList[,4]
# 	heightList <- max(heightList) - heightList + 1
#
# 	reOrder <- getDendrogramInfo$dataOrder
# 	for (i in 1:n) {
# 		reOrder[getDendrogramInfo$dataOrder[i]] <- i
# 	}
#
# 	dendr <- list(merge=getDendrogramInfo$mergeList[,2:3], order=reOrder,
#  		height=heightList*heightFactor)
#
# 	class(dendr) <- "hclust"
#
# 	dendr
#
# }
#
# getNestedClustering <- function(clustering) {
# 	# clustering: a list of un-nested clusterings
# 	# the function returns a list of nested clusterings
# 	# according to the appearing of pairs of points in the un-nested clusterings
#
# 	n <- dim(clustering)[1]
# 	m <- dim(clustering)[2]
#
# 	record <- array(1,dim = c(1,n))
# 	seed <- array(0,dim = c(1,n))
# 	result <- array(0,dim = c(n,m))
#
# 	for (i in 1:m) {
# 		j <- 1
# 		foundSeed <- 0
# 		while ((foundSeed !=1) && (j <= n)) {
# 			#find the first pair of seeds for the first component, update the record list
# 			if (record[j] == 1) {
# 				foundPair <- 0
# 				k <- 1
# 				while ((foundPair != 1) && (k <= n)){
# 					numAppearing <- sum((clustering[j,]!=0)&(clustering[k,]!=0)&(clustering[j,]==clustering[k,]))
# 					if ((k!=j)&&(numAppearing >= i)) {
# 						foundPair <- 1
# 						seedIndex <- 1
# 						seed[seedIndex] <- j
# 						record[j] <- 2
# 						seedIndex <- 2
# 						seed[seedIndex] <- k
# 						record[k] <- 2
# 						foundSeed <- 1
# 						componentNumber <- 1
# 					}
# 					else
# 						k <- k + 1
# 				}
# 				if (foundPair == 0) {
# 					record[j] <- 0
# 					j <- j + 1
# 				}
# 			}
# 			else
# 				j <- j + 1
# 		}
#
# 		#find all the components at the current level of co-appearing times for a pair of points
# 		while (seedIndex != 0) {
# 			#find current component
# 			while (seedIndex != 0) {
# 				firstPoint <- seed[seedIndex]
# 				seedIndex <- seedIndex - 1
# 				#update the result set
# 				result[firstPoint,i] <- componentNumber
# 				# find all the pair points of the first point
# 				j <- firstPoint
# 				k <- 1
# 				while (k <= n) {
# 					if (record[k] == 1) {
# 						numAppearing <- sum((clustering[j,]!=0)&(clustering[k,]!=0)&(clustering[j,]==clustering[k,]))
# 						if ((k!=j)&&(numAppearing >= i)) {
# 							seedIndex <- seedIndex + 1
# 							seed[seedIndex] <- k
# 							record[k] <- 2
# 						}
# 					}
# 					k <- k + 1
# 				}
# 			}
# 			componentNumber <- componentNumber + 1
#
# 			for (j in 1:n) {
# 				seed[j] <- 0
# 			}
#
# 			j <- 1
# 			foundSeed <- 0
# 			while ((foundSeed !=1) && (j <= n)) {
# 				#find the first pair of seeds for the next component, update the record list
# 				if (record[j] == 1) {
# 					foundPair <- 0
# 					k <- 1
# 					while ((foundPair != 1) && (k <= n)){
# 						numAppearing <- sum((clustering[j,]!=0)&(clustering[k,]!=0)&(clustering[j,]==clustering[k,]))
# 						if ((k!=j)&&(numAppearing >= i)) {
# 							foundPair <- 1
# 							seedIndex <- 1
# 							seed[seedIndex] <- j
# 							record[j] <- 2
# 							seedIndex <- 2
# 							seed[seedIndex] <- k
# 							record[k] <- 2
# 							foundSeed <- 1
# 						}
# 						else
# 							k <- k + 1
# 					}
# 					if (foundPair == 0) {
# 						record[j] <- 0
# 						j <- j + 1
# 					}
# 				}
# 				else
# 					j <- j + 1
# 			}
# 		}
#
# 		for (j in 1:n) {
# 			seed[j] <- 0
# 			if (record[j]!=0) {
# 				record[j] <- 1
# 			}
# 		}
# 	}
# 	for (i in 1:m) {
# 		if (max(result[,i]) > 0)
# 			m1 <- i
# 	}
# 	result[,1:m1]
# }
#
# getFamilyofBranchComponent <- function(clusteringFamily) {
# 	# generate the family of branch components from clusteringFamily
# 	# which is an n by m matrix with the value of <i,j> indicates the
# 	# component/cluster ID to which the point i of the jth member in the family belongs.
# 	# The result is an n by p matrix with <i,j> indicates the branch component ID
# 	# to which the point i of the jth member in the family belongs.
# 	n <- dim(clusteringFamily)[1]
# 	m <- dim(clusteringFamily)[2]
#
# 	branchComponentFamily <- array(0,dim=c(n,m))	# result of the generated family of branch components
# 	#familySize <- 1
# 	currentComponents <- array(0,dim=c(n,1)) 	# represents the components projecting down onto the current node level,
# 								# this structure is used for finding the immediate branches at the next node level.
# 	currentComponentsLevels <- array(0,dim=c(n,1))	# represents the current node level for each point which belongs to
# 									# a node
# 	nextComponents <- array(0,dim=c(n,1))	# the next level for the currentComponents
# 	componentNumList <- array(0,dim=c(m,1))	# the current total num of components in the current member
# 								# in the generated family of branch components
# 	newArray <- array(0,dim=c(n,1))		# used to reset an array
#
# 	branchComponentFamily[,1] <- clusteringFamily[,1]	# no change to the first member of the family of clustering/components
# 	currentComponents <- branchComponentFamily[,1]		# the first member of the family of clustering is also the
# 										# components projecting down onto the first node level
# 	componentNumList[1] <- max(currentComponents)
# 	for (i in 1:n) {
# 		if (currentComponents[i] > 0) {
# 			currentComponentsLevels[i] <- 1
# 		}
# 	}
#
# 	childComponents <- array(0,dim=c(n,1))	# the structure to store all the immediate branches/components of the working component
# 								# at the projected current node level
# 	childComponentsSize<-0
#
# 	for (g in 2:m) {
# 		#addFamilySize <- 0
# 		indexInNextComponents <- 0
# 		currentComponentNum <- max(currentComponents)	# get the num of components in the projected current node level
# 		#print(currentComponents)
# 		for (k in 1:currentComponentNum) {
# 			# in the next member of the given family, get the subset of points from the working component
# 			for (j in 1:n) {
# 				if (currentComponents[j] == k) {
# 					childComponents[j] <- clusteringFamily[j,g]
# 				}
# 			}
#
# 			# determine if the subset contains at least two disjoint sets
# 			childComponentIdList <- as.integer(levels(factor(childComponents)))
# 			listLength <- length(childComponentIdList)
# 			if (min(childComponentIdList) == 0) {	# remove the value zero from the list
# 				realListLength <- listLength - 1
# 				realChildComponentIdList <- array(0,dim=c(realListLength,1))
# 				index <- 1
# 				for (t in 1:listLength) {
# 					if  (childComponentIdList[t] > 0) {
# 						realChildComponentIdList[index] <- childComponentIdList[t]
# 						index <- index + 1
# 					}
# 				}
# 			}
# 			else {
# 				realListLength <- listLength
# 				realChildComponentIdList <- childComponentIdList
# 			}
# 			#print(g)
# 			#print(realChildComponentIdList )
#
# 			# after finding immediate branch components, add them to the result structure
# 			# and adjust the component IDs to be consecutive
# 			if ( realListLength > 1) {
# 				#if (addFamilySize == 0) {
# 				#	familySize <- familySize + 1
# 				#	addFamilySize <- 1
# 				#}
# 				getTheLevel <- 0
# 				for (j in 1:n) {
# 					if ( currentComponents[j] == k) { # for the branch components which are the children of the current working component
# 						# increase the node level for the points which are contained in the branch components
# 						currentComponentsLevels[j] <- currentComponentsLevels[j] + 1
# 						if (getTheLevel == 0) {
# 							getTheLevel <- currentComponentsLevels[j]
# 						}
# 						#branchComponentFamily[j,familySize] <- clusteringFamily[j,g]
# 						#branchComponentFamily[j,getTheLevel] <- clusteringFamily[j,g]
# 						for ( e in 1:realListLength) {
# 							#if (branchComponentFamily[j,familySize]==realChildComponentIdList[e]) {
# 							#	nextComponents[j] <- indexInNextComponents + e
# 							#}
# 							if (clusteringFamily[j,g]==realChildComponentIdList[e]) {
# 								nextComponents[j] <- indexInNextComponents + e
# 								branchComponentFamily[j,getTheLevel] <- componentNumList[getTheLevel] + e
# 							}
# 						}
# 					}
# 				}
# 				indexInNextComponents <- indexInNextComponents + realListLength # adjust the ID for the components at the next projected node level
# 				componentNumList[getTheLevel] <- componentNumList[getTheLevel] + realListLength
# 			}
#
# 			# project down the current working component to the next level if no branch components are found
# 			if ( realListLength == 1) {
# 				for (j in 1:n) {
# 					if ( currentComponents[j] == k) {
# 						for ( e in 1:realListLength) {
# 							if (clusteringFamily[j,g]==realChildComponentIdList[e]) {
# 								nextComponents[j] <- indexInNextComponents + e
# 							}
# 						}
# 					}
# 				}
# 				indexInNextComponents <- indexInNextComponents + realListLength
# 			}
# 			childComponents <- newArray
# 		}
# 		currentComponents <- nextComponents
# 		nextComponents <- newArray
# 	}
# 	# output the generated family from the first to the last node level
# 	if (max(currentComponentsLevels) > 1) {
# 		branchComponentFamily[,1:max(currentComponentsLevels)]
# 	}
# 	else {
# 		out <- array(branchComponentFamily[,1:1],dim=c(n,1))
# 		out
# 	}
# }

# getTreeStructure <- function(branchComponentFamily){
# 	# get the content for each node and the node relation (parent node and child node),
# 	# branchComponentFamily is an n by p matrix with <i,j> indicates the branch component ID
# 	# to which the point i of the jth member in the family belongs.
# 	# The result has two components: nodeContent and nodeMapping.
# 	# The first row in nodeContent is the node ID, the second row in nodeContent is the points
# 	# belonging to the node denoted by the ID in the first row at the same column; each row
# 	# in nodeMapping denotes a parent node and its child node
#
# 	n <- dim(branchComponentFamily)[1]
# 	m <- dim(branchComponentFamily)[2]
# 	notFound <- 1
# 	nodeMapping <- NULL
#
# 	totalComponentsNum <- 0
# 	for (i in 1:m) {
# 		componentsNum <- max(branchComponentFamily[,i])  # get the number of components in the current member of the family
# 		if (i > 1) {
# 			mappingCheck <- array(0,dim=c(componentsNum,1))  # check if the parent node of the current component is found
# 		}
# 		for (j in 1:n) {
# 			if (branchComponentFamily[j,i] > 0) { # the point belongs to a component at the current member of the family (level of the tree)
# 				if ((i > 1) && (mappingCheck[branchComponentFamily[j,i]] == 0)) {
# 					nodeMapping <- rbind(nodeMapping,c(max(nodeContent[1,nodeContent[2,]==j]),totalComponentsNum + branchComponentFamily[j,i]))
# 					mappingCheck[branchComponentFamily[j,i]] <- 1
# 				}
# 				if (notFound == 1) {
# 					nodeContent <- c(totalComponentsNum + branchComponentFamily[j,i],j)
# 					notFound <- 0
# 				}
# 				else {
# 					nodeContent <- cbind(nodeContent,c(totalComponentsNum + branchComponentFamily[j,i],j))
# 				}
# 			}
# 		}
# 		totalComponentsNum <- totalComponentsNum + componentsNum
# 	}
#
# 	#if (is.null(nodeMapping)==FALSE) {
# 	#	nodeMapping <- nodeMapping[order(nodeMapping[,1]),]
# 	#}
# 	#else {
# 	#	nodeMapping <- array(0,dim=c(1,1))
# 	#}
#
# 	addCommonRoot <- 0
# 	totalNodes <- max(nodeContent[1,])
# 	multiRoots <- array(0,dim=c(totalNodes,1))
# 	# if multiple roots, create a common root with the whole data set
# 	if (is.null(nodeMapping)==FALSE) {
# 		for (i in 1:totalNodes) {
# 			if (sum(nodeMapping == i) == 0) {
# 				multiRoots[i] <- 1
# 			}
# 		}
# 		for (i in 1:(dim(nodeMapping)[1])) {
# 			if (sum(nodeMapping[,2] == nodeMapping[i,1]) == 0) {
# 				multiRoots[nodeMapping[i,1]] <- 1
# 			}
# 		}
# 		if (sum(multiRoots) > 1) {
# 			for (i in 1:totalNodes) {
# 				if (multiRoots[i] == 1) {
# 					nodeMapping <- rbind(nodeMapping,c(0,i))
# 				}
# 			}
# 			for (i in 1:n) {
# 				nodeContent <- cbind(nodeContent,c(0,i))
# 			}
# 			addCommonRoot <- 1
# 		}
# 		else {
# 			nodeContent[1,] <- nodeContent[1,] - 1
# 			nodeMapping <- nodeMapping -1
# 		}
# 	}
# 	else {
# 		nodeMapping <- array(0,dim=c(totalNodes,2))
# 		for (i in 1:totalNodes) {
# 			nodeMapping[i,] <- c(0,i)
# 		}
# 		for (i in 1:n) {
# 			nodeContent <- cbind(nodeContent,c(0,i))
# 		}
# 		addCommonRoot <- 1
# 	}
#
# 	nodeContent <- nodeContent[,order(nodeContent[1,])]
# 	nodeMapping <- nodeMapping[order(nodeMapping[,1]),]
# 	treeStructure <- list(nodeContent = nodeContent,nodeMapping = nodeMapping,commonRoot = addCommonRoot)
# 	treeStructure
# }

# getFamilyofGraphsByBranchComponent <- function(branchComponentFamily) {
# 	# get the family of graphs from branchComponentFamily
# 	# the result is an n by m matrix with <i,j> indicates edge i appears
# 	# in the jth graph of the family
# 	n <- dim(branchComponentFamily)[1]
# 	m <- dim(branchComponentFamily)[2]
# 	len <- n*(n-1)/2
#
# 	graphFamily <- array(0, dim = c(len,m+2))
# 	familySize <- 1
# 	completeEdgeSet <- array(0, dim = c(n,n))
# 	emptyEdgeSet <- array(0, dim = c(n,n))
# 	currentCluster <- array(0, dim = c(n,1))
#
# 	for (g in 1:m) {
# 		componentIdList <- as.integer(levels(factor(branchComponentFamily[,g])))
#
# 		if (max(componentIdList) > 0) {
# 			for (k in 1:length(componentIdList)) {
# 				if (componentIdList[k] != 0) {
# 					nSize <- 0
# 					nCluster <- componentIdList[k]
# 					for (i in 1:n) {
# 						if (branchComponentFamily[i,g] == nCluster) {
# 							nSize <- nSize + 1
# 							currentCluster[nSize] <- i
# 						}
# 					}
# 					if ( nSize > 1) {
# 						for (i in 1:(nSize-1)) {
# 							for (j in (i+1):nSize) {
# 								completeEdgeSet[currentCluster[i],currentCluster[j]] <- 1
# 							}
# 						}
# 					}
# 				}
# 			}
# 			location <- 1
# 			for (i in 1:(n-1)) {
# 				for (j in (i+1):n) {
# 					graphFamily[location,familySize] <- completeEdgeSet[i,j]
# 					if (familySize == 1) {
# 						graphFamily[location,m+1] <- i
# 						graphFamily[location,m+2] <- j
# 					}
# 					location <- location + 1
# 				}
# 			}
# 			familySize <- familySize + 1
# 			completeEdgeSet <- emptyEdgeSet
# 		}
# 	}
# 	graphFamily
# }

# getGraphTot <- function(graphFamily) {
# 	# calculate the total of the family of graphs
# 	n <- dim(graphFamily)[1]
# 	m <- dim(graphFamily)[2]
#
# 	graphTot <- array(0,dim=c(n,1))
#
# 	for (j in 1:n) {
# 		for (i in 1:m) {
# 			graphTot[j] <- graphTot[j] + graphFamily[j,i]
# 		}
# 	}
# 	graphTot
# }

# generateTreebyBranchComponents <- function(branchComponentFamily) {
# 	# create an object of "componentTree" from the branch component family
#
# 	m <- dim(branchComponentFamily)[2]
# 	if (m == 1) {
# 		n <- dim(branchComponentFamily)[1]
# 		len <- choose(n,2)
# 		theMapping <- array(0,dim=c(1,2))
# 		theContent <- array(0,dim=c(2,n))
# 		theContent[2,] <- c(1:n)
# 		graphTot <- array(1,dim=c(len,1))
# 		out <- new("componentTree",components=branchComponentFamily,content = theContent,
# 			mapping = theMapping,total =graphTot)
# 	}
# 	else {
# 		theTreeStructure <- getTreeStructure(branchComponentFamily)
# 		if (theTreeStructure$commonRoot == 1) {
# 			branchComponentFamily <- cbind(rep(1,dim(branchComponentFamily)[1]),branchComponentFamily)
# 		}
# 		generatedGraphFamily <- getFamilyofGraphsByBranchComponent(branchComponentFamily)
# 		len <- dim(generatedGraphFamily)[2]-2  # not read the last two columns which are the vertices pairs
# 		if (len > 1) {
# 			graphTot <- getGraphTot(generatedGraphFamily[,1:len])
# 		}
# 		else {
# 			n <- dim(generatedGraphFamily)[1]
# 			tmp <- array(generatedGraphFamily[,1:1],dim=c(n,1))
# 			graphTot <- getGraphTot(tmp)
# 		}
# 		out <- new("componentTree",components=branchComponentFamily,content = theTreeStructure$nodeContent,
# 			mapping = theTreeStructure$nodeMapping,total =graphTot)
# 	}
# 	out
# }

# getTreeDistance <- function(T1,T2){
# 	# calculate the tree distance between T1 and T2
# 	if ((data.class(T1) != "componentTree")||(data.class(T2) != "componentTree")) {
# 		stop("Check if the argument is valid.")
# 	}
# 	class(T1) <- "componentTree"
# 	class(T2) <- "componentTree"
#
# 	w0 <- T1@total - 1
# 	if (sum(w0^2) > 0) {
# 		w0 <- w0/sqrt(sum(w0^2))
# 	}
#
# 	w1 <- T2@total - 1
# 	if (sum(w1^2) > 0) {
# 		w1 <- w1/sqrt(sum(w1^2))
# 	}
# 	sqrt(sum((w1-w0)^2))
# }
#
# getTreeAveraging <- function(tree_list) {
# 	# generate the tree averaging from a list of trees
# 	n <- length(tree_list)
# 	components_com <- tree_list[[1]]@components
# 	if (n>=2) {
# 		for (i in 2:n) {
# 			layer <- dim(tree_list[[i]]@components)[2]
# 			if (layer>1) {
# 				components_com <- cbind(components_com,tree_list[[i]]@components)
# 			}
# 		}
# 	}
# 	branchComponent2 <- getNestedClustering(components_com)
# 	components <- getFamilyofBranchComponent(branchComponent2)
# 	Tree <- generateTreebyBranchComponents(components)
#
# 	Tree
# }
