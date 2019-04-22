
#' plot a clusterTree object
#' @param x a clusterTree object
#' @param y NULL.  Will be ignored with a warning if non-NULL
#' @param labels labels for each data point
#' @param axes whether to plot axis on the left
#' @param frame.plot whether to plot frame for density plot
#' @param ann whether to annotate main, sub, xlab, ylab
#' @param main main title for the plot
#' @param sub subtitle for the plot
#' @param xlab label for x-axis
#' @param ylab label for y-axis
#' @param col color of rectangle
#' @param labels.cex a character expansion multiplier for the size of the object labels in the plot
#' @param labels.col colours of the object labels in the plot
#' @param ... remains to be processed
#' @return No return value (invisible()) 
#' @examples
#' data <- rbind(matrix(rnorm(100, mean = 10, sd = 2), nrow = 50),
#'               matrix(rnorm(100, mean = 0, sd = 1), nrow = 50),
#'               matrix(rnorm(100, mean = -10, sd = 3), nrow = 50)
#'               )
#' clustering1<-stats::hclust(dist(data),method='single')
#' clustering2<-kmeans(data,centers=3)
#' clustering3<-dbscan::dbscan(data,eps=.1)
#' res <- combineClusterings(clustering1,clustering2,clustering3)
#' plot(res)
#' @export 
plot.clusterTree <- function(x, y = NULL, labels = NULL, axes = TRUE, frame.plot = FALSE, ann = TRUE, 
                             main = "Cluster Tree Density Plot", sub = NULL, xlab = NULL, ylab = "Height",
                             col = NULL, labels.cex =1, labels.col = NULL, ...){
  clusterTreePlotInfo <- clusterTreeToClusterTreePlotInfo(x)
  graphics::plot(x = clusterTreePlotInfo, y = y, labels = labels, axes = axes, frame.plot = frame.plot, ann = ann, 
                             main = main, sub = sub, xlab = xlab, ylab = ylab,
                             col = col, labels.cex = labels.cex, labels.col = labels.col, ...)
}

#' transform a clusterTree object into a clusterTreePlotInfo object
#' clusterTreePlotInfo object is a unified data structure for recording information
#' related to plotting clusterTree object
#' @param x a clusterTree object
#' @param padding padding must have a value between 0 and 1, representing the proportion of space
#' @return a clusterTreePlotInfo object
#' @return clusterTreePlotInfo is a nested structure which record information for 
#' @return plotting clusterTree 
#' \item{tree}{a matrix, tree attribute of clusterTree object that remains to be plotted}
#' \item{rectangles}{a list of rectangles, each element record left,bottom,right,top coordinates of the rectangle}
#' \item{labels}{a list of list of labels' index, each element record labels' index for corresponding rectangle}
#' \item{lines}{a list of list of lines, each element record lines' start and end point}
#' \item{is.runts}{a list of boolean variable, each element record whether that rectangle is runts or not}
#' \item{subtree}{a list of subtree, each element represents the nested subtree of corresponding rectangle}
#' @examples
#' data <- rbind(matrix(rnorm(100, mean = 10, sd = 2), nrow = 50),
#'               matrix(rnorm(100, mean = 0, sd = 1), nrow = 50),
#'               matrix(rnorm(100, mean = -10, sd = 3), nrow = 50)
#'               )
#' clustering1<-stats::hclust(dist(data),method='single')
#' clustering2<-kmeans(data,centers=3)
#' clustering3<-dbscan::dbscan(data,eps=.1)
#' clusterTree <- combineClusterings(clustering1,clustering2,clustering3)
#' res <- clusterTreeToClusterTreePlotInfo(clusterTree)
#' @export 
clusterTreeToClusterTreePlotInfo <- function(x, padding = .2){
  if(class(x)[1]!='clusterTree')
    stop('input is not a clusterTree object!')
  if(padding <= 0 | padding >= 1)
    stop('padding must belong to (0,1)!')
  # reorder treeMatrix attribute so that data points within same 
  # cluster are close to each other
  order <- reOrderClusterTreeMatrix(x$treeMatrix)
  orderedTreeMatrix <- as.matrix(x$treeMatrix[order,])
  # add labels' index into plot info
  # this helps us identify the right index of labels for each label's position
  labels_index <- c(1:nrow(x$treeMatrix))
  labels_index <- labels_index[order]
  # call a helper function to help organize plot info for clusterTree in recursive fashion
  res <- recursivePlotInfo(treeMatrix = orderedTreeMatrix, leftBorderPlotRegion = 0, bottomBorderPlotRegion = 0, 
    rightBorderPlotRegion = 1, topBorderPlotRegion = 1, labels_index = labels_index, padding = padding)
  # class(res) <- c("clusterTreePlotInfo", class(res))
  res
}

#' a function to help organize plot info for clusterTree in a recursive fashion
#' First input is a ordered treeMatrix
#' My model assumes we need to plot all rectangles, lines, labels etc on (0,1)*(0,1) canvas
#' This function returns relative coordinates of rectangles, lines etc.
#' @param treeMatrix a ordered treeMatrix of clusterTree
#' @param leftBorderPlotRegion coordinate of the plot region's left border on (0,1)*(0,1) canvas
#' @param bottomBorderPlotRegion coordinate of the plot region's bottom border on (0,1)*(0,1) canvas
#' @param rightBorderPlotRegion coordinate of the plot region's right border on (0,1)*(0,1) canvas
#' @param topBorderPlotRegion coordinate of the plot region's top border on (0,1)*(0,1) canvas
#' @param labels_index index of labels to paste in each appropriate position
#' @param padding padding of rectangle in rectangle plot region
recursivePlotInfo <- function(treeMatrix, leftBorderPlotRegion, 
                              bottomBorderPlotRegion, rightBorderPlotRegion,
                              topBorderPlotRegion, labels_index, padding){
  # target of this function is to generate plot info of treeMatrix recursively
  # plot info includes treeMatrix, rectangles, labels, lines, is.runts, subtrees, numBranches
  # there is a correspondance relation between plot info, which is:
  # rectangles[[1]] <-> labels[[1]] <-> lines[[1]] <-> is.runts[[1]] <-> subtrees[[1]] <-> 'first branch'
  # rectangles[[2]] <-> labels[[2]] <-> lines[[2]] <-> is.runts[[2]] <-> subtrees[[2]] <-> 'second branch', etc...
  # function is implemented recursively
  n <- nrow(treeMatrix)
  # if ncol(treeMatrix) is 1, build plot info directly, and return
  # if ncol(treeMatrix) is not 1, build plot info for the first column, add
  # those info into return value
  # Then call function recursively to obtain plot info for treeMatrix[,2:n], 
  # add plot info of treeMatrix[,2:n] as subtree attribute of return value 
  if(ncol(treeMatrix)==1){
    # build plot info directly
    # use an row index iterator to iterate along first column
    rowIndexIter <- 1
    # left border of rectangle region, starting from left border of plot region
    leftBorderRectangleRegion <- leftBorderPlotRegion
    # initialize empty list for rectangles,labels,lines,is.runts and subtree
    rectangles <- list()
    subBranchLabelGroups <- list()
    subBranchLineGroups <- list()
    runtsFlag <- list()
    subtrees <- list()
    # start to iterate among rows
    while(rowIndexIter <= n){
        # find start row index and end row index of a cluster
        # record starting row index of a cluster
        clusterStartRowIndex <- rowIndexIter
        # initialize clusterEndRowIndex
        clusterEndRowIndex <- rowIndexIter
        # iterate if treeMatrix[clusterEndRowIndex,1] and treeMatrix[clusterStartRowIndex,1] have same value
        # have to deal with situation with NA values
        # in front of 'or' deal with non-NA situation 
        # after 'or' deal with NA situation
        while ((!is.na(treeMatrix[clusterStartRowIndex,1]) & 
                !is.na(treeMatrix[clusterEndRowIndex,1]) & 
                treeMatrix[clusterEndRowIndex,1]==treeMatrix[clusterStartRowIndex,1]) |
                (is.na(treeMatrix[clusterStartRowIndex,1]) & is.na(treeMatrix[clusterEndRowIndex,1]))) {
          if(clusterEndRowIndex<n){ clusterEndRowIndex <- clusterEndRowIndex + 1 }
          else{ break }
        }
        if(clusterEndRowIndex==n){ clusterEndRowIndex <- clusterEndRowIndex + 1 }
        # assign right border of rectangle region
        # length of right border - left border is proportional to (clusterEndRowIndex - clusterStartRowIndex)
        rightBorderRectangleRegion <- leftBorderRectangleRegion + 
                                      (clusterEndRowIndex-clusterStartRowIndex)/n*(rightBorderPlotRegion - leftBorderPlotRegion)
        # create a rectangle and add to rectangles
        rectangle <- createRectangleWithinRegion(leftBoundary = leftBorderRectangleRegion, bottomBoundary = bottomBorderPlotRegion,
                                                 rightBoundary = rightBorderRectangleRegion, topBoundary = topBorderPlotRegion, padding = padding)
        rectangles[[length(rectangles)+1]] <- rectangle
        # find x,y coordinates of labels
        labels_x <- seq(from = rectangle$left, to = rectangle$right, length.out = clusterEndRowIndex - clusterStartRowIndex) 
        labels_y <- (bottomBorderPlotRegion + rectangle$bottom)/2 
        indexOfLabels <- labels_index[clusterStartRowIndex:(clusterEndRowIndex-1)]
        labelsInfo <- list()
        for (j in c(1:(clusterEndRowIndex - clusterStartRowIndex))) {
          labelsInfo[[j]] <- list(x=labels_x[j], y=labels_y, index = indexOfLabels[j])
        }
        # add labelsInfo to labels
        subBranchLabelGroups[[length(subBranchLabelGroups)+1]] <- list(labels = labelsInfo)
        # empty lines info for this column
        subBranchLineGroups[[length(subBranchLineGroups)+1]] <- list(lines = list())
        # check whether it's a runt by checking whether it's a NA value
        runtsFlag[[length(runtsFlag)+1]] <- is.na(treeMatrix[clusterStartRowIndex,1])
        rowIndexIter <- clusterEndRowIndex
        leftBorderRectangleRegion <- rightBorderRectangleRegion
      }
    result <- list(treeMatrix = treeMatrix, plotLayerInfo = list(rectangles = rectangles, subBranchLabelGroups = subBranchLabelGroups,
                   subBranchLineGroups = subBranchLineGroups, runtsFlag = runtsFlag), subtrees = subtrees)
    class(result) <- c("clusterTreePlotInfo", class(result))
    result
  }
  else{
    # build plot info for the first column
    # use an row index iterator to iterate along first column
    rowIndexIter <- 1
    # relative position of leftBorderRectangleRegion to leftBorderPlotRegion
    # left border of rectangle region, starting from left border of plot region
    leftBorderRectangleRegion <- leftBorderPlotRegion
    # initialize empty list for rectangles,labels,lines,is.runts and subtree
    rectangles <- list()
    subBranchLabelGroups <- list()
    subBranchLineGroups <- list()
    runtsFlag <- list()
    subtrees <- list()
    # start to iterate among rows
    while(rowIndexIter <= n){
        # find start row index and end row index of a cluster
        # record starting row index of a cluster
        clusterStartRowIndex <- rowIndexIter
        # initialize clusterEndRowIndex
        clusterEndRowIndex <- rowIndexIter
        # iterate if treeMatrix[clusterEndRowIndex,1] and treeMatrix[clusterStartRowIndex,1] have same value
        # have to deal with situation with NA values
        # in front of 'or' deal with non-NA situation 
        # after 'or' deal with NA situation
        while ((!is.na(treeMatrix[clusterStartRowIndex,1]) & 
                !is.na(treeMatrix[clusterEndRowIndex,1]) & 
                treeMatrix[clusterEndRowIndex,1]==treeMatrix[clusterStartRowIndex,1]) |
                (is.na(treeMatrix[clusterStartRowIndex,1]) & is.na(treeMatrix[clusterEndRowIndex,1]))) {
          if(clusterEndRowIndex<n){ clusterEndRowIndex <- clusterEndRowIndex + 1 }
          else{ break }
        }
        if(clusterEndRowIndex==n){ clusterEndRowIndex <- clusterEndRowIndex + 1 }
        # assign right border of rectangle region
        # length of right border - left border is proportional to (clusterEndRowIndex - clusterStartRowIndex)
        rightBorderRectangleRegion <- leftBorderRectangleRegion + 
                                      (clusterEndRowIndex-clusterStartRowIndex)/n*(rightBorderPlotRegion - leftBorderPlotRegion)
        topBorderRectangleRegion <- bottomBorderPlotRegion + 1/ncol(treeMatrix)*(topBorderPlotRegion - bottomBorderPlotRegion)
        # create a rectangle and add to rectangles
        rectangle <- createRectangleWithinRegion(leftBoundary = leftBorderRectangleRegion, 
                                                 bottomBoundary = bottomBorderPlotRegion,
                                                 rightBoundary = rightBorderRectangleRegion,
                                                 topBoundary = topBorderRectangleRegion, padding = padding)
        rectangles[[length(rectangles)+1]] <- rectangle
        # find x,y coordinates of labels
        labels_x <- seq(from = rectangle$left, to = rectangle$right, length.out = clusterEndRowIndex - clusterStartRowIndex) 
        labels_y <- (bottomBorderPlotRegion + rectangle$bottom)/2 
        indexOfLabels <- labels_index[clusterStartRowIndex:(clusterEndRowIndex-1)]
        labelsInfo <- list()
        for (j in c(1:(clusterEndRowIndex - clusterStartRowIndex))) {
          labelsInfo[[j]] <- list(x=labels_x[j], y=labels_y, index = indexOfLabels[j])
        }
        # add labelsInfo to labels
        subBranchLabelGroups[[length(subBranchLabelGroups)+1]] <- list(labels = labelsInfo)
        # create subtree from recursive call
        subtree <- recursivePlotInfo(
          treeMatrix = as.matrix(treeMatrix[clusterStartRowIndex:(clusterEndRowIndex-1),2:ncol(treeMatrix)]),
          leftBorderPlotRegion = leftBorderRectangleRegion, 
          bottomBorderPlotRegion = topBorderRectangleRegion, 
          rightBorderPlotRegion = rightBorderRectangleRegion, 
          topBorderPlotRegion = topBorderPlotRegion, 
          labels_index = labels_index[clusterStartRowIndex:(clusterEndRowIndex-1)], padding = padding)
        # add subtree into subtrees list
        subtrees[[length(subtrees)+1]] <- subtree
        subtreeRectangles <- subtree$plotLayerInfo$rectangles
        # collect lines info for this particular subtree
        subtreeLines <- list()
        # for each rectangle in subtree, create a line
        for (j in 1:length(subtreeRectangles)) {
          # extract jth rectangle in subtree
          cur_rectangle <- subtreeRectangles[[j]]
          # add a line that links current rectangle to jth rectangle in subtree
          subtreeLines[[length(subtreeLines)+1]] <- list(start = list(x = (leftBorderRectangleRegion+rightBorderRectangleRegion)/2.0, y = rectangle$top), end = list(x = (cur_rectangle$left + cur_rectangle$right)/2.0, y = cur_rectangle$bottom))
        }
        # add this group of lines to list
        subBranchLineGroups[[length(subBranchLineGroups)+1]] <- list(lines = subtreeLines)
        # add is.runts information to list
        runtsFlag[[length(runtsFlag)+1]] <- is.na(treeMatrix[clusterStartRowIndex,1])
        # update rowIndexIter and leftBorderRectangleRegion
        rowIndexIter <- clusterEndRowIndex
        leftBorderRectangleRegion <- rightBorderRectangleRegion
      }
    result <- list(treeMatrix = treeMatrix, plotLayerInfo = list(rectangles = rectangles, subBranchLabelGroups = subBranchLabelGroups,
                   subBranchLineGroups = subBranchLineGroups, runtsFlag = runtsFlag), subtrees = subtrees)
    class(result) <- c("clusterTreePlotInfo", class(result))
    result
  }
}



#' create a rectangle in a region with padding parameter
#' @param leftBoundary left boundary of rectangle plot region
#' @param bottomBoundary bottom boundary of rectangle plot region
#' @param rightBoundary bottom boundary of rectangle plot region
#' @param topBoundary bottom boundary of rectangle plot region
#' @param padding padding of rectangle inside rectangle plot region
createRectangleWithinRegion <- function(leftBoundary, bottomBoundary, rightBoundary, topBoundary, padding){
  width <- rightBoundary - leftBoundary
  height <- topBoundary - bottomBoundary
  left <- leftBoundary + padding * width
  bottom <- bottomBoundary + padding * height
  right <- leftBoundary + (1-padding) * width
  top <- bottomBoundary + (1-padding)*height
  list(left = left, bottom = bottom, right = right, top = top)
}

####
# some function to understand structure of clusterTreePlotInfo

#' print number of branches of clusterTreePlotInfo object
#' @param x a clusterTreePlotInfo object
#' @return a integer, represents number of branches
#' @export
nBranches <- function(x){
  length(x$rectangles)
}

#' show treeMatrix of clusterTreePlotInfo object
#' @param x a clusterTreePlotInfo object
#' @return a matrix, represents treeMatrix
#' @export
showTreeMatrix <- function(x){
  x$treeMatrix
}

#' show coordinates of rectangles in clusterTreePlotInfo object
#' @param x a clusterTreePlotInfo object
#' @return a list of rectangle coordinates
#' @export
showRectangles <- function(x){
  x$rectangles
}

#' show coordinates of lines in clusterTreePlotInfo object
#' @param x a clusterTreePlotInfo object
#' @return a list of lines' coordinates
showLines <- function(x){
  x$lines
}

#' show coordinates of labels in clusterTreePlotInfo object
#' @param x a clusterTreePlotInfo object
#' @return a list of labels' coordinates
#' @export
showLabels <- function(x){
  x$labels
}

#' get plot information on a layer
#' @param x a clusterTreePlotInfo object
#' @param layer layer to check
#' @return a list which contains plot information 
#' @export
getBranchInfo <- function(x, layer = 1){
  if(layer <= 0)
    stop('layer must be greater than or equal to 1')
  if(layer == 1){
    res <- list(rectangles = x$rectangles, labels = x$labels, lines = x$lines, 
      is.runts = x$is.runts)
    res
  }
  else{
    num_branches <- length(x$subtree)
    rectangles <- list()
    labels <- list()
    lines <- list()
    is.runts <- list()
    if(num_branches > 0){
      for (ii in c(1:num_branches)) {
        info <- getBranchInfo(x$subtree[[ii]], layer = layer - 1)
        rectangles <- c(rectangles, info$rectangles)
        labels <- c(labels, info$labels)
        lines <- c(lines, info$lines)
        is.runts <- c(is.runts, info$is.runts)
      }
      res <- list(rectangles = rectangles, labels = labels, lines = lines, 
      is.runts = is.runts)
      res
    }
    else{
      return(list(numBranches = 0))
    }
  }
}
####

#' plot a clusterTreePlotInfo object
#' @param x a clusterTree object
#' @param y NULL.  Will be ignored with a warning if non-NULL
#' @param labels labels for each data point
#' @param axes whether to plot axis on the left
#' @param frame.plot whether to plot frame for density plot
#' @param ann whether to annotate main, sub, xlab, ylab
#' @param main main title for the plot
#' @param sub subtitle for the plot
#' @param xlab label for x-axis
#' @param ylab label for y-axis
#' @param col color of rectangle
#' @param labels.cex a character expansion multiplier for the size of the object labels in the plot
#' @param labels.col colours of the object labels in the plot
#' @param ... remains to be processed
#' @return No return value  
#' @examples
#' data <- rbind(matrix(rnorm(100, mean = 10, sd = 2), nrow = 50),
#'               matrix(rnorm(100, mean = 0, sd = 1), nrow = 50),
#'               matrix(rnorm(100, mean = -10, sd = 3), nrow = 50)
#'               )
#' clustering1<-stats::hclust(dist(data),method='single')
#' clustering2<-kmeans(data,centers=3)
#' clustering3<-dbscan::dbscan(data,eps=.8)
#' res <- combineClusterings(clustering1,clustering2,clustering3)
#' res2 <- clusterTreeToClusterTreePlotInfo(res)
#' plot(res2)
#' @export 
plot.clusterTreePlotInfo <- function(x, y = NULL, labels = NULL, axes = TRUE, frame.plot = FALSE, ann = TRUE, 
                             main = "Cluster Tree Density Plot", sub = NULL, xlab = NULL, ylab = "Height",
                             col = NULL, labels.cex =1, labels.col = NULL,  ...){
  if (!is.null(y)) warning("argument y is ignored")
  # gather what we need for clusterTree object and introduce to plotClusterTreeHelper function
  # variable checking
  if(is.null(col))
    col <- 'grey'
  if(is.null(labels.col))
    labels.col <- 'black'
  # process variable labels
  # processing accordingly with boolean input and vector input
  labels.plot <- FALSE
  if(is.logical(labels)){
    labels.plot <- labels
    labels <- paste("object", 1:nrow(x$treeMatrix))
  }else{
    if(is.null(labels)){
      labels.plot <- FALSE
      labels <- paste("object", 1:nrow(x$treeMatrix))
    }else{
      if(!is.vector(labels)){
        stop("labels must be a logical or vector!")
      }else{
        if(length(labels) != nrow(x$treeMatrix)){
          stop('length of labels and rows of clusterTree object does not match!')
        }else{
          labels.plot <- TRUE
        }
      }
    }
  }
  ###
  ###
  # plot clusterTreePlotInfo recursively
  graphics::plot.new()
  plotClusterTreePlotInfoRecursiveHelper(x = x, labels.plot = labels.plot, labels = labels, col = col, 
                                  labels.cex = 2*labels.cex/ncol(x$treeMatrix), labels.col = labels.col, ...)
  ###
  ###
  if(frame.plot){
    graphics::box()
  }
  if(ann){
    if(is.null(sub))
      sub <- "density plot"
    graphics::title(main = main, sub = sub, xlab = xlab, ylab = ylab)
  }
  if(axes)
    height <- as.double(ncol(x$treeMatrix))
    graphics::axis(2,(1/2/height)*seq(1,2*height,by = 2),labels = c(1:height))
}



#' helper function to plot a clusterTreePlotInfo object
#' @param x a clusterTreePlotInfo object
#' @param labels.plot whether to plot labels
#' @param labels labels of data points
#' @param col color of rectangles
#' @param labels.cex cex of labels
#' @param labels.col color of labels 
#' @param ... other parameters to pass to graphics::rect()
plotClusterTreePlotInfoRecursiveHelper <- function(x, labels.plot, labels, col, 
                                  labels.cex, labels.col, ...){
  # iterate through along each branches to plot rectangles, lines, etc
  for (i in c(1:length(x$plotLayerInfo$rectangles))) {
    # only plot if it is not a runt
    if(!x$plotLayerInfo$runtsFlag[[i]]){
      graphics::rect(x$plotLayerInfo$rectangles[[i]]$left, x$plotLayerInfo$rectangles[[i]]$bottom,
        x$plotLayerInfo$rectangles[[i]]$right, x$plotLayerInfo$rectangles[[i]]$top, col = col, ...)
      # if there exists ith group of lines to plot
      # note x$lines[i] corresponds to ith group of lines between ith rectangle and its bracnches,
      # not ith lines
      if(length(x$plotLayerInfo$subBranchLineGroups[[i]]$lines)>0){
        # plot jth lines in ith group of lines, note this line links ith rectangle and jth branch of ith rectangle
        for (j in c(1:length(x$plotLayerInfo$subBranchLineGroups[[i]]$lines))) {
          # only plot a line if jth branch of ith subtree is not a runt
          if(!x$subtree[[i]]$plotLayerInfo$runtsFlag[[j]]){
            graphics::segments(x0 = x$plotLayerInfo$subBranchLineGroups[[i]]$lines[[j]]$start$x,
              y0 = x$plotLayerInfo$subBranchLineGroups[[i]]$lines[[j]]$start$y,
              x1 = x$plotLayerInfo$subBranchLineGroups[[i]]$lines[[j]]$end$x,
              y1 = x$plotLayerInfo$subBranchLineGroups[[i]]$lines[[j]]$end$y)
          }
        }
      }
      # if we need to plot labels
      if(labels.plot){
        # if number of group of labels is greater than 0
        if(length(x$plotLayerInfo$subBranchLabelGroups)>0){
          # note x$labels[i][j] represents jth labels in ith group of labels
          # ith group of labels represent ith rectangle in this layer
          # jth labels represent jth labels in ith rectangle
          for (j in c(1:length(x$plotLayerInfo$subBranchLabelGroups[[i]]$labels))) {
            graphics::text(x = x$plotLayerInfo$subBranchLabelGroups[[i]]$labels[[j]]$x,
                           y = x$plotLayerInfo$subBranchLabelGroups[[i]]$labels[[j]]$y,
                           labels = labels[x$plotLayerInfo$subBranchLabelGroups[[i]]$labels[[j]]$index],
                           cex = labels.cex, 
                           col = labels.col[x$plotLayerInfo$subBranchLabelGroups[[i]]$labels[[j]]$index],
                           srt = 90)
          }
        }
      }
    }
  }
  # if subtrees exist
  if(length(x$subtrees)>0){
    for (i in c(1:length(x$subtrees))) {
      # call recursive function to plot
      plotClusterTreePlotInfoRecursiveHelper(x$subtree[[i]], labels.plot = labels.plot,
                                       labels = labels, col = col, 
                                       labels.cex = labels.cex, labels.col = labels.col)
    }
  }
}

