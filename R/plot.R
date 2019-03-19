
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
                             col = NULL, labels.cex =1, labels.col = NULL,  ...){
  if (!is.null(y)) warning("argument y is ignored")
  
  # gather what we need for clusterTree object and introduce to plotClusterTreeHelper function
  order <- reOrderClusterTreeMatrix(x$tree)
  orderedTree <- as.matrix(x$tree[order,])
  graphics::plot.new()
  height <- as.double(dim(x$tree)[2])
  
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
    labels <- x$labels
  }else{
    if(is.null(labels)){
      labels.plot <- FALSE
      labels <- x$labels
    }else{
      if(!is.vector(labels)){
        stop("labels must be a logical or vector!")
      }else{
        if(length(labels) != nrow(x$tree)){
          stop('length of labels and rows of clusterTree object does not match!')
        }else{
          labels.plot <- TRUE
        }
      }
    }
  }
  ###
  ###
  # plot clusterTree recursively
  plotClusterTreeHelper(x = orderedTree, xleft = 0, ybottom = 0, 
                        xright = 1, ytop = 1, labels.plot = labels.plot, 
                        labels = labels[order], col = col, labels.cex = labels.cex,
                        labels.col = if (length (labels.col) == length(labels)) {
                          labels.col[order]
                          } else {labels.col}
                        )
  ###
  ###
  if(frame.plot){
    graphics::box(...)
  }
  if(ann){
    if(is.null(sub))
      sub <- "density plot"
    graphics::title(main = main, sub = sub, xlab = xlab, ylab = ylab, ...)
  }
  if(axes)
    graphics::axis(2,(1/2/height)*seq(1,2*height,by = 2),labels = c(1:height))
  invisible()
}

#' a helper function of plot.clusterTree
#' this function plot clusterTree object recursively
#' @param x a matrix which is tree attribute of clusterTree object
#' @param xleft left x coordinate of rectangular plot region
#' @param ybottom lower y coordinate of rectangular plot region
#' @param xright left x coordinate of rectangular plot region
#' @param ytop upper top coordinate of rectangular plot region
#' @param labels.plot whether to plot labels
#' @param labels labels for the plot
#' @param col color of rectangle
#' @param labels.cex a character expansion multiplier for the size of the object labels in the plot
#' @param labels.col colours of the object labels in the plot
#' @return Invisibly returns the coordinates of the bottom level rectangles, each as
#'  (left bottom right top)
plotClusterTreeHelper <- function(x, xleft, ybottom, xright, ytop, 
                                  labels.plot, labels, col, 
                                  labels.cex = 1, labels.col){
  ###
  # This function plots clusterTree object recursively
  # for each recursive call, it plots rectangles, lines, labels in this layer
  # then this function calls itself recursively to plot next layer, etc
  ###
  n <- dim(x)[1]   
  if(dim(x)[2]==1){
    start <- 1
    left <- .0
    res <- c()
    while(start <= n){
      if(!is.na(x[start,1])){
        jj <- start
        while (!is.na(x[jj,1]) & x[jj,1]==x[start,1]) {
          if(jj<n){ jj <- jj+1 }
          else{ break }
        }
        if(jj==n){ jj <- jj+1 }
        right <- left + (jj-start)/n*(xright - xleft)
        lambda <- .2
        xl <- xleft+left+lambda*(right-left)
        yb <- ybottom+lambda*(ytop-ybottom)
        xr <- xleft+left+(1-lambda)*(right-left)
        yt <- ybottom+(1-lambda)*(ytop-ybottom)
        graphics::rect(xl,yb,xr,yt,col = col)
        if(labels.plot)
          graphics::text(x = seq(from = xl, to = xr, length.out = jj-start), 
                         y = rep((yb+ybottom)/2, jj-start), 
                         labels = labels[start:(jj-1)],
                         cex = labels.cex * 10*(yb-ybottom), 
                         col = labels.col[start:(jj-1)],
                         srt = 90)
        res <- c(res,c(xl,yb,xr,yt))
        start <- jj
        left <- right
      }
      else{
        break
      }
    }
    res
  }
  else{
    start <- 1
    left <- .0
    res <- c()
    while (start <= n) {
      if(!is.na(x[start,1])){
        jj <- start
        while(!is.na(x[jj,1]) & x[jj,1]==x[start,1]){
          if(jj<n){ jj <- jj+1 }
          else{ break }
        }
        if(jj==n){ jj <- jj+1 }
        right <- left + (jj-start)/n*(xright - xleft)
        ybottomupdate <- ybottom + (ytop-ybottom)/dim(x)[2]
        lambda <- .2
        xl <- xleft+left+lambda*(right-left)
        yb <- ybottom+lambda*(ybottomupdate-ybottom)
        xr <- xleft+left+(1-lambda)*(right-left)
        yt <- ybottom+(1-lambda)*(ybottomupdate-ybottom)
        graphics::rect(xl,yb,xr,yt,col = col)
        if(labels.plot)
          graphics::text(x = seq(from = xl, to = xr,length.out = jj-start), 
                         y = rep((yb+ybottom)/2, jj-start), 
                         labels = labels[start:(jj-1)],
                         cex = labels.cex * 10*(yb-ybottom), 
                         col = labels.col[start:(jj-1)],
                         srt = 90)
        res <- c(res, c(xl,yb,xr,yt))
        rectangles <- plotClusterTreeHelper(x = as.matrix(x[start:(jj-1),2:dim(x)[2]]),
                                            xleft = xleft+left, 
                                            ybottom = ybottomupdate, 
                                            xright = xleft+right, 
                                            ytop = ytop, 
                                            labels.plot = labels.plot, 
                                            labels = labels[start:(jj-1)],
                                            labels.cex = labels.cex,
                                            labels.col = labels.col[start:(jj-1)],
                                            col = col)
        if(length(rectangles)>=4){
          for (ii in 1:(length(rectangles)/4)) {
            graphics::segments(x0 = (xl+xr)/2,
                               y0 = yt,
                               x1 = (rectangles[ii*4-3]+rectangles[ii*4-1])/2,
                               y1 = rectangles[ii*4-2])
          }
        }
        start <- jj
        left <- right
      }
      else{
        break
      }
    }
    invisible(res)
  }
}


#' transform a clusterTree object into a clusterTreeDensityPlot object
#' clusterTreeDensityPlot object is a unified data structure for recording information
#' related to plotting clusterTree object
#' @param x a clusterTree object
#' @return a clusterTreeDensityPlot object
#' @return clusterTreeDensityPlot is a nested structure which record information for 
#' @return plotting clusterTree 
#' \item{tree}{a matrix, tree attribute of clusterTree object that remains to be plotted}
#' \item{rectangles}{a list of rectangles, each element record left,bottom,right,top coordinates of the rectangle}
#' \item{labels}{a list of list of labels' index, each element record labels' index for corresponding rectangle}
#' \item{lines}{a list of list of lines, each element record lines' start and end point}
#' \item{is.runts}{a list of boolean variable, each element record whether that rectangle is runts or not}
#' \item{subtree}{a list of subtrees, each element represents the nested subtree of corresponding rectangle}
#' @examples
#' data <- rbind(matrix(rnorm(100, mean = 10, sd = 2), nrow = 50),
#'               matrix(rnorm(100, mean = 0, sd = 1), nrow = 50),
#'               matrix(rnorm(100, mean = -10, sd = 3), nrow = 50)
#'               )
#' clustering1<-stats::hclust(dist(data),method='single')
#' clustering2<-kmeans(data,centers=3)
#' clustering3<-dbscan::dbscan(data,eps=.1)
#' clusterTree <- combineClusterings(clustering1,clustering2,clustering3)
#' res <- clusterTreeToClusterTreeDensityPlot(clusterTree)
#' @export 
clusterTreeToClusterTreeDensityPlot <- function(x){
  order <- reOrderClusterTreeMatrix(x$tree)
  orderedTree <- as.matrix(x$tree[order,])
  labels_index <- c(1:nrow(x$tree))
  labels_index <- labels_index[order]
  res <- clusterTreeToClusterTreeDensityPlotHelper(x = orderedTree, xleft = 0, ybottom = 0, 
    xright = 1, ytop = 1, labels_index = labels_index)
  class(res) <- c("clusterTreeDensityPlot", class(res))
  res
}


clusterTreeToClusterTreeDensityPlotHelper <- function(x, xleft, ybottom, xright, ytop, labels_index){
  n <- nrow(x)    
  if(ncol(x)==1){
    start <- 1
    left <- .0
    rectangles <- list()
    labels <- list()
    lines <- list()
    is.runts <- list()
    subtree <- list()
    while(start <= n){
        jj <- start
        while ((!is.na(x[start,1]) & !is.na(x[jj,1]) & x[jj,1]==x[start,1]) | (is.na(x[start,1]) & is.na(x[jj,1]))) {
          if(jj<n){ jj <- jj+1 }
          else{ break }
        }
        if(jj==n){ jj <- jj+1 }
        right <- left + (jj-start)/n*(xright - xleft)
        lambda <- .2
        xl <- xleft+left+lambda*(right-left)
        yb <- ybottom+lambda*(ytop-ybottom)
        xr <- xleft+left+(1-lambda)*(right-left)
        yt <- ybottom+(1-lambda)*(ytop-ybottom)
        rectangles[[length(rectangles)+1]] <- list(xleft = xl, ybottom = yb, xright = xr, ytop = yt)
        labels_x <- seq(from = xl, to = xr, length.out = jj-start) 
        labels_y <- (yb+ybottom)/2 
        labels_text <- labels_index[start:(jj-1)]
        labels_tmp <- list()
        for (kk in c(1:(jj-start))) {
          labels_tmp[[kk]] <- list(x=labels_x[kk], y=labels_y, text = labels_text[kk])
        }
        labels[[length(labels)+1]] <- labels_tmp
        lines[[length(lines)+1]] <- list()
        is.runts[[length(is.runts)+1]] <- is.na(x[start,1])
        start <- jj
        left <- right
      }
    list(tree = x, rectangles = rectangles, labels = labels, lines = lines, 
      is.runts = is.runts, subtree = subtree, n_branches = length(rectangles))
  }
  else{
    start <- 1
    left <- .0
    rectangles <- list()
    labels <- list()
    lines <- list()
    is.runts <- list()
    subtree <- list()
    while (start <= n) {
        jj <- start
        while((!is.na(x[start,1]) & !is.na(x[jj,1]) & x[jj,1]==x[start,1]) | (is.na(x[start,1]) & is.na(x[jj,1]))){
          if(jj<n){ jj <- jj+1 }
          else{ break }
        }
        if(jj==n){ jj <- jj+1 }
        right <- left + (jj-start)/n*(xright - xleft)
        ybottomupdate <- ybottom + (ytop-ybottom)/dim(x)[2]
        lambda <- .2
        xl <- xleft+left+lambda*(right-left)
        yb <- ybottom+lambda*(ybottomupdate-ybottom)
        xr <- xleft+left+(1-lambda)*(right-left)
        yt <- ybottom+(1-lambda)*(ybottomupdate-ybottom)
        rectangles[[length(rectangles)+1]] <- list(xleft = xl, ybottom = yb, xright = xr, ytop = yt)
        labels_x <- seq(from = xl, to = xr, length.out = jj-start) 
        labels_y <- (yb+ybottom)/2 
        labels_text <- labels_index[start:(jj-1)]
        labels_tmp <- list()
        for (kk in c(1:(jj-start))) {
          labels_tmp[[kk]] <- list(x=labels_x[kk], y=labels_y, text = labels_text[kk])
        }
        labels[[length(labels)+1]] <- labels_tmp
        subtree_tmp <- clusterTreeToClusterTreeDensityPlotHelper(x = as.matrix(x[start:(jj-1),2:dim(x)[2]]),
                                            xleft = xleft+left, 
                                            ybottom = ybottomupdate, 
                                            xright = xleft+right, 
                                            ytop = ytop, 
                                            labels_index = labels_index[start:(jj-1)])
        subtree[[length(subtree)+1]] <- subtree_tmp
        lines_tmp <- list()
        rectangle_tmp <- subtree_tmp$rectangles
        for (kk in 1:length(rectangle_tmp)) {
          cur_rectangle <- rectangle_tmp[[kk]]
          lines_tmp[[length(lines_tmp)+1]] <- list(start = list(x = (xl+xr)/2.0, y = yt), end = list(x = (cur_rectangle$xleft + cur_rectangle$xright)/2.0, y = cur_rectangle$ybottom))
        }
        lines[[length(lines)+1]] <- lines_tmp
        is.runts[[length(is.runts)+1]] <- is.na(x[start,1])
        start <- jj
        left <- right
      }
    list(tree = x, rectangles = rectangles, labels = labels, lines = lines, 
      is.runts = is.runts, subtree = subtree, n_branches = length(rectangles))
  }
}

#' plot a clusterTreeDensityPlot object
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
#' clustering3<-dbscan::dbscan(data,eps=.1)
#' res <- combineClusterings(clustering1,clustering2,clustering3)
#' res2 <- clusterTreeToClusterTreeDensityPlot(res)
#' plot(res2)
#' @export 
plot.clusterTreeDensityPlot <- function(x, y = NULL, labels = NULL, axes = TRUE, frame.plot = FALSE, ann = TRUE, 
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
    labels <- paste("object", 1:nrow(x$tree))
  }else{
    if(is.null(labels)){
      labels.plot <- FALSE
      labels <- paste("object", 1:nrow(x$tree))
    }else{
      if(!is.vector(labels)){
        stop("labels must be a logical or vector!")
      }else{
        if(length(labels) != nrow(x$tree)){
          stop('length of labels and rows of clusterTree object does not match!')
        }else{
          labels.plot <- TRUE
        }
      }
    }
  }
  ###
  ###
  # plot clusterTreeDensityPlot recursively
  graphics::plot.new()
  plotClusterTreeDensityPlotHelper(x = x, labels.plot = labels.plot, labels = labels, col = col, 
                                  labels.cex = 2*labels.cex/ncol(x$tree), labels.col = if (length (labels.col) == length(labels)){
                                                           labels.col[order]}else{labels.col})
  ###
  ###
  if(frame.plot){
    graphics::box(...)
  }
  if(ann){
    if(is.null(sub))
      sub <- "density plot"
    graphics::title(main = main, sub = sub, xlab = xlab, ylab = ylab, ...)
  }
  if(axes)
    height <- as.double(ncol(x$tree))
    graphics::axis(2,(1/2/height)*seq(1,2*height,by = 2),labels = c(1:height))
}

#' helper function to plot a clusterTreeDensityPlot object
#' @param x a clusterTreeDensityPlot object
#' @param labels.plot whether to plot labels
#' @param labels of data points
#' @param col color of rectangles
#' @param labels.cex cex of labels
#' @param labels.col color of labels 
plotClusterTreeDensityPlotHelper <- function(x, labels.plot, labels, col, 
                                  labels.cex, labels.col){
  for (ii in c(1:x$n_branches)) {
    if(!x$is.runts[[ii]]){
      graphics::rect(x$rectangles[[ii]]$xleft, x$rectangles[[ii]]$ybottom,
        x$rectangles[[ii]]$xright, x$rectangles[[ii]]$ytop, col = col)
      if(length(x$lines[[ii]])>0){
        for (jj in c(1:length(x$lines[[ii]]))) {
          if(!x$subtree[[ii]]$is.runts[[jj]]){
            graphics::segments(x0 = x$lines[[ii]][[jj]]$start$x,
              y0 = x$lines[[ii]][[jj]]$start$y,
              x1 = x$lines[[ii]][[jj]]$end$x,
              y1 = x$lines[[ii]][[jj]]$end$y)
          }
        }
      }
      if(labels.plot){
        if(length(x$labels[[ii]])>0){
          for (jj in c(1:length(x$labels[[ii]]))) {
            graphics::text(x = x$labels[[ii]][[jj]]$x,
                           y = x$labels[[ii]][[jj]]$y,
                           labels = labels[x$labels[[ii]][[jj]]$text],
                           cex = labels.cex, 
                           col = labels.col[x$labels[[ii]][[jj]]$text],
                           srt = 90)
          }
        }
      }
    }
  }
  if(length(x$subtree)>0){
    for (ii in c(1:length(x$subtree))) {
      plotClusterTreeDensityPlotHelper(x$subtree[[ii]], labels.plot = labels.plot,
                                       labels = labels, col = col, 
                                       labels.cex = labels.cex, labels.col = labels.col)
    }
  }
}

