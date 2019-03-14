
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
#' @export
plotClusterTreeHelper <- function(x, xleft, ybottom, xright, ytop, 
                                  labels.plot, labels, col, 
                                  labels.cex = 1, labels.col){
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
