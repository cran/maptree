# maptree package
#   for graphing and mapping of hierarchical clustering and
#   regression trees
# denis white, us epa, 6 June 2001, version 1.1-1
#
# function calls
#
############################################################
#
# draw.clust <- function (cluster, ps=par("ps"), size=10, 
#     col=NULL)
#
# draw.tree <- function (tree, ps=par("ps"), size=10, col=NULL, 
#     nodeinfo=FALSE, units="units", cases="obs", digits=0 )
#
# group.clust <- function (cluster, k=NULL, h=NULL)
#
# group.tree <- function (tree)
#
# map.groups <- function (pts, group, pch=19, size=5, col=NULL,
#     border=NULL)
#
# map.key <- function (x, y, lables=NULL, new=FALSE, size=5, 
#     ps=par("ps"), pch=19, head="", sep=1, col=NULL)
#
# ngon <- function (xydc, n=4, angle=0, type=1)
#
# prune.clust <- function (cluster, k=NULL, h=NULL)
#
# prune.Rpart <- function (tree, cp=NULL, best=NULL) 

############################################################

draw.clust <- function (cluster, ps=par("ps"), size=10, 
  col=NULL)
  # cluster is class hclust or twins
  # ps is par parameter, pointsize of text
  # size is in mm of colored square at observations
  # if col is NULL, use rainbow(),
  #   else if col="gray" or "grey", use gray()
  # returned value is col or generated colors
{
  if (class(cluster) == "hclust") clust <- cluster
  else if (inherits(cluster,"twins"))
    clust <- as.hclust (cluster)
  else 
    stop("draw: input not hclust or twins")
  merg <- clust$merge
  nmerg <- nrow(merg)
  if (nmerg<2) stop("draw: <3 clusters")
  hite <- clust$height
  cord <- order(clust$order)
  xmax <- nrow (merg) + 1
  ymax <- max (hite)
  pinx <- par ("pin")[1]
  piny <- par ("pin")[2]
  xmin <- 1
  box <- size/25.4
  xscale <- (xmax-xmin)/pinx
  xbh <- xscale*box/2
  tail <- 0.25
  yscale <- ymax/(piny - tail)
  ytail <- yscale*tail
  ybx <- yscale*box
  ymin <- - ytail
  xf <- 0.1 * (xmax-xmin)
  yf <- 0.1 * (ymax-ymin)
  x1 <- xmin - xf
  x2 <- xmax + xf
  y1 <- ymin - yf
  y2 <- ymax + yf
  plot (c(x1,x2),c(y1,y2),type="n",axes=FALSE,xlab="",ylab="")
  oldps <- par ("ps")
  par (ps=ps)
  if (is.null(col)) kol <- rainbow (xmax)
  else if (col=="gray" | col=="grey") kol <- 
    gray (seq(0.8,0.2,length=xmax))
  else kol <- col
  xmean <- rep (0, nmerg)
  i <- 1
  while (any(xmean == 0)) {
    if (xmean[i] == 0) {
      a <- merg[i,1]
      b <- merg[i,2]
      if (a<0) x1 <- cord[-a] else x1 <- xmean[a]
      if (b<0) x2 <- cord[-b] else x2 <- xmean[b]
      if (x1 != 0 && x2 != 0) xmean[i] <- mean(c(x1,x2))
      }
    i <- i + 1
    if (i > nmerg) i <- 1
    }
  for (i in 1:nmerg) {
    a <- merg[i,1]
    b <- merg[i,2]
    y2 <- hite[i]
    if (a > 0) {
      x1 <- xmean[a]
      y1a <- hite[a]
      }
    else {
      x1 <- cord[-a]
      y1a <- y2 - ytail
      px <- c(x1-xbh,x1+xbh,x1+xbh,x1-xbh,x1-xbh)
      py <- c(y1a-ybx,y1a-ybx,y1a,y1a,y1a-ybx)
      polygon(px,py,col=kol[x1],border=0)
      text.default (x1,y1a-(ybx/2),as.character(-a))
      }
    if (b > 0) {
      x2 <- xmean[b]
      y1b <- hite[b]
      }
    else {
      x2 <- cord[-b]
      y1b <- y2 - ytail
      px <- c(x2-xbh,x2+xbh,x2+xbh,x2-xbh,x2-xbh)
      py <- c(y1b-ybx,y1b-ybx,y1b,y1b,y1b-ybx)
      polygon(px,py,col=kol[x2],border=0)
      text.default (x2,y1b-(ybx/2),as.character(-b))
      }
    lines (c(x1,x2),c(y2,y2))
    lines (c(x1,x1),c(y1a,y2))
    lines (c(x2,x2),c(y1b,y2))
    }
  par (ps=oldps)
  invisible (kol)
}

############################################################

draw.tree <- function (tree, ps=par("ps"), size=10, col=NULL, 
  nodeinfo=FALSE, units="units", cases="obs", digits=0 )
  # tree is object of class tree
  # ps is par parameter, pointsize of text
  # if size=0, draw terminal symbol at leaves else a
  #   colored square with sides of length size in mm
  # if col is NULL, use rainbow(),
  #   else if col="gray" or "grey", use gray()
  # if nodeinfo=TRUE, add a line at each node with mean value
  #   of response, number of observations, and percent
  #   deviance explained (or classified correct)
  # units are for mean value of response (if regression tree)
  # cases are for observations
  # digits are rounding for response
  # returned value is col or generated colors
{
  rtree <- length (attr (tree, "ylevels")) == 0
  tframe <- tree$frame
  node <- as.numeric(row.names(tframe))
  depth <- tree.depth(node)
  maxdepth <- max(depth)
  x <-  - depth
  y <- x
  leaves <- tframe$var == "<leaf>"
  x[leaves] <- seq(sum(leaves))
  depth <- split(seq(node)[!leaves], depth[!leaves])
  parent <- match(node %/% 2, node)
  left.child <- match(node * 2, node)
  right.child <- match(node * 2 + 1, node)
  for(i in rev(depth))
    x[i] <- 0.5 * (x[left.child[i]] + x[right.child[i]])
  nleaves <- sum(leaves)
  nnodes <- length(node)
  if (rtree) {
    dev <- tframe$dev
    pcor <- rep (0, nnodes)
    for (i in 1:nnodes)
      if (! leaves[i]) {
        l <- dev[node == (node[i]*2)]
        r <- dev[node == (node[i]*2+1)]
        pcor[i] <- dev[i] - l - r
      }
    pcor <- round (pcor/dev[1],3)*100
    }
  else {
    crate <- rep (0, nnodes)
    trate <- 0
    for (i in 1:nnodes) {
      yval <- tframe$yval[i]
      string <- paste('tframe$yprob[,"',as.character(yval),'"]',
          sep="")
      crate[i] <- eval(parse(text=string))[i]
      if (leaves[i]) trate <- trate + tframe$n[i] * crate[i]
      }
    crate <- round (crate,3)*100
    trate <- round (trate/tframe$n[1],3)*100
    }
  if (is.null(col)) kol <- rainbow (nleaves)
  else if (col=="gray" | col=="grey") kol <- 
    gray (seq(0.8,0.2,length=nleaves))
  else kol <- col
  oldps <- par ("ps")
  par (ps=ps)
  xmax <- max(x)
  xmin <- min(x)
  ymax <- max(y)
  ymin <- min(y)
  pinx <- par ("pin")[1]
  piny <- par ("pin")[2]
  xscale <- (xmax - xmin)/pinx
  box <- size / 25.4
  if (box == 0) xbh <- xscale * 0.2
  else xbh <- xscale * box/2
  chr <- par("csi")
  tail <- box + chr
  yscale <- (ymax - ymin)/(piny - tail)
  ytail <- yscale * tail
  if (box == 0) ybx <- yscale * 0.2
  else ybx <- yscale * box
  ychr <- yscale * chr
  ymin <- ymin - ytail
  xf <- 0.1 * (xmax - xmin)
  yf <- 0.1 * (ymax - ymin)
  x1 <- xmin - xf
  x2 <- xmax + xf
  y1 <- ymin - yf
  y2 <- ymax + yf
  plot (c(x1,x2),c(y1,y2),type="n",axes=FALSE,xlab="",ylab="")
  string <- paste (as.character(tframe$var[1]), "<>", 
    substring(as.character(tframe$splits[1,1]),2))
  text.default (x[1], y[1], string)
  if (nodeinfo) {
    n <- tframe$n[1]
    if (rtree) {
      z <- round(tframe$yval[1], digits)
      r <- pcor[1]
      string <- paste (z," ",units,"; ",n," ",cases,"; ",
        r,"%",sep="")
      }
    else {
      z <- tframe$yval[1]
      r <- crate[1]
      string <- paste (z,"; ",n," ",cases,"; ",r,"%",
        sep="")
      }
    text.default (x[1], y[1]-ychr, string)
    }

  for (i in 2:nnodes) {
    ytop <- ychr * (as.integer(nodeinfo)+1)
    if (y[i] < y[i-1])
      lines(c(x[i-1], x[i]), c(y[i-1]-ytop, y[i]+ychr))
    else
      lines(c(x[parent[i]], x[i]), 
        c(y[parent[i]]-ytop, y[i]+ychr))
    if(! leaves[i]) {
      string <- paste (as.character(tframe$var[i]), "<>", 
        substring(as.character(tframe$splits[i,1]),2))
      text.default (x[i], y[i], string)
      if (nodeinfo) {
        n <- tframe$n[i]
        if (rtree) {
          z <- round(tframe$yval[i], digits)
          r <- pcor[i]
          string <- paste (z," ",units,"; ",n," ",cases,"; ",
            r,"%",sep="")
          }
        else {
          z <- tframe$yval[i]
          r <- crate[i]
          string <- paste (z,"; ",n," ",cases,"; ",r,"%",
            sep="")
          }
        text.default (x[i], y[i]-ychr, string)
        }
      }
    else {
      if (box == 0) {
        lines (c(x[i], x[i]), c(y[i], y[i]+ychr))
        lines (c(x[i]-xbh, x[i]+xbh), c(y[i], y[i]))
        }
      else {
        x1 <- x[i] - xbh
        x2 <- x[i] + xbh
        y1 <- y[i] - ybx + ychr
        y2 <- y[i] + ychr
        polygon(c(x1, x2, x2, x1, x1), c(y1, y1, y2, y2, y1), 
          col = kol[x[i]], border=0)
        # lines (c(x1, x2, x2, x1, x1), c(y1, y1, y2, y2, y1))
        }
      if (rtree) {
        z <- round(tframe$yval[i], digits)
        text.default(x[i], y[i]-ybx, paste(z,units,sep=" "))
        }
      else text.default(x[i], y[i]-ybx, tframe$yval[i])
      n <- tframe$n[i]
      text.default(x[i], y[i]-ybx-ychr, paste(n,cases,sep=" "))
      if (box != 0)
        text.default (x[i], y[i]-ychr/4, as.character(x[i]))
      }
    }
  if (nodeinfo) {
    if (rtree) string <- paste("Total deviance explained =",
      sum(pcor),"%")
    else string <- paste("Total classified correct =",trate,"%")
    par (ps=1.2*ps)
    if (box == 0) text.default (mean(x),ymin-3*ychr,string)
    else text.default (mean(x),ymin-ybx,string)
    }
  par (ps=oldps)
  invisible (kol)
}

############################################################

group.clust <- function (cluster, k=NULL, h=NULL)
  # alternative to cutree that orders groups from left to
  # right in draw order
  # cluster is class hclust or twins
  # k is desired number of groups
  # h is height at which to cut for grouping
  # needs k or h, k takes precedence
  # returns vector of membership
{
  if (is.null(h) && is.null(k)) return (cluster$order)
  if (!is.null(h) && h>max(cluster$height)) 
    stop("group: h>max(height)")
  if (!is.null(k) && (k==1 || k>length(cluster$height)))
    stop("group: k==1 || k=>nobs")
  if (class(cluster) == "hclust") clust <- cluster
  else if (inherits(cluster,"twins"))
    clust <- as.hclust (cluster)
  else
    stop("group: input not hclust or twins")
  merg <- clust$merge
  hite <- clust$height
  ordr <- clust$order
  nmerg <- nrow(merg)
  group <- rep (0, nmerg+1)
  if (!is.null(k)) keep <- rev(order(hite))[1:(k-1)]
  else keep <- seq(nmerg)[hite > h]

  mark.group <- function (node, grup)
  {
    a <- merg[node,1]
    b <- merg[node,2]
    if (a < 0) group[-a] <<- grup
    else mark.group (a, grup)
    if (b < 0) group[-b] <<- grup
    else mark.group (b, grup)
    invisible()
  }

  grup <- 0
  find.group <- function (node)
  {
    a <- merg[node,1]
    b <- merg[node,2]
    if (match(a,keep,0) != 0) find.group (a)
    else {
      grup <<- grup + 1
      if (a > 0) mark.group (a, grup)
      else group[-a] <<- grup
      }
    if (match(b,keep,0) != 0) find.group (b)
    else {
      grup <<- grup + 1
      if (b > 0) mark.group (b, grup)
      else group[-b] <<- grup
      }
    invisible()
  }

  find.group (length(hite))
  grup <- match(grup,unique(grup[clust$order]))
  invisible (group)
}

############################################################

group.tree <- function (tree)
  # alternative to tree$where that orders groups from left
  # to right in draw order
{
  group <- match (tree$where, sort(unique(tree$where)))
  names(group) <- names(tree$where)
  invisible (group)
}

############################################################

map.groups <- function (pts, group, pch=19, size=5, col=NULL,
  border=NULL)
  # pts must have components "x" and "y";
  # group is vector of length of number of cases (either
  #   polygons or points) and indexes colors;
  #   if pts are for polygons, then names(group) must match
  #   with pts$x[is.na(pts$y)]
  # if nrow(pts) != length(group) then map with polygon, 
  #   else if pch < 100 map with points, 
  #   else map with ngon (..., n=pch-100)
  # size is in mm, only for point symbol
  # if col is NULL, use rainbow(),
  #   else if col="gray" or "grey", use gray()
  # if border is NULL, use fill colors (col),
  #   else the specified color(s)
{
  n <- colnames (pts)
  if (! any (n == "x"))
    stop ("map.groups: pts has no $x")
  else if (! any (n == "y"))
    stop ("map.groups: pts has no $y")
  if (nrow(pts) != length(group)) 
    if (is.null (names (group)))
      stop ("map.groups: group has no names")
  dx <- diff(range(pts$x,na.rm=TRUE))
  dy <- diff(range(pts$y,na.rm=TRUE))
  px <- par("pin")[1]
  py <- par("pin")[2]
  if (dx/dy > px/py) py <- px * dy/dx else px <- py * dx/dy
  par(pin = c(px, py))
  plot(range(pts$x,na.rm=TRUE),range(pts$y,na.rm=TRUE),type="n",
    axes=FALSE,xlab="",ylab="")
  ncol <- length (unique (group))
  if (is.null(col)) fkol <- rainbow (ncol)
  else if (col=="gray" | col=="grey") fkol <- 
    gray (seq(0.8,0.2,length=ncol))
  else fkol <- rep (col, ncol)
  if (is.null(border)) bkol <- fkol
  else bkol <- rep (border, ncol)
  if (nrow(pts) != length(group)) {
    i <- pts$x[is.na(pts$y)]
    j <- match (i, as.integer(names(group)))
    polygon (pts$x, pts$y, lwd=0.1,
      col=fkol[group[j]], border=bkol[group[j]])
    }
  else if (pch < 100 | mode(pch) == "character") points (pts$x, 
    pts$y,col=fkol[group],pch=pch,cex=1.5*size/25.4/par("cin")[1])
  else apply(data.frame(x=pts$x,y=pts$y,d=size,c=I(fkol[group])),
    1,ngon,n=pch-100,type=1)
  invisible (fkol)
}

############################################################

map.key <- function (x, y, lables=NULL, new=FALSE, size=5,
  ps=par("ps"), pch=19, head="", sep=1, col=NULL)
  # x,y are lower left coordinates of key 
  #   in proportional units (0-1)
  # lables is vector of labels for classes, or if NULL,
  #   then integers 1:length(col), or "1"
  # if new=TRUE, call plot
  # size is (diameter) of polygon symbol in mm
  # ps is par parameter, pointsize of text
  # if pch < 100, use points for symbol, 
  #   else ngon (..., n=pch-100)
  # head is text heading for key
  # sep is separation in mm between adjacent symbols
  #   if sep=0 assume continuous scale and use gon=4
  #   and put lables at breaks between squares
  # if col is NULL, use rainbow(),
  #   else if col="gray" or "grey", use gray()
  # returned value is col or generated colors
{
  if (is.null(lables))
    if (is.null(col)) lables <- as.vector("1")
    else lables <- as.character(seq(length(col)-1))
  nsym <- length(lables)
  if (sep == 0) nsym <- nsym - 1
  if (is.null(col)) kol <- rainbow (nsym)
  else if (col=="gray" | col=="grey") kol <- 
    gray (seq(0.8,0.2,length=nsym))
  else kol <- col
  if (new)
    plot(c(0,1),c(0,1),type="n",axes=FALSE,xlab="",ylab="")
  oldps <- par("ps")
  par (ps=ps)
  oldadj <- par ("adj")
  par (adj=0)
  u <- par ("usr")
  ex <- par ("pin")[1]
  ey <- par ("pin")[2]
  uxr <- u[2] - u[1]
  uyr <- u[4] - u[3]
  halfx <- size/50.8
  xstep <- halfx + 0.05
  ystep <- (size + sep/2.0)/25.4
  px <- x * uxr + u[1]
  py <- y * uyr + u[3]
  hx <- halfx * uxr / ex
  dx <- xstep * uxr / ex
  dy <- ystep * uyr / ey
  qx <- px
  qy <- py - dy
  if (sep == 0) {
    for (i in 1:nsym) {
      qy <- qy + dy
      ngon (c(qx, qy, size=size, col=kol[i]), n=4, type=3)
      text (qx+dx, qy - dy/2, lables[i]) }
    text (qx+dx, qy + dy/2, lables[nsym+1])
    }
  else
    for (i in 1:nsym) {
      qy <- qy + dy
      if (pch < 100  | mode(pch) == "character") points (qx,
        qy, col=kol[i], pch=pch, cex=1.5*size/25.4/par("cin")[1])
      else ngon (c(qx, qy, size=size, col=kol[i]), n=pch-100, type=1)
      text (qx+dx, qy, lables[i]) }
  if (length(head)>0)
    {
    qy <- qy + (dy * length (grep ("$", head)))
    #qy <- qy + 0.2 * dy
    if (sep == 0) qy <- qy + 0.5 * dy
    text (qx-hx, qy, head)
    }
  par (ps=oldps, adj=oldadj)
  invisible (kol)
}

############################################################

ngon <- function (xydc, n=4, angle=0, type=1)
  # draw or fill regular polygon
  # xydc a four element vector with
  #   x and y of center, d diameter in mm, and c color
  # n number of sides of polygon, n>8 => circle
  #   if n odd, vertex at (0,d/2), else midpoint of side
  # angle is in degrees by which to rotate the figure
  # type=1 => interior filled, type=2 => edge
  # type=3 => both
{
      # scale factors for d based on n (ignoring angle)
      # n = 3, s = (2 + sqrt(3)) / 4     = 0.9330127
      # n = 4, s = 1 / sqrt(2)           = 0.7071068
      # n = 5, s = (1 + cos(.2*pi)) / 2  = 0.9045085
      # n = 6, s = sqrt(3) / 2           = 0.8660254
  u <- par("usr")
  p <- par("pin")
  d <- as.numeric(xydc[3])
  inch <- d/25.4
  s <- 1
  switch (n, stop ("ngon: n=1"), 
             stop ("ngon: n=2"),
             s <- 0.9330127,
             s <- 0.7071068,
             s <- 0.9045085,
             s <- 0.8660254)
  inch <- inch / s
  rad <- inch*((u[2]-u[1])/p[1])/2
  ys <- inch*((u[4]-u[3])/p[2])/2/rad
  if (n > 8) n <- d*4 + 1
  th <- pi*2/n
  costh <- cos (th)
  sinth <- sin (th)
  x <- y <- rep (0,n+1)
  if (n %% 2) {
    x0 <- 0
    y0 <- rad
    }
  else {
    x0 <- -rad*sin(th/2)
    y0 <-  rad*cos(th/2)
    }
  a <- pi*angle/180
  x[1] <- x0*cos(a) - y0*sin(a)
  y[1] <- x0*sin(a) + y0*cos(a)
  for (i in 2:(n+1)) {
    xl <- x[i-1]
    yl <- y[i-1]
    x[i] <- xl*costh - yl*sinth
    y[i] <- xl*sinth + yl*costh
    }
  x <- x + as.numeric(xydc[1])
  y <- y*ys + as.numeric(xydc[2])
  if (type %% 2) polygon (x,y,col=xydc[4],border=0)
  if (type %/% 2) lines (x,y,col=xydc[4])
  invisible ()
}

############################################################

prune.clust <- function (cluster, k=NULL, h=NULL)
  # analogous to prune.tree
  # cluster is class hclust or twins
  # k is desired number of groups
  # h is height at which to cut for grouping
  # needs k or h, k takes precedence
  # returns pruned cluster
{
  if (is.null(h) && is.null(k))
    stop ("prune: both h=NULL,k=NULL")
  if (!is.null(h) && h>max(cluster$height)) 
    stop("prune: h>max(height)")
  if (!is.null(k) && (k==1 || k>length(cluster$height)))
    stop("prune: k==1 || k=>nobs")
  if (class(cluster) == "hclust") clust <- cluster
  else if (inherits(cluster,"twins"))
    clust <- as.hclust (cluster)
  else
    stop("prune: input not hclust or twins")
  merg <- clust$merge
  hite <- clust$height
  nmerg <- nrow(merg)
  if (!is.null(k)) keep <- rev(order(hite))[1:(k-1)]
  else keep <- seq(nmerg)[hite > h]
  numerg <- matrix(0,nrow=length(keep),ncol=2)
  nuhite <- rep (0, length(keep))

  leaf <- 0
  node <- 0
  trim.clust <- function (oldnode)
  {
    a <- merg[oldnode,1]
    b <- merg[oldnode,2]
    if (match(a,keep,0) != 0) l <- trim.clust(a)
    else {
      leaf <<- leaf + 1
      l <- -leaf
      }
    if (match(b,keep,0) != 0) r <- trim.clust(b)
    else {
      leaf <<- leaf + 1
      r <- -leaf
      }
    node <<- node + 1
    numerg[node,] <<- c(l,r)
    nuhite[node] <<- hite[oldnode]
    return (node)
  }

  trim.clust(length(hite))
#  trim.clust(match(max(hite),hite))

  numerg <- matrix(as.integer(numerg),nrow=length(numerg)/2,ncol=2)
  nuhite <- as.double(nuhite - min(nuhite))
  nuordr <- as.double(seq(nrow(numerg)+1))
  nulabl <- as.character(nuordr)
  l <- list()
  l[[1]] <- numerg
  l[[2]] <- nuhite
  l[[3]] <- nuordr
  l[[4]] <- nulabl
  l[[5]] <- clust$method
  l[[6]] <- clust$call
  l[[7]] <- clust$dist.method
  names(l) <- names(clust)
  class(l) <- class(clust)
  l
}

############################################################

prune.Rpart <- function (tree, cp=NULL, best=NULL) 
  # modification to original prune.rpart to add best
{
  ff <- tree$frame
  id <- as.integer(row.names(ff))
  if (is.null (cp)) {
    cp <- ff$complexity
    cp <- cp[rev (order (cp))[best]]
    }
  toss <- id[ff$complexity <= cp & ff$var != "<leaf>"]
  if (length(toss) == 0) 
    return(tree)
  newx <- snip.rpart(tree, toss)
  temp <- pmax(tree$cptable[, 1], cp)
  keep <- match(unique(temp), temp)
  newx$cptable <- tree$cptable[keep, ]
  newx$cptable[max(keep), 1] <- cp
  newx
}
