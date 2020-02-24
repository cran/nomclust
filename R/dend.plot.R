#' Visualization of Cluster Hierarchy using a Dendrogram
#' 
#' @description The function \code{dend.plot()} visualizes the hierarchy of clusters using a dendrogram. The function also enables a user to mark the individual clusters with colors. 
#' The number of displayed clusters can be defined either by a user or by one of the five evaluation criteria.
#' 
#' @param x An output of the \code{nomclust()} or \code{nomprox()} functions containing the \code{dend} component.
#' 
#' @param clusters Either a \emph{numeric} value or a \emph{character} string with the name of the evaluation criterion expressing the number of displayed clusters in a dendrogram. The following evaluation criteria can be used: \code{"AIC"}, \code{"BIC"}, \code{"BK"}, \code{"PSFE"} and \code{"PSFM"}.
#' 
#' @param style A \emph{character} string or a \emph{vector} of colors defines a graphical style of the produced plots. There are two predefined styles in the \bold{nomclust} package, namely \code{"greys"} and \code{"dark"}, but a custom color scheme can be set by a user as a vector of a length four.
#'
#' @param colorful A \emph{logical} argument specifying if the output will be colorful or black and white.
#' 
#' @param clu.col An optional \emph{vector} of colors which allows a researcher to apply user-defined colors for displayed (marked) clusters in a dendrogram.
#' 
#' @param main A \emph{character} string with the chart title.
#' 
#' @param ac A \emph{logical} argument indicating if an agglomerative coefficient will be present in the output.
#' 
#' @param ... Other graphical arguments compatible with the generic \code{plot()} function.
#' 
#' @return The function returns a dendrogram describing the hierarchy of clusters that can help to identify the optimal number of clusters.
#' \cr
#'
#' @details The function can be applied to a \code{nomclust()} or \code{nomprox()} output containing the \code{dend} component. This component is not available when the optimization process is used.
#'
#' @seealso
#' \code{\link[nomclust]{eval.plot}}, \code{\link[nomclust]{nomclust}}, \code{\link[nomclust]{nomprox}}.
#' 
#' @author Jana Cibulkova and Zdenek Sulc. \cr Contact: \email{jana.cibulkova@@vse.cz}
#' 
#' @examples
#' # sample data
#' data(data20)
#' 
#' # creating an object with results of hierarchical clustering 
#' hca.object <- nomclust(data20, measure = "iof", eval = TRUE, opt = FALSE)
#' 
#' # a basic plot
#' dend.plot(hca.object)
#' 
#' # a dendrogram with color-coded clusters according to the BIC index
#' dend.plot(hca.object, clusters = "BIC", colorful = TRUE)
#' 
#' # using a dark style and specifying own colors in a solution with three clusters
#' dend.plot(hca.object, clusters = 3, style = "dark", clu.col = c("blue", "red", "green"))
#' 
#' # a black and white dendrogram
#' dend.plot(hca.object, clusters = 3, style = "dark", colorful = FALSE)
#' 
#' @export 


dend.plot <- function(x, clusters = "BIC", style = "greys", colorful = TRUE, clu.col = NA, main = "Dendrogram", ac = TRUE, ...) {
  #read input arguments
  if(style[1] == "dark"){style = rep("black",4)}
  if(style[1] == "greys"){    style = grey.colors(5)}
  
  if(typeof(x)=="list" & "dend" %in% names(x) & "opt" %in% names(x)){
    #create agnes object
    a.object=list()
    a.object$order=x$dend$order
    a.object$height=x$dend$height
    a.object$ac=x$dend$ac
    a.object$merge=x$dend$merge
    a.object$method=x$dend$method
    a.object$order.lab=x$dend$order.lab
    class(a.object) = c("agnes","twins")

    height=       sort(a.object$height)
    merge=        a.object$merge
    order=        a.object$order
    
    if(!is.numeric(clusters)){clusters=toupper(clusters)}
    
    if(!(clusters %in% c('BIC','AIC','PSFM','PSFE','BK') | if(is.numeric(clusters)){clusters%%1==0 & clusters>0}else{F})){
      stop('Choose the number of clusters to be visualized. It should be either a number or one of these criteria: AIC, BIC, BK, PSFM, PSFE.')
    }
    if(clusters=='BIC' ){clusters= x$opt$BIC}
    if(clusters=='AIC' ){clusters= x$opt$AIC}
    if(clusters=='PSFM'){clusters= x$opt$PSFM}
    if(clusters=='PSFE'){clusters= x$opt$PSFE}
    if(clusters=='BK'  ){clusters= x$opt$BK}
    
  }else{
    stop('Input argument x is missing or incorrect. Output from nomclust() or nomprox() functions is required.')

  }
  
  if(colorful==F){clu.col=rep('black',clusters)}
  if(colorful==T & length(clu.col)<=clusters){clu.col=rainbow(clusters)}

  #create plot
  pltree(a.object, hang=-1, lty=2, lwd=2, sub="", xlab='',axes=F, col.lab = style[3],col=style[1],main=main)
  mtext(paste(a.object$method,'linkage method'),col=style[1])
  if(ac == TRUE){
    mtext(paste("Agglomerative coefficient =",format(round(x$dend$ac,2),nsmall=2)),side=1,line=1)
  }
  h=seq(0,ceiling(max(height)*100)/100,length.out = 10)                   # seq for y axis
  axis(labels=F, side=2,las=2,col = style[3], at=h)                       # create y axis
  mtext(format(round((h*100)/100,2), nsmall = 2),side=2,las=2,col = style[1], at=h, line = 1) # add labels
  abline(h=h,lty=3,lwd=1,col=style[4])                                    # add horizontal lines
  
  #define helpful variables
  cutted=                    cutree(a.object,clusters) # cut the tree
  order.col=                 vector() # assign color for each observation
  for(i in 1:length(order)) {order.col[i]=cutted[order[i]]   }
  merge.backup=              merge    # for one-object clusters
  col.backup=                vector() # for one-object clusters
  starts.at.zero=           (merge<0) # check what vertical lines starts from 0
  cluster.height=            vector() # saves ancerstor's height for vertical lines to start from
  ancestors.col=             vector() # saves ancerstor's colors
  i=                         1        # start iterations
  
  
repeat{
  # get x coordinates
  step=merge[i,]
  x1=ifelse(starts.at.zero[i,1], which( order==(-step[1])),-merge[i,1] )
  x2=ifelse(starts.at.zero[i,2], which( order==(-step[2])),-merge[i,2] )
  
  #get y coordinates
  if(starts.at.zero[i,1]==T){
    y1.bottom=0
    y1.top=height[i]
    ancestors.col[x1]=order.col[x1]
  }else{
    y1.bottom=cluster.height[which(merge[1:i,1]==x1)] #start where the ancestor ended
    y1.top=height[i]
  }
  if(starts.at.zero[i,2]==T){
    y2.bottom=0
    y2.top=height[i]
    ancestors.col[x2]=order.col[x2]
  }else{
    y2.bottom=cluster.height[which(merge[1:i,1]==x2)] #start where the ancestor ended
    y2.top=height[i]
    
  }
  
  #update x coordinates of newly created cluster
  merge[which(merge==i)]=-mean(c(x1,x2))              #for future steps
  merge[i,]=c(mean(c(x1,x2)),mean(c(x1,x2)))          #keeping track of ancestors
  cluster.height[i]=max(c(y1.top,y2.top))
  
  #check for break
  if(ancestors.col[x1] != ancestors.col[x2]) break
  
  #draw segments
  col.backup=c(col.backup,ancestors.col[x1])
  segments(x1,y1.bottom,x1,y1.top,lwd=3, col=clu.col[ancestors.col[x1]])
  segments(x2,y2.bottom,x2,y2.top,lwd=3, col=clu.col[ancestors.col[x1]])
  segments(x1,y1.top,x2,y1.top,lwd=3, col=clu.col[ancestors.col[x1]])
  
  #move to next iteration
  i=i+1
  
}
  #work with one-object clusters
#  i=nrow(merge.backup)
#  i
  temp=as.vector(merge.backup[-c(1:(i-1)),])
  temp=-temp[which(temp<0)]
  temp=which(order %in% temp)
  
  if(length(temp)>0){
    col.temp=which(!(1:clusters %in% col.backup))
    col.temp
    clu.col
#    i=1
    for(i in 1:length(temp)){
      segments(temp[i],0,temp[i],y1.top,lwd=3, col=clu.col[col.temp[i]])
    }
  }
  

}

