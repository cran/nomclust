#' Visualization of Evaluation Criteria
#' 
#' @description The function visualizes the values of up to eight evaluation criteria for the range of cluster solutions defined by the user in the \bold{nomclust}, \bold{evalclust} or \bold{nomprox} functions.
#'  It also indicates the optimal number of clusters determined by these criteria. The charts for the evaluation criteria in the \bold{nomclust} package.
#' 
#' @param x An output of the "nomclust" object containing the \code{eval} and \code{opt} components.
#' 
#' @param criteria A character string or character vector specifying the criteria that are going to be visualized. It can be selected one particular criterion, a vector of criteria, or all the available criteria by typing \code{"all"}.
#' 
#' @param style A character string or a vector of colors defines the graphical style of the produced plots. There are two predefined styles in the \bold{nomclust} package, namely \code{"greys"} and \code{"dark"}, but a custom color scheme can be set by a user as a vector of a length four.
#' 
#' @param opt.col An argument specifying a color that is used for the optimal number of clusters identification.
#' 
#' @param main A character string with the chart title.
#' 
#' @param ... Other graphical arguments compatible with the generic \code{plot()} function.
#' 
#' @return The function returns a series of up to eight plots with evaluation criteria values and the graphical indication of the optimal numbers of clusters (for AIC, BIC, BK, PSFE, PSFM, SI).
#'
#' @details The function can display up to eight evaluation criteria. Namely, Within-cluster mutability coefficient (WCM), Within-cluster entropy coefficient (WCE),
#' Pseudo F Indices based on the mutability (PSFM) and the entropy (PSFE), Bayesian (BIC), and Akaike (AIC) information criteria for categorical data, the BK index, and the silhouette index (SI).
#' \cr
#'
#' @seealso
#' \code{\link[nomclust]{dend.plot}}, \code{\link[nomclust]{nomclust}}, \code{\link[nomclust]{evalclust}}, \code{\link[nomclust]{nomprox}}.
#' 
#' @author Jana Cibulkova and Zdenek Sulc. \cr Contact: \email{jana.cibulkova@@vse.cz}
#' 
#' @examples
#' # sample data
#' data(data20)
#' 
#' # creating an object with results of hierarchical clustering 
#' hca.object <- nomclust(data20, measure = "iof", eval = TRUE)
#' 
#' # a default series of plots
#' eval.plot(hca.object)
#' 
#' # changing the color indicating the optimum number of clusters
#' eval.plot(hca.object, opt.col= "darkorange")
#' 
#' # selecting only AIC and BIC criteria with the dark style
#' eval.plot(hca.object, criteria = c("AIC", "BIC"), style = "dark")
#' 
#' # selecting only SI
#' eval.plot(hca.object, criteria = "SI")
#' 
#' @export 


eval.plot = function(x, criteria = "all", style = "greys", opt.col = "red", main= "Cluster Evaluation", ...){
  if(style[1] == "dark"){style = rep("black",4)}
  if(style[1] == "greys"){    style = grey.colors(5)}
  
  if(typeof(x) == "list"  &  "eval" %in% names(x)  &  "opt" %in% names(x)){
    if(criteria[1] == "all"){
      par(ask=TRUE)
      x$eval=x$eval[-1]
    }
    if(criteria[1] != "all"){
      criteria=toupper(criteria)
      if(sum(!(criteria %in% names(x$eval)))>0){stop('Choose the number of clusters to be visualized. It should be either a number or one of these criteria: AIC, BIC, BK, PSFM, PSFE.')}
      x$eval = x$eval[which(names(x$eval) %in% criteria)]
      x$opt  = x$opt[ which(names(x$opt)  %in% criteria)]
    }
    if(length(criteria) > 1){par(ask=TRUE)}
    
      
    par(mgp=c(5,1,0))
    par(mar=c(7.1, 7.1, 4.1, 2.1))
    i=1
    for(i in 1:length(x$eval)){
      plot(unlist(x$eval[i]),type='n',yaxt='n',
           xlab="",ylab=paste(names(x$eval)[i],'values'),
           bty="n",col.lab = style[3],col=style[1],main=main)     # initialize plot
      mtext("Number of clusters",side=1,line=2,col=style[3])
      if(names(x$eval)[i] %in% c("WCM","WCE")){
        mtext(paste("Within-Cluster Variability based on:",names(x$eval)[i]),col= style[1])
      }else{
        mtext(paste("Optimal number of clusters based on:",names(x$eval)[i]),col= style[1])
      }
      h=seq(  floor(min(unlist(x$eval[i]),na.rm=T)*100)/100,
              ceiling(max(unlist(x$eval[i]),na.rm=T)*100)/100,length.out = 10)# seq for y axis
      axis(labels=F, side=2,las=2,col = style[3], at=h)                       # create y axis
      mtext(format(round((h*100)/100,2), nsmall = 2),side=2,las=2,col = style[1], at=h, line = 1) # add labels
      abline(h=h,lty=3,lwd=1,col=style[4])                                    # add horizontal lines
      lines(unlist(x$eval[i]),lty=1,col=style[1])                             # lines for eval. criterion
      points(unlist(x$eval[i]),pch=16)                                        # points for eval. criterion
      
      temp=which(names(x$opt)==names(x$eval)[i])
      abline(v=x$opt[temp],lty=2, lwd=1, col=opt.col)                                # optimal no. of clusters
      
    }
    
    
  }else{
    stop("An input argument x is missing or incorrect. Output from nomclust(), evalclust() or nomprox() with 'eval' and 'opt' components is required.")
    
  }
  par(ask=FALSE)
} 
