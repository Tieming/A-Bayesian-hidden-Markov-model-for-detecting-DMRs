#
# Raw Read Data Processing
#
# Tieming, 7/24/2017.
#
#

#' Read and Process Methylation Data into Bins
#'
#' This function reads in CpG methylation data and transforms observations into bins.
#' Users define starting position, end position, and bin size.
#'
#' @param start.pos Starting position for data analysis. Default to 1.
#' @param end.pos End position for data analysis. Default to the maximum position in the parameter 'pos'.
#' @param pos CpG positions
#' @param norm.m A matrix contains number of methylated reads of the normal group, each column is a biological sample
#' @param norm.um A matrix contains number of un-methylated reads of the normal group
#' @param los.m A matrix contains number of methylated reads of the abnormal group
#' @param los.um A matrix contains number of un-methylated reads of the abnormal group
#' @param bin.size User-defined bin size. Default to 40
#' @export
read.process <- function(start.pos=1, end.pos=max(pos), pos, norm.m, norm.um, los.m, los.um, bin.size){
	
	# apply Haldane-Anscombe correction.
	norm.m <- norm.m + 0.5
	norm.um <- norm.um + 0.5
	los.m <- los.m + 0.5
	los.um <- los.um + 0.5
	
	norm.c <- norm.m + norm.um
	los.c    <- los.m + los.um
	dist <- NULL
	dist.flag <- TRUE
	# when dist.flag==TRUE, "dist" starts to count from 1.
	# when dist.flag==FALSE, "dist" adds 1 from the last "dist".
    norm.p <- NULL
    los.p    <- NULL
    start.loc  <- NULL
    end.loc   <- NULL 
	for(i in 1:ceiling((end.pos-start.pos)/bin.size)){
		
		start <- start.pos + (i-1)*bin.size
		end  <- start.pos + i*bin.size -1
		temp.pos <- which(pos >= start & pos <= end)
	
	    if(length(temp.pos)>1){
	    	start.loc <- c(start.loc, start)
	    	end.loc  <- c(end.loc, end)
	        norm.p <- c(norm.p, sum(apply(as.matrix(norm.m[temp.pos,]), 2, sum))/sum(apply(as.matrix(norm.c[temp.pos,]), 2, sum)))
	        los.p    <- c(los.p, sum(apply(as.matrix(los.m[temp.pos,]), 2, sum))/sum(apply(as.matrix(los.c[temp.pos,]), 2, sum)))
	        if(dist.flag==FALSE){
	        	dist <- c(dist, temp.dist+1)
	        	dist.flag <- TRUE
	        }else if(dist.flag==TRUE){
	        	dist <- c(dist, 1)
	        	dist.flag <- TRUE
	        }
	    }else if(length(temp.pos)==1){
	    	start.loc <- c(start.loc, start)
	    	end.loc  <- c(end.loc, end)
	        norm.p <- c(norm.p, sum(as.matrix(norm.m[temp.pos,]))/sum(as.matrix(norm.c[temp.pos,])))
	        los.p     <- c(los.p, sum(as.matrix(los.m[temp.pos,]))/sum(as.matrix(los.c[temp.pos,])))
	        if(dist.flag==FALSE){
	        	dist <- c(dist, temp.dist+1)
	        	dist.flag <- TRUE
	        }else if(dist.flag==TRUE){
	        	dist <- c(dist, 1)
	        	dist.flag <- TRUE
	        }
	    }else if(length(temp.pos)==0 & dist.flag==FALSE){
	        temp.dist <- temp.dist + 1
	        dist.flag <- FALSE	
	    }else if(length(temp.pos)==0 & dist.flag==TRUE){
	    	temp.dist <- 1
	    	dist.flag <- FALSE
	    }
	}
	dist <- dist * bin.size
    norm.y <- log(norm.p/(1-norm.p))
    los.y     <- log(los.p/(1-los.p))
    o <- los.y - norm.y
    return(obs =cbind(o=o, dist=dist, los.y=los.y, norm.y=norm.y, los.p=los.p, norm.p=norm.p, start=start.loc, end=end.loc))
}



