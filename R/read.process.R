#' Read and Transform Observed Sequencing Data
#'
#' This function reads in methylated and unmethylated read count data,
#' and transform it into logarithm bin-wise data.
#'
#' @param pos CpG positions.
#' @param norm.m A matrix contains methylated read count data of normal group. Each column of the matrix represents a replicate and each row represents a CpG position.
#' @param norm.um A matrix contains unmethylated read count data of normal group.
#' @param abnorm.m A matrix contains methylated read count data of abnormal group.
#' @param abnorm.um A matrix contains unmethylated read count data of abnormal group.
#' @param bin.size Size of a bin. Default to 40.
#' @return obs A matrix contains transformed observations for each bin, distance between the current bin the bin ahead of it, average methylation rate of abnormal and normal groups, start and end position of each bin.
#' @export

read.process <- function(pos, norm.m, norm.um, abnorm.m, abnorm.um, bin.size=40){
	
	# apply Haldane-Anscombe correction.
	norm.m <- norm.m + 0.5
	norm.um <- norm.um + 0.5
	abnorm.m <- abnorm.m + 0.5
	abnorm.um <- abnorm.um + 0.5
	
	norm.c <- norm.m + norm.um
	abnorm.c    <- abnorm.m + abnorm.um
	dist <- NULL
	dist.flag <- TRUE
    norm.p <- NULL
    abnorm.p    <- NULL
    start.pos <- min(pos)
    end.pos  <- max(pos)
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
	        abnorm.p    <- c(abnorm.p, sum(apply(as.matrix(abnorm.m[temp.pos,]), 2, sum))/sum(apply(as.matrix(abnorm.c[temp.pos,]), 2, sum)))
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
	        abnorm.p <- c(abnorm.p, sum(as.matrix(abnorm.m[temp.pos,]))/sum(as.matrix(abnorm.c[temp.pos,])))
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
    abnorm.y <- log(abnorm.p/(1-abnorm.p))
    o <- abnorm.y - norm.y
    return(obs =cbind(o=o, dist=dist, abnorm.p=abnorm.p, norm.p=norm.p, start=start.loc, end=end.loc))
}



