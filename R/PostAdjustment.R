#' Post Adjustment
#'
#' This function adjusts HMM output such that each detected DMR has a minimum length and a minimum number of CpGs in each DMR.
#'
#' @param em.o Output from EM function.
#' @param CpG.pos Contains all CpG positions.
#' @param min.length Minimum length of a DMR. Default to 200.
#' @param min.CpGs Minimum number of CpGs contained in a DMR. Default to 3.
#' @return dmr.res A matrix contains transformed observations for each bin and the adjusted predictions for each bin (0-normal, 1-hyper, 2-hypo).
#' @return region A matrix contains detected regions satisfying user-defined length and number of CpGs.
#' @export
PostAdjustment <- function(em.o, CpG.pos, min.length=200, min.CpGs=3){
		direction <- em.o$res$direction
		flag <- FALSE
		region <- data.frame(NULL)
		region.cnt <- 0
		for(bin.cnt in 1:length(direction)){
				if(direction[bin.cnt]!=0 & flag==FALSE){
						flag <- TRUE
					    tmp.direction <- direction[bin.cnt]
					    start.pos <- em.o$res[bin.cnt, "start"]
					    start.bin <- bin.cnt
				}else if(direction[bin.cnt]!=0 & flag==TRUE){
						if(tmp.direction != direction[bin.cnt]){
								end.pos <- em.o$res[bin.cnt-1, "end"]
								end.bin <- bin.cnt-1
								DMR.length <- end.pos - start.pos
								num.CpGs <- sum(CpG.pos>=start.pos & CpG.pos<=end.pos)
								if(DMR.length < min.length | num.CpGs < min.CpGs){
										direction[start.bin:end.bin] <- 0
								}else{
										region.cnt <- region.cnt + 1
										if(tmp.direction==1){region.state="hyper"}
										if(tmp.direction==2){region.state="hypo"}
										region <- rbind(region, data.frame(region.cnt, 
										                 region.start=start.pos, region.end=end.pos, 
										                 region.state, num.CpGs))
								}
								start.pos <- em.o$res[bin.cnt, "start"]		
								start.bin <- bin.cnt
								tmp.direction <- direction[bin.cnt]
						}
						if(tmp.direction == direction[bin.cnt] & bin.cnt==length(direction)){
								end.pos <- em.o$res[bin.cnt, "end"]
								end.bin <- bin.cnt
								DMR.length <- end.pos - start.pos
								num.CpGs <- sum(CpG.pos >= start.pos & CpG.pos <= end.pos)
								if(DMR.length < min.length | num.CpGs < min.CpGs){
										direction[start.bin:end.bin] <- 0
								}else{
										region.cnt <- region.cnt + 1
										if(tmp.direction==1){region.state="hyper"}
										if(tmp.direction==2){region.state="hypo"}
										region <- rbind(region, data.frame(region.cnt, 
										                 region.start=start.pos, region.end=end.pos, 
										                 region.state, num.CpGs))
								}
								
						}
				}else if(direction[bin.cnt]==0 & flag==TRUE){
						flag <- FALSE
						end.pos <- em.o$res[bin.cnt-1, "end"]
						end.bin <- bin.cnt-1
						DMR.length <- end.pos - start.pos
						num.CpGs <- sum(CpG.pos>=start.pos & CpG.pos<=end.pos)
						if(DMR.length < min.length | num.CpGs < min.CpGs){
								direction[start.bin:end.bin] <- 0
						}else{
								region.cnt <- region.cnt + 1
								if(tmp.direction==1){region.state="hyper"}
								if(tmp.direction==2){region.state="hypo"}
								region <- rbind(region, data.frame(region.cnt, 
										         region.start=start.pos, region.end=end.pos, 
										         region.state, num.CpGs))

						}
				}		
		}
		region$length <- region$region.end - region$region.start
		em.o$res$direction <- direction
		dmr.res <- em.o$res
		return(list(dmr.res=dmr.res, region=region))
}




