#
# Post Adjustment
# Goal: (1) Require DMRs to have a minimum length, such as 500 bps.
#			(2) Require DMRs to contain a minimum number of CpGs, such as 3 CpGs.
#
# Tieming, 10/9/2017
#
#


#' Post Adjustment
#'
#' This function selects output from HMM to ensure the minimum length of DMRs and 
#' minimum number of CpGs in each DMR.
#'
#' @param obs Transformed observations for all bins. This is the output from read.process function
#' @param em.o Parameter estimation and output from HMM
#' @param CpG.pos CpG positions
#' @param min.length User-defined minimum DMR length. Default to 200
#' @param min.CpGs User-defined minimum number of CpGs in each DMR. Default to 3
#' @export
PostAdjustment <- function(obs, em.o, CpG.pos, min.length=200, min.CpGs=3){

		direction <- em.o$state
		flag <- FALSE
		region <- data.frame(NULL)
		region.cnt <- 0
		for(bin.cnt in 1:length(direction)){
				if(direction[bin.cnt]!=0 & flag==FALSE){
						flag <- TRUE
					    tmp.direction <- direction[bin.cnt]
					    start.pos <- obs[bin.cnt, "start"]
					    start.bin <- bin.cnt
				}else if(direction[bin.cnt]!=0 & flag==TRUE){
						if(tmp.direction != direction[bin.cnt]){
								end.pos <- obs[bin.cnt-1, "end"]
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
								start.pos <- obs[bin.cnt, "start"]		
								start.bin <- bin.cnt
								tmp.direction <- direction[bin.cnt]
						}
						if(tmp.direction == direction[bin.cnt] & bin.cnt==length(direction)){
								end.pos <- obs[bin.cnt, "end"]
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
						end.pos <- obs[bin.cnt-1, "end"]
						end.bin <- bin.cnt-1
						DMR.length <- end.pos - start.pos
						#cat("DMR.start is ", start.pos, fill=TRUE)
						#cat("DMR.end is ", end.pos, fill=TRUE)
						#cat("DMR.length is ", DMR.length, fill=TRUE)
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
		dmr.res <- data.frame(obs, direction)
		return(list(dmr.res=dmr.res, region=region))
}




