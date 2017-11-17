#
# Set up initial values for EM algorithm.
#
# Tieming Ji, 7/25/2017.
#
# InitialValue2.R: (1) requires equal absoluate values of mu+ and mu-; 
#	                      (2) requires equal values of sigma0^2, tao1^2, and tao2^2.
# 
#


#' Set up Initial Values for HMM
#'
#' This function sets up initial values for the EM algorithm to
#' estimate parameters of the hidden Markov model.
#'
#' @param obs Transformed observations for bins. This is output from read.process function
#' @export
initial.value <- function(obs){
### set up initial values for EM algorithm.
    sd <- sd(obs[,1])/2
    s1 <- which(obs[,1] > 2*sd)
    s2 <- which(obs[,1] < -2*sd)
    initial.state <- rep(0, length(obs[,1]))
    initial.state[s1] <- 1
    initial.state[s2] <- 2
    
    # emission parameters.
    if(sum(initial.state==1)!=0 & sum(initial.state==1)!=1){
        tmp.mu.pos <- mean(obs[initial.state==1,1])
    }else{
    	tmp.mu.pos=3*sd
    }
    if(sum(initial.state==2)!=0 & sum(initial.state==2)!=1){
    	tmp.mu.neg <- mean(obs[initial.state==2,1])
    }else{
    	tmp.mu.neg= -3*sd
    }
    mu.pos <- (tmp.mu.pos+abs(tmp.mu.neg))/2
    mu.neg <- (-abs(tmp.mu.pos)+tmp.mu.neg)/2
    sd0 <- sd
    tao1 <- sd
    tao2 <- sd
    
    p0=p1=p2=1/3

    rho <- 0.5
    tran.p <- array(0.1, dim=c(3,3))

	
	return(initial.value=list(p0=p0, p1=p1, p2=p2,
	           mu.pos=mu.pos, mu.neg=mu.neg, 
	           sd0=sd0, tao1=tao1, tao2=tao2, 
	           rho=rho, tran.p=tran.p))
	
}







