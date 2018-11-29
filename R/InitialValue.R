#' Set up Initial Value for EM algorithm
#'
#' This function sets up initial values for estimating hidden Markov model 
#' parameters by EM algorithm
#'
#' @param obs Output from read.process function
#' @return initial.value Initial values for EM algorithm
#' @export
initial.value <- function(obs){
    sd <- sd(obs[,1])
    s1 <- which(obs[,1] > sd)
    s2 <- which(obs[,1] < -sd)
    initial.state <- rep(0, length(obs[,1]))
    initial.state[s1] <- 1
    initial.state[s2] <- 2
    
    # emission parameters.
    if(sum(initial.state==1)!=0 & sum(initial.state==1)!=1){
        tmp.mu.pos <- mean(obs[initial.state==1,1])
    }else{
    	tmp.mu.pos=3/2*sd
    }
    if(sum(initial.state==2)!=0 & sum(initial.state==2)!=1){
    	tmp.mu.neg <- mean(obs[initial.state==2,1])
    }else{
    	tmp.mu.neg=-3/2*sd
    }
    mu.pos <- (tmp.mu.pos+abs(tmp.mu.neg))/2
    mu.neg <- (-abs(tmp.mu.pos)+tmp.mu.neg)/2
    sd0 <- sd
    tao1 <- sd
    tao2 <- sd
    
    # transition parameters.
    p0=p1=p2=1/3
    rho <- 0.05
    tran.p <- array(0.1, dim=c(3,3))
	
	return(initial.value=list(p0=p0, p1=p1, p2=p2,
	           mu.pos=mu.pos, mu.neg=mu.neg, 
	           sd0=sd0, tao1=tao1, tao2=tao2, 
	           rho=rho, tran.p=tran.p))
	
}







