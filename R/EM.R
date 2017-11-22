#' EM Algorithm for Fitting the Hidden Markov Model
#'
#' This function implements the EM algorithm to estimate HMM parameters,
#' and find the best sequence of hidden states based on model fitting.
#'
#' @param initial.value Output from initial.value function.
#' @param obs Output from read.process function.
#' @return res A matrix contains tansformed observations for each bin and the predicted direction for each bin (0-normal, 1-hyper, 2-hypo).
#' @return para Fitted parameters for HMM.
#' @export
EM <- function(initial.value, obs){
    cat("EM algorithm starts:", fill=TRUE)
    iter <- 1
    para.current <- initial.value
    
    repeat{
    	cat("Iteration ", iter, fill=TRUE)
	    e.loglikhood <- e.step(obs, para.current)	
	    para.new <- m.step(obs, para.current, e.loglikhood)	 	    
	    tol <- 0.001
	    if(abs(para.current$mu.pos-para.new$mu.pos) < tol &
	       abs(para.current$sd0 - para.new$sd0) < tol &
	       abs(para.current$rho - para.new$rho) < tol &
	       abs(sum(para.current$tran.p)-sum(para.new$tran.p)) < 10*tol
	       ){
	       	   # Find the sequence of best states.
	       	   state <- apply(e.loglikhood$P_k, 1, which.max)-1
	       	   para.current <- para.new
	           break
	     }else{
	           para.current <- para.new 
	           iter <- iter + 1
	     }
	}
	cat("EM algorithm converges.", fill=TRUE)
	res <- data.frame(cbind(obs, direction=as.numeric(state))) 
	return(list(res=res, para=para.current))	
}

m.step <- function(obs, para.current, e.loglikhood) {
    P_k     <- e.loglikhood$P_k  # P_k(b)
    P_kl    <- e.loglikhood$P_kl  # P_{kl}(b)
    rho     <- para.current$rho
    tran.p <- para.current$tran.p
    
    ## Update p_kl where (l !=k).   
    hks <- rep(NA, 3)
    for (k in 1:3) {
        l <- k.column(k)
        hk.init <- mean(apply(P_kl[, l], 2, sum) / as.vector(t(tran.p))[l])
        o <- optimize(f=update.p, interval=c(sum(P_kl[,l]), 100*sum(P_kl[,l])), 
                               maximum=TRUE, tol=0.01,
                               dist = obs[-1,2], P_kl=P_kl, rho=rho, k = k)
        hks[k]  <- o$maximum
    }
    tran.p.new <- matrix(apply(P_kl, 2, sum) / rep(hks, each = 3), 
                                    nrow=3, ncol=3, byrow=TRUE)
    
    ## Using grid search to find a new value for rho.    
    o <- lapply(seq(0.01, 10, 0.01), FUN=update.rho, dist=obs[-1,2], P_kl=P_kl, tran.p.new=tran.p.new)
    rho.new <- seq(0.01, 10, 0.01)[which.max(unlist(o))]
  
    ## Update mu+, mu-, sd0, tao1, and tao2.
    tmp.mu.pos.new <- sum(obs[,1]*P_k[,2])/sum(P_k[,2])
    tmp.mu.neg.new <- sum(obs[,1]*P_k[,3])/sum(P_k[,3])
    mu.pos.new <- (tmp.mu.pos.new + abs(tmp.mu.neg.new))/2
    mu.neg.new <- (-abs(tmp.mu.pos.new) + tmp.mu.neg.new)/2 
    sd0.new <- sqrt(sum(obs[,1]^2*P_k[,1])/sum(P_k[,1]))
    tao1.new <- sqrt(sum((obs[,1]-mu.pos.new)^2*P_k[,2])/sum(P_k[,2]))
    tao2.new <- sqrt(sum((obs[,1]-mu.neg.new)^2*P_k[,3])/sum(P_k[,3]))
    tmp.sigma.new <- (sd0.new + tao1.new + tao2.new)/3
    sd0.new =tao1.new=tao2.new=tmp.sigma.new
    
    ## Update p0, p1, and p2.
    p0.new <- P_k[1,1]
    p1.new <- P_k[1,2]
    p2.new <- P_k[1,3]
    
    return(list(p0=p0.new, p1=p1.new, p2=p2.new, 
                    mu.pos=mu.pos.new, mu.neg=mu.neg.new,
                    sd0=sd0.new, tao1=tao1.new, tao2=tao2.new,
                    rho=rho.new, tran.p=tran.p.new))       
}

update.rho <- function(rho.init, dist, P_kl, tran.p.new){
	tran.p.new[tran.p.new==0] <- 0.0001
	G2 <- sum(P_kl[,1]*log(1-sum(tran.p.new[1,c(2,3)])*(1-exp(-rho.init*dist)))+
	                 P_kl[,5]*log(1-sum(tran.p.new[2,c(1,3)])*(1-exp(-rho.init*dist)))+
	                 P_kl[,9]*log(1-sum(tran.p.new[3,c(1,2)])*(1-exp(-rho.init*dist))))
	G3 <- sum(P_kl[,2]*log(tran.p.new[1,2]*(1-exp(-rho.init*dist)))+
	                 P_kl[,3]*log(tran.p.new[1,3]*(1-exp(-rho.init*dist)))+
	                 P_kl[,4]*log(tran.p.new[2,1]*(1-exp(-rho.init*dist)))+
	                 P_kl[,6]*log(tran.p.new[2,3]*(1-exp(-rho.init*dist)))+
	                 P_kl[,7]*log(tran.p.new[3,1]*(1-exp(-rho.init*dist)))+
	                 P_kl[,8]*log(tran.p.new[3,2]*(1-exp(-rho.init*dist))))
	return(G2+G3)
}

update.p <- function(hk.init, dist, P_kl, rho, k) {
  
    if (k == 1) {
        tmp1 <- sum(P_kl[, c(2,3)])
        term1 <- sum(P_kl[, 1] * log(1-tmp1*(1-exp(-rho*dist))/hk.init))
        term2 <- sum(P_kl[, 2] * log(sum(P_kl[, 2])*(1-exp(-rho*dist)) / hk.init) +
                              P_kl[, 3] * log(sum(P_kl[, 3])*(1-exp(-rho*dist)) / hk.init)) 
    }
    if (k == 2) {
        tmp1 <- sum(P_kl[, c(4,6)])
        term1 <- sum(P_kl[, 5] * log(1-tmp1*(1-exp(-rho*dist))/hk.init))
        term2 <- sum(P_kl[, 4] * log(sum(P_kl[, 4])*(1-exp(-rho*dist)) / hk.init) +
                              P_kl[, 6] * log(sum(P_kl[, 6])*(1-exp(-rho*dist)) / hk.init)) 
    }
    if (k == 3) {
        tmp1 <- sum(P_kl[, c(7,8)])
        term1 <- sum(P_kl[, 9] * log(1-tmp1*(1-exp(-rho*dist))/hk.init))
        term2 <- sum(P_kl[, 7] * log(sum(P_kl[, 7])*(1-exp(-rho*dist)) / hk.init) +
                              P_kl[, 8] * log(sum(P_kl[, 8])*(1-exp(-rho*dist)) / hk.init)) 
    }

    return(term1 + term2)
}


k.column <- function(k) {
    if (k == 1) 
        out <- c(2, 3)
    if (k == 2)
        out <- c(4, 6)
    if (k == 3)
        out <- c(7, 8)
    return(out)
}    
    

e.step <- function(obs, para){
    B <- dim(obs)[1]  # B: total number of bins.
    
    ## forward algorithm.
    f <- NULL    
    mu.pos <- para$mu.pos  
    mu.neg <- para$mu.neg    
    sd0 <- para$sd0   
    tao1 <- para$tao1   
    tao2 <- para$tao2   
    bin.cnt <- 1
    f <- rbind(f, c(para$p0*dnorm(obs[bin.cnt,1], 0, sd=sd0),
                         para$p1*dnorm(obs[bin.cnt,1], mu.pos, sd=tao1), 
                         para$p2*dnorm(obs[bin.cnt,1], mu.neg, sd=tao2)))
    repeat{
          bin.cnt <- bin.cnt + 1
    	  f.bin.b <- forward(f[bin.cnt-1,], para, obs[bin.cnt,c(1,2)])
    	  f <- rbind(f, f.bin.b)
    	  if(bin.cnt>=B){break}
    }
    
    ## backward algorithm.
    bin.cnt <- B
    b <- NULL  
	b <- rbind(b, c(1,1,1))
	repeat{  
		    bin.cnt <- bin.cnt - 1
		    b.bin.b <- backward(b[1,], para, obs[bin.cnt+1,c(1,2)])
		    b <- rbind(b.bin.b, b)
		    if(bin.cnt <= 1){break}
	  }
	
  	P_k <- f*b/apply(f*b, 1, sum)	
  	tran.p <- para$tran.p
  	rho <- para$rho
  	sum_P_kl_loga_kl <- 0
  	P_kl_all <- NULL
   	sum_P_k_logOb <- 0
  	for(bin.cnt in 1:(B-1)){
  		dist.rate <- 1-exp(-rho*obs[bin.cnt+1,2])
  	    tran.a <- tran.p*dist.rate
  	    tran.a[1,1] <- 1-sum(tran.a[1,c(2,3)])
  	    tran.a[2,2] <- 1-sum(tran.a[2,c(1,3)])
  	    tran.a[3,3] <- 1-sum(tran.a[3,c(1,2)])
        
        P_kl_temp <- NULL
        for(i in 1:3){
               P_kl_temp <- c(P_kl_temp, 
                                      f[bin.cnt,i]*tran.a[i,]*dnorm(obs[bin.cnt+1,1],
                                      mean=c(0,mu.pos,mu.neg),
  	                                  sd=c(sd0,tao1,tao2))*b[bin.cnt+1,])
  	    }  
  	    tran.a[which(tran.a==0)] <- 0.0001
  	    P_kl_all <- rbind(P_kl_all, P_kl_temp/sum(P_kl_temp))
  	    sum_P_kl_loga_kl <- sum_P_kl_loga_kl + 
  	                                    sum(P_kl_temp/sum(P_kl_temp)*log(as.vector(t(tran.a))))
  	  
  	    Ob <- dnorm(obs[bin.cnt,1],c(0, mu.pos, mu.neg),
  	                           sd=c(sd0, tao1, tao2))
  	    Ob[Ob==0] <- 0.0001
  	    sum_P_k_logOb <- sum_P_k_logOb + sum(P_k[bin.cnt, ]*log(Ob))
  	    if(bin.cnt == (B-1)){
  	        Ob <- dnorm(obs[B,1],c(0, mu.pos, mu.neg),
  	                     sd=c(sd0, tao1, tao2))
  	        Ob[Ob==0] <- 0.0001
  	        sum_P_k_logOb <- sum_P_k_logOb + sum(P_k[B, ]*log(Ob))
  	    }  
  	}
    if(para$p0==0){para$p0=0.0001}
    if(para$p1==0){para$p1=0.0001}
    if(para$p2==0){para$p2=0.0001}
    e.loglik <- sum(P_k[1,]*log(c(para$p0, para$p1, para$p2))) + 
                       sum_P_kl_loga_kl +
                       sum_P_k_logOb                    
    return(list(e.loglik=e.loglik, P_k = P_k, P_kl = P_kl_all))
	
}


forward <- function(f, para, obs){
    tran.p <- para$tran.p
    rho <- para$rho
    o.b <- obs[1]
    dist.b <- obs[2]

    dist.rate <- 1-exp(-rho*dist.b)
	tran.a <- tran.p*dist.rate
	tran.a[1,1] <- 1-sum(tran.a[1,c(2,3)])
	tran.a[2,2] <- 1-sum(tran.a[2,c(1,3)])
	tran.a[3,3] <- 1-sum(tran.a[3,c(1,2)])

    f.bin.b <- rbind(f[1]*tran.a[1,], f[2]*tran.a[2,], f[3]*tran.a[3,])
                              
    mu.pos <- para$mu.pos  
    mu.neg <- para$mu.neg   
    sd0 <- para$sd0  
    tao1 <- para$tao1  
    tao2 <- para$tao2  
     f.bin.b <- cbind(f.bin.b[,1]*dnorm(o.b, 0, sd=sd0), 
                               f.bin.b[,2]*dnorm(o.b, mu.pos, sd=tao1),
                               f.bin.b[,3]*dnorm(o.b, mu.neg, sd=tao2)) 
    f.bin.b <- apply(f.bin.b, 2, sum)
    f.bin.b <- f.bin.b*(1/f.bin.b)[which.min(1/f.bin.b)]
    return(f.bin.b)
}


backward <- function(b, para, obs){
	tran.p <- para$tran.p
	rho <- para$rho
	o.b <- obs[1]
	dist.b <- obs[2]
	
	mu.pos <- para$mu.pos  
    mu.neg <- para$mu.neg   
    sd0 <- para$sd0  
    tao1 <- para$tao1  
    tao2 <- para$tao2  

	dist.rate <- 1-exp(-rho*dist.b)
	tran.a <- tran.p*dist.rate
	tran.a[1,1] <- 1-sum(tran.a[1,c(2,3)])
	tran.a[2,2] <- 1-sum(tran.a[2,c(1,3)])
	tran.a[3,3] <- 1-sum(tran.a[3,c(1,2)])
	
	b.bin.b <- rbind(tran.a[,1]*dnorm(o.b, 0, sd=sd0),
	                                tran.a[,2]*dnorm(o.b, mu.pos, sd=tao1),
                                    tran.a[,3]*dnorm(o.b, mu.neg, sd=tao2))
	b.bin.b <- c(sum(b.bin.b[,1]*b), sum(b.bin.b[,2]*b), sum(b.bin.b[,3]*b))
	b.bin.b <- b.bin.b*(1/b.bin.b)[which.min(1/b.bin.b)]
	return(b.bin.b)
}








