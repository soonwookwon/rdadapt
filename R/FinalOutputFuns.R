which_proc = function(Cgrid,C,gamma,X,tind,kplus_ind,procmat,ProcNames,Nsim_proc,
                      sigma_t, sigma_c, alpha,lower){
  
	# Cgrid: vector of grids for C
	# C: list of C's to adapt
	# gamma: list of gamma's to adapt (=1 for now)
	# X: (n * d) data matrix
	# tind: n-dim treatment indicator vector
	# kplus_ind: monotone variable indicator 
	# procmat: adaptation indicator matrix
	# ProcNames: names for adaptive procedures
	# Nsim: number of simulations to use when calculating minimax procedure
	# sig_fun: (estimated) variance function
  # alpha: significance level
  # lower: True if we are considering lower CI
	

	gamgrid <- rep(1,length(Cgrid))
	num_grid <- length(Cgrid)
	gamma <- rep(1,length(C))   # gamma = 1
	N_proc <- length(ProcNames)

	Xt <- X[tind,,drop=FALSE]   
	Xc <- X[!tind,,drop=FALSE]
	nt <- sum(tind)   # Number of treated sample
	nc <- n - nt   # Number of control sample
	k <- length(X[1,])
	
	truef <- function (gamma, C) {
	  
	  f <- function(x) {
	    
	    x <- matrix(x,1,k)
	    
	    fc <- C*Norm(Vplus(x, kplus_ind))^gamma - C*Norm(Vminus(x, kplus_ind))^gamma   
	    res <- fc
	    
	    return(res)
	    
	    
	  }
	  return(f)
	  
	}
	
  
	# fX_wc is a (n * num_grid)-dim matrix where f(X)[,j] = (f_j(x_1), ..., f_j(x_n))
	# with f_j denoting the worst-case function under (gamma_j, C_j)

	fX_wc <- matrix(NA, nrow = n, ncol = num_grid)
	for(j in 1:num_grid) {
	  fX_wc[,j] <- apply(X, 1, function(x) truef(gamgrid[j],Cgrid[j])(matrix(x,nrow=1)))
	}

	fXt_wc <- fX_wc[tind == 1,]
	fXc_wc <- fX_wc[tind == 0,]
	
	# Calculating moduli of continuity and related quantities

	proc_bmat = matrix(0,sum(procmat),4)
	proc_dmat = matrix(0,sum(procmat),4)
	proc_alphavec = numeric(N_proc)

	beg_ind <- 1

	# This for-loop part sometimes runs without error but it doesn't other times. What's the problem?
	for(i in 1:N_proc){

		print(paste("proc",i))

		adpt_ind <- procmat[i,]
		adpt_len <- sum(procmat[i,])
		end_ind <- beg_ind + adpt_len -1
		modres <- mod_del_cal(gamma[adpt_ind],C[adpt_ind],Xt,Xc,kplus_ind,sigma_t,sigma_c)
		proc_bmat[beg_ind:end_ind,] <- modres$b_mat
		proc_dmat[beg_ind:(beg_ind + adpt_len -1),] <- modres$d_mat
		
		if(lower == T){
		  proc_alphavec[i] <- modres$alnewL
		}else{
		  proc_alphavec[i] <- modres$alnewU
		}
		

		beg_ind <- beg_ind + adpt_len

	}

	# Oracle
	bmat_orc = matrix(0,num_grid,4)
	dmat_orc = matrix(0,num_grid,4)

	for(j in 1:num_grid){
		res_mod_orc = mod_del_cal_orc(gamgrid[j],Cgrid[j],max(gamgrid),max(Cgrid),
		                              Xt,Xc,kplus_ind,sigma_t,sigma_c)
		bmat_orc[j,] = res_mod_orc$b_mat
		dmat_orc[j,] = res_mod_orc$d_mat
	}

	# Generating y_i's and computating CIs

	CIs <- array(0,c(num_grid,Nsim_proc,N_proc))
	CIs_orc <- matrix(0,Nsim_proc,num_grid)


	for (i in 1:Nsim_proc) {

		ut <- rnorm(nt, sd = sigma_t)   # generate u_i's for the treated
		Yt <- fXt_wc + ut
		uc <- rnorm(nc, sd = sigma_c)   # generate u_i's for the contol
		Yc <- fXc_wc + uc

		for (j in 1:num_grid){

		beg_ind <- 1

  		for (j2 in 1:N_proc){
  		  
  		  adpt_ind <- procmat[j2,]
  		  adpt_len <- sum(procmat[j2,])
  		  end_ind <- beg_ind + adpt_len -1
  		  
  		  CIs[j,i,j2] <- CI_one_sd_RD(proc_bmat[beg_ind:end_ind,,drop=F],
  										proc_dmat[beg_ind:end_ind,,drop=F],
  										gamma[adpt_ind],C[adpt_ind],Xt,Xc,
  										kplus_ind,Yt[,j],Yc[,j],lower,
  										al = proc_alphavec[j2],Dir=F)
  		  
  		  beg_ind <- beg_ind + adpt_len
  	  
  		}
		
		CIs_orc[i,j] <- CI_one_sd_RD(bmat_orc[j,,drop=F],dmat_orc[j,,drop=F],
									   gamgrid[j],Cgrid[j],
									   Xt,Xc,kplus_ind,Yt[,j],Yc[,j],lower,
									   al = alpha,Dir=F)
		}


		if(i%%50 == 0) print(i)

	}
	
	CIlen <- apply(0 - CIs,MARGIN = c(1,3),FUN=mean,na.rm=T)
	CIlen_orc <- apply(0 - CIs_orc,MARGIN = 2,FUN=mean,na.rm=T)
	
	if(lower == F){
	  CIlen <- -CIlen
	  CIlen_orc <- -CIlen_orc
	}
	
	regret_mat <- rhofun(CIlen,
	                     matrix(rep(CIlen_orc,N_proc),num_grid,N_proc,byrow = F))
	
	maxregret <- apply(regret_mat,2,max)
	names(maxregret) <- ProcNames
	
	minMRind <- which.min(maxregret)
	bestProc <- C[procmat[minMRind,]]
	bestProcName <- ProcNames[minMRind]
	
	print(paste("C",seq(1:length(bestProc)),":",bestProc),sep="")
	
	res <- list(CIlen = CIlen,CIlen_orc = CIlen_orc,regret_mat = regret_mat,
	            maxregret = maxregret,bestProc = bestProc,bestProcName = bestProcName,
	            proc_bmat = proc_bmat,proc_dmat = proc_dmat,
	            proc_alphavec = proc_alphavec)
	return(res)
  
}



which_proc_fig = function (pdfname,Cgrid,proc_res_list,ProcNames,
                           leg_loc = "bottomright",disp_mmreg = F){
  
  N_proc <- length(ProcNames)
  
  bestProcName <- proc_res_list$bestProcName
  
  CIlen <- proc_res_list$CIlen
  CIlen_orc <- proc_res_list$CIlen_orc
  
  pdf(pdfname)
  
  if(disp_mmreg == T){
    
    regret_mat <- proc_res_list$regret_mat
    par(mfrow=c(1,2))
    
  } 
  
  miny = min(CIlen)*0.9
  maxy = max(CIlen)*1.1
  plot(Cgrid,CIlen[,1],type="l",ylim=c(miny,maxy),
       ylab="Average excess length of lower CI",xlab="Bound on f'(x)")
  
  if(N_proc > 1){
    
    for(i in 2:N_proc){
      lines(Cgrid,CIlen[,i],lty=i)
    }
    
  }
  
  lines(Cgrid,CIlen_orc,col="red")
  
  legend(leg_loc,c(ProcNames,"oracle"),
         lty=c(1:N_proc,1),col=c(rep("black",N_proc),"red"),bg="white",cex=0.75)
  
  if(disp_mmreg == T){
    
    miny = min(regret_mat)*0.9
    maxy = max(regret_mat)*1.1
    plot(Cgrid,regret_mat[,1],type="l",ylim=c(miny,maxy),xlab="Bound on f'(x)",
         ylab="Distance to the oracle")
    
    if(N_proc > 1){
      
      for(i in 2:N_proc){
        lines(Cgrid,regret_mat[,i],type="l",lty=i)
      }
      
    }
    
    
    legend("topright",ProcNames,lty=1:N_proc,bg="white",cex=0.75)
    
  }
  
  titlename = paste("Comparison of CI Lengths (n=",
                    n,", C in [",min_C,", ",max_C,"])",sep="")
  mtext(titlename,line = -2.5, outer = T)
  mtext(paste("DGP: ",DGPname,"//","Chosen Proc:",bestProcName),line = -3.5, outer = T, cex = 0.8)
  dev.off()
  
  return(0)
  
}

CI_gen <- function(Y,X,tind,C,kplus_ind,lower,modres = NULL){
  
  Xt <- X[tind,,drop=FALSE]   
  Xc <- X[!tind,,drop=FALSE]
  
  Yt <- Y[tind]
  Yc <- Y[!tind]
  
  gamma <- length(C)
  
  if(is.null(modres)){
    modres <- mod_del_cal(gamma,C,Xt,Xc,kplus_ind,sigma_t,sigma_c)
  }
 
  bmat <- modres$b_mat
  dmat <- modres$d_mat
  
  if(lower == T){
    al_adj <- modres$alnewL
  }else{
    al_adj <- modres$alnewU
  }
  
  
  res <- CI_one_sd_RD(bmat,
                      dmat,
                      gamma,C,Xt,Xc,
                      kplus_ind,Yt,Yc,lower,
                      al = al_adj,Dir=F)
  
  return(res)
  
  
}
