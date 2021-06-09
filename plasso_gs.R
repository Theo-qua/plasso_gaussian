


S_func <- function(x, a) {  # Soft Thresholding Operator
  
  
  return( pmax(abs(x) - a,0) * sign(x))
}

#library(RGCCA)# soft threshold
#_func=soft.threshold






error.bars <-function(x, upper, lower, width = 0.02, ...) {
  xlim <- range(x)
  barw <- diff(xlim) * width
  segments(x, upper, x, lower, ...)
  segments(x - barw, upper, x + barw, upper, ...)
  segments(x - barw, lower, x + barw, lower, ...)
  range(upper, lower)
}



reg<-function(r,Z){
  my_one<-matrix(1,nrow(Z))
  my_w=data.frame(Z,my_one)
  my_w<-as.matrix(my_w)
  my_inv<-solve(t(my_w)%*%my_w)
  my_res<-my_inv%*%(t(my_w)%*%r)
  # new<- lm(r~1,na.action=na.exclude)
  beta0<-matrix(my_res[(K+1)])
  # new1<- lm(r~Z,singular.ok = TRUE)
  
  theta0<- matrix(my_res[c(1:(K))])
  
  return(list(beta0,theta0))
}








gradient_j<-function(beta,theta,U,U2,U3,y,X,W,r,alpha,lambda,K,j,N){
  
  theta_j<- matrix(as.numeric(theta[j,]),nrow = 1)
  
  dJ_dtheta_j<-matrix(0,nrow = 1,ncol = K)
  R<-c()
  
  if (isTRUE(any(c((beta[j]), theta_j)!=0))==T) {
    u<-as.numeric(beta[j])/norm(matrix(c(beta[j], theta_j)),type= "F")
    U<-c(U,u)
    
  } else{
    if(length(U)>=1){
      for (z in 1:length(U)) {
        if (norm(matrix(U[z]),type = "F")<=1) {
          R<-c(R,U[z])
          next(z)
        } else{
          next(z)
        }
        
        
      }
    }else{R=c(R,0)}
    u<-sample (c(R), size=1)
    U<-c(U,u)
  }
  
  
  
  
  dJ_dbeta_j=-(t(matrix(X[,j] )) %*% matrix((matrix(r)) ) )/N+(1-alpha)*lambda*u
  #DJ_BETA[j]<-dJ_dbeta_j
  
  
  
  
  
  #for (n in 1:K) {
  
  R<-c()
  if (isTRUE(any(c((beta[j]), theta_j)!=0))==T) {
    u2<-as.numeric(theta_j)/norm(matrix(c(as.numeric(beta[j]), theta_j)),type= "F")
    U2<-c(U2,u2)
    
  } else{
    if(length(U2)>=1){
      for (z in 1:length(U2)) {
        if (norm(matrix(U2[z]),type = "F")<=1) {
          R<-c(R,U2[z])
          next(z)
        } else{
          next(z)
        }
        
        
      }
      
    }else{R=c(R,0)}
    u2<-sample (c(R), size=1)
    U2<-c(U2,u2)
  }
  
  
  d<-c()
  if (isTRUE(any(c(theta_j)!=0))==T) {
    u3<-as.numeric(theta_j)/norm(matrix(c(theta_j)),type= "F")
    U3<-c(U3,u3)
    
  } else{
    if(length(U3)>=1){
      for (z in 1:length(U3)) {
        if (norm(matrix(U3[z]),type = "F")<=1) {
          d<-c(d,U3[z])
          next(z)
          
        } else{
          next(z)
        }
        
        
      }
    }else{d=c(d,0)}
    u3<-sample (c(d), size=1)
    U3<-c(U3,u3)
  }
  
  
  
  
  v<-sign(as.numeric(theta_j))
  
  
  
  
  
  
  
  
  dJ_dtheta_j= -((t(data.frame(W[j]) )%*% matrix(matrix((r) ) ) ))/N+(1-alpha)*lambda*(u2+u3)+alpha*lambda*v
  #}
  
  
  #DJ_THETA[j,]<-dJ_dtheta_j
  
  
  
  
  
  
  L1<- dJ_dbeta_j
  L2<- dJ_dtheta_j
  
  
  return(list(L1,L2,U,U2,U3))
}


concat_beta_theta<-function(beta, theta){
  p=length(beta); K =ncol(theta) 
  my_matrix = matrix(0,p, (K+1))
  my_matrix[, c(1:K)] = theta
  my_matrix[, (K+1)] = beta
  return (my_matrix)
  
}


count_nonzero<-function(x){
  n = 0
  for (i in 1: length(x)){
    if (x[i] != 0){
      n=n + 1
    }
  }
  return (n)
  
}



# 
# compute_w_j<-function(x_j, z,theta_t){
#   xz_theta<- (matrix(x_j, nrow = length(x_j), ncol = ncol(z)) * z)
#   
#   
#   w_j <- Reduce(f = '+', x = xz_theta)
#   
#   return(w_j)
#   
# }


compute_w<-function(X, Z){
  p=ncol(X)
  K=ncol(Z)
  N=nrow(X)
  W <- lapply(seq_len(p),
              function(j) (matrix(X[, j], nrow = N, ncol = K) * Z))
  return(W)
}



compute_pliable<-function(X, Z, theta){
  p=ncol(X)
  N=nrow(X)
  K=ncol(Z)
  
  xz_theta <- lapply(seq_len(p),
                     function(j) (matrix(X[, j], nrow = N, ncol = K) * Z) %*% t(theta)[, j])
  xz_term<- (Reduce(f = '+', x = xz_theta))
  
  return(xz_term)
  
  
}


model<-function(beta0, theta0, beta, theta, X, Z){
  p=ncol(X)
  N=nrow(X)
  K=ncol(Z)
  #The pliable lasso model described in the paper
  #y ~ f(X)
  
  #formulated as
  
  #y ~ b_0 + Z theta_0 + X b + \sum( w_j theta_ji )
  
  
  
  intercepts = as.numeric(beta0)+Z%*%(matrix(theta0))
  shared_model = X%*%matrix(beta)
  pliable = compute_pliable(X, Z, theta)
  return( intercepts + shared_model + pliable)
}



model_min_j<-function(beta0, theta0, beta, theta, X, Z, j, W ){
  
  #y ~ f(x) with X_j removed from the model
  
  
  beta[j] = 0.0
  theta[j, ] = 0.0
  return (model(beta0, theta0, beta, theta, X, Z))
  
}



model_j<-function(beta_j, theta_j, x_j, W, j, Z){
  
  
  #r_j ~ beta_j * X_j + W_j @ theta_j
  
  #Only a the residual fit on a single predictor is made
  
  
  # Caching is disabled
  #zz=matrix(compute_w_j(x_j, Z,matrix(theta_j,1)))
  w_j = as.matrix(data.frame(W[j]))
  zz<- w_j%*%theta_j
  return (beta_j * x_j + zz)
}


objective_ADMM<-function(beta0,theta0,beta,theta,X,Z,y,W,alpha,lambda,p,K,N){
  
  # theta_transpose<-t(theta)
  # 
  # xz_theta <- lapply(seq_len(p),
  #                    function(jj) (matrix(X[, jj], nrow = N, ncol = K) * Z) %*% theta_transpose[, jj])
  # XZ_term <- Reduce(f = '+', x = xz_theta)
  # n_l <-   as.numeric(beta0)+Z%*%(theta0) + X %*% beta  + XZ_term
  # n_l <- unlist(model(beta0, theta0, beta, theta, X, Z))
  mse = (1 / (2 * N)) * sum((y -model(beta0, theta0, beta, theta, X, Z))^2)
  
  

  coef_matrix = concat_beta_theta(beta, theta)

  #Numba does not support the axis kwarg on la.norm() [axis 1 means operates across rows]
  penalty_1=0; penalty_2 = 0
  for (jj in 1: p){
    penalty_1= penalty_1 + twonorm(matrix(coef_matrix[jj, ]))
    penalty_2=penalty_2 + twonorm(matrix(theta[jj, ]))
  }
  penalty_3 =sum(abs(theta))
  #
  objective_l = mse + (1-alpha) * lambda * (penalty_1 + penalty_2) + alpha * lambda * penalty_3

  
  
  
  
  # y_hat_l=n_l
  
  # pr_l<-matrix(sign(n_l))
  
  
  #pr[x]<-list(matrix(as.numeric(1/( exp(-n_i_b)+1),N)))
  
  
  
  
  
  
  
  
  #r_l<-hinge(y*y_hat_l)
  
  
  
  
  
  # norm_1_l=   lapply(seq_len(p),
   #                   function(g)(  lambda*(1-alpha)*(norm(matrix(c(beta[g],theta[g,])),type = "F") +norm(matrix(c(theta[g,])),type = "F") ) + lambda*alpha* sum(abs(theta[g,]))      ))
  # 
  
  
  
 # objective_l <- mse +sum(unlist(norm_1_l))
  
  
  #pr<-(sign(y_hat_l))
  
  
  # prob_d <-matrix((pr))
  #objective_l=apply(apply(prob_d, 2, FUN="!=", y), 2, sum)/length(y)
  
  return( objective_l)
}





penalties_min_j<-function(beta_0, theta_0, beta, theta, X, Z, y, W, j){
  # Compute the MSE, penalty 1 and penalty 2 when the jth predictor is not in the model.
  N=nrow(X)
  p=ncol(X)
  
  
  y_hat_min_j = model_min_j(beta_0, theta_0, beta, theta, X, Z, j, W)
  mse = (1/(2*N))/ sum((y-y_hat_min_j)^2)
  
  coef_matrix = concat_beta_theta(beta, theta)
  
  # Ignore the jth modifier from the model
  coef_matrix[j, ] = 0.0
  theta[j, ] = 0.0
  
  # Compute penalties
  penalty_1=0; penalty_2 = 0.0
  
  for (jj in 1: p){
    penalty_1= penalty_1 + twonorm(matrix(coef_matrix[jj, ]))
    penalty_2=penalty_2 + twonorm(matrix(theta[jj, ]))
  }
  penalty_3 = sum(abs(theta))
  
  
  return (list(mse, penalty_1, penalty_2, penalty_3))
  
}




objective_j<-function(beta0,theta0,beta,theta,X,Z,y,W,alpha,lambda,p,K,N,j,r_min_j){
  beta_j<-beta[j]
  theta_j<-theta[j,]
  theta_transpose<-t(theta_j)
  K=length(theta_j)
  # print(j)
  
  precomputed_penalties_minus_j = penalties_min_j(
    beta0, theta0, beta, theta, X, Z, y, W, j)
  # 
  # 
  mse<-as.numeric(unlist(precomputed_penalties_minus_j[1]))
  penalty_1<-as.numeric(unlist(precomputed_penalties_minus_j[2]))
  penalty_2<-as.numeric(unlist(precomputed_penalties_minus_j[3]))
  penalty_3<-as.numeric(unlist(precomputed_penalties_minus_j[4]))
  
  
  #r_min_j = ((model_min_j(beta0, theta0, beta, theta, X, Z, j, W )))
  
  # print(theta_j)
  # zz<- as.matrix(data.frame(W[j]))%*%theta_j
  #print(zz)
  r_hat = model_j(beta[j], theta[j,], X[, j], W, j, Z)
  mse_1 = (1 / (2*N)) * ( sum(  (r_min_j - r_hat )^2 )) +mse
  
  
  # Penalty 1
  coef_vector = matrix(0,(K+1))
  coef_vector[c(1:K)] = as.numeric(theta_j)
  coef_vector[(K+1)] = beta_j
  penalty_1 =penalty_1+ twonorm(matrix(coef_vector) )
  
  # Penalty 2
  penalty_2 =penalty_2+ twonorm(matrix(theta_j))
  
  # Penalty 3
  penalty_3 =penalty_3+ sum(abs(theta_j))
  
  objective_l = mse_1 + (1-alpha) * lambda * (penalty_1 + penalty_2) + alpha * lambda * penalty_3
  
  
  
  
  
  
  
  
  # norm_1_l=     lambda*(1-alpha)*(norm(matrix(c(beta_j,theta_j)),type = "F") +norm(matrix(c(theta_j)),type = "F") ) + lambda*alpha* sum(abs(theta_j))      
  
  #objective_l=apply(apply(prob_d, 2, FUN="!=", y), 2, sum)/length(y)
  
  
  
  return( objective_l)
}




###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
####################TRAIN TRUE WORK#################################

for_v=5;sv=0;fq=50;st=20;mv=50;ms=50;tol=1e-2;max_it=50;bb_it=100# default values

coordinate_descent<-function(y, X, Z, nlambda, alpha,new_t,my_mbeta,max_interaction_terms,maxgrid,lambda_min,my_lambda=NULL,tt=NULL,for_v,sv,fq,st,mv,ms,max_iter,bb_it,tol){
  
  
  X=as.matrix(X)
  Z=as.matrix(Z)
  
  N <- length(y)
  
  p <- ncol(X)
  K <- ncol(Z)
  
  rat=lambda_min
 
  
  
  
  
  
  beta0 = 0.0#estimates$Beta0
  theta0 = matrix(0,K)
  beta = matrix(0,p)
  theta = matrix(0,p,K)
  
  # Precomputed variables
  W<-compute_w(X, Z)
  
  if(is.null(my_lambda)){
    
   
    
    #print(r1)
    
    O_1=max(abs(t(X)%*%matrix(y))/(N*(1-alpha)) )
    #print(O_1)
    #lambda<-   max(O_1)
    big_lambda<-max(O_1)
    Lambda_min<- rat*big_lambda
    
    lambda_i<- exp(seq(log(big_lambda),log(big_lambda*rat),length=maxgrid))
    
    lambda_i[1]<-big_lambda;lambda_i[maxgrid]<-Lambda_min
    
    
    print(lambda_i)
  }else{
    lambda_i<-my_lambda
    # r_current = hinge(y*model(beta_0, theta_0, beta, theta, X, Z))
    # b = reg(r_current,Z)  # Analytic solution how no sample lower bound (Z.T @ Z + cI)^-1 @ (Z.T @ r)
    # beta_0<-matrix(unlist(b[1]))
    # theta_0<-matrix(unlist(b[2])) 
    # 
    # 
    # n_i_b<-model(beta_0, theta_0, beta, theta, X, Z)
    # r1<-(hinge(x=y*n_i_b))
    # print()
    
  }
  nlambda=length(lambda_i)
  
  U<-c(); U2<-c(); U3<-c()
  
  
  
  # Solve ABG Parameters
  #tt=1/max(eigen((t(X)%*%X)/N)$values)
  
  #tt = 0.1 / mean(X**2)
  # Lists
  lam_list = c()
  beta_0_list = c()
  theta_0_list = list()
  beta_list = list()
  theta_list = list()
  DEV=c()
  n_main_terms=c()
  non_zero_theta=c()
  
  loss<-matrix(0,N,nlambda)
  
  
  
  active_set1<-c()
  active_set2<-c()
  
  strong_set<-matrix(0,p,1)
  strong_set1<-matrix(0,p,1)
  my_q=1
  my_v<-sv
  my_V<-for_v
  
  q=1
  tolerance = tol
  my_ok = 0
  
 
  
  
  
  beta0 = 0.0#estimates$Beta0
  theta0 = matrix(0,K)
  beta = matrix(0,p)
  theta = matrix(0,p,K)
  
  #i=q
  while (isTRUE(my_ok==0)==T  ){
    
    
    
    if( (q<=1 | isTRUE(q<=fq)==T)==T  ){
      
      i=q
      lam=lambda_i[i]
      my_v<-my_v
      my_V<-my_V
      
      r_current = (y-model(beta0, theta0, beta, theta, X, Z))
      #b = reg(r_current,Z)  # Analytic solution how no sample lower bound (Z.T @ Z + cI)^-1 @ (Z.T @ r)
      #beta_0<-matrix(unlist(b[1]))
      # theta_0<-matrix(unlist(b[2]))
      # beta_1<-beta; theta_1<-theta
      # Iterate through all p features
      
      if(i>40 & max_it>1000){max_iter=1000}
      beta_sup<-matrix(NA,p)
      for( iii in 1: max_iter){
        # iter_prev_score = objective(
        #   beta0, theta0, beta, theta,
        #   X, Z, y,
        #   alpha, lam, W
        # )
        
        iter_prev_score= objective_ADMM(
          beta0,theta0,beta,theta,
          X,Z,y,W,alpha,lambda=lam,p,K,N)
        
        beta01=beta0;theta01=theta0;beta1=beta;theta1=theta
        # Compute beta_0 and theta_0 from the least square regression of the current residual on Z
        # Z + 1s = W matrix where Z is for theta_0 and 1s is for beta_0
        
        r_current = (y-model(beta0, theta0, beta, theta, X, Z))
        b = reg(r_current,Z)  # Analytic solution how no sample lower bound (Z.T @ Z + cI)^-1 @ (Z.T @ r)
        beta0<-matrix(unlist(b[1]))
        theta0<-matrix(unlist(b[2]))
        
        
        r = y-model(beta0, theta0, beta, theta, X, Z)
        
        #rr[,1]=r1;rr[,2]=r
        #print(as.data.frame(matrix(r1,N,1),matrix(r,N,1)))
       
        
        my_p<-sample(c(1:p))
        my_p<-c(my_p)
        
        for (jj in 1: p){
          j=jj#my_p[jj]
          # Strong screening rule (https://statweb.stanford.edu/~tibs/ftp/strong.pdf)
          if(i<=1){checkk=2*lam-lambda_i[1] }else{checkk=( 2*lam-lambda_i[1] )}
          
          
          
          # print(c(i,as.numeric(abs( (t(matrix(X[,j]*c(y) ))%*% matrix((r1) ) )  )/(N*alpha) ), checkk ))
          
          if(as.numeric(abs( (t(matrix( (X[,j]) ) )%*% matrix(((r))   ) )  )/(N*alpha) )< checkk  ){
            next(jj)
          }else{ 
            
            
            
            w_j = ( data.frame(W[j]) )
            r_min_j =( matrix((r + model_j(beta[j], theta[j, ], X[, j], W, j, Z))  ) )
           
            
            #r_min_j = hinge(matrix(y*model_min_j(beta0, theta0, beta, theta, X, Z, j, W )))
            # Check if beta_j == 0 and theta_j == 0
            cond_17a = abs(t(matrix( (X[,j]) ) )%*% r_min_j / N)
            #print( cond_17a)
            cond_17b = as.numeric(twonorm(matrix(S_func((t( (w_j) )%*% r_min_j   )/N, alpha*(lam)) )))
            
            
            
            if (cond_17a<= (1-alpha)*lam){
              
              strong_set[j]<-  0
            } else{
              strong_set[j]<-1
            }
            
            if (cond_17b<= 2*(1-alpha)*(lam)){
              
              strong_set1[j]<-0
            } else{
              strong_set1[j]<-1
            }
            
            
            
            
            
            if(isTRUE(strong_set[j]==0 & strong_set1[j]==0)==T ){
              
              #beta[j]<-0;theta[j,]<-0
              next(jj)
              
            } else{
              # print(c(iii,j))
              beta_j_hat = (N/sum( ( matrix( ( (X[,j]^2 ) )  ) ) ))*  S_func((t( matrix( (X[,j]) ))%*%r_min_j )/N,(1-alpha)*lam ) 
              #print(c(beta_j_hat,cond_17a,(1-alpha)*lam))
              
              # r_min_j = model_min_j(beta0, theta0, beta, theta, X, Z, j, W )
              
              cond_19 = twonorm(matrix(S_func((t( (data.frame(W[j]) ) )%*%matrix( (( (r_min_j-matrix(X[,j])*as.numeric(beta_j_hat) ) ) ) )  )/N,alpha*lam)))
              
              
                
                if (isTRUE(cond_19<= (1-alpha)*lam)==T  ) {
                  # beta_j != 0 and theta_j == 0
                  beta[j]<-as.numeric(beta_j_hat)
                  #theta[j,]<-as.numeric(0)
                  active_set1=c(active_set1,j)
                  
                  
                  next(jj)
                  
                  
                }else{ 
                  
                  
                  # r_min_j =( matrix( (r - model_j(beta[j], theta[j, ], orig_X[, j], orig_W, j, orig_Z))) )
                 
                  
                 # grad_j<-gradient_j(beta,theta,U,U2,U3,y,X,W,r,alpha,lambda=lam,K,j,N)
                 # L1=matrix(unlist(grad_j[1]));L2<-matrix(unlist(grad_j[2]))
                  
                #  U=c(unlist(grad_j[3]));U2=c(unlist(grad_j[4]));U3=c(unlist(grad_j[5]))
                  
                  value1<-   quadratic(beta,theta, alpha,lambda=lam,beta0,theta0,j,b=1,W,X=X,Z=Z,y=y,N,n_i=n_i_b,big_delta_1=NULL,big_delta_2=NULL,b_1=NULL,b_2=NULL,t=t,bb_it = bb_it,r_min_j,tol,U,U2,U3)
                  #value1<-     quadratic(L1, L2, beta_j=as.numeric(beta[j]),theta_j=as.numeric(theta[j,]), alpha,lambda,t,mbeta,v_beta[j],v_theta[j,],my_n=(i-1))
                  
                  beta[j]<-unlist(value1[[1]]);theta[j,]<-unlist(value1[[2]]);t=unlist(value1[[3]]);U<-unlist(value1[[4]]);U2<-unlist(value1[[5]]);U3<-unlist(value1[[6]])
                  active_set2=c(active_set2,j)
                  active_set1=c(active_set1,j)
                  # 
                  # t=1
                  # for (ff in 1:bb_it){  # Max steps
                  #   r = (r_min_j + model_j(beta[j], theta[j, ], X[, j], W, j, Z) )
                  # 
                  # 
                  #
                  
                }####beta_j and theta_j not zero
                
              ####beta_j not zero
            }####first decision cond17
          }
        }####j
        
        #beta0=beta01;theta0=theta01;beta=beta1;theta=theta1
      
        if(i>new_t){
          # iter_current_score = objective(
          #   beta0, theta0, beta, theta,
          #   X, Z, y=orig_y,
          #   alpha, lam, W
          # )
          
          iter_current_score= objective_ADMM(beta0,theta0,
                                        beta,theta,X,Z,
                                        y,W,alpha,lambda=lam,p,K,N)
          #print(c(iter_prev_score,iter_current_score,beta))
          if(isTRUE(abs(iter_prev_score - iter_current_score) <tolerance)==T ){
            # print(c(iter_prev_score,iter_current_score))
            break  # Converged on lam_i
            
          }else{next(iii)}
          
        }else{break}
        
        
        
        
      }#max_it
      
      
      n_i_b<-(model(beta0, theta0, beta, theta, X, Z))
      
      
      
      
      
      
      
      #pr[max(y)]<- list(matrix(1- rowSums(my_pr),N))
      
      
      #v1_d<-1*(v1==d)
     
      
      # }
      Dev<-  sum((y-n_i_b)^2)/N
      #Dev<-apply(apply(prob_d, 2, FUN="!=", orig_y), 2, sum)/length(orig_y)
      
      
      #(y,yhat,w=rep(1,length(y)))
      # }
      
      
      
      
      
      
      
      
      
      # Check maximum interaction terms reached. If so early stop just like Tibs.
      n_interaction_terms = count_nonzero(c(theta))
      
      if( n_interaction_terms >= max_interaction_terms){
        print('Maximum Interaction Terms reached.')
        my_ok = 1
        break
      }
      # Save coefficients
      
      DEV<-(c(DEV,Dev))
      err=1e-5
      
      
      
      
      
      
      n_main_terms = (c(n_main_terms,count_nonzero(c(beta)) ) )
      
      
      
      non_zero_theta<- (c(non_zero_theta,count_nonzero(c(theta))))
      print(c(i,iii,(n_main_terms[i]),non_zero_theta[i] , DEV[i]))
      
      
      
      
      lam_list<-(c(lam_list,lam))
      beta_0_list<-(c(beta_0_list,beta0))
      theta_0_list[[i]]<-list(matrix(theta0,1,K))
      beta_list[[i]]<-beta
      theta_list[[i]]<-list(as.matrix(theta,p,K))
      my_q<-my_q+1
      q=q+1
      #active set 
    }else{
      i=q
      lam=lambda_i[i]
      my_v<-my_v
      my_V<-my_V
      
      # 
      r_current = (y-model(beta0, theta0, beta, theta, X, Z))
      # b = reg(r_current,Z)  # Analytic solution how no sample lower bound (Z.T @ Z + cI)^-1 @ (Z.T @ r)
      # beta_0<-matrix(unlist(b[1]))
      # theta_0<-matrix(unlist(b[2]))
      # beta_1<-beta; theta_1<-theta
      # Iterate through all p features
      
      if(i>40 & max_it>1000){max_iter=1000}
      
      #for( iii in 1: max_iter){
      # iter_prev_score = objective(
      #   beta0, theta0, beta, theta,
      #   X, Z, y,
      #   alpha, lam, W
      # )
      
      iter_prev_score= objective_ADMM(
        beta0,theta0,beta,theta,
        X,Z,y,W,alpha,lambda=lam,p,K,N)
      
      beta01=beta0;theta01=theta0;beta1=beta;theta1=theta
      # Compute beta_0 and theta_0 from the least square regression of the current residual on Z
      # Z + 1s = W matrix where Z is for theta_0 and 1s is for beta_0
      r_current = (y-model(beta0, theta0, beta, theta, X, Z))
      b = reg(r_current,Z)  # Analytic solution how no sample lower bound (Z.T @ Z + cI)^-1 @ (Z.T @ r)
      beta0<-matrix(unlist(b[1]))
      theta0<-matrix(unlist(b[2]))
      
      
      r = y-model(beta0, theta0, beta, theta, X, Z)
   
   
      my_p<-sample(c(1:p))
      my_p<-c(my_p)
      
      for (jj in 1: p){
        j=jj#my_p[jj]
        # SAFE screening rule (https://statweb.stanford.edu/~tibs/ftp/strong.pdf)
        if(i<=1){checkk=2*lam-lambda_i[1] }else{checkk=( 2*lam-lambda_i[i-1] )}
        
        
        
        # print(c(i,as.numeric(abs( (t(matrix(X[,j]*c(y) ))%*% matrix((r1) ) )  )/(N*alpha) ), checkk ))
        
        if(as.numeric(abs( (t(matrix(X[,j] ))%*% matrix((r) ) )  )/(N*alpha) )< checkk  ){
          next(jj)
        }else{ 
          
          
          
          w_j = ( data.frame(W[j]) )
          r_min_j =( matrix((r + model_j(beta[j], theta[j, ], X[, j], W, j, Z))) )
         
          #r_min_j = hinge(matrix(y*model_min_j(beta0, theta0, beta, theta, X, Z, j, W )))
          # Check if beta_j == 0 and theta_j == 0
          cond_17a = abs(t(matrix(X[,j] ) )%*% r_min_j / N)
          #print( cond_17a)
          cond_17b = as.numeric(twonorm(matrix(S_func((t( (w_j) )%*% r_min_j   )/N, alpha*(lam)) )))
          
          
          
          if (cond_17a<= (1-alpha)*lam){
            
            strong_set[j]<-  0
          } else{
            strong_set[j]<-1
          }
          
          if (cond_17b<= 2*(1-alpha)*(lam)){
            
            strong_set1[j]<-0
          } else{
            strong_set1[j]<-1
          }
          
          
          
          
          
          if(isTRUE(strong_set[j]==0 & strong_set1[j]==0)==T ){
            
            #beta[j]<-0;theta[j,]<-0
            next(jj)
            
          } else{
            # print(c(iii,j))
            beta_j_hat = (N/sum( ( matrix( ( (X[,j]^2 ) )  ) ) ))*  S_func((t( matrix( (X[,j]) ))%*%r_min_j )/N,(1-alpha)*lam ) 
            #print(c(beta_j_hat,cond_17a,(1-alpha)*lam))
            
            # r_min_j = model_min_j(beta0, theta0, beta, theta, X, Z, j, W )
            
            
            
            
            
            
            
              
              #r=1-y*model(beta0, theta0, beta, theta, X, Z)
             # r_min_j =( matrix( (r + model_j(beta2[j], theta2[j, ], X[, j], W, j, Z)) ) )
              #r_min_j = model_min_j(beta0, theta0, beta, theta, X, Z, j, W )
              cond_19 = twonorm(matrix(S_func((t( (w_j) )%*%matrix( ((r_min_j -matrix(X[,j])*as.numeric(beta_j_hat) ) ) )  /N),alpha*lam)))
              
              
              
              
              if (isTRUE(cond_19<= (1-alpha)*lam)==T  ) {
                # beta_j != 0 and theta_j == 0
                # beta[j]<-as.numeric(beta_j_hat)
                #theta[j,]<-as.numeric(0)
                active_set1=c(active_set1,j)
                
                
                next(jj)
                
                
              }else{ 
                
                
               
                active_set2=c(active_set2,j)
                active_set1=c(active_set1,j)
                # 
                # t=1
                # for (ff in 1:bb_it){  # Max steps
                #   r = (r_min_j + model_j(beta[j], theta[j, ], X[, j], W, j, Z) )
                # 
                # 
                #
                
              }####beta_j and theta_j not zero
              
            ####beta_j not zero
          }####first decision cond17
        }
      }####j
      
      
      
      while (my_v<=my_V) {
        if(isTRUE(q>mv)==T){
          my_V=ms}else{my_V=my_V}
        print(c(my_v,my_V))
        i=q
        lam<-lambda_i[q]
        #print(c(q,lambda))
        
        my_v<-my_v
        my_V<-my_V
        
        # 
        r_current = (y-model(beta0, theta0, beta, theta, X, Z))
        # b = reg(r_current,Z)  # Analytic solution how no sample lower bound (Z.T @ Z + cI)^-1 @ (Z.T @ r)
        # beta_0<-matrix(unlist(b[1]))
        # theta_0<-matrix(unlist(b[2]))
        # beta_1<-beta; theta_1<-theta
        # Iterate through all p features
       
        if(i>40 & max_it>1000){max_iter=1000}
        beta_sup<-matrix(NA,p)
        for( iii in 1: max_iter){
          # iter_prev_score = objective(
          #   beta0, theta0, beta, theta,
          #   X, Z, y,
          #   alpha, lam, W
          # )
          
          iter_prev_score= objective_ADMM(
            beta0,theta0,beta,theta,
            X,Z,y,W,alpha,lambda=lam,p,K,N)
          
          beta01=beta0;theta01=theta0;beta1=beta;theta1=theta
          # Compute beta_0 and theta_0 from the least square regression of the current residual on Z
          # Z + 1s = W matrix where Z is for theta_0 and 1s is for beta_0
          r_current = (y-model(beta0, theta0, beta, theta, X, Z))
          b = reg(r_current,Z)  # Analytic solution how no sample lower bound (Z.T @ Z + cI)^-1 @ (Z.T @ r)
          beta0<-matrix(unlist(b[1]))
          theta0<-matrix(unlist(b[2]))
          
          
          r =y- model(beta0, theta0, beta, theta, X, Z)
          
         
          
          # W <- lapply(seq_len(p),
          #            function(j) (matrix(X[, j], nrow = nrow(X), ncol = K) * Z))
          # # # # ##
          #   
          #r=1-y*model(beta0, theta0, beta, theta, X, Z)
          
          
          
          
          my_active_set1= as.matrix(as.numeric(levels(factor(active_set1))))
          my_active_set2= as.matrix(as.numeric(levels(factor(active_set2))))
          if(length(my_active_set1)>=1){
            for (jj in 1: length(my_active_set1)){
              j=my_active_set1[jj]
              # print(c(iii,j))
              
              
              
              w_j = ( data.frame(W[j]) )
              r_min_j =( matrix((r + model_j(beta[j], theta[j, ], X[, j], W, j, Z))) )
              
              
              # r_min_j = hinge(matrix(y*model_min_j(beta0, theta0, beta, theta, X, Z, j, W )))
              # Check if beta_j == 0 and theta_j == 0
              cond_17a = abs(t(matrix(X[,j] ) )%*% r_min_j / N)
              #print( cond_17a)
              cond_17b = as.numeric(twonorm(matrix(S_func((t( (w_j) )%*% r_min_j   )/N, alpha*(lam)) )))
              
              
              
              if (cond_17a<= (1-alpha)*lam){
                
                strong_set[j]<-  0
              } else{
                strong_set[j]<-1
              }
              
              if (cond_17b<= 2*(1-alpha)*(lam)){
                
                strong_set1[j]<-0
              } else{
                strong_set1[j]<-1
              }
              
              
              
              
              
              if(isTRUE(strong_set[j]==0 & strong_set1[j]==0)==T ){
                
                #beta[j]<-0;theta[j,]<-0
                next(jj)
                
              } else{
                # print(c(iii,j))
                beta_j_hat = (N/sum( ( (matrix(X[,j]^2))) ) )*  S_func((t(matrix(X[,j] ))%*%r_min_j )/N,(1-alpha)*lam ) 
                #print(c(beta_j_hat,cond_17a,(1-alpha)*lam))
                
                # r_min_j = model_min_j(beta0, theta0, beta, theta, X, Z, j, W )
              
                
                
             
                  
                 
                  cond_19 = twonorm(matrix(S_func((t( (w_j) )%*%matrix( ((r_min_j -matrix(X[,j])*as.numeric( beta_j_hat) ) ) )  /N),alpha*lam)))
                  
                  
                  
                  
                  if (isTRUE(cond_19<= (1-alpha)*lam)==T  ) {
                    # beta_j != 0 and theta_j == 0
                    beta[j]<-as.numeric(beta_j_hat)
                    #theta[j,]<-as.numeric(0)
                    active_set1=c(active_set1,j)
                    
                    
                    
                    
                    
                    next(jj) }else{ 
                      
                    
                      # 
                    #  grad_j<-gradient_j(beta,theta,U,U2,U3,y,X,W,r,alpha,lambda=lam,K,j,N)
                      # 
                     # L1=matrix(unlist(grad_j[1]));L2<-matrix(unlist(grad_j[2]))
                      # 
                    #  U=c(unlist(grad_j[3]));U2=c(unlist(grad_j[4]));U3=c(unlist(grad_j[5]))
                      # 
                      value1<-   quadratic(beta,theta, alpha,lambda=lam,beta0,theta0,j,b=1,W,X=X,Z=Z,y=y,N,n_i=n_i_b,big_delta_1=NULL,big_delta_2=NULL,b_1=NULL,b_2=NULL,t=t,bb_it = bb_it,r_min_j,tol,U,U2,U3)
                      #value1<-     quadratic(L1, L2, beta_j=as.numeric(beta[j]),theta_j=as.numeric(theta[j,]), alpha,lambda,t,mbeta,v_beta[j],v_theta[j,],my_n=(i-1))
                      
                      beta[j]<-unlist(value1[[1]]);theta[j,]<-unlist(value1[[2]]);t=unlist(value1[[3]]);U<-unlist(value1[[4]]);U2<-unlist(value1[[5]]);U3<-unlist(value1[[6]])
                      active_set2=c(active_set2,j)
                      active_set1=c(active_set1,j)
                      
                      
                      next(jj)}####beta_j and theta_j not zero
                  
                ####beta_j not zero
              }####first decision cond17
            }####j
            # print(c(iii,j,beta[j]))
            
          }else{
            my_v=my_V
            iii=max_it
            break
          }
        
          if(i>2){
            # iter_current_score = objective(
            #   beta0, theta0, beta, theta,
            #   X, Z, y=orig_y,
            #   alpha, lam, W
            # )
            
            iter_current_score= objective_ADMM(beta0,theta0,
                                          beta,theta,X,Z,
                                          y,W,alpha,lambda=lam,p,K,N)
            if(isTRUE((iter_prev_score - iter_current_score) <tolerance)==T  ){
              # print(c(iter_prev_score,iter_current_score))
              break  # Converged on lam_i
              
            }else{next(iii)}
            
          }else{break}
          
          
          
          
        }#max_it
        
        
        n_i_b<-(model(beta0, theta0, beta, theta, X, Z))
        
        
        Dev<-  sum((y-n_i_b)^2)/N
        #Dev<-apply(apply(prob_d, 2, FUN="!=", orig_y), 2, sum)/length(orig_y)
        
        
        #(y,yhat,w=rep(1,length(y)))
        # }
        
        
        
        
        
        
        
        
        
        # Check maximum interaction terms reached. If so early stop just like Tibs.
        n_interaction_terms = count_nonzero(c(theta))
        
        if( n_interaction_terms >= max_interaction_terms){
          print('Maximum Interaction Terms reached.')
          my_ok = 1
          break
        }
        # Save coefficients
        
        DEV<-(c(DEV,Dev))
        err=1e-5
        
        
        
        
        
        
        
        n_main_terms = (c(n_main_terms,count_nonzero(c(beta)) ) )
        
        
        
        non_zero_theta<- (c(non_zero_theta,count_nonzero(c(theta))))
        print(c(i,iii,(n_main_terms[i]),non_zero_theta[i] , DEV[i]))
        
        
        
        
        lam_list<-(c(lam_list,lam))
        beta_0_list<-(c(beta_0_list,beta0))
        theta_0_list[[i]]<-list(matrix(theta0,1,K))
        beta_list[[i]]<-beta
        theta_list[[i]]<-list(as.matrix(theta,p,K))
        
        
        
        
        
        
        
        my_v<-my_v+1
        q=q+1
        print(q)
        if( isTRUE(q>nlambda)==T){
          #my_ok<-1
          print(q)
          break
        }
        
      }#for while loop 
      
      
    }### end of active set 
    
    
    
    my_v<-sv
    my_V<-for_v
    q=q
    i=q
    if(isTRUE(q>nlambda)==T | n_interaction_terms >= max_interaction_terms  ){
      my_ok=1}else{my_ok=0}
  }#while
  pred<-data.frame(Lambda=matrix(lam_list,ncol = 1),nzero=n_main_terms,nzero_inter=non_zero_theta,DEV=matrix(DEV,ncol=1))
  
  # Return results
  return (list(beta0=beta_0_list,beta=beta_list,theta0=theta_0_list,theta=theta_list,path=pred,Lambdas=lam_list,non_zero=n_main_terms))
  
  
  
  
}













pliable_predict_lasso<-function(object ,X,Z,y,lambda=NULL){
  lambda.arg=lambda
  if(is.null(lambda.arg))
  { lambda=object$Lambdas;isel=1:length(lambda)}
  
  if(!is.null(lambda.arg)){
    
    isel=as.numeric(knn1(matrix(object$Lambdas,ncol=1),matrix(lambda.arg,ncol=1),1:length(object$Lambdas)))
    
  }
  
  # print(c(isel,length(isel)))
  
  N <- nrow(X)
  
  p <- ncol(X)
  #print(Z)
  K <- ncol(as.matrix(Z))
  #print(K)
  Z=matrix(as.numeric(Z),N,K)
  
  yh=array(0,length(isel))
  DEV=matrix(NA,length(isel))
  my_theta<-array(0,c(ncol(X),ncol(Z), length(isel)))
  
  
  #pBETA0<-matrix(0,nrow = length(isel)); pBETA<-matrix(0,nrow = length(isel)); pTHETA0<-matrix(0,nrow = length(isel)); pTHETA<-matrix(0,nrow =length(isel))
  
  
  pBETA0<-lapply(seq_len(1),
                 function(j)(matrix(0,nrow = length(isel))))
  
  pBETA<-lapply(seq_len(1),
                function(j)(matrix(0,nrow=p,ncol=length(isel))))
  
  pTHETA0<-lapply(seq_len(1),
                  function(j)(matrix(0,nrow=K,ncol=length(isel))))
  
  pTHETA<-lapply(seq_len(1),
                 function(j)(array(0,c(p,K,length(isel)))))
  
  
  
  
  ii=0
  for(m in isel){
    ii=ii+1
    
    # pred_beta<-lapply(seq_len(max(y)),
    #                  function(j)(matrix(0,nrow = 1,ncol = ncol(X))))
    #pred_theta<-lapply(seq_len(max(y)),
    #                    function(j)(matrix(0,nrow = p,ncol = K)))
    #  pred_theta0<-lapply(seq_len(max(y)),
    #                    function(j)(matrix(0,nrow = 1,ncol = K)))
    
    #  pred_beta0<-lapply(seq_len(max(y)),
    #                  function(j)(0))
    
    
    z=m
    #BETA0<-matrix(unlist(object$beta0[z]),1,max(y))
    
    
    # BETA<-as.data.frame(object$beta[z])
    # THETA0<-as.data.frame(object$theta0[z])
    # THETA<-object$theta[z]
    
    
    
    
    n_i<-lapply(seq_len(1),
                function(j)(matrix(0,nrow = N)))
    pr<-lapply(seq_len(1),
               function(j)(matrix(0,nrow = N)))
    
    
    
    for (x in 1:1) {
      #theta<- matrix(unlist(THETA[[1]][x]),p,K)
      
      
      beta0<-as.numeric(object$beta0[z])
      beta <- as.numeric(unlist(object$beta[z]))
      
      theta <- matrix(unlist(object$theta[[z]]),p,K) ; theta_transpose <- t(theta)
      
      
      theta0 <- unlist(object$theta0[[z]])
      
      
      pBETA0[[x]][ii] <-beta0
      pBETA[[x]][,ii] <-beta
      pTHETA[[x]][,,ii] <-theta
      pTHETA0[[x]][,ii] <-theta0
      
      
      
      
      
      xz_theta <- lapply(seq_len(p),
                         function(j) (matrix(X[, j], nrow = N, ncol = K) * Z) %*% t(theta)[, j])
      xz_term<- (Reduce(f = '+', x = xz_theta))
      
      # pred_theta[x]<-list(as.matrix(theta,p,K))
      # pred_beta0[x]<-beta0
      # pred_theta0[x]<-theta0
      # pred_beta[x]<-beta
      # beta=matrix(unlist(BETA[x]),p)
      #n_i<-beta0*as.numeric(beta0)*matrix(1,nrow = N,ncol = 1)+Z%*%(theta0)+X%*%matrix(beta) + XZ_term
      n_i[x]<-list(as.numeric((beta0))+Z%*% ((theta0))+X %*%(beta) + xz_term)
      
    }
    
    
    
    # pr[max(y)]<- list(matrix(1- rowSums(my_pr),N))
    
    #v1_d<-1*(v1==d)
    
    
    
    #pr[max(y)]<- list(matrix(1- rowSums(my_pr),N))
    
    
    #v1_d<-1*(v1==d)
    
    # n_i_d<-matrix(unlist(n_i[d]),N)
    
    # deviance_y1[d]<-y_d*n_i_d[l]
    # deviance_y2[d]<-exp(n_i_d[l])
    
    # }
    Dev<-  sum((y-unlist(n_i[1]) )^2)/N
    #Dev<-apply(apply(prob_d, 2, FUN="!=", y), 2, sum)/length(y)
    
    
    #Dev1[l]<-  sum(deviance_y2)
    
    
    #(y,yhat,w=rep(1,length(y)))
    #}
    #DEV1[i]<-((2)*sum(Dev1))/length(y)
    DEV[ii]<-(Dev)
    
    # DEV[ii]<- ((2)*sum(Dev))#+
    yh[ii]<- list(matrix(unlist(n_i[1]),N,1))
    
  }
  return(list(y_hat=yh,beta0=pBETA0,beta=pBETA,theta0=pTHETA0,theta=pTHETA,deviance=DEV))
}

