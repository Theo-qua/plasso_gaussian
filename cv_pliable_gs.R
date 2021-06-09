cv_pliable<-function(fit,nfolds,X,Z,y,alpha,nlambda=length(fit$Lambdas),new_t,maxgrid,tol,lambda_min,for_v,sv,fq,st,mv,ms,max_it=1,bb_it,foldid=NULL){
  BIG=10e9
  no<-nrow(X)
  ni<-ncol(X)
  nz<-ncol(Z)
  ggg=vector("list",nfolds)
  
  yhat=array(NA,c(no,length(fit$Lambdas)))
  
  # yhat=matrix(0,nfolds,length(result$Lambdas))
  my_nzero<-matrix(0,nfolds,length(fit$Lambdas))
  
  
  if(is.null(foldid)) foldid = sample(rep(1:nfolds, ceiling(no/nfolds)), no, replace=FALSE)  #foldid = sample(rep(seq(nfolds), length = no))
  
  nfolds=length(table(foldid))
  
  status.in=NULL
  
  for(ii in 1:nfolds){
    print(c("fold,", ii))
    oo=foldid==ii
    
   
    
    ggg[[ii]]<-    coordinate_descent(y=y[!oo], X=X[!oo,,drop=F], Z=Z[!oo,,drop=F], nlambda=length(fit$Lambdas), alpha = alpha,new_t=new_t,my_mbeta=.09,max_interaction_terms=500,maxgrid=50,lambda_min=.001,my_lambda=fit$Lambdas,tt=NULL,for_v,sv,fq,st,mv,ms,max_iter=max_it,bb_it=bb_it,tol=tol)
    
    cv_p<-pliable_predict_lasso(ggg[[ii]] ,X=X[oo,,drop=F],Z=Z[oo,],y=y[oo])
    
    
    #  Coeff<-lapply(seq_len(max(y)),
    #               function(j)(matrix(0,ncol(X),nlambda)))
    
    
    #for (i in 1:nlambda) {
    #  f<-Matrix(unlist(ggg[[ii]]$beta[i]),ncol(X),max(y),sparse = T)
    # for (j in 1:max(y)) {
    #  Coeff[[j]][,i]<-f[,j]
    #}
    
    # }
    
    #nnn_zero<-matrix(0,max(y),nlambda)
    
    #for (i in 1:max(y)) {
    #  q<-matrix(unlist(Coeff[[i]]),ncol(X),nlambda)
    # nnn_zero[i,]<-colSums(q!=0)
    
    
    # }
    #non_zero<-matrix(0,nlambda)
    #for (i in 1:nlambda) {
    # non_zero[i]<-max(nnn_zero[,i])
    #}
    
    #print(cv_p$deviance)
    # my_nzero[ii,]<-non_zero
    yhat[oo,]<- unlist(cv_p$y_hat)
    #(result ,X,Z,y,lambda=NULL)
  }
  
  #print(yhat)
  ym=array(y,dim(yhat))
  #print(ym)
  err=( errfun.gaussian(y=ym,yhat=yhat))
  
  
  
  non_zero<-matrix(0,nlambda)
  
  non_zero<-c(fit$path$nzero)
  
  cvm=apply(err,2,mean,na.rm=T)
  nn=apply(!is.na(err),2,sum,na.rm=T)
  cvsd=sqrt(apply(err,2,var,na.rm=T)/nn)
  cvm.nz=cvm; cvm.nz[non_zero==0]=BIG
  imin=which.min(cvm.nz)
  imin.1se=which(cvm< cvm[imin]+cvsd[imin])[1]
  
  out=list(lambda=fit$Lambdas,cvm=cvm,cvsd=cvsd,cvup = cvm +
             cvsd, cvlo = cvm - cvsd, nz=c(fit$path$nzero),lambda.min=fit$Lambdas[imin],lambda.1se=fit$Lambdas[imin.1se])
  
  
  return(out) 
}



quad_solution<-function(u, v, w){
  temp = ((v^2) - (4 * u * w))^0.5
  root1 = (-v + temp) / (2 * u)
  root2 = (-v - temp) / (2 * u)
  roots<-list(root1, root2)
  return (roots)
}

##############
####### iterations 



quadratic=function(beta ,theta, alpha,lambda,beta0,theta0,j,b,W,X,Z,y,N,n_i,big_delta_1=NULL,big_delta_2=NULL,b_1=NULL,b_2=NULL,t=NULL,bb_it,r_min_j,tol,U,U2,U3){
  # print(c(b,j,t))
  
  big=10e9; eps = 1e-5
  
 # U<-matrix(0); U2<-matrix(0); U3<-matrix(0)
  
  
  t=1;nes=1; beta_u=beta;theta_u=theta;beta_new=beta;theta_new=theta;n_l=n_i;n_r=n_i; v_beta=beta;v_theta=theta;beta_uu=beta;theta_uu=theta;old_beta=beta;old_theta=theta
  j=j
  
  # U<-matrix(0); U2<-matrix(0); U3<-matrix(0)
  # N=(nrow(X))
  
  okay=0
  #  teta_k=2/(nes+1)
  #  beta_new=(1-teta_k)*old_beta+teta_k*beta_uu
  # theta_new=(1-teta_k)*old_theta+teta_k*theta_uu
  #
  # #
  # old_beta_u=beta_u; old_theta_u=theta_u
  # teta_k<-2/(nes+1)
  # #
  #  beta_new<-(1-teta_k)*old_beta_u+teta_k*v_beta
  #  theta_new<-(1-teta_k)*old_theta_u+teta_k*v_theta
  
  
  # 
  # 
  # beta_j=as.numeric(beta_new[j]);theta_j=as.numeric(theta_new[j,])
  # 
  # a<- as.numeric(norm(matrix(beta_j), type="F")) ; bb<-as.numeric(twonorm(matrix(theta_j))) ; rho_2<- as.numeric(sqrt(a^2+bb^2))
  # 
  # 
  # ##### after norming 
  # c<- as.numeric(t*(1-alpha)*lambda );  g_1<- as.numeric( abs(beta_j-t*as.numeric(L1)))
  # v<-matrix(0,nrow = 1,ncol = K)
  # 
  # # for (b in 1:K) {
  # v<-S_func(theta_j-t*L2,t*alpha*lambda)
  # # }
  # 
  # g_2<-as.numeric(twonorm(matrix(v)))
  # 
  # #posroot<- -c+sqrt(c^2+2*c*g_2+g_1^2+g_2^2)
  # 
  # #negroot<- -c-sqrt(c^2+2*c*g_2+g_1^2+g_2^2)
  # root = quad_solution(1, 2 * c, 2 * c * g_2 - g_1^ 2 - g_2^ 2)
  # root1 <- unlist(root[1])
  # root2<- unlist(root[2])
  # # posroot <- -c+sqrt(c^2-2*c*g_2+g_1^2+g_2^2)
  # # negroot<- -c-sqrt(c^2-2*c*g_2+g_1^2+g_2^2)
  # 
  # # a: norm of beta, b: norm of theta
  # # Hence, we choose the largest value of a and b to take positive value
  # 
  # a = c(
  #   g_1 * root1 / (c + root1),
  #   g_1 * root2 / (c + root2),
  #   g_1 * root1 / (c + root2),
  #   g_1 * root2 / (c + root1)
  # )
  # 
  # bb = c(
  #   g_1*root1 * (c - g_2) / (c + root1),
  #   g_1* root2 * (c - g_2) / (c + root2),
  #   g_1*root1 * (c - g_2) / (c + root2),
  #   g_1* root2 * (c - g_2) / (c + root1)
  # )
  # 
  # 
  # 
  # 
  # x_min = big
  # 
  # j_hat=0; k_hat = 0
  # for (jjj in 1:4){
  #   for (kkk in 1:4){
  #     denominator = (a[jjj]^2 + bb[kkk]^ 2)^ 0.5  # l2 norm
  #     if (isTRUE(denominator > 0)==T){
  #       val1 = (1 + (c / denominator)) * a[jjj] - g_1
  #       val2 = (1 + c * (1 / bb[kkk] + 1 / denominator)) * bb[kkk] - g_2
  #       
  #       temp = abs(val1) + abs(val2)  # l1 norm
  #       if (isTRUE(temp < x_min)==T){
  #         j_hat=jjj; k_hat = kkk
  #         x_min = temp
  #       }
  #     }
  #     
  #     
  #     
  #   }
  #   
  # }  
  # 
  # xnorm = (a[j_hat]^2 + bb[k_hat]^2)^0.5  # l2 norm
  # 
  # 
  # 
  # 
  # #beta_hat<-((beta_j-h*L1))/c_1
  # #new_beta_hat<- beta_hat+ my_b*(beta_hat-old_beta_hat)
  # new_v_beta<-((beta_j-t*L1))/(1 + c / xnorm)
  # # print(L2)
  # # print((theta_j-t*L2))
  # # print(t*alpha*lambda)
  # # print(S_func(as.matrix(theta_j-t*L2,0)))
  # new_v_theta<- (matrix(S_func(theta_j-t*L2,t*alpha*lambda)))/(1 + c * ( (1 / xnorm) + (1 /abs(bb[k_hat])) ) ) 
  # 
  # 
  # beta_u[j]<-new_v_beta;theta_u[j,]<-new_v_theta
  # 
  # theta_transpose_l<-t(theta_u);theta_transpose_r<-t(theta)
  # 
  # 
  # objective_l<- objective_j(beta0,theta0,beta_u,theta_u,X,Z,y,W,alpha,lambda,p,K,N,j,r_min_j)
  # 
  # objective_r<- objective_j(beta0,theta0,beta_new,theta_new,X,Z,y,W,alpha,lambda,p,K,N,j,r_min_j)
  # 
  # 
  # rhs<-objective_r+ matrix(c(L1,L2),nrow = 1)%*%matrix((c(beta_u[j],theta_u[j,])-c(beta_new[j],theta_new[j,])),ncol = 1)+(1/(2*t))*(norm(matrix((c(beta_u[j],theta_u[j,])-c(beta_new[j],theta_new[j,])),ncol = 1),"F"))^2
  # 
  # 
  # 
  # tol=1e-5
  # improvement = rhs - objective_l
  # if (objective_l<=rhs ){
  #   # if(abs(beta_u[j])<abs(beta[j])  ){
  #   #   beta_u[j]=beta[j]
  #   #   
  #   # }
  #   # print(c(b,j,beta_u[j],beta[j],t,nes))
  #   #print(c(b,j,objective_l,rhs,t,nes))
  #   #beta_new[j]<-new_v_beta;theta_new[j,]<-new_v_theta
  #   okay=1
  #   
  # }  else{
  #   # nes=nes+1
  #   # teta_k=2/(nes+1)
  #   # old_beta=v_beta;old_theta=v_theta
  #   # v_beta=beta_u
  #   # v_theta=theta_u
  #   # 
  #   # beta_uu=old_beta+(1/(teta_k))*(v_beta-old_beta)
  #   # theta_uu=old_theta+(1/(teta_k))*(v_theta-old_theta)
  #   # 
  #   # beta_new=(1-teta_k)*old_beta+teta_k*beta_uu
  #   # theta_new=(1-teta_k)*old_theta+teta_k*theta_uu
  #   
  #   
  #   
  #   # 
  #   old_beta=v_beta;old_theta=v_theta
  #   v_beta=beta_u
  #   v_theta=theta_u
  #   beta_new=old_beta+(nes/(nes+3))*(v_beta-old_beta)
  #   theta_new=old_theta+(nes/(nes+3))*(v_theta-old_theta)
  #   
  #   
  #   # print(beta_new)
  #   #print(theta_new)
  #   nes=nes+1
  #   t=.8*t
  # }
  # 
  # 
  # 
  
  
  while (okay<1) {
    #print(beta[j])
    # 
    # old_beta_u=beta_u; old_theta_u=theta_u
    # teta_k<-2/(nes+1)
    # 
    # beta_new<-(1-teta_k)*old_beta_u+teta_k*v_beta
    # theta_new<-(1-teta_k)*old_theta_u+teta_k*v_theta
    # 
    
    theta_transpose_l<-t(theta_u);theta_transpose_r<-t(theta_new)
    
    
    
    
    # xz_theta_r <- lapply(seq_len(p),
    #                      function(jj) (matrix(X[, jj], nrow = length(y), ncol = K) * Z) %*% theta_transpose_r[, jj])
    # XZ_term <- Reduce(f = '+', x = xz_theta_r)
    
    n_r <- unlist(model(beta0, theta0, beta_new, theta_new, X, Z))
    # n_r <- (as.numeric(beta0) + X %*% beta_new + Z %*% theta0 + XZ_term)
    
    
    
    
    n_r_b<-matrix((n_r),length(y))
    y_hat_r<-n_r_b
    
    
    
    
    
    
    
    
    
    r_r<-(y-y_hat_r)
    
    
    grad_j<- gradient_j(beta=beta_new,theta=theta_new,U,U2,U3,y,X,W,r=r_r,alpha,lambda,K,j,N)
    
   
    L1=matrix(unlist(grad_j[1]));L2<-matrix(unlist(grad_j[2]))
    
    U=c(unlist(grad_j[3]));U2=c(unlist(grad_j[4]));U3=c(unlist(grad_j[5]))
    
    # 
    # theta_j<- matrix(as.numeric(theta_new[j,]),nrow = 1)
    # 
    # 
    # dJ_dtheta_j<-matrix(0,nrow = 1,ncol = K)
    # R<-matrix(0,nrow = length(U))
    # 
    # if (isTRUE(any(c((beta_new[j]), theta_j)!=0))==T) {
    #   u<-as.numeric(beta_new[j])/norm(matrix(c(beta_new[j], theta_j)),type= "F")
    #   U<-c(U,u)
    #   
    # } else{
    #   for (z in 1:length(U)) {
    #     if (norm(matrix(U[z]),type = "F")<=1) {
    #       R[z]<-U[z]
    #       next(z)
    #     } else{
    #       next(z)
    #     }
    #     
    #     
    #   }
    #   u<-sample (c(R), size=1)
    #   U<-c(U,u)
    # }
    # 
    # 
    # 
    # dJ_dbeta_j=-(t(matrix(X[,j])) %*% (matrix((y-y_hat_r)) ) )/N+(1-alpha)*lambda*u
    # 
    # 
    # #DJ_BETA[j]<-dJ_dbeta_j
    # 
    # #for (n in 1:K) {
    # R<-matrix(0,nrow = length(U2))
    # if (isTRUE(any(c((beta_new[j]), theta_j)!=0))==T) {
    #   u2<-as.numeric(theta_j)/norm(matrix(c(as.numeric(beta_new[j]), theta_j)),type= "F")
    #   U2<-c(U2,u2)
    #   
    # } else{
    #   for (z in 1:length(U2)) {
    #     if (norm(matrix(U2[z]),type = "F")<=1) {
    #       R[z]<-U2[z]
    #       next(z)
    #     } else{
    #       next(z)
    #     }
    #     
    #     
    #   }
    #   u2<-sample (c(R), size=1)
    #   U2<-c(U2,u2)
    # }
    # 
    # 
    # d<-matrix(0,nrow = length(U3))
    # if (isTRUE(any(c(theta_j)!=0))==T) {
    #   u3<-as.numeric(theta_j)/norm(matrix(c(theta_j)),type= "F")
    #   U3<-c(U3,u3)
    #   
    # } else{
    #   for (z in 1:length(U3)) {
    #     if (norm(matrix(U3[z]),type = "F")<=1) {
    #       d[z]<-U3[z]
    #       next(z)
    #       
    #     } else{
    #       next(z)
    #     }
    #     
    #     
    #   }
    #   u3<-sample (c(R), size=1)
    #   U3<-c(U3,u3)
    # }
    # 
    # 
    # 
    # 
    # v<-sign(as.numeric(theta_j))
    # 
    # 
    # 
    # 
    # dJ_dtheta_j= -((t(data.frame(W[j]))%*% (matrix((y-y_hat_r)) ) ) )/N+(1-alpha)*lambda*(u2+u3)+alpha*lambda*v
    # 
    # # }
    # 
    # 
    # 
    # L1<- dJ_dbeta_j## gradient 
    # L2<- dJ_dtheta_j
    # 
    # #print(L2)
    
    
    
    
    
    
    beta_j=as.numeric(beta_new[j]);theta_j=as.numeric(theta_new[j,])
    
    a<- as.numeric(norm(matrix(beta_j), type="F")) ; bb<-as.numeric(twonorm(matrix(theta_j))) ; rho_2<- as.numeric(sqrt(a^2+bb^2))
    
    
    ##### after norming 
    c<- as.numeric(t*(1-alpha)*lambda );  g_1<- as.numeric( abs(beta_j-t*as.numeric(L1)))
    
    
    # for (b in 1:K) {
    v<-S_func(theta_j-t*L2,t*alpha*lambda)
    # }
    
    g_2<-as.numeric(twonorm(matrix(v)))
    
    #posroot<- -c+sqrt(c^2+2*c*g_2+g_1^2+g_2^2)
    
    #negroot<- -c-sqrt(c^2+2*c*g_2+g_1^2+g_2^2)
    root = quad_solution(1, 2 * c, 2 * c * g_2 - g_1^ 2 - g_2^ 2)
    root1 <- unlist(root[1])
    root2<- unlist(root[2])
    # posroot <- -c+sqrt(c^2-2*c*g_2+g_1^2+g_2^2)
    # negroot<- -c-sqrt(c^2-2*c*g_2+g_1^2+g_2^2)
    
    # a: norm of beta, b: norm of theta
    # Hence, we choose the largest value of a and b to take positive value
    
    a = c(
      g_1 * root1 / (c + root1),
      g_1 * root2 / (c + root2),
      g_1 * root1 / (c + root2),
      g_1 * root2 / (c + root1)
    )
    
    bb = c(
      g_1*root1 * (c - g_2) / (c + root1),
      g_1*root2 * (c - g_2) / (c + root2),
      g_1*root1 * (c - g_2) / (c + root2),
      g_1*root2 * (c - g_2) / (c + root1)
    )
    
    
    
    
    x_min = big
    
    j_hat=0; k_hat = 0
    for (jjj in 1:4){
      for (kkk in 1:4){
        denominator = (a[jjj]^2 + bb[kkk]^ 2)^ 0.5  # l2 norm
        if (isTRUE(denominator > 0)==T){
          val1 = (1 + (c / denominator)) * a[jjj] - g_1
          val2 = (1 + c * (1 / bb[kkk] + 1 / denominator)) * bb[kkk] - g_2
          
          temp = abs(val1) + abs(val2)  # l1 norm
          if (isTRUE(temp < x_min)==T){
            j_hat=jjj; k_hat = kkk
            x_min = temp
          }
        }
        
        
        
      }
      
    }  
    
    xnorm = (a[j_hat]^2 + bb[k_hat]^2)^0.5  # l2 norm
    
    
    
    
    #beta_hat<-((beta_j-h*L1))/c_1
    #new_beta_hat<- beta_hat+ my_b*(beta_hat-old_beta_hat)
    new_v_beta<-((beta_j-t*L1))/(1 + c / xnorm)
    #print(new_v_beta)
    
    new_v_theta<- (matrix(S_func(theta_j-t*L2,t*alpha*lambda)))/(1 + c * ( (1 / xnorm) + (1 /abs(bb[k_hat])) ) ) 
    #		new_v_theta[n]<- norm(matrix(S_func(theta_j[n]-t*L2[n],t*alpha*lambda)),type="F")/c_2
    #theta_hat[n]<-new_v_theta[n]+(my_n/(my_n+3))*(new_v_theta[n]-theta_j[n])
    # }
    
    beta_u[j]<-new_v_beta;theta_u[j,]<-new_v_theta
    
    theta_transpose_l<-t(theta_u);theta_transpose_r<-t(theta_new)
    
    
    objective_l<- objective_j(beta0,theta0,beta_u,theta_u,X,Z,y,W,alpha,lambda,p,K,N,j,r_min_j)
    
    objective_r<- objective_j(beta0,theta0,beta_new,theta_new,X,Z,y,W,alpha,lambda,p,K,N,j,r_min_j)
    
    
    rhs<-objective_r+ matrix(c(L1,L2),nrow = 1)%*%matrix((c(beta_u[j],theta_u[j,])-c(beta_new[j],theta_new[j,])),ncol = 1)+(1/(2*t))*(norm(matrix((c(beta_u[j],theta_u[j,])-c(beta_new[j],theta_new[j,])),ncol = 1),"F"))^2
    
    
    
    
    
    
    # 
    # if(isTRUE(objective_l<=rhs)==T & abs(beta_new[j])>= beta[j]){
    #   
    #   beta_new[j]<-beta_u[j]
    #   print(c(b,j,beta_new[j],beta[j],t,nes))
    #   #print(c(b,j,objective_l,rhs,t,nes))
    #   #beta_new[j]<-new_v_beta;theta_new[j,]<-new_v_theta
    #   okay=1
    # }else 
    tol=1e-5
    improvement = rhs - objective_l
    if (objective_l<=rhs  ){
      # if(abs(beta_u[j])<abs(beta[j])  ){
      #   beta_u[j]=beta[j]
      # }
      #print(c(b,j,beta_u[j],beta[j],t,nes))
      #print(c(b,j,objective_l,rhs,t,nes))
      #beta_new[j]<-new_v_beta;theta_new[j,]<-new_v_theta
      okay=1
      
    }  else{
      
      # 
      # nes=nes+1
      # teta_k=2/(nes+1)
      # old_beta=v_beta;old_theta=v_theta
      # v_beta=beta_u
      # v_theta=theta_u
      # 
      # beta_uu=old_beta+(1/(teta_k))*(v_beta-old_beta)
      # theta_uu=old_theta+(1/(teta_k))*(v_theta-old_theta)
      # 
      # beta_new=(1-teta_k)*old_beta+teta_k*beta_uu
      # theta_new=(1-teta_k)*old_theta+teta_k*theta_uu
      
      
      
      old_beta=v_beta;old_theta=v_theta
      v_beta=beta_u
      v_theta=theta_u
      beta_new=old_beta+(nes/(nes+3))*(v_beta-old_beta)
      theta_new=old_theta+(nes/(nes+3))*(v_theta-old_theta)
      
      
      # print(beta_new)
      #print(theta_new)
      nes=nes+1
      t=.8*t
    }
    
    if(isTRUE(nes>bb_it)==T){
      # beta_u[j]=beta[j]
      #print(c(b,j,beta_u[j],beta[j],t,nes))
      okay=1
    }
    
  }#while
  
  beta1<-beta_u[j];theta1<-theta_u[j,];t=t
   #print(c(nes,j,beta[j],beta1))
  #while
  
  # print(c(b,j))
  #print(c(b,j))
  #print(beta1)
  #print(theta1)
  return(list(beta1,theta1,t,U,U2,U3))
}



