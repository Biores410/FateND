
#' @name  get_fit_HSC
#' @title fit the three lineages of time series HSC data to the trend like switch genes and transient genes
#' @description  function to fit time series data by nls
#' @param lineage1 refer to time series data of lineage1 in HSC
#' @param lineage2 refer to time series data of lineage2 in HSC
#' @param lineage3 refer to time series data of lineage3 in HSC
#' @return  g1 & g2 & fitdata of lineage1 lineage2  and lineage3
#' @export
get_fit_HSC <-function(lineage1,lineage2,lineage3){
  
  #################################
  ######logis fit function-fitting expression of switch genes
  #################################
  logis_Fit=function(gene_exp){     
    
    d=data.frame(cbind(1:length(gene_exp),as.numeric(gene_exp)))
    colnames(d)=c("x","y")
    nlc <- nls.control(maxiter = 1000) 
    
    fm1DNase1 <- try(nls(y~ SSlogis(x, Asym, xmid, scal),control = nlc, d),silent = T)
    ###Only genes with switch trends will be fitted successfully, so the results need to be screened
    if(!'try-error' %in% class(fm1DNase1))            
    {
      pred=predict(fm1DNase1)
      goodness_of_fit=cor(as.numeric(gene_exp),pred)  
      coefficient=coef(fm1DNase1)
      
    }else{
      pred=rep(NA,length(as.numeric(gene_exp)))
      goodness_of_fit=NA
      coefficient=rep(NA,3)
    }
    #Returns the results of the fit, which are the coefficient: Asym,xmid,scal; goodness of fit;fit value.
    return(c(coefficient,goodness_of_fit,pred))
  }
  
  
  
  #################################
  ######norm fit function - fitting expression of transient genes
  #################################
  
  norm_nls_Fit1=function(gene_exp){ 
    
    gene_exp=as.numeric(gene_exp)
    d=data.frame(cbind(1:length(gene_exp),gene_exp))
    colnames(d)=c("x","y")
    nlc <- nls.control(maxiter = 10000)
    
    sigma_model<- try(nls(y ~ k*exp(-σ*((x-μ)^2)),control=nlc,
                          d,start = list(k=4,μ=1000,σ=0.00001)),silent = T)
    ###Only genes with transient trends will be fitted successfully, so the results need to be screened
    if(!'try-error' %in% class(sigma_model))            
    {
      
      k=environment(environment(sigma_model[["m"]][["fitted"]])[["getPars"]])[["env"]][["k"]]
      μ=environment(environment(sigma_model[["m"]][["fitted"]])[["getPars"]])[["env"]][["μ"]]
      σ=environment(environment(sigma_model[["m"]][["fitted"]])[["getPars"]])[["env"]][["σ"]]
      
      x=1:length(gene_exp)
      y=k*exp(-σ*((x-μ)^2))
      goodness_of_fit=cor(gene_exp,y)  
      coefficient=coef(sigma_model)
      
    }else{
      
      y=rep(NA,length(gene_exp))
      goodness_of_fit=NA
      coefficient=rep(NA,3)
    }
    
    #Returns the results of the fit, which are the coefficient: k,μ,σ; goodness of fit;fit value.
    return(c(coefficient,goodness_of_fit,y))
  }
  
  
  norm_nls_Fit2=function(gene_exp){ 
    
    gene_exp=as.numeric(gene_exp)
    d=data.frame(cbind(1:length(gene_exp),gene_exp))
    colnames(d)=c("x","y")
    nlc <- nls.control(maxiter = 10000)
    
    sigma_model<- try(nls(y ~ k*exp(-σ*((x-μ)^2)),control=nlc,
                          d,start = list(k=4,μ=1200,σ=0.00001)),silent = T)
    
    if(!'try-error' %in% class(sigma_model))            # 判断当前循环的try语句中的表达式是否运行正确
    {
      
      k=environment(environment(sigma_model[["m"]][["fitted"]])[["getPars"]])[["env"]][["k"]]
      μ=environment(environment(sigma_model[["m"]][["fitted"]])[["getPars"]])[["env"]][["μ"]]
      σ=environment(environment(sigma_model[["m"]][["fitted"]])[["getPars"]])[["env"]][["σ"]]
      
      x=1:length(gene_exp)
      y=k*exp(-σ*((x-μ)^2))
      goodness_of_fit=cor(gene_exp,y)  
      coefficient=coef(sigma_model)
      
    }else{
      
      y=rep(NA,length(gene_exp))
      goodness_of_fit=NA
      coefficient=rep(NA,3)
    }
    
    
    return(c(coefficient,goodness_of_fit,y))
  }
  
  
  
  ###lineage1 
  increade_lineage1_F=t(apply(lineage1,1,logis_Fit))  ##卡scale(+)    ##c("Asym","xmid","scal","goodness")
  decrease_lineage1_F=t(apply(lineage1,1,norm_nls_Fit1))  ###c("k","μ","σ","goodness") (k=4,μ=1000,σ=0.00001)
  
  increade_lineage1_S=increade_lineage1_F[intersect(which(increade_lineage1_F[,3]>0),which(increade_lineage1_F[,1]>0)),]
  increade_lineage1_S=increade_lineage1_S[order(increade_lineage1_S[,4],decreasing=T),]
  
  t1=which(decrease_lineage1_F[,1]>0)
  t2=which((decrease_lineage1_F[,2]<=(0.95*dim(lineage1)[2]))&(decrease_lineage1_F[,2]>=(0.05*dim(lineage1)[2])))
  t3=which((decrease_lineage1_F[,3]<=0.001)&(decrease_lineage1_F[,3]>=0.000001))
  decreade_lineage1_S=decrease_lineage1_F[intersect(intersect(t1,t2),t3),]
  decreade_lineage1_S=decreade_lineage1_S[order(decreade_lineage1_S[,4],decreasing=T),]
  
  
  ###lineage2
  increade_lineage2_F=t(apply(lineage2,1,logis_Fit))      ##卡scale(+)    ##c("Asym","xmid","scal","goodness") 
  decrease_lineage2_F=t(apply(lineage2,1,norm_nls_Fit2))     ##c("k","μ","σ","goodness") k=4,μ=1200,σ=0.00001)
  
  increade_lineage2_S=increade_lineage2_F[intersect(which(increade_lineage2_F[,3]>0),which(increade_lineage2_F[,1]>0)),]
  increade_lineage2_S=increade_lineage2_S[order(increade_lineage2_S[,4],decreasing=T),]
  
  t1=which(decrease_lineage2_F[,1]>0)
  t2=which((decrease_lineage2_F[,2]<=(0.95*dim(lineage2)[2]))&(decrease_lineage2_F[,2]>=(0.05*dim(lineage2)[2])))
  t3=which((decrease_lineage2_F[,3]<=0.001)&(decrease_lineage2_F[,3]>=0.000001))
  decreade_lineage2_S=decrease_lineage2_F[intersect(intersect(t1,t2),t3),]
  decreade_lineage2_S=decreade_lineage2_S[order(decreade_lineage2_S[,4],decreasing=T),]
  
  
  ###lineage3
  increade_lineage3_F=t(apply(lineage3,1,logis_Fit))     ##卡scale(+)    ##c("Asym","xmid","scal","goodness") 
  decrease_lineage3_F=t(apply(lineage3,1,norm_nls_Fit2))   ##c("k","μ","σ","goodness")  k=4,μ=1200,σ=0.00001
  
  increade_lineage3_S=increade_lineage3_F[intersect(which(increade_lineage3_F[,3]>0),which(increade_lineage3_F[,1]>0.05)),]
  increade_lineage3_S=increade_lineage3_S[order(increade_lineage3_S[,4],decreasing=T),]
  
  t1=which(decrease_lineage3_F[,1]>0)
  t2=which((decrease_lineage3_F[,2]<=(0.95*dim(lineage3)[2]))&(decrease_lineage3_F[,2]>=(0.05*dim(lineage3)[2])))
  t3=which((decrease_lineage3_F[,3]<=0.001)&(decrease_lineage3_F[,3]>=0.000001))
  decreade_lineage3_S=decrease_lineage3_F[intersect(intersect(t1,t2),t3),]
  decreade_lineage3_S=decreade_lineage3_S[order(decreade_lineage3_S[,4],decreasing=T),]
  
  
  A=rownames(increade_lineage1_S)
  
  B=rownames(increade_lineage2_S)
  
  C=rownames(increade_lineage3_S)
  
  m=intersect(intersect(A,B),C)
  AA=c()
  BB=c()
  CC=c()
  for( i in 1:length(m)){
    if(max(increade_lineage1_S[m[i],4],increade_lineage2_S[m[i],4],increade_lineage3_S[m[i],4])==increade_lineage1_S[m[i],4])
      AA=c(AA,m[i])else if(max(increade_lineage1_S[m[i],4],increade_lineage2_S[m[i],4],increade_lineage3_S[m[i],4])==increade_lineage2_S[m[i],4])
        BB=c(BB,m[i])else if(max(increade_lineage1_S[m[i],4],increade_lineage2_S[m[i],4],increade_lineage3_S[m[i],4])==increade_lineage3_S[m[i],4])
          CC=c(CC,m[i])
  }
  
  A1=c(setdiff(A,m),AA)
  B1=c(setdiff(B,m),BB)
  C1=c(setdiff(C,m),CC)
  
  
  m=intersect(A1,B1)
  AA=c()
  BB=c()
  for( i in 1:length(m)){
    if(max(increade_lineage1_S[m[i],4],increade_lineage2_S[m[i],4])==increade_lineage1_S[m[i],4])
      AA=c(AA,m[i])else if(max(increade_lineage1_S[m[i],4],increade_lineage2_S[m[i],4])==increade_lineage2_S[m[i],4])
        BB=c(BB,m[i])
  }
  
  A2=c(setdiff(A1,m),AA)
  B2=c(setdiff(B1,m),BB)
  C1=c(setdiff(C,m),CC)
  
  m=intersect(A2,C1)
  AA=c()
  CC=c()
  for( i in 1:length(m)){
    if(max(increade_lineage1_S[m[i],4],increade_lineage3_S[m[i],4])==increade_lineage1_S[m[i],4])
      AA=c(AA,m[i])else if(max(increade_lineage1_S[m[i],4],increade_lineage3_S[m[i],4])==increade_lineage3_S[m[i],4])
        CC=c(CC,m[i])
  }
  
  A3=c(setdiff(A2,m),AA)
  B2=c(setdiff(B1,m),BB)
  C2=c(setdiff(C1,m),CC)
  
  m=intersect(B2,C2)
  BB=c()
  CC=c()
  for( i in 1:length(m)){
    if(max(increade_lineage2_S[m[i],4],increade_lineage3_S[m[i],4])==increade_lineage2_S[m[i],4])
      BB=c(BB,m[i])else if(max(increade_lineage2_S[m[i],4],increade_lineage3_S[m[i],4])==increade_lineage3_S[m[i],4])
        CC=c(CC,m[i])
  }
  A3=c(setdiff(A2,m),AA)
  B3=c(setdiff(B2,m),BB)
  C3=c(setdiff(C2,m),CC)
  
  lineage1_fit=rbind(increade_lineage1_S[A3,-c(1:4)],decreade_lineage1_S[setdiff(rownames(decreade_lineage1_S),A3),-c(1:4)])
  
  lineage2_fit=rbind(increade_lineage2_S[B3,-c(1:4)],decreade_lineage2_S[setdiff(rownames(decreade_lineage2_S),B3),-c(1:4)])
  
  lineage3_fit=rbind(increade_lineage3_S[C3,-c(1:4)],decreade_lineage3_S[setdiff(rownames(decreade_lineage3_S),C3),-c(1:4)])
  
  a=intersect(rownames(lineage1_fit),rownames(lineage2_fit))
  U=intersect(a,rownames(lineage3_fit))
  
  new_A=intersect(U,A3)
  new_A <- new_A[!grepl("ZNF", new_A)]
  
  new_B=intersect(U,B3)
  new_B <- new_B[!grepl("ZNF", new_B)]
  
  new_C=intersect(U,C3)
  new_C <- new_C[!grepl("ZNF", new_C)]
  
  results=list(new_A, union(new_B,new_C), lineage1_fit,lineage2_fit,lineage3_fit)
  return(results)
  
}
