
#' @name  get_fit_committed
#' @title fit the time series committed data to the trend like switch genes and transient genes
#' @description  function to fit time series data by nls
#' @param lineage1 refer to time series of committed data
#' @return  g1 & g2 & fitdata 
#' @export
get_fit_committed <-function(lineage1){
  
  #################################
  ######logis fit function- fitting expression of switch genes  
  #################################
  logis_Fit=function(gene_exp){    
    
    gene_exp=as.numeric(gene_exp)
    d=data.frame(cbind(1:length(gene_exp),gene_exp))
    colnames(d)=c("x","y")
    
    nlc <- nls.control(maxiter = 1000,minFactor = 1/1024) 
    fm1DNase1 <- try(nls(y~ SSlogis(x, Asym, xmid, scal),control = nlc, d),silent = T)
    ###Only genes with switch trends will be fitted successfully, so the results need to be screened
    if(!'try-error' %in% class(fm1DNase1))            
    {
      pred=predict(fm1DNase1)
      goodness_of_fit=cor(gene_exp,pred)  
      coefficient=coef(fm1DNase1)
      
    }else{
      pred=rep(NA,length(gene_exp))
      goodness_of_fit=NA
      coefficient=rep(NA,3)
    }
    #Returns the results of the fit, which are the coefficient: Asym,xmid,scal; goodness of fit;fit value.
    return(c(coefficient,goodness_of_fit,pred))
  }
  
  #################################
  ######norm fit function - fitting expression of transient genes
  #################################
  norm_nls_Fit=function(gene_exp){     
    
    gene_exp=as.numeric(gene_exp)
    d=data.frame(cbind(1:length(gene_exp),gene_exp))
    colnames(d)=c("x","y")
    nlc <- nls.control(maxiter = 10000)
    sigma_model<- try(nls(y ~ k*exp(-σ*((x-μ)^2)),control=nlc,
                          d,start = list(k=10,μ=400,σ=0.0001)),silent = T)
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
  
  ###lineage1 fit 
  increade_lineage1_F=t(apply(lineage1,1,logis_Fit))    #c("Asym","xmid","scal","goodness")
  decrease_lineage1_F=t(apply(lineage1,1,norm_nls_Fit))  #c("k","μ","σ","goodness")
  
  #for increade_lineage1_F, the first column is Asym;the second column is xmid;the third column is scal,the fourth column is goodness.
  increade_lineage1_S=increade_lineage1_S=increade_lineage1_F[intersect(intersect(which((increade_lineage1_F[,3]>5)&(increade_lineage1_F[,3]<90)),which((increade_lineage1_F[,2]>0.15*(dim(lineage1)[2]))&(increade_lineage1_F[,2]<0.85*(dim(lineage1)[2])))),which(increade_lineage1_F[,1]>2)),]
  increade_lineage1_S=increade_lineage1_S[order(increade_lineage1_S[,4],decreasing=T),]
  increade_lineage1_S=increade_lineage1_S[c(1:26),]
  
  #for decrease_lineage1_F,the first column is k;the second column is μ;the third column is σ,the fourth column is goodness.
  t1=which(decrease_lineage1_F[,1]>2)
  t2=which((decrease_lineage1_F[,2]<(0.8*dim(lineage1)[2]))&(decrease_lineage1_F[,2]>(0.2*dim(lineage1)[2])))
  t3=which((decrease_lineage1_F[,3]<0.001)&(decrease_lineage1_F[,3]>0.00001))
  decreade_lineage1_S=decrease_lineage1_F[intersect(intersect(t1,t2),t3),]
  decreade_lineage1_S=decreade_lineage1_S[order(decreade_lineage1_S[,4],decreasing=T),]
  decreade_lineage1_S=decreade_lineage1_S[c(1:27),]
  
  
  A=rownames(increade_lineage1_S)
  B=rownames(decreade_lineage1_S)
  
  ####If A and B are duplicated, they are reclassified by goodness of fit.
  AA=c()
  BB=c()
  for (i in 1:length(intersect(A,B))){
    
    if((increade_lineage1_S[intersect(A,B)[i],4]>decreade_lineage1_S[intersect(A,B)[i],4])){
      
      AA=c(AA,intersect(A,B)[i])
    }else {
      
      BB=c(BB,intersect(A,B)[i])
    }
    
    
  }
  
  m=intersect(A,B)
  
  A=c(setdiff(A,m),AA)
  B=c(setdiff(B,m),BB)
  
  lineage1_fit=rbind(increade_lineage1_S[A,-c(1:4)],decreade_lineage1_S[B,-c(1:4)])
  ####A is g1,B is g2
  results=list(A, B, lineage1_fit)
  return(results)
}
