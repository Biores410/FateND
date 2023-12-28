
#' @name  get_fit_embryos
#' @title fit the two lineages of time series embryos data to the trend like switch genes and transient genes
#' @description  function to fit time series data by nls
#' @param lineage1 refer to time series data of lineage1 in embryos
#' @param lineage2 refer to time series data of lineage2 in embryos
#' @return  g1 & g2 & fitdata of lineage1 and lineage2 
#' @export
get_fit_embryos <-function(lineage1,lineage2){
  
  #################################
  ######logis fit function-fitting expression of switch genes
  #################################
  logis_Fit=function(gene_exp){     
    
    d=data.frame(cbind(1:length(gene_exp),as.numeric(gene_exp)))
    colnames(d)=c("x","y")
    nlc <- nls.control(maxiter = 1000,minFactor = 1/2048) 
    
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
    
    
    d=data.frame(cbind(1:length(gene_exp),as.numeric(gene_exp)))
    colnames(d)=c("x","y")
    
    nlc <- nls.control(maxiter = 1000)
    
    sigma_model<- try(nls(y ~ k*(1/(σ*sqrt(2*pi)))*exp((-1/2)*(x-μ)^2/σ^2),control=nlc,
                          d,start = list(k=100,μ=150,σ=sqrt(150))),silent = T)
    ###Only genes with transient trends will be fitted successfully, so the results need to be screened
    if(!'try-error' %in% class(sigma_model))           
    {
      
      k=environment(environment(sigma_model[["m"]][["fitted"]])[["getPars"]])[["env"]][["k"]]
      μ=environment(environment(sigma_model[["m"]][["fitted"]])[["getPars"]])[["env"]][["μ"]]
      σ=environment(environment(sigma_model[["m"]][["fitted"]])[["getPars"]])[["env"]][["σ"]]
      
      x=1:length(gene_exp)
      y=k*(1/(σ*sqrt(2*pi)))*exp((-1/2)*(x-μ)^2/σ^2)
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
  
  ###lineage1
  increade_lineage1_F=t(apply(lineage1,1,logis_Fit))   ##c("Asym","xmid","scal","goodness")
  decrease_lineage1_F=t(apply(lineage1,1,norm_nls_Fit))  ##c("k","μ","σ","goodness")
  
  #for increade_lineage1_F, the first column is Asym;the second column is xmid;the third column is scal,the fourth column is goodness.
  increade_lineage1_S=increade_lineage1_F[intersect(which((increade_lineage1_F[,3]>0)&(increade_lineage1_F[,3]<(dim(lineage1)[2]/7.6))),which((increade_lineage1_F[,2]>0.2*(dim(lineage1)[2]))&(increade_lineage1_F[,2]<0.8*(dim(lineage1)[2])))),]
  increade_lineage1_S=increade_lineage1_S[order(increade_lineage1_S[,4],decreasing=T),]
  #for decrease_lineage1_F,the first column is k;the second column is μ;the third column is σ,the fourth column is goodness.
  decreade_lineage1_S=decrease_lineage1_F[intersect(which((decrease_lineage1_F[,3]<(0.2*dim(lineage1)[2]))&(decrease_lineage1_F[,3]>(0.01*dim(lineage1)[2]))),which((decrease_lineage1_F[,2]<(0.85*dim(lineage1)[2]))&(decrease_lineage1_F[,2]>(0.15*dim(lineage1)[2])))),]
  decreade_lineage1_S=decreade_lineage1_S[order(decreade_lineage1_S[,4],decreasing=T),]
  
  
  
  ###lineage2
  increade_lineage2_F=t(apply(lineage2,1,logis_Fit))  ## c("Asym","xmid","scal","goodness")
  decrease_lineage2_F=t(apply(lineage2,1,norm_nls_Fit))  ## c("k","μ","σ","goodness")
  
  #for increade_lineage2_F, the first column is Asym;the second column is xmid;the third column is scal,the fourth column is goodness.
  increade_lineage2_S=increade_lineage2_F[intersect(which((increade_lineage2_F[,3]>0)&(increade_lineage2_F[,3]<(dim(lineage2)[2]/7.6))),which((increade_lineage2_F[,2]>0.2*(dim(lineage2)[2]))&(increade_lineage2_F[,2]<0.8*(dim(lineage2)[2])))),]
  increade_lineage2_S=increade_lineage2_S[order(increade_lineage2_S[,4],decreasing=T),]
  #for decrease_lineage2_F,the first column is k;the second column is μ;the third column is σ,the fourth column is goodness.
  decreade_lineage2_S=decrease_lineage2_F[intersect(which((decrease_lineage2_F[,3]<(0.2*dim(lineage2)[2]))&(decrease_lineage2_F[,3]>(0.01*dim(lineage2)[2]))),which((decrease_lineage2_F[,2]<(0.85*dim(lineage2)[2]))&(decrease_lineage2_F[,2]>(0.15*dim(lineage2)[2])))),]
  decreade_lineage2_S=decreade_lineage2_S[order(decreade_lineage2_S[,4],decreasing=T),]
  
  
  A=intersect(rownames(increade_lineage1_S),rownames(decreade_lineage2_S))
  B=intersect(rownames(increade_lineage2_S),rownames(decreade_lineage1_S))
  
  ####If A and B are duplicated, they are reclassified by goodness of fit.
  AA=c()
  BB=c()
  for (i in 1:length(intersect(A,B))){
    
    if((increade_lineage1_S[intersect(A,B)[i],4]>increade_lineage2_S[intersect(A,B)[i],4])){
      
      AA=c(AA,intersect(A,B)[i])
    }else {
      
      BB=c(BB,intersect(A,B)[i])
    }
    
    
  }
  
  m=intersect(A,B)
  
  A=c(setdiff(A,m),AA)
  B=c(setdiff(B,m),BB)
  
  
  lineage1_fit=rbind(increade_lineage1_S[A,-c(1:4)],decreade_lineage1_S[B,-c(1:4)])
  lineage2_fit=rbind(increade_lineage2_S[B,-c(1:4)],decreade_lineage2_S[A,-c(1:4)])
  
  ####A is g1,B is g2
  results=list(A, B, lineage1_fit,lineage2_fit)
  return(results)
  
}

