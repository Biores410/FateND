
#' @name FateND_HSC
#' @title  by inputing the fit data of three lineages and genepairs and all genes of HSC,we can get the cell fate determinates
#' @description  function to get cell fate determinate TFs
#' @param lineage1_fit refer to the fitting  time series data of lineage1 in HSC
#' @param lineage2_fit refer to the fitting  time series data of lineage2 in HSC
#' @param lineage3_fit refer to the fitting  time series data of lineage3 in HSC
#' @param overall_genepair refer to all gene pairs
#' @param U refer to all genes
#' @return  cell fate determinate TFs
#' @export
FateND_HSC <- function(lineage1_fit,lineage2_fit,lineage3_fit,overall_genepair,U){
  
  #####prepare diff_term,Autocat_term,Inhibit_term,degrad_term for SDE
  F_Pre_process_data=function(branch_1,branch_2,branch_3,gene_pair){
    
    diff_branch_1=branch_1[gene_pair,-1]-branch_1[gene_pair,-ncol(branch_1)]
    diff_branch_2=branch_2[gene_pair,-1]-branch_2[gene_pair,-ncol(branch_2)]
    diff_branch_3=branch_3[gene_pair,-1]-branch_3[gene_pair,-ncol(branch_3)]
    diff_term=cbind(diff_branch_1,diff_branch_2,diff_branch_3)     
    
    degrad_term=cbind(branch_1[gene_pair,-1],branch_2[gene_pair,-1],branch_3[gene_pair,-1])  
    
    #S
    mean_each_branch_pair=rbind(apply(branch_1[gene_pair,],1,mean),apply(branch_2[gene_pair,],1,mean),apply(branch_3[gene_pair,],1,mean))
    mean_all_branch=apply(mean_each_branch_pair,2,mean)
    S=sqrt(mean_all_branch[1]*mean_all_branch[2])         
    
    Autocat_term=degrad_term^4/(S^4+degrad_term^4)       
    Inhibit_term=S^4/(S^4+degrad_term^4)                
    
    output=list(diff_term,Autocat_term,Inhibit_term,degrad_term)
    
    names(output)=c("diff_term","Autocat_term","Inhibit_term","degrad_term")
    return(output)
  }
  
  ####lm regression and compute F-value 
  F_Gene_interaction_prediction=function(Pre_gene_pair_data){
    
    
    data1=data.frame(cbind(Pre_gene_pair_data[["diff_term"]][1,],Pre_gene_pair_data[["Autocat_term"]][1,],Pre_gene_pair_data[["Inhibit_term"]][2,],Pre_gene_pair_data[["degrad_term"]][1,]))
    colnames(data1)=c("P","M","k","x")
    data2=data.frame(cbind(Pre_gene_pair_data[["diff_term"]][2,],Pre_gene_pair_data[["Autocat_term"]][2,],Pre_gene_pair_data[["Inhibit_term"]][1,],Pre_gene_pair_data[["degrad_term"]][2,]))
    colnames(data2)=c("Q","N","h","y")
    
    N=dim(data1)[1]
    
    ####The residual sum of squares of the first equation in H1
    lm11=lm(P~0+M+k+x,data1)
    b1=lm11$coefficients["k"]
    SSE11=sum((lm11$residuals)^2)
    
    ####The residual sum of squares of the second equation in H1
    lm12=lm(Q~0+N+h+y,data2)
    b2=lm12$coefficients["h"]
    SSE12=sum((lm12$residuals)^2)
    
    ####The residual sum of squares of the first equation in H0
    lm01=lm(P~0+M+x,data1)
    SSE01=sum((lm01$residuals)^2)
    
    ####The residual sum of squares of the second equation in H0
    lm02=lm(Q~0+N+y,data2)
    SSE02=sum((lm02$residuals)^2)
    
    Fvalue1=(SSE01-SSE11)*(N-1)/SSE11
    Fvalue2=(SSE02-SSE12)*(N-1)/SSE12
    # in order to balance the effection of Fvalue1 and Fvalue2
    Fvalue=(Fvalue1*Fvalue2)/(Fvalue1+Fvalue2)
    
    ###return Estimate_inhibit and F_value
    results_lm=rbind(cbind(b1,Fvalue),cbind(b2,Fvalue))
    colnames(results_lm)=c("Estimate_inhibit","F-value")
    results_lm=as.matrix(results_lm)
    return(results_lm)
  }
  
  
  ####Cycle to count all gene pairs
  F_Overall_gene_pairs_prediction=function(branch_1,branch_2,branch_3){
    
    Interaction_prediction=list()
    for(i in 1:nrow(overall_genepair)){
      
      gene_pair=overall_genepair[i,]
      pre_gene_pair_data=F_Pre_process_data(branch_1,branch_2,branch_3,gene_pair)
      interaction_prediction=F_Gene_interaction_prediction(pre_gene_pair_data)
      rownames(interaction_prediction)=gene_pair
      Interaction_prediction[[i]]=interaction_prediction
      
    }
    names(Interaction_prediction)=apply(overall_genepair,1,function(x) paste(x,collapse  = "_"))
    return(Interaction_prediction)
  }
  
  
  ####sifting the gene pairs with positive b1(Estimate_inhibit) 
  F_significance_pair=function(Overall_gene_pairs_prediction){
    
    significant_inhibit_P_value=list()
    
    for(i in 1:length(Overall_gene_pairs_prediction))   {
      interaction_pred=Overall_gene_pairs_prediction[[i]]
      if (interaction_pred[1,"Estimate_inhibit"]>0& interaction_pred[2,"Estimate_inhibit"]>0){
        significant_inhibit_P_value[[i]]=interaction_pred
      } else {significant_inhibit_P_value[[i]]=NULL}
      
    }
    
    significant_inhibit_P_value_out=significant_inhibit_P_value[!unlist(lapply(significant_inhibit_P_value,is.null))]
    
    
    return(significant_inhibit_P_value_out)
    
  } 
  
  ###Rank according to the F-value
  F.ordered_significant_pair=function(significance_pair){
    
    Multip_F_value=c()
    for(i in 1:length(significance_pair)){
      significance_pair_i=significance_pair[[i]]
      multip_F_value=significance_pair_i[1,"F-value"]
      Multip_F_value=c(Multip_F_value,multip_F_value)
    }
    
    return(significance_pair[order(Multip_F_value,decreasing = T)])
  }
  
  #make results as a matrix
  ourDE_pair_output=function(NA_results){
    
    results=c()
    for(i in 1:length(NA_results)){
      aa=mean(NA_results[[i]][,2])
      aa=as.matrix(aa)
      gene=c(rownames(NA_results[[i]])[1],rownames(NA_results[[i]])[2])
      rownames(aa)=paste(gene[1],gene[2], sep = "_")
      results=rbind(results,aa)
      
    }
    
    results=as.matrix(results)
    colnames(results)=c("F_value")
    return(results)
  }
  
  
  Overall_gene_pairs_prediction=F_Overall_gene_pairs_prediction(lineage1_fit[U,],lineage2_fit[U,],lineage3_fit[U,])
  significance_pair=F_significance_pair(Overall_gene_pairs_prediction)
  results=F.ordered_significant_pair(significance_pair)
  our_pair_Fvalue_results=ourDE_pair_output(results)
  
  
  #show the results of top 20
  par(mfrow=c(1,1))
  our_pair_top=our_pair_Fvalue_results[c(1:20),]
  our_pair_top=as.matrix(our_pair_top)
  our_pair_top=cbind(c(1:20),our_pair_top)
  our_pair_top_decrease=our_pair_top[order(our_pair_top[,1],decreasing=TRUE),]
  plot(our_pair_top_decrease[,2],our_pair_top_decrease[,1],xlim=c(600,0),ylim=c(20,1),cex = 1,col="red",xlab=c("F-Value"),ylab=c("rank"),pch=2,)
  grid(lty = "dotted")
  text(our_pair_top_decrease[,2]+80, our_pair_top_decrease[,1], labels = rownames(our_pair_top_decrease),srt=0,cex=0.8)
  
  ##plot three TFs expression (high F-value and low F-value)
  par(mfrow=c(3,3),mar=c(2,2,1,1),oma=c(2,2,2,2))
  red<- brewer.pal(9,"Reds") 
  blue<-brewer.pal(9,"Blues") 
  gene_pair=c(c("LEF1"),c("SATB2"))
  plot(lineage1_fit[gene_pair[1],],col=c("white"),ylim=c(0,0.2))
  points(as.numeric(Ery_data[gene_pair[1],])/7,col=red[5],cex=0.5,pch=16)
  points(as.numeric(Ery_data[gene_pair[2],])/7,col=blue[5],cex=0.5,pch=16)
  lines(lineage1_fit[gene_pair[1],],col="darkred",lwd=3)
  lines(lineage1_fit[gene_pair[2],],col="darkblue",lwd=3)
  legend("topleft", legend = c(gene_pair[1],gene_pair[2],c("lineage_Ery"),c("lineage_DC"),c("lineage_Mono")), col=c("darkred","darkblue",c("black"),c("black"),c("black")),cex = 0.8,lty=c(1,1,6,6,2,2), bty = "n")
  
  plot(lineage2_fit[gene_pair[2],],col=c("white"),ylim=c(0,0.1))
  points(as.numeric(DC_data[gene_pair[1],])/10,col=red[5],cex=0.5,pch=16)
  points(as.numeric(DC_data[gene_pair[2],])/10,col=blue[5],cex=0.5,pch=16)
  lines(lineage2_fit[gene_pair[2],],col="darkblue",lwd=3,lty=c("twodash"))
  lines(lineage2_fit[gene_pair[1],],col="darkred",lwd=3,lty=c("twodash"))
  
  plot(lineage3_fit[gene_pair[2],],col=c("white"),ylim=c(0,0.1))
  points(as.numeric(Mono_data[gene_pair[1],])/10,col=red[5],cex=0.5,pch=16)
  points(as.numeric(Mono_data[gene_pair[2],])/10,col=blue[5],cex=0.5,pch=16)
  lines(lineage3_fit[gene_pair[2],],col="darkblue",lwd=3,lty=c("dashed"))
  lines(lineage3_fit[gene_pair[1],],col="darkred",lwd=3,lty=c("dashed"))
  
  text(700,0.08,labels = c("F-value = 458.48"),srt=0,cex=1)
  
  
  gene_pair=c(c("ZBTB16"),c("IKZF2"))
  plot(lineage1_fit[gene_pair[1],],col=c("white"),ylim=c(0,0.3))
  points(as.numeric(Ery_data[gene_pair[1],])/7,col=red[5],cex=0.5,pch=16)
  points(as.numeric(Ery_data[gene_pair[2],])/7,col=blue[5],cex=0.5,pch=16)
  lines(lineage1_fit[gene_pair[1],],col="darkred",lwd=3)
  lines(lineage1_fit[gene_pair[2],],col="darkblue",lwd=3)
  legend("topleft", legend = c(gene_pair[1],gene_pair[2],c("lineage_Ery"),c("lineage_DC"),c("lineage_Mono")), col=c("darkred","darkblue",c("black"),c("black"),c("black")),cex = 0.8,lty=c(1,1,6,6,2,2), bty = "n")
  
  plot(lineage2_fit[gene_pair[2],],col=c("white"),ylim=c(0,0.3))
  points(as.numeric(DC_data[gene_pair[1],])/10,col=red[5],cex=0.5,pch=16)
  points(as.numeric(DC_data[gene_pair[2],])/10,col=blue[5],cex=0.5,pch=16)
  lines(lineage2_fit[gene_pair[2],],col="darkblue",lwd=3,lty=c("twodash"))
  lines(lineage2_fit[gene_pair[1],],col="darkred",lwd=3,lty=c("twodash"))
  
  plot(lineage3_fit[gene_pair[2],],col=c("white"),ylim=c(0,0.3))
  points(as.numeric(Mono_data[gene_pair[1],])/10,col=red[5],cex=0.5,pch=16)
  points(as.numeric(Mono_data[gene_pair[2],])/10,col=blue[5],cex=0.5,pch=16)
  lines(lineage3_fit[gene_pair[2],],col="darkblue",lwd=3,lty=c("dashed"))
  lines(lineage3_fit[gene_pair[1],],col="darkred",lwd=3,lty=c("dashed"))
  
  text(700, 0.25, labels = c("F-value = 79.33"),srt=0,cex=1)
  
  
  
  gene_pair=c(c("KLF1"),c("CEBPA"))
  plot(lineage1_fit[gene_pair[1],],col=c("white"),ylim=c(0,0.8))
  points(as.numeric(Ery_data[gene_pair[1],])/7,col=red[5],cex=0.5,pch=16)
  points(as.numeric(Ery_data[gene_pair[2],])/7,col=blue[5],cex=0.5,pch=16)
  lines(lineage1_fit[gene_pair[1],],col="darkred",lwd=3)
  lines(lineage1_fit[gene_pair[2],],col="darkblue",lwd=3)
  legend("topleft", legend = c(gene_pair[1],gene_pair[2],c("lineage_Ery"),c("lineage_DC"),c("lineage_Mono")), col=c("darkred","darkblue",c("black"),c("black"),c("black")),cex = 0.8,lty=c(1,1,6,6,2,2), bty = "n")
  
  plot(lineage2_fit[gene_pair[2],],col=c("white"),ylim=c(0,0.4))
  points(as.numeric(DC_data[gene_pair[1],])/10,col=red[5],cex=0.5,pch=16)
  points(as.numeric(DC_data[gene_pair[2],])/10,col=blue[5],cex=0.5,pch=16)
  lines(lineage2_fit[gene_pair[2],],col="darkblue",lwd=3,lty=c("twodash"))
  lines(lineage2_fit[gene_pair[1],],col="darkred",lwd=3,lty=c("twodash"))
  
  plot(lineage3_fit[gene_pair[2],],col=c("white"),ylim=c(0,0.6))
  points(as.numeric(Mono_data[gene_pair[1],])/10,col=red[5],cex=0.5,pch=16)
  points(as.numeric(Mono_data[gene_pair[2],])/10,col=blue[5],cex=0.5,pch=16)
  lines(lineage3_fit[gene_pair[2],],col="darkblue",lwd=3,lty=c("dashed"))
  lines(lineage3_fit[gene_pair[1],],col="darkred",lwd=3,lty=c("dashed"))
  
  text(700, 0.5, labels = c("F-value = 8.18"),srt=0,cex=1)
  
  return(results)
}


