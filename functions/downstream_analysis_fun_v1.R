library(ggrepel)
library(msigdbr)
library(fgsea)
library(clusterProfiler)



gene_comparison_tmaker_adj<-function(genes,ids0,ids1){
  cor.list=c()
  for (i in gg) {
    cor.tcell=cor(sp.counts.norm[,i],rowSums(sp.counts.norm[,tcell_markers]))
  }
  names(cor.list)=gg
  genes.tadj=names(cor.list[which(abs(cor.list)>0.2)])
  genes.noadj=genes[!genes %in% genes.tadj]
  

  #Genes that do not need to adjust t cell marker, just adjust for UMI
  group_0 <- sp.counts.norm[ids0, genes.noadj]
  group_1 <- sp.counts.norm[ids1, genes.noadj]
  
  # Calculate Fold Change and P-values
  results <- data.frame(
    gene = genes.noadj,
    logFC = numeric(length(genes.noadj)),
    p_value = numeric(length(genes.noadj))
  )
  
  for (j in 1:ncol(dat)) {
    # Get gene expression for this gene
    gene_0 <- group_0[, j]
    gene_1 <- group_1[, j]
    print(j)
    
    # Perform t-test
    ttest <- t.test(gene_0,gene_1)
    m0=t.test$estimates[1]
    m1=t.test$estimates[2]
    
    if (mean(gene_0)!=0) {logFC <- log2(mean(gene_1) / mean(gene_0))}

    # Store results
    results$logFC[j] <- logFC
    results$p_value[j] <- ttest$p.value
  }
  results1=results
  
  
  #Genes that associated with t cell marker, just adjust for sum of t cell markers
  
  results <- data.frame(
    gene = genes.adj,
    logFC = numeric(length(genes.adj)),
    p_value = numeric(length(genes.adj))
  )
  rownames(results)=genes.adj
  xx=cbind(c(rep(0,length(ids0)),rep(1,length(ids1))),rowSums(sp.counts.r[c(ids0,ids1),tcell_markers]))
  for (i in genes.adj) {
    YY=sp.counts.r[c(ids0,ids1),i]
    tmp.dat = data.frame(trr = xx[,1],tcell_sum=xx[,2], y = YY)
    mm = glm.nb(y ~ tcell_sum+trr, data = tmp.dat)
    beta=summary(mm)$coef['trr','Estimate']
    logFC=log2(exp(beta))
    p_value=summary(mm)$coef['trr','Pr(>|z|)']
    
    # Store results
    results[i,"logFC"] <- logFC
    results[i,"p_value"] <- p_value
  }
  
  results2=results
  return(list(results1,results2))
}




gene_comparison_tmaker_adj_all<-function(genes,ids0,ids1,umi.adj=TRUE){
  
  results <- data.frame(
    gene = genes,
    logFC = numeric(length(genes)),
    p_value = numeric(length(genes))
  )
  rownames(results)=genes
  xx=cbind(c(rep(0,length(ids0)),rep(1,length(ids1))),meta_dat[c(ids0,ids1),'nCount_Spatial'],rowSums(sp.counts.r[c(ids0,ids1),tcell_markers]))
  
  for (i in genes) {
    YY=sp.counts.r[c(ids0,ids1),i]
    tmp.dat = data.frame(trr = xx[,1],umi=xx[,2],tcell_sum=xx[,3], y = YY)
    if (umi.adj==TRUE) {mm = glm.nb(y ~ tcell_sum+umi+trr, data = tmp.dat)} else {mm = glm.nb(y ~ tcell_sum+trr, data = tmp.dat)}

    beta=summary(mm)$coef['trr','Estimate']
    logFC=log2(exp(beta))
    p_value=summary(mm)$coef['trr','Pr(>|z|)']
    
    # Store results
    results[i,"logFC"] <- logFC
    results[i,"p_value"] <- p_value
  }
  
  results$adj_p_value <- p.adjust(results$p_value, method = "BH")
  results$log_adj_p_value <- -log10(results$adj_p_value)
  
  results$adj_p_value=round(results$adj_p_value,4)
  return(results)
}


gene_comparison_umi_adj<-function(genes,ids0,ids1){
  #Genes that do not need to adjust t cell marker, just adjust for UMI
  group_0 <- sp.counts.norm[ids0, genes]
  group_1 <- sp.counts.norm[ids1, genes]
  
  # Calculate Fold Change and P-values
  results <- data.frame(
    gene = genes,
    logFC = numeric(length(genes)),
    p_value = numeric(length(genes))
  )
  
  for (j in 1:length(genes)) {
    # Get gene expression for this gene
    gene_0 <- group_0[, j]
    gene_1 <- group_1[, j]
    print(j)
    
    # Perform t-test
    ttest <- t.test(gene_0,gene_1)
    m0=ttest$estimate[1]
    m1=ttest$estimate[2]
    logFC <- log2(m1/m0)
    # Store results
    results$logFC[j] <- logFC
    results$p_value[j] <- ttest$p.value
  }
  
  results$adj_p_value <- p.adjust(results$p_value, method = "BH")
  results$log_adj_p_value <- -log10(results$adj_p_value)
  
  results$adj_p_value_round=round(results$adj_p_value,4)
  return(results)
  
}




gene_comparison_tpm<-function(genes,ids0,ids1){
  #Genes that do not need to adjust t cell marker, just adjust for UMI
  group_0 <- sp.counts.tpm[ids0, genes]
  group_1 <- sp.counts.tpm[ids1, genes]
  
  # Calculate Fold Change and P-values
  results <- data.frame(
    gene = genes,
    logFC = numeric(length(genes)),
    p_value = numeric(length(genes))
  )
  
  for (j in 1:length(genes)) {
    # Get gene expression for this gene
    gene_0 <- group_0[, j]
    gene_1 <- group_1[, j]
    print(j)
    
    # Perform t-test
    ttest <- t.test(gene_0,gene_1)
    m0=ttest$estimate[1]
    m1=ttest$estimate[2]
    logFC <- log2(m1/m0)
    # Store results
    results$logFC[j] <- logFC
    results$p_value[j] <- ttest$p.value
  }
  
  results$adj_p_value <- p.adjust(results$p_value, method = "BH")
  results$log_adj_p_value <- -log10(results$adj_p_value)
  
  results$adj_p_value_round=round(results$adj_p_value,4)
  return(results)
  
}



compare_with_tumor<-function(sp.dat=sp.counts.norm,tumor.loc=tumor.loc, cut.off= -1){
  dat=sp.dat[rownames(tumor.loc),]
  TRR=rep(0,nrow(dat))
  dat=cbind(dat,TRR)
  dat[which(rownames(tumor.loc) %in% ids1),'TRR']=1
  
  tmp.mat=c()
  for (i in 1:ncol(sp.dat)) {
    ttest=t.test(dat[which(rownames(tumor.loc) %in% ids1),colnames(sp.dat)[i]],
                 dat[which(!rownames(tumor.loc) %in% ids1),colnames(sp.dat)[i]])
    m1=ttest$estimate[1]
    m0=ttest$estimate[2]
    out=c(m1,m0,ttest$p.value)
    tmp.mat=rbind(tmp.mat,out)
  }
  colnames(tmp.mat)=c('mean_TRR','mean_other',"pvalue")
  rownames(tmp.mat)=colnames(sp.dat)
  tmp.mat=tmp.mat[order(tmp.mat[,'pvalue']),]
  tmp.mat=data.frame(tmp.mat)
  tmp.mat$adj_p_value=p.adjust(tmp.mat$pvalue,method='BH')
  tmp.mat=data.frame(tmp.mat)
  tmp.mat$logFC=log2(tmp.mat$mean_TRR/tmp.mat$mean_other)
  
  low.mat=tmp.mat[which(tmp.mat$logFC< cut.off & tmp.mat$adj_p_value<0.05),]
  
  return(low.mat)
}
