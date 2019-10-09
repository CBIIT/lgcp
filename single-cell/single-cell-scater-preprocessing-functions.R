cor_dist <- function(som_codes, method){
  corr_mat <- cor(t(som_codes),
                  method = method)
  return(1 - corr_mat)
}

binomial_deviance <- function(x,p,n){
  term1<-sum(x*log(x/(n*p)), na.rm=TRUE)
  nx<-n-x
  term2<-sum(nx*log(nx/(n*(1-p))), na.rm=TRUE)
  2*(term1+term2)
}

compute_gene_info <- function(m,gmeta=NULL,mod=c("binomial","multinomial","poisson","geometric")){
  #m a data matrix with genes=rows
  #gmeta a pre-existing data frame with gene-level metadata
  mod<-match.arg(mod)
  if(!is.null(gmeta)){ stopifnot(nrow(m)==nrow(gmeta)) }
  gnz<-Matrix::rowSums(m>0)
  sz<-compute_size_factors(m,mod)
  gof<-function(g){ gof_func(m[g,],sz,mod) }
  gof<-as.data.frame(t(vapply(1:nrow(m),gof,FUN.VALUE=rep(0.0,2))))
  #colnames(gof)<-c("deviance","pval")
  gmu<-Matrix::rowMeans(m)
  gvar<-apply(m,1,var)
  gfano<-ifelse(gvar>0 & gmu>0, gvar/gmu, 0)
  res<-cbind(nzsum=gnz,fano=gfano,gof)
  res$pval_fdr<-p.adjust(res$pval,"BH")
  if(is.null(gmeta)){ return(res) } else { return(cbind(gmeta,res)) }
}

compute_size_factors <- function(m,mod=c("binomial","multinomial","poisson","geometric")){
  #given matrix m with samples in the columns
  #compute size factors suitable for the discrete model in 'mod'
  mod<-match.arg(mod)
  sz<-Matrix::colSums(m) #base case, multinomial or binomial
  if(mod %in% c("multinomial","binomial")){ return(sz) }
  sz<-log(sz)
  sz<-sz - mean(sz) #make geometric mean of sz be 1 for poisson, geometric
  if(mod=="poisson"){ return(exp(sz)) }
  sz #geometric, use log scale size factors
}

filterDev <- function(sce,nkeep=nrow(sce),dev=c("binomial","poisson","geometric"),ret=c("sce","ranks")){
  dev<-match.arg(dev)
  ret<-match.arg(ret)
  gm<-compute_gene_info(counts(sce),gmeta=rowData(sce),mod=dev)
  o<-order(gm$deviance,decreasing=TRUE,na.last=FALSE)[1:nkeep]
  #NA deviance => badly fitting null model=> highly variable gene
  if(ret=="sce"){
    res<-sce[o,]
    return(res[,colSums(counts(res))>0])
  } else {
    return(rownames(sce)[o])
  }
}

filterExpr <- function(sce, nkeep=nrow(sce),ret=c("sce","ranks")){
  #modified from https://github.com/markrobinsonuzh/scRNAseq_clustering_comparison/blob/master/Rscripts/filtering/filterExpr.R
  ret<-match.arg(ret)
  exprsn<-rowMeans(logcounts(sce))
  o<-order(exprsn,decreasing=TRUE)[1:nkeep]
  if(ret=="sce"){
    res<-sce[o,]
    return(res[,colSums(counts(res))>0])
  } else {
    return(rownames(sce)[o])
  }
}

filterHVG <- function(sce, nkeep=nrow(sce), total_umi="nUMI", ret=c("sce","ranks")){
  #modified from https://github.com/markrobinsonuzh/scRNAseq_clustering_comparison/blob/master/Rscripts/filtering/filterHVG.R
  ret<-match.arg(ret)
  if(!(total_umi %in% colnames(colData(sce)))){ stop("total_umi must be in coldata") }
  seu<-CreateSeuratObject(counts(sce),meta.data=as.data.frame(colData(sce)),min.cells=0,min.features=0,project="scRNAseq")
  seu<-NormalizeData(seu,verbose=FALSE)
  seu<-ScaleData(seu,vars.to.regress=total_umi,verbose=FALSE)
  seu<-FindVariableFeatures(seu,selection.method="dispersion",nfeatures=nkeep,verbose=FALSE)
  vf<-head(VariableFeatures(seu), nkeep)
  #sometimes Seurat coerces rownames, so use numeric IDs instead
  o<-match(vf,rownames(seu))
  if(ret=="sce"){
    res<-sce[o, ]
    return(res[,colSums(counts(res))>0])
  } else {
    return(rownames(sce)[o])
  }
}

gof_func <- function(x,sz,mod=c("binomial","multinomial","poisson","geometric")){
  #Let n=colSums(original matrix where x is a row)
  #if binomial, assumes sz=n, required! So sz>0 for whole vector
  #if poisson, assumes sz=n/geometric_mean(n), so again all of sz>0
  #if geometric, assumes sz=log(n/geometric_mean(n)) which helps numerical stability. Here sz can be <>0
  #note sum(x)/sum(sz) is the (scalar) MLE for "mu" in Poisson and "p" in Binomial
  mod<-match.arg(mod)
  fit<-list(deviance=0,df.residual=length(x)-1,converged=TRUE)
  if(mod=="multinomial"){
    fit$deviance<-multinomial_deviance(x,sum(x)/sum(sz))
  } else if(mod=="binomial"){
    fit$deviance<-binomial_deviance(x,sum(x)/sum(sz),sz)
  } else if(mod=="poisson"){
    fit$deviance<-poisson_deviance(x,sum(x)/sum(sz),sz)
  } else if(mod=="geometric"){
    if(any(x>0)) {
      fit<-glm(x~offset(sz),family=MASS::negative.binomial(theta=1))
    }
  } else { stop("invalid model") }
  if(fit$converged){
    dev<-fit$deviance
    df<-fit$df.residual #length(x)-1
    pval<-pchisq(dev,df,lower.tail=FALSE)
    res<-c(dev,pval)
  } else {
    res<-rep(NA,2)
  }
  names(res)<-c("deviance","pval")
  res
}

multinomial_deviance <- function(x,p){
  -2*sum(x*log(p))
}

poisson_deviance <- function(x,mu,sz){
  #assumes log link and size factor sz on the same scale as x (not logged)
  #stopifnot(all(x>=0 & sz>0))
  2*sum(x*log(x/(sz*mu)),na.rm=TRUE)-2*sum(x-sz*mu)
}

rank_all_genes <- function(sce, total_umi="nUMI"){
  #sce=SingleCellExperiment object with assays "counts" and "logcounts"
  #returns a dataframe with same rownames as sce
  #columns: dev,hvg,expr
  #each column is the rank order of genes by each criteria
  #low rank=highly informative gene
  gg<-rownames(sce)
  gfs<-list()
  gfs$dev<-filterDev(sce,ret="ranks")
  gfs$expr<-filterExpr(sce,ret="ranks")
  gfs$hvg<-filterHVG(sce,total_umi=total_umi,ret="ranks")
  res<-list()
  for(i in names(gfs)){
    rk<-seq_along(gfs[[i]])
    names(rk)<-gfs[[i]]
    res[[i]]<-rk[gg]
  }
  res<-as.data.frame(res)
  rownames(res)<-gg
  res
}

tidy_sparse_matrix <- function (x) 
{
  s <- Matrix::summary(x)
  row <- s$i
  if (!is.null(rownames(x))) {
    row <- rownames(x)[row]
  }
  col <- s$j
  if (!is.null(colnames(x))) {
    col <- colnames(x)[col]
  }
  ret <- data.frame(row = row, column = col, value = s$x, stringsAsFactors = FALSE)
  ret
}
