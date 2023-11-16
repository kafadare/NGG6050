#function to fix ID format between different files
id_format_fix <- function(v) {
  v <- gsub("^(\\d)", "X\\1", v)
  out <- gsub("-", ".", v)
  return(out)
}

get_gene_info <- function(genid){
  library(biomaRt)
  ensembl <- useEnsembl(biomart = "genes", 
                        dataset = "hsapiens_gene_ensembl")
  bm_gene = getBM(attributes = c('ensembl_gene_id',
                                 'chromosome_name',
                                 'start_position',
                                 'end_position',
                                 'strand'),
                  filters = c('ensembl_gene_id'),
                  mart = ensembl,
                  values = genid)
  colnames(bm_gene) = c('genid','Chr','start','end','strand')
  return(bm_gene)
}

#prep data
prep_datExp <- function(.data, gtf){
  gtf <- gtf %>% filter(!(Chr %in% c("X","Y","M")))
  .data <- .data %>% filter(genid %in% gtf$genid)
  gtf <- gtf %>% filter(genid %in% .data$genid)
  rownames(.data) <- .data$genid
  .data <- .data[,!colnames(.data) %in% "genid"]
  .data <- round(.data)
  return(.data)
}

#function to plot gene density
plot_gen_density <- function(.data, offset = 0.1, title_str = 'Scaled Gene Counts ', xlim = c(-10,30), ylim = c(0,1), log = TRUE){
  if(log == TRUE){
    # View the distribution of expression for each sample.
    # box plot, looking for big differences in read depth (raw counts), symmetry in distribution across samples
    par(mfrow=c(1,2))
    boxplot(.data, range = 0, main = paste(title_str), xlab = 'Samples', xaxt = "n")
    boxplot(log2(offset+.data), range = 0, main = paste('log2(counts+',offset,')'), xlab = 'Samples', xaxt = "n")
    # Histogram/density plot
    # Look for: how well do the distributions line up, outlier samples, zero counts
    par(mfrow=c(1,1))
    i <- 1
    plot <- plot(density(log2(offset+.data[,i])), main = title_str, xlab = paste('log2(counts+',offset,')'), 
                 xlim = xlim, ylim = ylim)
    for(i in 1:ncol(.data)){
      lines(density(log2(.1+.data[,i])), col = i)
    }
  }else{
    # View the distribution of expression for each sample.
    # box plot, looking for big differences in read depth (raw counts), symmetry in distribution across samples
    par(mfrow=c(1,2))
    boxplot(.data, range = 0, main = paste(title_str), xlab = 'Samples', xaxt = "n")
    boxplot(.data, range = 0, main = paste('No Log Counts'), xlab = 'Samples', xaxt = "n")
    # Histogram/density plot
    # Look for: how well do the distributions line up, outlier samples, zero counts
    par(mfrow=c(1,1))
    i <- 1
    plot <- plot(density(.data[,i]), main = title_str,   xlab = 'No Log Counts', 
                 xlim = xlim, ylim = ylim)
    for(i in 1:ncol(.data)){
      lines(density(.data[,i]), col = i)
    }
  }
  plot
  return(plot)
}

#function to remove lowly expressed genes
rm_low <- function(.data,datExpr.tpm, gtf, cutoff = 0.1, percent = 0.25){
  #remove unwanted Chr
  #remove lowly expressed
  # cutoff default 0.1, default 25% subjects
  datExpr.tpm <- datExpr.tpm %>% filter(genid %in% gtf$genid)
  rownames(datExpr.tpm) <- datExpr.tpm$genid
  datExpr.tpm <- datExpr.tpm[,-1]
  # cutoff=0.1, 25% subjects
  keep <- (rowSums(datExpr.tpm > cutoff)) > percent*ncol(datExpr.tpm)
  # # Count filter, use cpm, not raw or scaled from TPM counts!
  # keep <- filterByExpr(datExpr, min.count = (cutoff*100), min.prop=percent)
  datExpr <- .data[keep,]
  return(datExpr)
}

#function to normalize and batch correct
norm_batch <- function(.data, meta, output_dir){
  datExpr <- .data
  # Normalize
  datExpr.vst <- varianceStabilizingTransformation(as.matrix(datExpr), blind = TRUE)
  datExpr.vst <- as.data.frame(datExpr.vst)
  #### 4
  # Remove outliers
  normadj <- adjacency(datExpr.vst,type = 'signed',corFnc = 'bicor')   #Calculate network adjacency
  netsummary <- fundamentalNetworkConcepts(normadj)
  C <- netsummary$Connectivity   #Extract connectivity of each sample
  Z.C <- (C-mean(C))/sqrt(var(C))   #Covert to Z-score
  outliers <- (Z.C < -3)
  datExpr.final <- datExpr.vst
  datExpr.final <- datExpr.final[,!outliers]
  write.table(datExpr.final, paste0(output_dir,"counts.norm.noComBat.tsv"), col.names = T, row.names = T, quote = F, sep = "\t")
  exprMat <- as.matrix(datExpr.final)
  # ComBat
  data.batch <- c()
  ##need to get data.batch
  for (i in 1:ncol(datExpr.final)) {
    sample <- colnames(datExpr.final)[i]
    #match the sample to id in meta file and find corresponding "study"
    data.batch[i] <- meta[which(meta$Subject %in% sample), "study"]
  }
  combat_expr <- ComBat(dat = exprMat, batch = data.batch, mod = NULL, par.prior = TRUE, prior.plots = FALSE)
  combat_expr <- as.data.frame(combat_expr)
  
  write.table(combat_expr, paste0(output_dir, "counts.batch.processed.tsv"), col.names = T, row.names = T, quote = F, sep = "\t")
  
  outlier.df <- data.frame(c(which(outliers)))
  outlier.df$c.which.outliers.. <- rownames(outlier.df)
  write.table(outlier.df,paste0(output_dir, "outlier.txt"), sep="\t", quote=F, col.names = F, row.names = F)
  out <- list(noCombat = datExpr.final, combat = combat_expr, outlier = outlier.df)
  return(out)
}

standardize<- function(X)
{
  X = as.matrix(X)
  # n = dim(X)[1]
  # p = dim(X)[2]
  
  X = scale(X, center = TRUE, scale = TRUE)
  # X = scale(X,center=FALSE, scale=sqrt(apply(X^2,2,sum)))
  
  # m = apply(X,2,mean)
  # st = sqrt(apply(X^2,2,sum));
  # st_mat = matrix(st, nrow = length(st), ncol = dim(X)[2], byrow=FALSE)
  # X2 = X / st_mat
  return (X)
}

